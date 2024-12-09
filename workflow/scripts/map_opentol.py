import sqlite3
import requests
import json
import util
import os
import argparse
import pandas as pd
import numpy as np


"""
This script is responsible for mapping taxonomic names from a SQLite database to the Open Tree of Life (OpenTOL) using 
the TNRS API.

The script performs the following steps:
1. Connects to the SQLite database.
2. Loads all unmatched records from the 'taxon' table into a DataFrame.
3. Iterates over the DataFrame in chunks and checks all names within each chunk against the TNRS API.
4. Processes the results from the TNRS API. If there is a single match that is not a synonym and whose score is high 
   enough, the script updates the 'taxon' table with the OpenTOL ID.
5. If the initial exact match fails, the script attempts to recover through a fuzzy match.
6. After all matches have been processed, the script creates an index on the 'opentol_id' column in the 'taxon' table 
   for optimized queries.
7. Commits the changes to the database and closes the connection.

The script uses command line arguments for the database file to enrich, marker name, output status file, and log level.
The script is invoked by the Snakefile as a shell command with the required arguments in the rule 'map_opentol'.
"""

# Set the Open Tree of Life API endpoint
# https://github.com/OpenTreeOfLife/germinator/wiki/TNRS-API-v3#match_names
endpoint = "https://api.opentreeoflife.org/v3/tnrs/match_names"

# TODO compare matches counts b/w this and previous implementation


def tnrs_request(kingdom, names, fuzzy=False):
    """
    Performs a taxonomic name resolution (TNRS) request on the OpenToL API for a list of names . The service does this
    within a given taxonomic context, to specify which code applies (botany or zoology).
    :param kingdom: Taxonomic context, 'Animals' or 'Plants'
    :param names: Taxon name to match
    :param fuzzy: Boolean, if True does approximate matching
    :return: A result object with the matches, parsed from JSON output
    """

    tnrs_params = {
        "names": names,
        "do_approximate_matching": fuzzy,
        "verbose": False,
        "context": kingdom
    }
    json_response = requests.post(endpoint, json=tnrs_params)
    return json.loads(json_response.text)


def tnrs_response(matches, score=1.0):
    """
    Processes a TNRS response. Accepts responses if and only if there is a single match that is not a synonym whose
    score is greater than or equal to the score parameter. Returns the OTT ID if successful, None otherwise.
    :param matches: The data structure resulting from the parsed JSON returned by the TNRS service
    :param score: A value between 0 and 1 (inclusive)
    :return: OTT ID or None
    """

    # There are results
    if len(matches) != 0:

        # There should be only one result that is not a synonym and whose score is high enough
        filtered = [d for d in matches if
                    not d['is_synonym'] and d['score'] >= score and d['taxon']['rank'] == 'species']
        if len(filtered) == 1:
            return filtered[0]['taxon']['ott_id']
        else:
            logger.debug("No single match: %s", filtered)

    # Either there was no result, or there were too many ambiguous results
    return None


def match_opentol(kingdom, chunksize, fuzzy):
    """
    Matches names against the OpenTOL. Does this for all records in the taxon table that don't have an opentol_id.
    Attempts first an exact match, then tries to recover through a fuzzy match if initially there was no match.
    Updates the database records with matches.
    :param kingdom: The applicable kingdom to search within.
    :param chunksize: How many names to check in batch. Literal<=10000, Fuzzy<=250
    :param fuzzy: If True do fuzzy matching, otherwise literal matching
    :return:
    """

    # Load all unmatched records into df, iterate over it in chunks
    df = pd.read_sql("SELECT * FROM taxon WHERE opentol_id IS NULL", conn)
    for _, chunk_df in df.groupby(np.arange(len(df)) // chunksize):

        # Check all names within the chunk as a list against TNRS
        tnrs_data = tnrs_request(kingdom, chunk_df['species'].to_list(), fuzzy=fuzzy)

        # Process the results, update temp table
        for i in range(len(tnrs_data['results'])):

            # Prepare record to insert, coerce types
            taxon_id = int(chunk_df['taxon_id'].iloc[i])

            # Attempt to get OTT ID
            ott_id = tnrs_response(tnrs_data['results'][i]['matches'])
            if ott_id is not None:
                opentol_id = int(ott_id)
                cursor.execute(f'UPDATE taxon SET opentol_id = {opentol_id} WHERE taxon_id = {taxon_id}')


def postprocess_db():
    # create indexes on opentol_id
    logger.info('Going to index opentol_id and analyze the database for optimized queries')
    cursor.execute('CREATE INDEX IF NOT EXISTS taxon_opentol_id_idx ON taxon (opentol_id)')
    cursor.execute('ANALYZE')


if __name__ == '__main__':

    # Define and process command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-d', '--database', required=True, help='Database file to enrich')
    parser.add_argument('-m', '--marker', required=True, help='Marker name (e.g. COI-5P)')
    parser.add_argument('-o', '--output', required=True, help='0-byte status file')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Configure logger
    logger = util.get_formatted_logger('map_opentol', args.verbosity)

    # Connect to database
    logger.info(f'Going to connect to database {args.database}')
    conn = sqlite3.connect(args.database)
    cursor = conn.cursor()
    cursor.execute('pragma journal_mode=OFF')
    cursor.execute('PRAGMA synchronous=OFF')
    cursor.execute('PRAGMA cache_size=100000')
    cursor.execute('PRAGMA temp_store = MEMORY')

    # Infer taxonomic context from marker name
    if args.marker == "COI-5P":
        focal_kingdom = 'Animals'
    else:
        focal_kingdom = "Plants"

    # Literal matching in big steps
    logger.info("Going to match literally in chunks of 10,000 names")
    match_opentol(focal_kingdom, chunksize=10000, fuzzy=False)

    # Fuzzy matching in small steps
    logger.info("Going to match fuzzily in chunks of 250 names")
    match_opentol(focal_kingdom, chunksize=250, fuzzy=True)

    # Compute indexes and write to final file
    postprocess_db()
    conn.commit()
    conn.close()

    # Touch the file
    filename = args.output
    with open(filename, 'a'):
        os.utime(filename, None)
