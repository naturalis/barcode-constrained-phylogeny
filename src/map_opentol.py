import sqlite3
import requests
import json
import logging
import os
import pandas as pd
import numpy as np

# Set the Open Tree of Life API endpoint
# https://github.com/OpenTreeOfLife/germinator/wiki/TNRS-API-v3#match_names
endpoint = "https://api.opentreeoflife.org/v3/tnrs/match_names"

# Instantiate logger
logging.basicConfig(level=snakemake.params.log_level)  # noqa: F821
logger = logging.getLogger(__name__)


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


def match_opentol(database, kingdom, chunksize, fuzzy):
    """
    Matches names against the OpenTOL. Does this for all records in the taxon table that don't have an opentol_id.
    Attempts first an exact match, then tries to recover through a fuzzy match if initially there was no match.
    Updates the database records with matches.
    :param database: Name of SQLite database
    :param kingdom: The applicable kingdom to search within.
    :param chunksize: How many names to check in batch. Literal<=10000, Fuzzy<=250
    :param fuzzy: If True do fuzzy matching, otherwise literal matching
    :return:
    """

    # Try all unmatched names in database in chunks
    conn = sqlite3.connect(database, isolation_level=None)
    cursor = conn.cursor()

    # Create temporary table
    cursor.execute("""CREATE TABLE IF NOT EXISTS taxon_temp (
        taxon_id INTEGER PRIMARY KEY,
        taxon TEXT,
        kingdom TEXT NOT NULL,
        class TEXT NOT NULL,
        family TEXT NOT NULL,
        genus TEXT NOT NULL,
        opentol_id INTEGER)""")

    # Load all unmatched records into df, iterate over it in chunks
    df = pd.read_sql("SELECT * FROM taxon WHERE opentol_id IS NULL", conn)
    for _, chunk_df in df.groupby(np.arange(len(df)) // chunksize):

        # Check all names within the chunk as a list against TNRS
        tnrs_data = tnrs_request(kingdom, chunk_df['taxon'].to_list(), fuzzy=fuzzy)

        # Process the results, update temp table
        for i in range(len(tnrs_data['results'])):

            # Prepare record to insert, coerce types
            taxon_id = int(chunk_df['taxon_id'].iloc[i])
            taxon = str(chunk_df['taxon'].iloc[i])
            kingdom = str(chunk_df['kingdom'].iloc[i])
            classe = str(chunk_df['class'].iloc[i])
            family = str(chunk_df['family'].iloc[i])
            genus = str(chunk_df['genus'].iloc[i])

            # Attempt to get OTT ID
            ott_id = tnrs_response(tnrs_data['results'][i]['matches'])
            if ott_id is not None:
                opentol_id = int(ott_id)
                data_to_insert = (taxon_id, taxon, kingdom, classe, family, genus, opentol_id)
                sql_command = """
                    INSERT INTO taxon_temp (taxon_id, taxon, kingdom, class, family, genus, opentol_id)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                """
                cursor.execute(sql_command, data_to_insert)

            else:
                data_to_insert = (taxon_id, taxon, kingdom, classe, family, genus)
                sql_command = """
                    INSERT INTO taxon_temp (taxon_id, taxon, kingdom, class, family, genus)
                    VALUES (?, ?, ?, ?, ?, ?)
                """
                cursor.execute(sql_command, data_to_insert)
            conn.commit()
            logger.info(data_to_insert)

    # Copy over the already matched records
    cursor.execute("""INSERT INTO taxon_temp SELECT * FROM taxon WHERE opentol_id IS NOT NULL""")

    # Remove old taxon table
    cursor.execute("""DROP TABLE taxon""")

    # Rename new taxon table from temp to taxon
    cursor.execute("""ALTER TABLE taxon_temp RENAME TO taxon""")

    conn.close()


def postprocess_db(database):
    # create indexes on family, genus, opentol_id, taxon_id country and nucraw
    conn = sqlite3.connect(database, isolation_level=None)
    cursor = conn.cursor()
    cursor.execute("""CREATE INDEX IF NOT EXISTS class_idx ON taxon (class)""")
    cursor.execute("""CREATE INDEX IF NOT EXISTS family_idx ON taxon (family)""")
    cursor.execute("""CREATE INDEX IF NOT EXISTS genus_idx ON taxon (genus)""")
    cursor.execute("""CREATE INDEX IF NOT EXISTS opentol_id_idx ON taxon (opentol_id)""")
    cursor.execute("""CREATE INDEX IF NOT EXISTS taxon_id_idx ON barcode (taxon_id)""")
    cursor.execute("""CREATE INDEX IF NOT EXISTS country_idx ON barcode (country)""")
    cursor.execute("""CREATE INDEX IF NOT EXISTS nucraw_idx ON barcode (nucraw)""")
    conn.commit()
    conn.close()


if __name__ == '__main__':
    temp_database_name = snakemake.input[0]  # noqa: F821
    marker = snakemake.params.marker  # noqa: F821
    database_file = snakemake.output[0]  # noqa: F821

    # Infer taxonomic context from marker name
    if marker == "COI-5P":
        kingdom = 'Animals'
    else:
        kingdom = "Plants"

    # Literal matching in big steps
    logger.info("Going to match literally in chunks of 10,000 names")
    match_opentol(temp_database_name, kingdom, chunksize=10000, fuzzy=False)

    # Fuzzy matching in small steps
    logger.info("Going to match fuzzily in chunks of 250 names")
    match_opentol(temp_database_name, kingdom, chunksize=250, fuzzy=True)

    # Compute indexes and write to final file
    postprocess_db(temp_database_name)
    os.rename(temp_database_name, database_file)
