import argparse
import os
import numpy as np
import pandas as pd
import sqlite3

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()

parser.add_argument('-db', default="BOLD_COI-5P_barcodes.db",
                    help="Name of the the database file: {file_name}.db")
args = parser.parse_args()

def map_checklistbank(conn, cursor):
    #TODO Option: Use TNRS of OpenTOL to map the rest of the taxon names
    #TODO When it is working, drop taxon table and rename table opentol_temp to taxon
    """
    This function maps BOLD taxon names from the custom database to an
    Open Tree of Life ID. These IDs are accessed with a request to the
    Open tree of life API, which needs a list of taxon names as paramater. The output is saved as a temporary
    table in the database.
    :param conn: Connection to SQLite database
    :param cursor: Cursor for SQLite database
    """
    # # Change column name from taxon to scientificName
    df = pd.read_sql("SELECT * FROM taxon", conn)
    import requests
    import json

    # Set the Open Tree of Life API endpoint
    endpoint = "https://api.opentreeoflife.org/v3"

    # Split df into chunks
    list_df = np.array_split(df, 50)

    ott_ids = []
    names = []
    for taxons in list_df:
        tnrs_url = f"{endpoint}/tnrs/match_names"
        # TODO instead of Animals check if its animalia or plantea kingdom and put appropiate context
        tnrs_params = {
                        "names": taxons['taxon'].to_list(),
                        "do_approximate_matching": False,
                        "verbose": False,
                        "context": 'Animals'
                    }
        tnrs_response = requests.post(tnrs_url, json=tnrs_params)
        tnrs_data = json.loads(tnrs_response.text)

        # Get the first OTT ID returned from the search
        for result in tnrs_data['results']:
            if len(result['matches']) != 0:
                # Extract the ott_id and name values from the 'taxon' dictionary
                ott_id = result['matches'][0]['taxon']['ott_id']
                name = result['matches'][0]['taxon']['unique_name']
                # Append the values to the corresponding lists
                ott_ids.append('ott' + str(ott_id))
                names.append(name)

    # Create a Pandas DataFrame from the lists of values
    df_ot = pd.DataFrame({'ott_id': ott_ids, 'taxon': names})

    # Left join old table with new table containing opentol IDs
    df_merged = pd.merge(df, df_ot, left_on='taxon', right_on='taxon', how='left').drop(columns=['opentol_id'])\
        .rename(columns={'ott_id': ' opentol_id'})

    # Drop duplicate rows
    df_merged.drop_duplicates(ignore_index=True, inplace=True)

    # Safe df to a table in the db
    df_merged.to_sql('opentol_temp', conn, if_exists='replace', index=False)

    conn.commit()


def alter_db(conn, cursor):
    """
    Creates a new tables named taxonomy and populates it with the columns of two
    tables combined. It takes the taxon table and the opentol_temp table and
    joins it so that all opentol_id's are put on taxons that was returned after
    the CheckListBank request. The tables opentol_temp is deleted and so is
    taxon, the new taxonomy table is then renamed to taxon.
    :param conn: Connection to SQLite database
    :param cursor: Cursor for SQLite database
    """
    # Create new table
    cursor.execute("""CREATE TABLE taxonomy(
                      taxon_id INTEGER,
                      taxon TEXT,
                      family TEXT,
                      kingdom TEXT,
                      opentol_id TEXT)""")

    # Left join taxon and opentol_temp
    cursor.execute("""INSERT INTO taxonomy SELECT taxon.taxon_id, 
                      taxon.taxon, taxon.family, taxon.kingdom, opentol_temp.ID
                      FROM taxon LEFT JOIN opentol_temp ON
                      opentol_temp.inputName = taxon.taxon""")

    # Remove temporary table
    cursor.execute("""DROP TABLE opentol_temp""")

    # Remove old taxon table
    cursor.execute("""DROP TABLE taxon""")

    # Rename new taxon table from taxonomy to taxon
    cursor.execute("""ALTER TABLE taxonomy RENAME TO taxon""")

    conn.commit()


if __name__ == '__main__':
    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect(args.db)
    # Create a cursor
    cursor = conn.cursor()

    # Make request to checklistbank and map BOLD taxons to opentol id's
    map_checklistbank(conn, cursor)

    # Alter new table and remove old tables
    # alter_db(conn, cursor)

    # Close the connection
    conn.close()
