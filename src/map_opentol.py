import argparse
import os
import numpy as np
import pandas as pd
import sqlite3

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()

parser.add_argument('-db', default="BOLD_COI_barcodes.db",
                    help="Name of the the database file: {file_name}.db")
args = parser.parse_args()
import requests
import json


def map_checklistbank(conn):
    # Change column name from taxon to scientificName
    df = pd.read_sql("SELECT * FROM taxon", conn)

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
        .rename(columns={'ott_id': 'opentol_id'})

    # Drop duplicate rows
    df_merged.drop_duplicates(ignore_index=True, inplace=True)

    # Safe df to a table in the db
    df_merged.to_sql('opentol_temp', conn, if_exists='replace', index=False)

    conn.commit()


def map_checklistbank_fuzzy():
    df = pd.read_sql("SELECT * FROM opentol_temp WHERE opentol_id IS NULL",
                     conn)

    # Set the Open Tree of Life API endpoint
    endpoint = "https://api.opentreeoflife.org/v3"

    # Split df into chunks
    list_df = np.array_split(df, 50)

    ott_ids = []
    names = []
    new_taxons = []

    # Loop through chunks from the df
    for taxons in list_df:
        tnrs_url = f"{endpoint}/tnrs/match_names"
        # TODO instead of Animals check if its animalia or plantea kingdom and put appropiate context
        tnrs_params = {
                        "names": taxons['taxon'].to_list(),
                        "do_approximate_matching": True,
                        "verbose": False,
                        "context": 'Animals'
                    }
        tnrs_response = requests.post(tnrs_url, json=tnrs_params)
        tnrs_data = json.loads(tnrs_response.text)

        # Loop trough response results
        for result in tnrs_data['results']:
            # Check if there was a match found
            if len(result['matches']) != 0:
                # Extract the ott_id, bold name and opentol ID  values from the response
                ott_id = result['matches'][0]['taxon']['ott_id']
                name = result['matches'][0]['matched_name']
                new_taxon = result['matches'][0]['taxon']['unique_name']

                # Shows a lsit of flags like incertae_sedis, not used now but might be later
                flags = result['matches'][0]['taxon']['flags']

                # Append the values to the corresponding lists
                ott_ids.append('ott' + str(ott_id))
                names.append(name)
                new_taxons.append(new_taxon)

    # Create a Pandas DataFrame from the lists of values
    df_ot = pd.DataFrame({'ott_id': ott_ids, 'taxon': names, 'new_taxon': new_taxons})

    # Left join old table with new table containing opentol IDs
    df_old = pd.read_sql("SELECT * FROM opentol_temp", conn)

    df_merged = pd.merge(df_old, df_ot, left_on='taxon', right_on='taxon',
                         how='left')

    # If opentol_id is empty, replace with ott_id
    df_merged['opentol_id'] = df_merged['opentol_id'].fillna(df_merged['ott_id'])

    # If opentol_id and ott_id are the same, safe mask
    mask = df_merged['opentol_id'] == df_merged['ott_id']

    # Replace values in column taxon with the new_taxon value, where mask is True
    df_merged.loc[mask, 'taxon'] = df_merged.loc[mask, 'new_taxon']

    df_merged.drop(columns=['new_taxon', 'ott_id'], inplace=True)

    df_merged.to_sql('temp', conn, if_exists='replace', index=False)
    conn.commit()


def alter_db(conn, cursor):
    # Remove temporary table
    cursor.execute("""DROP TABLE opentol_temp""")

    # Remove old taxon table
    cursor.execute("""DROP TABLE taxon""")

    # Rename new taxon table from temp to taxon
    cursor.execute("""ALTER TABLE temp RENAME TO taxon""")

    # Commit changes
    conn.commit()

if __name__ == '__main__':
    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect(args.db)
    # Create a cursor
    cursor = conn.cursor()

    print("Looking for exact matches between BOLD and OpenTOL taxon names...")
    map_checklistbank(conn, cursor)
    print("Fuzzy matching remainder taxons...")
    map_checklistbank_fuzzy()

    # Alter new table and remove old tables
    alter_db(conn, cursor)

    # Close the connection
    conn.close()
