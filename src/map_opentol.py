import argparse
import os
from io import StringIO
import numpy as np
import pandas as pd
import sqlite3
import requests

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()

parser.add_argument('-db', default="BOLD_COI_barcodes.db",
                    help="Name of the the database file: {file_name}.db")
args = parser.parse_args()


def map_checklistbank(conn, cursor):
    # TODO Fix Checklistbank request error (405,HTTP 405 Method Not Allowed)
    """
    This function maps BOLD taxon names from the custom database to an
    Open Tree of Life ID. These IDs are accessed with a request to the
    CheckListBank API, which needs a csv containing (atleast) a column named
    scientificName which hold taxon names. The output is saved as a temporary
    table in the database.
    :param conn: Connection to SQLite database
    :param cursor: Cursor for SQLite database
    """
    results = cursor.execute("""SELECT taxon FROM taxon""").fetchall()

    # Change column name from taxon to scientificName
    df = pd.DataFrame(results, columns=['scientificName'])

    # Split dataframe into batches and put them in a list
    list_df = np.array_split(df, 50)

    # CheckListBank request url
    url = "https://api.checklistbank.org/dataset/201891/nameusage/match"

    # Make headers
    headers = {
        'Accept': 'text,csv',
        'Content-Type': 'text/tsv',
    }
    ## BACK END FETCH FAILED eventhough this worked a week before, try later
    ## ALSO DOESNT WORK IN COMMANDLINE

    # Loop trough list of df batches and make a
    # for i in list_df:
    #     data = i.to_csv(index=False)
    #     response = requests.post(url, headers=headers,
    #                          data=data)
    #     content = response.content.decode('utf-8')
    #     print(pd.read_csv(StringIO(content)))

    #     # Put response in dataframe format
    # df_batch = pd.read_csv(StringIO(content), sep=',',
    #                        usecols=['ID', 'inputName'])
    # Rename columns names
    # df_batch.rename(columns={'ID':'opentol_id'})

    # Append rows to sql table 'opentol_temp', create if not exists yet
    # df_batch.to_sql('opentol_temp', conn, if_exists='append')

    ## REMOVE WHEN REQUEST WORKS, uses old response from request (around 10.000
    # taxons as test data)
    pd.read_csv('match_test.csv', sep='\t').to_sql('opentol_temp', conn,
                                                       if_exists='append')
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
    alter_db(conn, cursor)

    # Close the connection
    conn.close()
