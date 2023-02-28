from io import StringIO
import numpy as np
import pandas as pd
import sqlite3
import requests

def map_checklistbank(conn, cursor):
    count = cursor.execute("""SELECT COUNT(taxon) FROM taxon""").fetchall()[0]
    results = cursor.execute("""SELECT taxon FROM taxon""").fetchall()
    df = pd.DataFrame(results, columns=['scientificName'])
    list_df = np.array_split(df, 50)
    url = "https://api.checklistbank.org/dataset/201891/nameusage/match"
    headers = {
        'Accept': 'text,csv',
        'Content-Type': 'text/tsv',
    }
    ## BACK END FETCH FAILED eenthough this worked a week before, try later
    ## ALSO DOESNT WORK IN COMMANDLINE
    # for i in list_df:
    #     data = i.to_csv(index=False)
    #     response = requests.post(url, headers=headers,
    #                          data=data)
    #     content = response.content.decode('utf-8')
    #     print(pd.read_csv(StringIO(content)))
    # f = pd.read_csv(StringIO(content), sep=',', usecols=['ID', 'inputName'])
    # f.rename(columns={'ID':'opentol_id'})
    # f.to_sql('opentol_temp', conn, if_exists='append')

    ## REMOVE WHEN REQUEST WORKS
    df = pd.read_csv('match_test.csv', sep='\t').to_sql('opentol_temp', conn,
                                                       if_exists='append')
    conn.commit()

def alter_db(conn, cursor):
    cursor.execute("""CREATE TABLE taxonomy(
                      taxon_id INTEGER,
                      taxon TEXT,
                      family TEXT,
                      kingdom TEXT,
                      opentol_id TEXT)""")
    cursor.execute("""INSERT INTO taxonomy SELECT taxon.taxon_id, 
                      taxon.taxon, taxon.family, taxon.kingdom, opentol_temp.ID
                      FROM taxon LEFT JOIN opentol_temp ON
                      opentol_temp.inputName = taxon.taxon""")
    cursor.execute("""DROP TABLE opentol_temp""")
    cursor.execute("""DROP TABLE taxon""")
    cursor.execute("""ALTER TABLE taxonomy RENAME TO taxon""")
    conn.commit()


if __name__ == '__main__':
    # Connect to the database (creates a new file if it doesn't exist)
    # TODO replace db name with user argument
    conn = sqlite3.connect("BOLD_COI_barcodes.db")
    # Create a cursor
    cursor = conn.cursor()

    # Make request to checklistbank and map BOLD taxons to opentol id's
    map_checklistbank(conn, cursor)

    # Alter new table and remove old ones
    alter_db(conn, cursor)

    # Close the connection
    conn.close()
