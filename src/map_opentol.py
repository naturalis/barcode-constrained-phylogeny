import json
from io import StringIO

import numpy as np
# import requests
import pandas as pd
import sqlite3

import requests


if __name__ == '__main__':
    # Connect to the database (creates a new file if it doesn't exist)
    # TODO replace db name with user argument
    conn = sqlite3.connect("BOLD_COI_barcodes.db")
    # Create a cursor
    cursor = conn.cursor()
    count = cursor.execute("""SELECT COUNT(taxon) FROM taxon""").fetchall()[0]
    results = cursor.execute("""SELECT taxon FROM taxon""").fetchall()
    df = pd.DataFrame(results, columns=['scientificName'])
    list_df = np.array_split(df, 50)
    url = "https://api.checklistbank.org/dataset/201891/nameusage/match"
    headers = {
        'Accept': 'text,csv',
        'Content-Type': 'text/tsv'
    }
    # cursor.execute("""DROP TABLE taxonomy""")
    # cursor.execute("""DROP TABLE opentol_temp""")
    for i in list_df:
        data = i.to_csv(index=False)
        response = requests.post(url, headers=headers,
                             data=data)
        content = response.content.decode('utf-8')
        f = pd.read_csv(StringIO(content), sep=',', usecols=['ID', 'inputName'])
        f.rename(columns={'ID':'opentol_id'})
        f.to_sql('opentol_temp', conn, if_exists='append')

    query_insert = """CREATE TABLE taxonomy AS
     SELECT taxon.taxon_id, taxon.taxon, taxon.family, taxon.kingdom, opentol_temp.id
     FROM taxon INNER JOIN opentol_temp ON opentol_temp.inputName = taxon.taxon"""
    result = cursor.execute(query_insert)
    conn.commit()
    # Close the connection
    conn.close()