# Directory
import argparse
import csv
import os
import sqlite3
import pandas as pd

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()
parser.add_argument('-marker', default="matK",
                    help="Which barcode marker(s) to select: COI-5P, MatK, R...")
parser.add_argument('-kingdom', default="Plantae",
                    help="Which kingdom to filter barcodes from: Animalia or Plantae")
parser.add_argument('-indir', default=par_path+"/data/mnt/bold_public/datapackages/recent-data/BOLD_Public.30-Dec-2022.tsv",
                    help="BOLD tsv file location")
args = parser.parse_args()

def make_db(conn, cursor):
    cursor.execute("""CREATE TABLE IF NOT EXISTS taxon (
        taxon TEXT,
        kingdom TEXT NOT NULL,
        family TEXT NOT NULL,
        opentol_id TEXT)
    """)
    # Create a table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS barcode (
        barcode_id INTEGER PRIMARY KEY,
        processid TEXT,
        marker_code TEXT,
        nucraw TEXT,
        country TEXT,
        taxon_id INTEGER,
        FOREIGN KEY (taxon_id) REFERENCES taxonomy (taxon_id)
    )""")
    # Commit the changes
    conn.commit()

def extract_bold(conn, cursor):
    # pd.set_option('display.max_columns', None)
    for chunk in pd.read_csv(args.indir,quoting=csv.QUOTE_NONE,
                             low_memory=False, sep="\t", chunksize=10000):
        df = chunk.loc[(chunk['marker_code'] == args.marker) &
                       (chunk["kingdom"] == args.kingdom)]
        df_temp = df[['taxon', 'kingdom', 'family']].dropna()

        df_temp.to_sql('taxon', conn, if_exists='append',
                            index=False)
        df_temp = df[['processid', 'marker_code', 'nucraw', 'country']]
        df_temp.to_sql('barcode', conn, if_exists='append', index=False)
        conn.commit()
    conn.close()


if __name__ == '__main__':
    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect("bold_test.db")

    # Create a cursor
    cursor = conn.cursor()

    # Make SQLite tables
    make_db(conn, cursor)

    # Dump BOLD data into DB
    extract_bold(conn, cursor)

    # Close the connection
    conn.close()