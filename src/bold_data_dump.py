# Directory
import argparse
import csv
import os
import sqlite3
import pandas as pd

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()
# FIXME How to do this so that user can give multiple markers?
parser.add_argument('-marker', default="COI-5P",
                    help="Which barcode marker(s) to select: COI-5P, MatK,"
                         " RbcL")
parser.add_argument('-kingdom', default="Animalia",
                    help="Which kingdom to filter barcodes from: Animalia"
                         " or Plantae")
parser.add_argument('-indir', default=par_path+"/data/mnt/bold_public/"
                                               "datapackages/recent-data/"
                                               "BOLD_Public.30-Dec-2022.tsv",
                    help="BOLD tsv file location")
parser.add_argument('-db', default="BOLD_COI_barcodes.db",
                    help="Name of the the database file: {file_name}.db")
args = parser.parse_args()


def extract_bold(conn):
    for chunk in pd.read_csv(args.indir,quoting=csv.QUOTE_NONE,
                             low_memory=False, sep="\t", chunksize=10000):
        #  FIXME How to do this so that user can give multiple markers?
        # Keep rows that match userarguments
        df = chunk.loc[
            (chunk['marker_code'] == args.marker) &
            (chunk["kingdom"] == args.kingdom)]

        # Keep stated columns, do not keep rows where NAs are present
        df_temp = df[['taxon', 'kingdom', 'family']].dropna()

        # Add rows to SQLite table (makes table if not exitst yet)
        df_temp.to_sql('taxon_temp', conn, if_exists='append',
                            index=False)

        # Keep stated columns
        df_temp = df[['processid', 'marker_code', 'nucraw', 'country', 'taxon']]

        # Add rows to SQLite table (makes table if not exitst yet)
        df_temp.to_sql('barcode_temp', conn, if_exists='append', index=False)

        conn.commit()
    conn.close()


if __name__ == '__main__':
    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect(args.db)

    # Create a cursor
    cursor = conn.cursor()

    # Dump BOLD data into DB
    extract_bold(conn)

    # Close the connection
    conn.close()