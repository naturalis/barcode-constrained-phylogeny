import argparse
import csv
import os
import sqlite3
import pandas as pd
import numpy as np

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()

def select_marker():
    """Ask user to enter their desired marker.
    Based on their answer add argument (type marker and kingdom).
    :return: arguments and desired marker.
    """
    print("Select your desired marker. \n"
          "1) For COI-5P \n"
          "2) For matK and rbcL")
    marker = input("Selected number: ")

    # Use userinput to determine desired marker
    match int(marker):
        case 1:
            desired_marker = "COI-5P"

            # Add argument based on desired marker
            parser.add_argument('-marker', default=desired_marker,
                                help="Which barcode marker(s) to select: COI-5P, MatK,"
                                     " RbcL")
            kingdom = "Animalia"
            add_args(kingdom, desired_marker)
            args = parser.parse_args()
            return args, desired_marker
        case 2:
            kingdom = "Plantae"
            desired_marker = "matK_rbcL"

            # Add argument based on desired marker
            parser.add_argument('-marker', action='append', nargs=2, metavar=('matK', 'rbcL'),
                                help="Which barcode marker(s) to select: COI-5P, MatK,"
                                     " RbcL")
            add_args(kingdom, desired_marker)
            args = parser.parse_args('-marker matK rbcL'.split())
            return args, desired_marker
        case _:
            return "Something went wrong"


def add_args(kingdom, desired_marker):
    """Add argument kingdom based on user's input, select input directory and create database file.
    :param kingdom: Chosen kingdom based on the user's input.
    :param desired_marker: Chosen marker based on the user's input.
    """
    parser.add_argument('-kingdom', default=kingdom,
                        help="Which kingdom to filter barcodes from: Animalia"
                             " or Plantae")
    parser.add_argument('-indir', default=par_path + "/data/mnt/bold_public/"
                                                     "datapackages/recent-data/"
                                                     "BOLD_Public.30-Dec-2022.tsv",
                        help="BOLD tsv file location")

    # Add to file containing the marker
    parser.add_argument('-db', default=par_path + "/data/databases/BOLD_{}_barcodes.db".format(desired_marker),
                        help="Name of the the database file: {file_name}.db")


def extract_bold(conn, args, desired_marker):
    """
     It reads a TSV file with a snapshot of the BOLD database in chunks, selects rows that match the user's
    arguments for the desired marker and kingdom, and writes them to two tables in the
    SQLite database (taxon_temp and barcode_temp) with a subset of their columns.
    :param conn: Connection object to the database.
    :param args: Command line arguments
    :param desired_marker: String representing the desired marker
    """
    for chunk in pd.read_csv(args.indir, quoting=csv.QUOTE_NONE,
                             low_memory=False, sep="\t", chunksize=10000):
        # Keep rows that match user arguments
        if desired_marker == "COI-5P":
            df = chunk.loc[
                (chunk['marker_code'] == args.marker) & (chunk["kingdom"] == args.kingdom)]

        # args.marker[0][0] is matK and args.marker[0][1] is rbcL
        else:
            df = chunk.loc[
                ((chunk['marker_code'] == args.marker[0][0]) & (chunk["kingdom"] == args.kingdom)) |
                (chunk['marker_code'] == args.marker[0][1]) & (chunk["kingdom"] == args.kingdom)
                ]

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
    # Set args
    args, desired_marker = select_marker()

    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect(args.db)

    # Create a cursor
    cursor = conn.cursor()

    # Dump BOLD data into DB
    extract_bold(conn, args, desired_marker)

    # Close the connection
    conn.close()