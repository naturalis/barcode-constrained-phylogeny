import sqlite3
import argparse
import os
import numpy as np
import pandas as pd

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()

# TODO NOAH: change default to how you named your DB or call with commandline
parser.add_argument('-db', default="BOLD_COI-5P_barcodes.db",
                    help="Name of the the database file: {file_name}.db")

args = parser.parse_args()


def divide_fastafiles(conn):
    """
    Takes the barcodes from the SQLite database and divides them into their
    taxonomic family names. For every family a fasta file is made named
    'fasta/family/{family name}.fasta'. All barcodes in that family are
    put in the fasta file with the barcode id from the SQLite db and their
    nucleotide sequence in FASTA format.
    :param conn: Connection to SQLite database
    """
    # Make directory to put FASTA files in
    os.makedirs('fasta/family', exist_ok=True)

    # Put needed data from db into a dataframe
    df = pd.read_sql_query("SELECT barcode.barcode_id, taxon.family, "
                           "barcode.nucraw, taxon.opentol_id FROM barcode LEFT JOIN taxon ON "
                           "barcode.taxon_id = taxon.taxon_id", conn)

    # Dropping rows where there is not an opentol_id
    df = df.dropna(subset=['opentol_id'])
    # Loop through unique family names
    for family in set(df['family']):
        # Takes about 25 min for COI (7,5 million barcodes)
        print('Making FASTA file for sequences from the family %s...' % family)

        # Grab records from specific family and put them in temp dataframe
        df_family = df[df['family'] == family]
        # Continue if there are more than one barcode from family
        if len(df_family) > 1:
            # Ignore warning
            pd.options.mode.chained_assignment = None

            # Make column with the FASTA header for every barcode
            df_family['fasta'] = df_family['barcode_id'].apply(
                lambda barcode_id: '>' + str(barcode_id) + '\n')

            # Add nucraw sequence with the header
            fasta_out = df_family['fasta'] + df_family["nucraw"]

            # Write to FASTA file and name it as their respective family name
            np.savetxt("fasta/family/%s.fasta" % family, fasta_out.values,
                               fmt="%s")


if __name__ == '__main__':
    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect(args.db)

    # Create a cursor
    cursor = conn.cursor()

    # Write barcodes to FASTA in family groups
    divide_fastafiles(conn)

    # Close the connection
    conn.close()
