import sqlite3
import argparse
import os

import numpy as np
import pandas as pd

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()
# Add to file containing the marker
parser.add_argument('-db', default="BOLD_COI_barcodes.db",
                    help="Name of the the database file: {file_name}.db")
args = parser.parse_args()


def divide_fastafiles(conn, cursor):
    os.makedirs('fasta/family', exist_ok=True)
    df = pd.read_sql_query("SELECT barcode.barcode_id, taxon.family, barcode.nucraw FROM barcode LEFT JOIN taxon ON barcode.taxon_id = taxon.taxon_id LIMIT 10000", conn)

    # Split dataframe into batches and put them in a list
    for family in set(df['family']):
        print('Making FASTA file for sequences from the family %s...' % family)
        # Grab records from specific family
        df_family = df[df['family'] == family]

        pd.options.mode.chained_assignment = None

        df_family['fasta'] = df_family['barcode_id'].apply(lambda x: '>' + str(x) + '\n')
        fasta_out = df_family['fasta'] + df_family["nucraw"]
        # Save BOLD sequences and their respective header in a fastafile
        np.savetxt("fasta/family/%s.fasta" % family, fasta_out.values,
                           fmt="%s")


if __name__ == '__main__':
    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect(args.db)

    # Create a cursor
    cursor = conn.cursor()
    divide_fastafiles(conn,cursor)

    # Close the connection
    conn.close()