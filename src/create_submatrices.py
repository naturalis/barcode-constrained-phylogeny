"""Create file containing all distances matrices.

This script uses a bestTree file from RAXML to create a distance matrix.

PREREQUISITES: The otol database has been downloaded.
Megatree tool is installed and used to create database from otol file.

The otol version used for this project was 13.4
Get otol database from https://files.opentreeoflife.org/synthesis/opentree13.4/opentree13.4tree.tgz
The database used is based on the /opentree13.4_tree/labelled_supertree/labelled_supertree_ottnames.tre

Install Megatree using the GITHUB doc.
https://github.com/rvosa/bio-phylo-forest-dbtree.git

creating the otol database with megatree loader:
megatree-loader -i labelled_supertree.tre -d {outputfile.db}
"""


import sqlite3
import os
import pandas as pd
import logging
#logging.basicConfig(snakemake.params.log_level) # noqa: F821
logger = logging.getLogger(__name__)


def loop_over_families(cursor, input_file, outputfile):
    #family = input_file.rstrip("fasta")
    #family = "fasta/alignment/Aethionema_hmmer.fasta"
    logger.info(f"Creating altered submatrix {input_file}")
    df = pd.read_csv(input_file,delimiter="\t")
    print(df)
    names = get_otts_in_df(df)
    print(names)
    # Write barcodes to FASTA in family groups
    blacklist = get_blacklist(names, cursor)
    print(blacklist)
    alter_matrix(df, blacklist, outputfile)

def get_otts_in_df(df):
    names = []
    for c in df.columns:
        names.append(c)
    names.pop(0)
    return names


def get_blacklist(names, cursor):
    """Create matrix containing only the ott found in the database"""
    blacklist = []
    for name in names:
        ott = name.split("|")
        result = cursor.execute("SELECT COUNT() FROM node WHERE name = '{}';".format(ott[1]))
        results = result.fetchall()
        print(results)
        found = True if results[0][0] == 1 else False
        if not found:
            blacklist.append(name)
    return blacklist


def alter_matrix(dist_df, blacklist, outputfile):
    df = dist_df[dist_df.columns.drop(blacklist)]   # Working
    #df.drop(blacklist, inplace=True, axis=0) # -> Why not working, it worked last week :(
    df.to_csv(outputfile, sep='\t', header=True)
    logger.info(f"Wrote distance matrix to file {outputfile}")
    return df


if __name__ == "__main__":
    #path2 = os.getcwd()  # Get current working directory
    #path = snakemake.input[0]  # noqa: F821
    #total_path = os.path.join(path2, path)
    #family_list = os.listdir(os.path.join(path2, path))
    #print(family_list)
    #db_file = snakemake.input[1]    # noqa: F821
   # db_file="../data/databases/outfile.db"
    input_file = snakemake.input[0] # noqa: F821
    db_file = snakemake.input[1] # noqa: F821
    output_file = snakemake.output[0]   # noqa: F821
    conn = sqlite3.connect(db_file)
    # Create a cursor
    cursor = conn.cursor()
    loop_over_families(cursor, input_file, output_file)
    # Close the connection
    conn.close()
