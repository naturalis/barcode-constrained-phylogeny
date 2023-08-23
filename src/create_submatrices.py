"""Create file containing all distances matrices.

This script uses a bestTree file from RAXML to create a distance matrix.
"""

import sqlite3
import os
import pandas as pd
import logging
logger = logging.getLogger(__name__)


def loop_over_families(cursor, input_file, outputfile):
    """Call all the functions that need to be done"""
    logger.info(f"Creating altered submatrix {input_file}")
    df = pd.read_csv(input_file,delimiter="\t")
    names = get_otts_in_df(df)
    blacklist = get_blacklist(names, cursor)
    alter_matrix(df, blacklist, outputfile)

def get_otts_in_df(df):
    """Get all the ott_names in the matrix.
    :return: list of names in the matrix."""
    names = []
    for c in df.columns:
        names.append(c)
    names.pop(0)
    return names


def get_blacklist(names, cursor):
    """Create matrix containing only the ott found in the database.
    Iterate over the names list. 
    Check if the name in names list can be found in the open tree of life database.
    :return: list of names that can not  be  found in the opent tree of life database."""
    blacklist = []
    for name in names:
        ott = name.split("|")
        result = cursor.execute("SELECT COUNT() FROM node WHERE name = '{}';".format(ott[1]))
        results = result.fetchall()
        found = True if results[0][0] == 1 else False
        if not found:
            blacklist.append(name)
    return blacklist


def alter_matrix(dist_df, blacklist, outputfile):
    """Alter the matrix and write to file.
    Drop the names that can not be found in the open tree of life database. 
    Write the altered dataframe as csv file to {outputfile}.
    :return: dataframe containing entries that can be found in the open tree of life."""
    df = dist_df[dist_df.columns.drop(blacklist)]   
    df.to_csv(outputfile, sep='\t', header=True)
    logger.info(f"Wrote distance matrix to file {outputfile}")
    return df


if __name__ == "__main__":
    input_file = snakemake.input[0] # noqa: F821
    db_file = snakemake.input[1] # noqa: F821
    output_file = snakemake.output[0]   # noqa: F821
    conn = sqlite3.connect(db_file)
    # Create a cursor
    cursor = conn.cursor()
    loop_over_families(cursor, input_file, output_file)
    # Close the connection
    conn.close()
