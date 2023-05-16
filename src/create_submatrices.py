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

par_path = os.path.abspath(os.path.join(os.pardir))
# User arguments
family_list = os.listdir(par_path + "/src/test/matK/")
file = "test/matK/matrix/dist_matK.xlsx"
par_path = os.path.abspath(os.path.join(os.pardir))


def create_fam_list():
    par_path = os.path.abspath(os.path.join(os.pardir))
    family_list = os.listdir(par_path + "/src/test/matK/trees/")
    flist = []
    for family in family_list:
        family = family.split(".")
        flist.append(family[0])
    return flist

def loop_over_df(flist, cursor, df):
    for family in flist:
        df1 = df.get(family)
        names = get_otts_in_df(df1)
        # Write barcodes to FASTA in family groups
        mode = get_mode_excelwriter(family, flist)
        blacklist = get_blacklist(names, cursor)
        alter_matrix(df1, blacklist, mode, outputfile, family)

def get_otts_in_df(df):
    names = []
    for c in df.columns:
        names.append(c)
    names.pop(0)
    return names

def get_mode_excelwriter(family, flist):
    if family == flist[0]:
        mode = "w"
    else:
        mode = "a"
    return mode
def read_xlsx(excel_file,family):
    df = pd.read_excel(excel_file, sheet_name=family, index_col=[0])
    return df

def get_blacklist(names, cursor):
    """Create matrix containing only the ott found in the database"""
    blacklist = []
    for name in names:
        ott = name.split("_")
        result = cursor.execute("SELECT COUNT() FROM node WHERE name = '{}';".format(ott[0]))
        results = result.fetchall()
        found = True if results[0][0] == 1 else False
        if not found:
            blacklist.append(name)
    return blacklist


def alter_matrix(dist_df, blacklist, mode, outputfile, family):
    #TODO returns warning but no error, seems to work
    df = dist_df[dist_df.columns.drop(blacklist)]   # Working
    df.drop(blacklist, inplace=True, axis=0)
    path = r"test/matK/{}".format(outputfile)
    with pd.ExcelWriter(path, engine='openpyxl', mode=mode) as writer:
        df.to_excel(writer, sheet_name=family)
    return df


if __name__ == "__main__":
    # TODO change with snakemake
    # TODO get rid of calling db like this
    db_file = snakemake.input[0]
    excel_file = "test/matK/matrix/dist_matK.xlsx"
    #outputfile = snakemake.output[0]
    #conn = sqlite3.connect(db_file)
    # Create a cursor
    #cursor = conn.cursor()
    flist = create_fam_list()
    df = read_xlsx(excel_file, flist)
    print(df.keys())
    #loop_over_df(flist, cursor, df)
    # Close the connection
    #conn.close()