import numpy as np
from ete3 import Tree       # To use for distance matrix
import pandas as pd
from Bio import SeqIO
import sqlite3
import argparse
import os
import xlsxwriter

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()


parser.add_argument('-db', default="data/databases/outfile.db",
                    help="Name of the the database file: {file_name}.db")

args = parser.parse_args()


def create_matrix():
    #TODO set tree file as input variable
    """Keep this just writing to file
    :return:
    """
    tree = Tree("test/Abyl/Abyl.bestTree",
                format=1)  # Format indicates newick structure in file
    nodes = list(tree.get_leaves())  # Get all leaves
    names = []
    dmat = np.zeros((len(nodes), len(nodes)))
    for num in range(len(nodes)):
        names.append(nodes[num].name)
        for num2 in range(num, len(nodes)):
            d = tree.get_distance(nodes[num],
                                  nodes[num2],
                                  topology_only=False)
            dmat[num, num2] = dmat[num2, num] = round(d, 6)  # dmat makes matrix
    dist_df = pd.DataFrame(data=dmat, index=names, columns=names)  # Make into dataframe
    return dist_df, names




def get_blacklist(names):
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


def alter_matrix(dist_df, blacklist):
    #TODO returns warning but no error, seems to work
    df = dist_df[dist_df.columns.drop(blacklist)]
    df.drop(blacklist, inplace=True)
    df.to_csv("data/control_abyl.txt",
                   sep='\t',
                   header=True)
    return df


def get_highest(df):
    highest = df.max()  # <- cant use
    print(max(highest))  # Get highest values
    result_list = search_pos(df, {max(highest)})
    print(result_list)
    # Get row and column name
    representatives = get_corresponding_ott(df, result_list[0][0], result_list[0][1])
    return representatives



def search_pos(df_data: pd.DataFrame, search: set) -> list:
    nda_values = df_data.values
    tuple_index = np.where(np.isin(nda_values, [e for e in search]))
    return [(row, col, nda_values[row][col]) for row, col in zip(tuple_index[0], tuple_index[1])]


def get_corresponding_ott(df, col_pos, row_pos):
    representatives = []
    print("pos",col_pos)
    colname = df.columns[col_pos]
    print(colname)
    print("pos",row_pos)
    rowname = df.index[row_pos]
    print("index",df.index[row_pos])
    representatives.append(colname)
    representatives.append(rowname)
    return representatives


def representatives_to_file(representatives):
    """Important to get the normal sequences in this file"""

    f = "test/Abyl/Abylidae_with_ott.fasta"
    #with open("{}/{}".format(dir, filename), "r") as input:
    with open(f, "r") as input:
        with open("representatives.fasta", "a") as outputs:
            for record in SeqIO.parse(input, "fasta"):
                for seq in representatives:
                    if record.id == seq:
                        SeqIO.write(record, outputs, 'fasta')


def get_species(representatives):
    with open("representatives.txt", "a") as output:
        for header in representatives:
            header = header.split("_")
            output.write(header[0]+"\n")


if __name__ == "__main__":

    #TODO change with snakemake
    #conn = sqlite3.connect(args.db)
    # Create a cursor
    #cursor = conn.cursor()
    # Write barcodes to FASTA in family groups
    dist_df, names = create_matrix()
    # Check to add to sheet instead of new file
    #dist_df.to_csv("data/tree_distancematrix2.txt",
    #               sep='\t',
    #               header=True)  # Write distance matrix to txt file

    path = r"data/tree_distance_matrices.xlsx"
    writer = pd.ExcelWriter(path, engine='xlsxwriter')
    dist_df.to_excel(writer, sheet_name='Abyl')
    writer.close()
    #blacklist = get_blacklist(names)
    #altered_matrix = alter_matrix(dist_df, blacklist)
    #representatives = get_highest(altered_matrix)
    #representatives_to_file(representatives)
    #get_species(representatives)
    # Close the connection
    #conn.close()







