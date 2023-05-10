"""Create file containing all distances matrices.

This script uses a bestTree file from RAXML to create a distance matrix.

PREREQUISITES:

"""
import numpy as np
from ete3 import Tree       # To use for distance matrix
import pandas as pd
import argparse
import os

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()


parser.add_argument('-db', default="data/databases/outfile.db",
                    help="Name of the the database file: {file_name}.db")

args = parser.parse_args()


def create_matrix(family):
    #TODO set tree file as input variable
    """Keep this just writing to file
    :return:
    """
    tree = Tree("test/Abyl/{}.bestTree".format(family),
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

def write_matrix_to_file(dist_df, family):
    """Append the matrix to file.
    Path contains the path to the excel file.
    The matrix will be appended to the path file.
    The matrix will on his own sheet with the corresponding family name.
    """
    path = r"data/tree_distance_matrices.xlsx"
    with pd.ExcelWriter(path, engine='openpyxl', mode='a') as writer:
        dist_df.to_excel(writer, sheet_name=family)


if __name__ == "__main__":
    #TODO change with snakemake
    #TODO get rid of calling db like this
    family="Abronicidae"
    dist_df, names = create_matrix(family)
    write_matrix_to_file(dist_df, family)







