import numpy as np
from ete3 import Tree       # To use for distance matrix
import pandas as pd
import os


def loop_over_fam():
    for family in family_list:
        mode = get_mode_excelwriter(family, family_list)
        family = family.split(".")  # Expect a X.fasta.raxml.bestTree
        dist_df = create_matrix(family[0])
        write_matrix_to_file(dist_df, family[0], outputfile, mode)


def get_mode_excelwriter(family, flist):
    if family == flist[0]:
        mode = "w"
    else:
        mode = "a"
    return mode

def create_matrix(family):
    """Get newick from RAXML output.
    Get leave info from the newick.
    Create a matrix based on the amount of leaves
    Get all names from the leaves and append to list.
    Use the get distance from ete3, to get the distance between two nodes.
    Use that distance to add to the dataframe.
    :return: dataframe containing distances.
    """
    tree = Tree("test/matK/{}.fasta.raxml.bestTree".format(family),
                format=1)  # Format indicates newick structure in file
    leave = list(tree.get_leaves())  # Get all leaves
    names = []
    dmat = np.zeros((len(leave), len(leave)))
    for num in range(len(leave)):
        names.append(leave[num].name)
        for num2 in range(num, len(leave)):
            d = tree.get_distance(leave[num],
                                  leave[num2],
                                  topology_only=False)
            dmat[num, num2] = dmat[num2, num] = round(d, 6)  # dmat makes matrix
    dist_df = pd.DataFrame(data=dmat, index=names, columns=names)  # Make into dataframe
    return dist_df

def write_matrix_to_file(dist_df, family, outputfile, mode):
    """Append the matrix to file.
    Path contains the path to the excel file.
    The matrix will be appended to the path file.
    The matrix will on his own sheet with the corresponding family name.
    """
    path = r"test/matK/{}".format(outputfile)
    with pd.ExcelWriter(path, engine='openpyxl', mode=mode) as writer:
        dist_df.to_excel(writer, sheet_name=family)


if __name__ == "__main__":
    par_path = os.path.abspath(os.path.join(os.pardir))
    path = snakemake.input[0]
    outputfile = snakemake.output[0]
    family_list = os.listdir(par_path + path)
    loop_over_fam()

