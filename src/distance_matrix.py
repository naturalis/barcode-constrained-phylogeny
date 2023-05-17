import numpy as np
from ete3 import Tree       # To use for distance matrix
import pandas as pd
import os

# Getting rid of test with excel and changing it back to csv file
# Loop can maybe be done in snakemake
# From data/fasta/family get  per family dir
#

def loop_over_fam(path):
    for family in family_list:
        outputfile = ""
        #mode = get_mode_excelwriter(family, family_list)
        family = family.split(".")  # Expect a X.fasta.raxml.bestTree
        dist_df = create_matrix(family[0])

        #write_matrix_to_file(dist_df, family[0], outputfile, mode) -> write_matrix_to_csv


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
    dist_df.to_csv("data/control_abyl.txt",sep='\t',header=True)


if __name__ == "__main__":
    path2 = os.getcwd()     # Get current working directory
    path = snakemake.input[0]   # noqa: F821
    path = os.path.normpath(path)
    total_path = os.path.join(path2, path)
    print(total_path)
    family_list = os.listdir(os.path.join(path2, path))
    print(family_list)
    #loop_over_fam(path)

