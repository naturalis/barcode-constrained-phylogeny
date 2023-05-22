import numpy as np
from ete3 import Tree       # To use for distance matrix
import pandas as pd
import os

# Getting rid of test with excel and changing it back to csv file
# Loop can maybe be done in snakemake
# From data/fasta/family get  per family dir
#

def loop_over_fam(path):
    for family in path:
        if family.__contains__("bestTree"):
            print("fam",family)
            family = family.split(".")  # Expect a X.fasta.raxml.bestTree
            outputfile = "data/fasta/{}/dist_matrix_{}.txt".format(family[0],family[0])
            dist_df = create_matrix(family[0])
            write_matrix_to_file(dist_df, outputfile) # write matrix to csv


def create_matrix(family):
    """Get newick from RAXML output.
    Get leave info from the newick.
    Create a matrix based on the amount of leaves
    Get all names from the leaves and append to list.
    Use the get distance from ete3, to get the distance between two nodes.
    Use that distance to add to the dataframe.
    :return: dataframe containing distances.
    """
    print("data/fasta/{}/{}.fasta.raxml.bestTree".format(family,family))
    tree = Tree("data/fasta/{}/{}.fasta.raxml.bestTree".format(family,family),
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

def write_matrix_to_file(dist_df, outputfile):
    """Append the matrix to file.
    Path contains the path to the excel file.
    The matrix will be appended to the path file.
    The matrix will on his own sheet with the corresponding family name.
    """
    dist_df.to_csv(outputfile,sep='\t',header=True)


if __name__ == "__main__":
    path2 = os.getcwd()     # Get current working directory
    path = snakemake.input[0]   # noqa: F821
    total_path = os.path.join(path2, path)
    family_list = os.listdir(os.path.join(path2, path))
    print(family_list)
    loop_over_fam(family_list)

