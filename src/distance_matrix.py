import numpy as np
from ete3 import Tree       # To use for distance matrix
import pandas as pd


def create_matrix(family):
    """Get newick from RAXML output.
    Get leave info from the newick.
    Create a matrix based on the amount of leaves
    Get all names from the leaves and append to list.
    Use the get distance from ete3, to get the distance between two nodes.
    Use that distance to add to the dataframe.
    :return: dataframe containing distances.
    """
    tree = Tree("fasta/{}/{}.bestTree".format(family, family),
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

def write_matrix_to_file(dist_df, family, outputfile):
    """Append the matrix to file.
    Path contains the path to the excel file.
    The matrix will be appended to the path file.
    The matrix will on his own sheet with the corresponding family name.
    """
    path = r"data/{}".format(outputfile)
    with pd.ExcelWriter(path, engine='openpyxl', mode='a') as writer:
        dist_df.to_excel(writer, sheet_name=family)


if __name__ == "__main__":
    family = snakemake.input[0]
    outputfile = snakemake.output[1]
    dist_df = create_matrix(family)
    write_matrix_to_file(dist_df, family, outputfile)







