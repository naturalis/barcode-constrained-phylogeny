import numpy as np
from ete3 import Tree       # To use for distance matrix
import pandas as pd
import os
import logging


#logging.basicConfig(snakemake.params.log_level) # noqa: F821
logger = logging.getLogger(__name__)

def loop_over_fam(family, outfile):
    # Does not work anymore with the current directory structure
    #for family in path:
        #if family != "backbone":
    logger.info(f"Creating distance matrix for {family}")
    #outputfile = "distance_matrix.txt"
    dist_df = create_matrix(family)
    write_matrix_to_file(dist_df, outfile) # write matrix to csv
    logger.info(f"Distance matrix for {family} is written to outputfile")


def create_matrix(family):
    """Get newick from RAXML output.
    Get leave info from the newick.
    Create a matrix based on the amount of leaves
    Get all names from the leaves and append to list.
    Use the get distance from ete3, to get the distance between two nodes.
    Use that distance to add to the dataframe.
    :return: dataframe containing distances.
    """
    #print("fasta/alignment/{}/{}.fasta.raxml.bestTree".format(family,family))
    # Add a check for empty file
    #try:
    tree = Tree(family, format=1)  # Format indicates newick structure in file
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
    #except:
        #logger.warning("Something went wrong")
    #dist_df = pd.DataFrame()    # Returns an empty dataframe
    return dist_df

def write_matrix_to_file(dist_df, outputfile):
    """Append the matrix to file.
    Path contains the path to the excel file.
    The matrix will be appended to the path file.
    The matrix will on his own sheet with the corresponding family name.
    """
    dist_df.to_csv(outputfile,sep='\t',header=True)


if __name__ == "__main__":
    input_file = snakemake.input[0]   # noqa: F821
    output_file = snakemake.output[0]    # noqa: F821
    loop_over_fam(input_file, output_file)

