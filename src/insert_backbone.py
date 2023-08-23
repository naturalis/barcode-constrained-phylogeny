# Get backbone tree
# insert
from ete3 import Tree       # To use for distance matrix
import os
import re
import logging


#logging.basicConfig(snakemake.params.log_level) # noqa: F821
logger = logging.getLogger(__name__)



def read_backbone(backbone_file):
    """Read the backbone tree file and return it.
    :return: newick string"""
    with open(backbone_file, "r") as output:
        backbone = output.read()
    return backbone

def alter_backbone(backbone_file, subtree_file):
    """Get newick from RAXML output.
    Get leave info from the newick.
    Create a matrix based on the amount of leaves
    Get all names from the leaves and append to list.
    Use the get distance from ete3, to get the distance between two nodes.
    Use that distance to add to the dataframe.
    :return: dataframe containing distances.
    """

    tree = Tree(backbone_file,format=1)  # Format indicates newick structure in file
    leaves = list(tree.get_leaves())  # Get all leaves
    for i in range(len(leaves)):
        if (i % 2) == 0:
            leave_to_replace, tree2 = evaluate(leaves[i], subtree_file)
            print("replace", leave_to_replace)
            if leave_to_replace != "":
                leaves[i].add_child(tree2)
                leaves[i].delete()
        else:
            leaves[i].delete()
    return tree

def evaluate(leave,subtree):
    """Compare the newick subtree with the names from the leave.
    Open the subtree file.
    Get the newick from the subtree.
    If the newick is empty the tree is empty.
    If the newick is not empty perform operations.
    Transform the newick string to a Tree object. 
    Get the leave(from backbone), use [3:] so you get the name.
    Change the format of the leave.
    Iterate over the leaves and check:
        Check if leave from subtree is the same as leave from backbone.
    :return: {found} empty if no similarities are found, else return the leave"""
    
    found = ""
    with open(subtree, "r") as sub:
        subtree = sub.read()
    if subtree == "":
        tree2 = ""
    else:
        tree2 = Tree(subtree, format=1)
        leave2 = str(leave)[3:]
        leave2 = leave2.replace(" ", "")
        leave2 = leave2.replace("'", "")
        found = tree2.get_leaves_by_name(name=f"{leave2}")
    return found, tree2

def call_functions(subtree, backbone, altered_file):
    """Call the functions that need to be performed"""
    backbone = read_backbone(backbone)
    tree = alter_backbone(backbone, subtree)
    logger.info("Added subtree to backbone")
    logger.info("Writing tree to file")
    tree.write(format=1, outfile=altered_file)



if __name__ == "__main__":
    subtree = snakemake.input[0]    # noqa: F821
    backbone = snakemake.input[1]   # noqa: F821
    altered_backbone = snakemake.output[0]  # noqa: F821
    call_functions(subtree, backbone, altered_backbone)
