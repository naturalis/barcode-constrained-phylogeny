# Get backbone tree
# insert
from ete3 import Tree       # To use for distance matrix
import os
import re
import logging


#logging.basicConfig(snakemake.params.log_level) # noqa: F821
logger = logging.getLogger(__name__)



def read_backbone(backbone_file):
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
            if leave_to_replace != "":
                leaves[i].add_child(tree2)
                leaves[i].delete()
        else:
            leaves[i].delete()
    return tree

def evaluate(leave,subtree):
    found = ""
    with open(subtree, "r") as sub:
        subtree = sub.read()
    if subtree == "":
        tree2 = ""
    else:
        tree2 = Tree(subtree, format=1)
        leave2 = str(leave)[3:]
        found = tree2.search_nodes(name=leave2)
    return found, tree2

def loop_over_fam(path):
    tree = ""
    #for family in path:
    #    logger.info("Looping over family folders")
    #if family != "backbone":
    #subtree_file = "../data/fasta/{}/{}.fasta.raxml.bestTree".format(family,family)
    subtree_file = "fasta/filles/{}/{}1.bestTree".format(family, family)
    #subtree_file=family
    print(subtree_file)
    #backbone = read_backbone(backbone_file="../data/fasta/backbone/rep.bestTree")
    backbone = read_backbone("rep.bestTree")
    tree = alter_backbone(backbone, subtree_file)
    logger.info("Added subtree to backbone")
    logger.info("Writing tree to file")
    #tree.write(format=1, outfile="../data/fasta/backbone/filled_backbone.tree")
    tree.write(format=1, outfile="filled.tree")



if __name__ == "__main__":
    path2 = os.getcwd()  # Get current working directory
    path = "fasta/filles/"  # noqa: F821
    #outputfile = snakemake.output[0]    # noqa: F821
    total_path = os.path.join(path2, path)
    family_list = os.listdir(os.path.join(path2, path))
    loop_over_fam(family_list)