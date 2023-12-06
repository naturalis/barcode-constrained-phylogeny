import subprocess

import ete3
from ete3 import Tree


def check_constraint_length(tree_file):
    """Check length of the constraint tree. If less than 3 leaves ignore the constraint tree.
    :return:
    """
    try:
        tree = Tree(tree_file, format=1)  # Format indicates newick structure in file
        leaves = list(tree.get_leaves())  # Get all leaves
    except ete3.parser.newick.NewickError:
        leaves = []
    return len(leaves)


def raxml(wildcard, alignment_file, constraint_tree):
    if check_constraint_length(constraint_tree) > 3:
        subprocess.run("mkdir -p ../results/raxml/{} | raxml-ng --msa {} --model GTR+G --tree-constraint {} --prefix ../results/raxml/{}/{}".format(wildcard, alignment_file, constraint_tree, wildcard, wildcard), shell=True)
    else:
        subprocess.run("mkdir -p ../results/raxml/{} | raxml-ng --msa {} --model GTR+G  --prefix ../results/raxml/{}/{}".format(wildcard, alignment_file, wildcard, wildcard), shell=True)


if __name__ == '__main__':
    wildcard = snakemake.wildcards # noqa: F821
    alignment_file = snakemake.input[0] # noqa: F821
    constraint_tree = snakemake.input[1] # noqa: F821
    raxml_output = snakemake.output[0] # noqa: F821

    raxml(wildcard, alignment_file, constraint_tree)