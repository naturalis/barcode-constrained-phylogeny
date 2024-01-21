import argparse
import dendropy
import logging
import os

from Bio import SeqIO


def read_tree(filename, rooting='default-rooted', schema='newick'):
    """
    Reads the provided tree file (in the format specified as `schema`) using Dendropy with the
    provided rooting policy.
    :param filename: file location
    :param rooting: rooting policy
    :param schema: file format, e.g. newick
    :return: a dendropy.Tree object
    """
    logger.info(f'Going to read {filename} as {rooting} {schema}')
    return dendropy.Tree.get(
        path=filename,
        schema=schema,
        rooting=rooting
    )


def preprocess_subtree(tree_file, fasta_file):
    """
    Preprocesses the input tree file by pruning all the tips that are not also in the input FASTA file
    :param tree_file: Newick tree file
    :param fasta_file: FASTA sequence file
    :return: pruned dendropy tree
    """
    logger.info(f'Going to preprocess tree {tree_file}')
    logger.info(f'Using alignment as ingroup {fasta_file}')

    # Read the tree
    tree = read_tree(tree_file)

    # Read the alignment
    ingroup = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        ingroup.add(record.id)

    # Iterate over the leaves
    outgroup = set()
    for leaf in tree.leaf_nodes():
        if leaf.taxon.label not in ingroup:
            logger.info(f'Going to prune outgroup sequence {leaf.taxon.label}')
            outgroup.add(leaf.taxon)

    # Prune the outgroup
    tree.prune_taxa(outgroup)

    return tree


if __name__ == '__main__':

    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-t', '--tree', required=True, help='Input backbone Newick tree')
    parser.add_argument('-f', '--folder', required=True, help='Location of folder with subtree folders')
    parser.add_argument('-o', '--out', required=True, help="Output grafted newick")
    parser.add_argument('-n', '--nfamilies', required=True, help='Number of families')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Configure logger
    logging.basicConfig()
    logger = logging.getLogger('graft_clades')
    logger.setLevel(args.verbosity)

    # Read the backbone
    backbone = read_tree(args.tree)

    # Iterate over folders
    base_folder = os.path.abspath(args.folder)
    for i in range(1, int(args.nfamilies) + 1):
        logger.info(f'Processing subtree {i}')

        # Peprocess the focal family tree
        subfolder = f'{i}-of-{args.nfamilies}'
        alignment_file = os.path.join(base_folder, subfolder, 'aligned.fa')
        subtree_file = os.path.join(base_folder, subfolder, 'raxml-ready.fa.raxml.bestTree.rooted')
        subtree = preprocess_subtree(subtree_file, alignment_file)

        # Intersect the subtree labels and the backbone set
        subtree_leaf_labels = set([leaf.taxon.label for leaf in subtree.leaf_nodes()])
        backbone_leaf_labels = set([leaf.taxon.label for leaf in backbone.leaf_nodes()])
        intersection = set()
        for label in subtree_leaf_labels:
            if label in backbone_leaf_labels:
                intersection.add(label)

        # Graft the subtree
        mrca = backbone.mrca(taxon_labels=intersection)
        mrca.clear_child_nodes()
        for child in subtree.seed_node.child_nodes():
            mrca.add_child(child)

    backbone.write(path=args.out, schema="newick")

