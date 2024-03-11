import argparse
import dendropy
import util


def read_tree(filename, rooting='default-rooted', schema='newick'):
    """
    Reads the provided tree file (in the format specified as `schema`) using Dendropy with the
    provided rooting policy.
    :param filename: file location
    :param rooting: rooting policy
    :param schema: file format, e.g. newick
    :return: a dendropy.Tree object
    """
    logger.info(f'Going to read {filename} as {schema} with rooting {rooting}')
    return dendropy.Tree.get(
        path=filename,
        schema=schema,
        rooting=rooting
    )


def get_pseudo_outgroup(rooted_tree):
    """
    Given the provided rooted input tree, finds the smallest clade subtended by the root and returns for it
    the set of leaf labels of its descendants.
    :param rooted_tree: a rooted dendropy.Tree object
    :return: a set of leaf labels
    """
    logger.info(f"Going to find pseudo outgroup in input tree")
    root = rooted_tree.seed_node

    # Check if tree is bifurcating at the root
    if len(root.child_nodes()) != 2:
        logger.error("The tree is not bifurcating at the root.")

    # Find the smallest clade immediately subtended by the root
    smallest_clade = min(root.child_nodes(), key=lambda x: len(x.leaf_nodes()))

    # Return the leaves subtended by this clade
    return set([leaf.taxon.label for leaf in smallest_clade.leaf_nodes()])


def find_set_bipartition(tree, query_set):
    """
    Given an input set of leaf labels, finds the smallest edge bipartition of which the input
    is a subset.
    :param tree: an unrooted dendropy.Tree object
    :param query_set: a set of leaf labels
    :return:
    """
    logger.info(f"Going to find bipartition for provided input set")
    logger.debug(query_set)
    tree.update_bipartitions()

    # Identify the smallest split that monophylizes tip_set
    smallest_split = None
    for edge in tree.postorder_edge_iter():
        if edge.bipartition.leafset_bitmask is not None:
            leaves = set([taxon.label for taxon in edge.bipartition.leafset_taxa(tree.taxon_namespace)])
            logger.debug(f'Processing bipartition {leaves}')
            if query_set.issubset(leaves):
                logger.debug(f'Found smallest bipartition that has monophyletic {query_set}')
                smallest_split = edge
                break
    return smallest_split


if __name__ == '__main__':
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-i', '--intree', required=True, help='Input unrooted backbone tree file')
    parser.add_argument('-c', '--constraint', required=True, help="Input rooted backbone constraint")
    parser.add_argument('-o', '--outtree', required=True, help="Output rooted tree")
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Configure logging
    logger = util.get_formatted_logger('backbone_constraint', args.verbosity)

    # Read trees
    backbone = read_tree(args.intree)
    try:
        constraint = read_tree(args.constraint)

        # Get the pseudo-outgroup, i.e. the tips of the smallest clade subtended by the root in the constraint tree
        outgroup = get_pseudo_outgroup(constraint)

        # Find the bipartition that most closely resembles the pseudo-outgroup in the unrooted backbone
        bipartition = find_set_bipartition(backbone, outgroup)

        # Midpoint root on the bipartition
        new_length = bipartition.length / 2
        backbone.reroot_at_edge(bipartition, length1=new_length, length2=new_length, update_bipartitions=True)
    except Exception as e:
        logger.warning(f'Reading the constraint tree failed: {e}')

    # Write to file
    logger.info(f'Going to write rerooted tree to {args.outtree}')
    newick_string = backbone.as_string(schema="newick")
    with open(args.outtree, 'w') as file:
        file.write(newick_string)



