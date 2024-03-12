import argparse
import dendropy
import os
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
    logger.info(f'Going to read {filename} as {rooting} {schema}')
    return dendropy.Tree.get(
        path=filename,
        schema=schema,
        rooting=rooting
    )


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
    logger = util.get_formatted_logger('graft_clades', args.verbosity)

    # Read the backbone tree as a dendropy object, calculate distances to root, and get its leaves
    backbone = read_tree(args.tree)
    backbone.calc_node_root_distances()
    backbone_leaf_labels = set([leaf.taxon.label for leaf in backbone.leaf_nodes()])

    # Iterate over folders
    base_folder = os.path.abspath(args.folder)
    for i in range(1, int(args.nfamilies) + 1):
        logger.info(f'Processing subtree {i}')

        # Peprocess the focal family tree
        subfolder = f'{i}-of-{args.nfamilies}'
        subtree_file = os.path.join(base_folder, subfolder, 'aligned.fa.raxml.bestTree.rooted')
        subtree = read_tree(subtree_file)
        subtree.calc_node_root_distances()

        # Intersect the subtree labels and the backbone set if subtree >= 3 tips
        subtree_leaf_labels = set([leaf.taxon.label for leaf in subtree.leaf_nodes()])
        if len(subtree_leaf_labels) >= 3:
            intersection = set()
            for label in subtree_leaf_labels:
                if label in backbone_leaf_labels:
                    intersection.add(label)

            # Find the mrca, calculate distance between exemplars
            logger.info(f'Intersection {intersection}')
            mrca = backbone.mrca(taxon_labels=intersection)
            bbdist = 0
            stdist = 0
            for label in intersection:
                bbleaf = backbone.find_node_with_taxon_label(label)
                bbdist += (bbleaf.root_distance - mrca.root_distance)
                stleaf = subtree.find_node_with_taxon_label(label)
                stdist += stleaf.root_distance

            # Adjust subtree
            rescale = bbdist / stdist
            logger.info(f'Rescaling by factor {rescale}')
            for node in subtree.preorder_node_iter():
                node.edge_length = node.edge_length * rescale

            # Graft subtree
            mrca.clear_child_nodes()
            for child in subtree.seed_node.child_nodes():
                mrca.add_child(child)

    backbone.write(path=args.out, schema="newick")
