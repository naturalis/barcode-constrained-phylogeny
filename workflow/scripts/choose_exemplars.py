import util
import argparse
import dendropy
import statistics

from Bio import SeqIO


"""
This script, `choose_exemplars.py`, is responsible for selecting exemplar tips from a phylogenetic tree based on a 
specified strategy and writing the corresponding sequences from an input alignment to an output file.

The script performs the following steps:
1. Reads the input Newick tree.
2. Picks two exemplar tips from the tree based on the specified strategy. The strategies can be the two tallest, the 
   two shallowest, or the two closest to the median.
3. Reads the input FASTA alignment file.
4. Writes the sequences corresponding to the selected exemplar tips to an output file.

The script uses command line arguments for the input Newick tree file, input FASTA alignment file, output FASTA 
alignment file, tip picking strategy, and log level. The script is invoked by the Snakefile as a shell command with the
required arguments in the rule `choose_exemplars`.
"""


def pick_tips(tree, strategy):
    """
    Picks two exemplar tips, either the two tallest, the two shallowest, or the two closest to the median
    :param tree: Dendropy tree
    :param strategy: How to pick exemplars
    :return:
    """
    logger.info(f'Going to pick exemplars from {tree} following {strategy} strategy')

    # List of exemplars to return
    representatives = []

    # Get root and its immediate children
    root = tree.seed_node
    children = root.child_nodes()
    if len(children) == 2:
        for child in children:
            dists = []

            # Calculate distances to root from all leaves
            for tip in child.leaf_nodes():
                dist = tip.distance_from_root()
                dists.append({'dist': dist, 'id': tip.taxon.label})

            # Get shallowest tip
            dists = sorted(dists, key=lambda x: x['dist'])
            if str(strategy).capitalize().startswith('S'):
                representatives.append(dists[0]['id'])

            # Get tallest tip
            elif str(strategy).capitalize().startswith('T'):
                representatives.append(dists[-1]['id'])

            # Get median tip
            else:
                median_dist = statistics.median([d['dist'] for d in dists])
                closest = min(dists, key=lambda x: abs(x['dist'] - median_dist))
                representatives.append(closest['id'])

    else:
        logger.warning(f'Ingroup root in {tree} is not bifurcating, will approximate rooting')
        return None
    logger.debug(representatives)
    return representatives


def write_sequences(inaln, outaln, subset):
    """
    Appends the specified subset of sequences from the input alignment to the output
    :param inaln:
    :param outaln:
    :param subset:
    :return:
    """
    if subset == None:
        subset = []
    logger.info(f'Going to write {len(subset)} sequences to {outaln}')
    with open(outaln, 'a') as outfh:
        for record in SeqIO.parse(inaln, "fasta"):
            if record.id in subset:
                outfh.write(f'>{record.id}\n')
                outfh.write(f'{record.seq}\n')


if __name__ == "__main__":

    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-t', '--tree', required=True, help='Input Newick tree')
    parser.add_argument('-i', '--inaln', required=True, help='Input aligned FASTA file')
    parser.add_argument('-o', '--outaln', required=True, help="Output FASTA alignment")
    parser.add_argument('-s', '--strategy', required=True, help='Tip picking strategy: [t]all, [s]hallow, [m]edian')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Configure logger
    logger = util.get_formatted_logger('choose_exemplars', args.verbosity)

    # Read newick, get list of tips
    tree = dendropy.Tree.get(
        path=args.tree,
        schema="newick",
        rooting='default-rooted'
    )
    tips = []
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            tips.append(node.taxon.label)

    # Include the whole file if fewer than 3 ingroup tips
    if len(tips) < 3:
        logger.info(f'Infile {args.tree} has fewer than 3 tips, will include all')
        write_sequences(args.inaln, args.outaln, tips)

    # Pick exemplars by strategy
    else:
        exemplars = pick_tips(tree, args.strategy)
        write_sequences(args.inaln, args.outaln, exemplars)
