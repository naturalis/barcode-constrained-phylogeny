import logging
import argparse
import sqlite3

from Bio import SeqIO
from Bio.Phylo import read as read_newick

logging.basicConfig()
logger = logging.getLogger('choose_exemplars')


def get_ingroup_labels(input_file):
    """
    Reads input FASTA, returns list of IDs (first words in defline)
    :param input_file: FASTA file
    :return:
    """
    labels = []
    for record in SeqIO.parse(input_file, "fasta"):
        labels.append(record.id)
    return labels


def pick_shallowest_tips(tree_file, ingroup):
    """
    Picks two tips nearest to the ingroup root on either side
    :param tree_file: Newick tree file
    :param ingroup: List of labels to consider
    :return:
    """
    logger.info(f'Going to pick exemplars from {tree_file}')
    tree = read_newick(tree_file, 'newick')
    leaves = [tree.find_any(name=label) for label in ingroup]
    mrca = tree.common_ancestor(leaves)
    children = mrca.clades
    representatives = []
    if len(children) == 2:
        for child in children:
            tips = child.get_terminals()
            dists = []
            for tip in tips:

                # Skip if focal tip is aberrant outgroup
                if tip.id not in ingroup:
                    logger.warning(f'Have aberrant outgroup {tip.id} in ingroup of tree {tree_file}')
                    continue

                # Calculate distance to mrca
                dist = tree.distance(tip, mrca)
                dists.append({'dist': dist, 'id': tip.id})

            # Get shallowest tip
            dists = sorted(dists, key=lambda x: x['dist'])
            representatives.append(dists[0]['id'])
    else:
        logger.error(f'Ingroup root in {tree_file} is not bifurcating')
    return representatives


def write_sequences(inaln, outaln, subset):
    """
    Appends the specified subset of sequences from the input alignment to the output
    :param inaln:
    :param outaln:
    :param subset:
    :return:
    """
    logger.info(f'Going to write {len(subset)} sequences to {outaln}')
    with open(outaln, 'a') as outfh:
        for record in SeqIO.parse(inaln, "fasta"):
            if record.id in subset:
                outfh.write(f'>{record.id}\n')
                outfh.write(f'{record.seq}\n')


if __name__ == "__main__":

    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-d', '--database', required=True, help='SQLite database file')
    parser.add_argument('-t', '--tree', required=True, help='Input Newick tree')
    parser.add_argument('-i', '--inaln', required=True, help='Input aligned FASTA file')
    parser.add_argument('-o', '--outaln', required=True, help="Output FASTA alignment")
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Configure logger
    logger.setLevel(args.verbosity)

    # Connect to the database (creates a new file if it doesn't exist)
    logger.info(f"Going to connect to database {args.db}")
    connection = sqlite3.connect(args.db)

    # Read FASTA, get list of ingroup tips
    seq_labels = get_ingroup_labels(args.inaln)

    # Maybe just include the whole file
    if len(seq_labels) < 3:
        logger.info(f'Infile {args.inaln} has fewer than 3 sequences, will include all')
        write_sequences(args.inaln, args.outaln, seq_labels)
    else:

        # ... or just the exemplars
        exemplars = pick_shallowest_tips(args.tree, seq_labels)
        write_sequences(args.inaln, args.outaln, exemplars)

