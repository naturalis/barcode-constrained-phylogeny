import logging
import argparse
import sqlite3
import dendropy
import io

from Bio import SeqIO
from Bio.Phylo import read as read_newick

logging.basicConfig()
logger = logging.getLogger('choose_exemplars')

"""
In this script, we attempt to select the exemplar(s) that represent the focal family in the construction of the backbone
topology. Barring unforeseen errors, there are three possibilities:

1. The focal family is monotypic. Example: the Aye aye, which is the sole member of the family Daubentoniidae. In this
   case, there is nothing to choose, and the one sequence for the Aye Aye is used in the backbone.
2. The focal family has two species. Here, too, we are forced to use everything we got, i.e. the two sequences for the
   two species in the family. However, there will be some repercussions downstream, because we have to have these two
   sequences as a monophyletic group in the constraint tree.
3. There are three or more species/sequences. We choose as exemplars two species, which must straddle the root. This so
   that we can graft the subtree we have for this family in lieu of that split and the depths will be somewhat ok-ish.
   As exemplars we choose the shallowest tips, so that as little saturation as possible may have accumulated in the 
   sequences.
"""

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


def pick_shallowest_tips(tree_file, tree, ingroup):
    """
    Picks two tips nearest to the ingroup root on either side
    :param tree_file: Newick tree file
    :param tree: Parsed tree
    :param ingroup: List of labels to consider
    :return:
    """
    logger.info(f'Going to pick exemplars from {tree_file}')

    # First, get all leaves (including None), then filter out None values
    all_leaves = [tree.find_any(name=label) for label in ingroup]
    leaves = [leaf for leaf in all_leaves if leaf is not None]
    logger.debug(f'Leaves: {leaves}')

    # Get MRCA of ingroup leaves, and its immediate children
    mrca = tree.common_ancestor(leaves)
    children = mrca.clades
    representatives = []
    if len(children) == 2:
        for child in children:
            tips = child.get_terminals()
            dists = []
            for tip in tips:

                # Skip if focal tip is aberrant outgroup
                if tip.name not in ingroup:
                    logger.warning(f'Have aberrant outgroup {tip.name} in ingroup of tree {tree_file}')
                    continue

                # Calculate distance to mrca
                dist = tree.distance(tip, mrca)
                dists.append({'dist': dist, 'id': tip.name})

            # Get shallowest tip
            dists = sorted(dists, key=lambda x: x['dist'])
            representatives.append(dists[0]['id'])
    else:
        logger.warning(f'Ingroup root in {tree_file} is not bifurcating, will approximate rooting')
        return None
    logger.debug(representatives)
    return representatives


def reroot_on_split(tree_file, ingroup):
    """
    This function is visited in situations when pick_shallowest_tips encounters a basal trichotomy. This happens when
    raxml hasn't actually rooted the tree, which in turn happens when the ostensible outgroup is not actually
    monophyletic with respect to the ingroup. In the primates I saw this once, in the Neotropical monkeys squirrel
    monkeys. In that case, one of the genera within the family had tall tips and a long basal branch, while the
    selected outgroup (howlers) is shallow and therefore became enmeshed in the ingroup. The solution we implement here
    is to look for the internal edge that has the smallest bipartition that encloses all outgroup taxa and root on the
    midpoint of that edge.
    :param tree_file: a Newick tree file
    :param ingroup: List of labels to consider
    :return:
    """
    # Load Newick tree as unrooted with DendroPy, which computes split bitmasks when updating the bipartitions of a tree
    tree = dendropy.Tree.get(path=tree_file, schema="newick", rooting="force-unrooted")
    tree.update_bipartitions()

    # Make the set of outgroup nodes
    tip_set = set()
    for tip in tree.leaf_nodes():
        if tip.taxon.label not in ingroup:
            tip_set.add(tip.taxon)
    logger.debug(f'Constructed outgroup set {tip_set}')

    # Identify the smallest split that monophylizes tip_set
    smallest_split = None
    for edge in tree.postorder_edge_iter():
        if edge.bipartition.leafset_bitmask is not None:
            logger.debug(f'Processing bipartition {edge.bipartition.leafset_taxa(tree.taxon_namespace)}')
            if tip_set.issubset(edge.bipartition.leafset_taxa(tree.taxon_namespace)):
                logger.debug(f'Found smallest bipartition that has monophyletic {tip_set}')
                smallest_split = edge
                break

    # Found the smallest split, rerooting
    if smallest_split:
        new_length = smallest_split.length / 2
        tree.reroot_at_edge(smallest_split, length1=new_length, length2=new_length, update_bipartitions=True)

        # TODO overwrite *.bestTree with this tree for grafting into backbone
        tree.prune_taxa(tip_set)

        # Convert the DendroPy tree to a Newick string
        newick_string = tree.as_string(schema="newick")
        biopython_tree = read_newick(io.StringIO(newick_string), "newick")
        logger.debug(newick_string)
        logger.info('Rerooted, attempting to pick exemplars again.')
        return pick_shallowest_tips(tree_file, biopython_tree, ingroup)
    else:
        logger.error('Something bad is happening for which we have no solution')


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
    logger.info(f"Going to connect to database {args.database}")
    connection = sqlite3.connect(args.database)

    # Read FASTA, get list of ingroup tips
    seq_labels = get_ingroup_labels(args.inaln)

    # Maybe just include the whole file
    if len(seq_labels) < 3:
        logger.info(f'Infile {args.inaln} has fewer than 3 sequences, will include all')
        write_sequences(args.inaln, args.outaln, seq_labels)
    else:

        # ... or just the exemplars
        tree = read_newick(args.tree, 'newick')
        exemplars = pick_shallowest_tips(args.tree, tree, seq_labels)
        if exemplars is None:
            exemplars = reroot_on_split(args.tree, seq_labels)
        write_sequences(args.inaln, args.outaln, exemplars)

