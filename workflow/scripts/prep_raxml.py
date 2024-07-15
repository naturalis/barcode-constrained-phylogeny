import util
import argparse
import sqlite3
import os

from Bio.AlignIO import read as read_alignment
from Bio.Phylo import read as read_newick, write as write_newick


"""
This script, `prep_raxml.py`, is responsible for preparing the input for RAxML, a tool for phylogenetic analysis.

The script performs the following steps:
1. Connects to the SQLite database.
2. Reads the input FASTA alignment file.
3. Creates a mapping between Open Tree of Life (OpenTOL) IDs and process IDs from the alignment.
4. Expands the input Newick tree based on the mapping. If a tip in the tree corresponds to multiple process IDs, it 
   splits the tip into multiple tips, each labeled with a process ID. If a tip does not correspond to any process ID, 
   it removes the tip.
5. Writes the expanded tree to an output file in Newick format.

The script uses command line arguments for the input Newick tree file, input FASTA alignment file, output Newick tree 
file, SQLite database file, and log level. The script is invoked by the Snakefile as a shell command with the required
arguments in the rule `prep_raxml`.
"""


def make_constraint(intree, outtree, processmap):
    """
    Makes a constraint tree compatible with the opentol input tree, but expanded
    to all process IDs that correspond with the longest sequence in each distinct
    BIN for that opentol species.
    :param intree:
    :param outtree:
    :param processmap:
    :return:
    """
    logger.info(f"Going to create constraint tree from {intree} to {outtree}")
    try:
        tree = read_newick(intree, 'newick')
    except ValueError:
        logger.info("No trees made for this family.")
        return

    # Map opentol_id to process_id, possibly adding tips if there are multiple process_ids
    # for this opentol_id (which means there are multiple BINs in this species)
    have_tree = True
    for tip in tree.get_terminals():
        if tip.name is not None:
            processes = processmap[tip.name]

            # Add and label tips if needed
            if len(processes) > 1:
                tip.split(n=len(processes))
                for child, process in zip(tip.get_terminals(), processes):
                    child.name = process
                    logger.info(f'Grafted child {process} to {tip.name}')
            else:
                tip.name = processes[0]
        else:
            logger.warning(f'Encountered None tip in {intree}, probably empty tree')
            have_tree = False

    # Write without branch lengths
    if have_tree:
        write_newick(tree, outtree, 'newick', plain=True)
    else:
        with open(outtree, 'a'):
            pass


def make_mapping(aln, conn):
    """
    Creates a one-to-many mapping between opentol_id (the keys in the dictionary) and a list
    of processids belonging to that opentol_id. Is used to expand the constraint tree, which
    has opentol leaves, so that it holds the distinct BINs as exemplified by the process ID
    producing the longest barcode sequence as leaves instead.
    :param aln:
    :param conn:
    :return: dict
    """
    logger.info('Looking up OpenToL IDs for process IDs in the alignment')
    map_dict = {}
    for seq in aln:
        process_id = seq.id

        # Because we are querying on the basis of the alignment, we may encounter cases
        # where there is a process_id without an opentol_id. However, this is not going
        # to be a problem later on because we only need to remap the tree, which is
        # based on opentol_ids.
        query = """
            SELECT t.opentol_id
            FROM barcode b, taxon t
            WHERE b.processid = ?
            AND b.taxon_id = t.taxon_id
            AND t.opentol_id IS NOT NULL        
        """
        cursor = conn.execute(query, (process_id,))
        record = cursor.fetchone()

        # Check if record is not empty
        if record is not None:
            opentol_id = f'ott{record[0]}'  # tree has ott prefixes
            if opentol_id not in map_dict:
                map_dict[opentol_id] = []
            map_dict[opentol_id].append(process_id)

    return map_dict


if __name__ == '__main__':

    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    parser.add_argument('-t', '--intree', required=True, help='Input Newick tree')
    parser.add_argument('-a', '--inaln', required=True, help='Input FASTA alignment')
    parser.add_argument('-o', '--outtree', required=True, help='Output Newick tree')
    parser.add_argument('-d', '--db', required=True, help="SQLite database")
    args = parser.parse_args()

    # Configure logger
    logger = util.get_formatted_logger('prep_raxml', args.verbosity)

    # Connect to the database (creates a new file if it doesn't exist)
    logger.info(f"Going to connect to database {args.db}")
    connection = sqlite3.connect(args.db)

    # Read input file, exit if it is empty
    logger.info(f'Going to read FASTA file {args.inaln}')
    infile = os.path.realpath(os.path.abspath(args.inaln))
    try:
        alignment = read_alignment(infile, 'fasta')
    except:
        logger.info("No records in the family.")
        alignment = []
        open(args.outtree, 'a')

    # Write Newick tree using process_id as labels, grafting subtended split BINS and outgroups
    logger.info(f'Going to expand OpenToL constraint to subtended process IDs')
    mapping = make_mapping(alignment, connection)
    make_constraint(args.intree, args.outtree, mapping)
