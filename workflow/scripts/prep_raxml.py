import logging
import argparse
import sqlite3
import os

from Bio.AlignIO import read as read_alignment
from Bio.Phylo import read as read_newick, write as write_newick

logging.basicConfig()
logger = logging.getLogger('prep_raxml')


def concat_seqfile(input_file, output_file, top_sequences):
    """
    Concatenates the sequences from input file and top list into the output file
    :param input_file:
    :param output_file:
    :param top_sequences:
    :return:
    """
    logger.info(f"Going to write sequences to {output_file}")

    # Copy input to output
    with open(input_file, 'r') as i, open(output_file, 'w') as o:
        for line in i:
            o.write(line)

    # Append sequences
    with open(output_file, 'a') as handle:
        for record in top_sequences:
            logger.info(f'Distance {record["dist"]} - Sequence {record["seq"].id}')
            handle.write(f'>{record["seq"].id}\n')
            handle.write(f'{record["seq"].seq}\n')


def add_sequence(sequence, distance, top_sequences, max_length):
    """
    Adds the candidate sequence object to the running tally of
    top_sequences with the
    :param sequence:
    :param distance:
    :param top_sequences:
    :param max_length:
    :return:
    """

    # Add, sort, truncate
    seq_record = {
        'seq': sequence,
        'dist': distance
    }
    top_sequences.append(seq_record)
    top_sequences = sorted(top_sequences, key=lambda x: x['dist'])
    top_sequences = top_sequences[:max_length]
    return top_sequences


def calc_dist(seq1, seq2):
    """
    Calculates the edit distance between two sequences. Skips gaps '-'.
    :param seq1:
    :param seq2:
    :return:
    """

    # Calculate mismatches excluding gaps
    mismatches = 0
    non_gap_positions = 0
    for char1, char2 in zip(seq1.seq, seq2.seq):
        if char1 != '-' and char2 != '-':
            non_gap_positions += 1
            if char1 != char2:
                mismatches += 1

    # Calculate distance only for non-gap positions
    if non_gap_positions > 0:
        distance = mismatches / non_gap_positions
    else:
        distance = None

    # Return results
    logger.debug(f'Pairwise distance {seq1.id} <=> {seq2.id} = {distance}')
    return distance


def is_valid_folder(folder_name):
    """
    Checks if folder_name is a folder with naming pattern n-of-m
    :param folder_name:
    :return:
    """
    if '-' in folder_name and 'of' in folder_name:
        parts = folder_name.split('-of-')
        n, m = int(parts[0]), int(parts[1])
        if n > m:
            return False
        else:
            return True
    return False


def find_parallel_files(reference):
    """
    Finds raxml-ready files
    :param reference:
    :return:
    """

    # Get the base directory
    base_dir = os.path.dirname(os.path.abspath(reference))

    # Get the parent directory of the base directory
    parent_dir = os.path.dirname(base_dir)

    # List to store file paths
    file_paths = []

    # Iterate over all directories in the parent directory
    for folder in os.listdir(parent_dir):
        folder_path = os.path.join(parent_dir, folder)
        if os.path.isdir(folder_path) and is_valid_folder(folder):
            file_name = os.path.join(folder_path, 'aligned.fa')
            if os.path.exists(file_name):
                file_paths.append(file_name)

    return file_paths


def find_nearest_neighbors(reference, aln, max_length):
    """
    Finds max_length nearest sequences in the sister files relative to the reference alignment file
    :param reference:
    :param aln:
    :param max_length:
    :return:
    """
    nearest_neighbours = []

    # Iterate over alignments
    for f in find_parallel_files(reference):

        # Don't compare infile with itself
        if os.path.realpath(os.path.abspath(f)) == reference:
            continue

        # Iterate over sequences in candidate file
        logger.debug(f'Going to compare sequences from {f}')
        candidates = read_alignment(f, 'fasta')
        for candidate in candidates:

            # Calculate average distance to ingroup
            total_dist = 0
            counts = 0
            for seq in aln:
                total_dist += calc_dist(candidate, seq)
                counts += 1
            mean_dist = total_dist / counts
            nearest_neighbours = add_sequence(candidate, mean_dist, nearest_neighbours, max_length)
    return nearest_neighbours


def make_constraint(intree, outtree, processmap, top_sequences, aln):
    """
    Makes a constraint tree compatible with the opentol input tree, but expanded
    to all process IDs that correspond with the longest sequence in each distinct
    BIN for that opentol species.
    :param intree:
    :param outtree:
    :param processmap:
    :param top_sequences:
    :param aln:
    :return:
    """
    logger.info(f"Going to create constraint tree from {intree} to {outtree}")
    tree = read_newick(intree, 'newick')

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
                    logger.info(f'Added child {process} to {tip.name}')
            else:
                tip.name = processes[0]
        else:
            logger.warning(f'Encountered None tip in {intree}, probably empty tree')
            have_tree = False

    # Add ingroup sequences not in mapping to the root?
    # Add outgroups, needs new node below root?

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
    parser.add_argument('-a', '--inaln', required=True, help='Input aligned PHYLIP file')
    parser.add_argument('-o', '--outtree', required=True, help='Output Newick tree')
    parser.add_argument('-f', '--outaln', required=True, help="Output FASTA alignment")
    parser.add_argument('-n', '--num_outgroups', required=True, type=int, help='Number of outgroups to add')
    parser.add_argument('-d', '--db', required=True, help="SQLite database")
    args = parser.parse_args()

    # Configure logger
    logger.setLevel(args.verbosity)

    # Connect to the database (creates a new file if it doesn't exist)
    logger.info(f"Going to connect to database {args.db}")
    connection = sqlite3.connect(args.db)

    # Read input file, exit if it is empty
    logger.info(f'Going to read FASTA file {args.inaln}')
    infile = os.path.realpath(os.path.abspath(args.inaln))
    alignment = read_alignment(infile, 'fasta')

    # This will return the top num_outgroups sequences by shortest mean distance to ingroup
    nn = find_nearest_neighbors(infile, alignment, args.num_outgroups)

    # Write FASTA output using process_id is label, including outgroups
    logger.info(f'Have {len(nn)} nearest neighbor sequences')
    concat_seqfile(infile, args.outaln, nn)

    # Write Newick tree using process_id as labels, grafting subtended split BINS and outgroups
    logger.info(f'Going to expand OpenToL constraint to subtended process IDs')
    mapping = make_mapping(alignment, connection)
    make_constraint(args.intree, args.outtree, mapping, nn, alignment)
