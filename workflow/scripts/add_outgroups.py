import argparse
import os
import logging
from Bio.AlignIO import read as read_alignment

logging.basicConfig()
logger = logging.getLogger('add_outgroups')


def add_sequence(sequence, distance, top_sequences, max_length=5):
    """
    Adds the candidate sequence object to the running tally of
    top_sequences with the
    :param sequence:
    :param distance:
    :param top_sequences
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
            file_name = os.path.join(folder_path, 'raxml-ready.fa')
            if os.path.exists(file_name):
                file_paths.append(file_name)

    return file_paths


def concat_seqfile(input_file, output_file, top_sequences):
    """
    Appends the sequences from the top list to infile
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
    with open(input_file, 'a') as handle:
        for record in top_sequences:
            logger.info(f'Distance {record.dist} - Sequence {record.seq.id}')
            handle.write(f'>{record.seq.id}\n')
            handle.write(f'{record.seq.seq}\n')


if __name__ == '__main__':

    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-i', '--infile', required=True, help='Input alignment FASTA file')
    parser.add_argument('-i', '--outfile', required=True, help='Output file name')
    parser.add_argument('-n', '--num_outgroups', required=True, help='Number of outgroups to add')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Read input file
    logger.info(f'Going to read FASTA file {args.infile}')
    infile = os.path.realpath(os.path.abspath(args.infile))
    alignment = read_alignment(infile, 'fasta')

    # This will store the top num_outgroups sequences
    # by shortest mean distance to ingroup
    nearest_neighbours = []
    for f in find_parallel_files(infile):

        # Don't compare infile with itself
        if os.path.realpath(os.path.abspath(f)) == infile:
            continue

        # Iterate over sequences in candidate file
        logger.info(f'Going to compare sequences from {f}')
        candidates = read_alignment(f, 'fasta')
        for candidate in candidates:

            # Calculate average distance to ingroup
            total_dist = 0
            counts = 0
            for seq in alignment:
                total_dist += calc_dist(candidate, seq)
                counts += 1
            mean_dist = total_dist / counts
            nearest_neighbours = add_sequence(candidate, mean_dist, nearest_neighbours, args.num_outgroups)

    # Write output
    logger.info(f'Have {len(nearest_neighbours)} nearest neighbor sequences')
    concat_seqfile(args.infile, args.outfile, nearest_neighbours)
