import logging
import tempfile
import argparse
import sqlite3
import os
import util

from Bio import SeqIO
from Bio.AlignIO import read as read_alignment
from subprocess import run


"""
This script, `msa_hmm.py`, is responsible for aligning sequences using a Hidden Markov Model (HMM) and writing the 
aligned sequences to an output file.

The script performs the following steps:
1. Merges the provided input files and maps the process IDs to barcode IDs, which are shorter. This is useful because 
   the alignment steps truncate IDs. At a later stage we map them back to process IDs. In both cases this is done via 
   the database.
2. Corrects any sequences that are reverse complemented. It does this by aligning the sequence with and without reverse 
   complementing the sequence, and keeps the sequence that gives the highest alignment score.
3. Aligns the reverse complement-corrected sequences. Remaps barcode_id to process_id using a database lookup. Writes 
   output to outfile.

The script uses command line arguments for the HMM model file, log level, input unaligned ingroup FASTA file, unaligned 
outgroup FASTA file, output aligned FASTA file, and SQLite database. The script is invoked by the Snakefile as a shell
command with the required arguments in the rule `msa_hmm`.
"""


def correct_revcom(hmmfile, inseqs):
    """
    Read a list of SeqIO records using hmmer with and without reverse complementing the sequence,
    check which one to keep(highest "*" count returned from with align_score()),
    and return all in the right orientation
    :param hmmfile: HMM model file
    :param inseqs: list of sequences
    :return:
    """

    # Read each sequence record from fasta file
    sequences = []
    log_reverse = 0
    for record in inseqs:

        # Run align_score with original sequence
        count_1, alignment_original = align_score(record, hmmfile)
        logger.debug(f'Forward alignment score {count_1} for {record.id}')

        # Reverse complement the sequence and align
        record.seq = record.seq.reverse_complement()
        count_2, alignment_reverse = align_score(record, hmmfile)
        logger.debug(f'Reverse complemented alignment score {count_2} for {record.id}')

        # Check if the first count is higher
        if count_1 >= count_2:

            # Reverse the reverse complement
            record.seq = record.seq.reverse_complement()
            sequences.append(record)
        else:

            # Add alignment to list
            sequences.append(record)

            # Log counter
            log_reverse += 1

    logger.info(f'Corrected {log_reverse} reverse complemented sequences out of {len(sequences)}')
    return sequences


def align_write(sequences, outfile, conn):
    """
    Aligns the reverse complement-corrected sequences. Remaps barcode_id to process_id using a database
    lookup by way of the provided conn. Writes output to outfile.
    :param sequences:
    :param outfile:
    :param conn:
    :return:
    """

    # Save the sequence to the temporary file
    with open(f'{outfile}.tmp1', mode='w+') as temp_fasta:
        for seq in sequences:
            temp_fasta.write(f'>{seq.id}\n')
            temp_fasta.write(f'{seq.seq}\n')

    # Align with hmmalign, capture and parse output as phylip
    run(['hmmalign', '--trim', '-o', f'{outfile}.tmp2', '--outformat', 'phylip', hmmfile, f'{outfile}.tmp1'])
    try:
        aligned = read_alignment(f'{outfile}.tmp2', 'phylip')
    except ValueError:
        logger.info("No sequences found.")
        aligned = []
    os.remove(f'{outfile}.tmp1')
    os.remove(f'{outfile}.tmp2')

    # Map barcode_id to process_id and write to outfile
    with open(f'{outfile}', mode='a') as output:
        for seq in aligned:

            # Need to lookup the process IDs in the database
            bid = seq.id
            process_id = conn.execute(f'SELECT processid FROM barcode WHERE barcode_id={bid}').fetchone()
            output.write(f'>{process_id[0]}\n')
            output.write(f'{seq.seq}\n')

    logger.info(f'Wrote/appended aligned sequences to {outfile}')


def align_score(record, hmmfile):
    """
    Uses a Hidden Markov Model to align a sequence using hmmalign. This gives an stockholm file as output.
    Returns the number of "*" in the stockholm file. Reads the output file and saves its contents in the
    variable 'alignment'.
    :param record: a sequence from the fasta input
    :param hmmfile: HMM model file
    :return: int number of "*" in a stockholm file, content of the hmmalignment output file
    """
    # Open two temporary files to store a fasta sequence and a stockholm sequence
    with tempfile.NamedTemporaryFile(mode='w+') as temp_fasta, tempfile.NamedTemporaryFile(mode='w+') as temp_stockholm:

        # Save the sequence to the temporary file
        SeqIO.write(record, temp_fasta.name, 'fasta')

        # Run hmm align, read the aligned sequence
        run(['hmmalign', '--trim', '-o', temp_stockholm.name, hmmfile, temp_fasta.name])
        alignment = read_alignment(temp_stockholm.name, "stockholm")
        seq = alignment[0]

    # Return probability colum of the stockholm file
    quality_string = alignment.column_annotations['posterior_probability']

    count = 0

    # Count . and * characters
    dot_count = quality_string.count('.')
    star_count = quality_string.count('*')

    # Give value 0 to . and value 10 to *
    count += dot_count * 0
    count += star_count * 10

    # Add all numbers
    digit_sum = sum(int(char) for char in quality_string if char.isdigit())

    # Add numbers to the values calculated from . and *
    count += digit_sum

    # Calculate average count to account for gaps (0's)
    average_count = count / len(quality_string)
    return average_count, seq


def merge_files(ingroup, outgroup, conn):
    """
    Merges the provided input files and maps the process IDs to barcode IDs, which are shorter.
    This is useful because the alignment steps truncate IDs. At a later stage we map them back
    to process IDs. In both cases this is done via the database, hence the connection object needs
    to be provided here. In the case of the ingroup file, the process ID is assumed to be the 3rd
    item in a pipe-separated list on the definition line. In the outgroup file, it is the entire
    first word of the definition line, i.e. the record ID according to SeqIO.
    :param ingroup: FASTA file
    :param outgroup: FASTA file
    :param conn: Database connection
    """
    merged = []

    # ingroup has defline where process ID is 3rd element in pipe-separated list
    for record in SeqIO.parse(ingroup, 'fasta'):
        procid = record.id.split('|')[2]
        record.id = procid
        merged.append(record)

    # outgroup already has this, just needs appending
    for record in SeqIO.parse(outgroup, 'fasta'):
        merged.append(record)

    # map all to barcode ID
    for record in merged:
        procid = record.id
        query = f'select barcode_id from barcode where processid="{procid}";'
        record.id = str(conn.execute(query).fetchone()[0])

    # return result
    return merged


if __name__ == '__main__':

    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-m', '--model', required=True, help='Location of HMM file, should match the marker')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    parser.add_argument('-i', '--ingroup', required=True, help='Input unaligned ingroup FASTA file')
    parser.add_argument('-g', '--outgroup', required=True, help='Unaligned outgroup FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output aligned FASTA file')
    parser.add_argument('-d', '--db', required=True, help="SQLite database")
    args = parser.parse_args()

    # Configure logger
    logger = util.get_formatted_logger('msa_hmm', args.verbosity)

    # Connect to the database (creates a new file if it doesn't exist)
    logger.info(f"Going to connect to database {args.db}")
    connection = sqlite3.connect(args.db)

    # Merge input files with deflines mapped to barcode_id
    logger.info(f"Going to merge input files {args.ingroup} and {args.outgroup}")
    uncorrected = merge_files(args.ingroup, args.outgroup, connection)

    # Report the location of the HMM file
    hmmfile = args.model
    logger.info(f"Going to use Hidden Markov Model from {hmmfile}")

    # Do the reverse complement correction
    logger.info('Going to correct any reverse-complemented sequences')
    corrected = correct_revcom(hmmfile, uncorrected)

    # Write alignment
    align_write(corrected, args.output, connection)
