import logging
import tempfile
import argparse

from Bio import SeqIO
from Bio.AlignIO import read as read_alignment
from subprocess import run

logging.basicConfig()
logger = logging.getLogger('msa_hmm')


def write_alignments(hmmfile, seqfile, outfile):
    """
    Read a fasta file per record, align record using hmmer with and without reverse complementing the sequence,
    check which one to keep(highest "*" count returned from with align_score()),
    and write all sequences together to a FASTA file.
    :param hmmfile: HMM model file
    :param seqfile: name of FASTA file input
    :param outfile: name of alignment fasta file output
    :return:
    """
    logger.info(f'Aligning sequences in FASTA file {seqfile}')

    # Read each sequence record from fasta file
    alignments = []
    log_reverse = 0
    for record in SeqIO.parse(seqfile, "fasta"):

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
            alignments.append(record)
        else:

            # Add alignment to list
            alignments.append(record)

            # Log counter
            log_reverse += 1

    logger.info(f'Corrected {log_reverse} reverse complemented sequences out of {len(alignments)}')

    # Save the sequence to the temporary file
    with open(f'{outfile}-revcomfix.fa', mode='a+') as temp_fasta:
        for seq in alignments:
            temp_fasta.write(f'>{seq.description}\n')
            temp_fasta.write(f'{seq.seq}\n')
    run(['hmmalign', '--trim', '--outformat', 'phylip', '-o', outfile, hmmfile, f'{outfile}-revcomfix.fa'])

    logger.info(f'Wrote all aligned sequences to {outfile}')


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


if __name__ == '__main__':

    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-m', '--model', required=True, help='Location of HMM file, should match the marker')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    parser.add_argument('-i', '--input', required=True, help='Input unaligned FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output aligned FASTA file')
    args = parser.parse_args()

    # Configure logger
    logger.setLevel(args.verbosity)

    # Check the HMM file
    hmmfile = args.model
    logger.info(f"Going to use Hidden Markov Model from {hmmfile}")

    # Check the input file
    in_file = args.input
    logger.info(f"Going to align sequences from input file {in_file}")

    # Announce the output file
    out_file = args.output
    logger.info(f"Will write alignment to output file {out_file}")

    write_alignments(hmmfile, in_file, out_file)
