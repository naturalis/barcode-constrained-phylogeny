import logging
import tempfile

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO import read as read_alignment
from subprocess import run

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=snakemake.params.log_level)  # noqa: F821
logger = logging.getLogger(__name__)


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
    logger.info("Aligning sequences in FASTA file %s" % in_file)
    alignments = []
    log_reverse = 0
    # Read each sequence record from fasta file
    for record in SeqIO.parse(seqfile, "fasta"):
        # Run align_score with original sequence
        count_1, alignment_original = align_score(record, hmmfile)
        # Reverse complement the sequence
        record.seq = record.seq.reverse_complement()
        # run align_score with reverse complemented sequence
        count_2, alignment_reverse = align_score(record, hmmfile)
        logger.debug('Record %s original alignment score:%.2f\nRecord %s reverse complement alignment score:%.2f'
                     % (record.id, count_1, record.id, count_2))
        # Check if the first count is higher
        if count_1 >= count_2:
            # Add alignment to list
            record.seq = record.seq.reverse_complement()
            alignments.append(record)
        else:
            logger.debug('Record %s is in reverse complement.' % record.id)
            # Add alignment to list
            alignments.append(record)
            # Log counter
            log_reverse += 1
    logger.info("There were %i sequences that were reverse complemented before alignment." % log_reverse)
    # Rewrite stockholm file to a fasta alignment file
    # msa = MultipleSeqAlignment(alignments)
    # AlignIO.write(alignments, output, "fasta")
    with tempfile.NamedTemporaryFile(mode='w+') as temp_fasta:
        # Save the sequence to the temporary file
        SeqIO.write(alignments, temp_fasta.name, 'fasta')
        logger.info("Alignn all family barcodes in one file.")
        # Run hmm align (arguments: output file, model file, input file)
        run(['hmmalign','-o', outfile, hmmfile, temp_fasta.name])


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
        # Run hmm align (arguments: output file, model file, input file)
        run(['hmmalign', '-o',temp_stockholm.name, hmmfile, temp_fasta.name])
        # Read the stockholm alignment
        alignment = read_alignment(temp_stockholm.name, "stockholm")
        test = SeqIO.read(temp_stockholm.name, "stockholm")
    # Return probability colum of the stockholm file
    quality_string = alignment.column_annotations['posterior_probability']

    count = 0
    # Count . and * characters
    dot_count = quality_string .count('.')
    star_count = quality_string .count('*')
    # Give value 0 to . and value 10 to *
    count += dot_count * 0
    count += star_count * 10
    # Add all numbers
    digit_sum = sum(int(char) for char in quality_string if char.isdigit())
    # Add numbers to the values calculated from . and *
    count += digit_sum
    # Calculate average count to account for gaps (0's)
    average_count = count/len(quality_string)
    return average_count, test




if __name__ == '__main__':
    hmmfile = snakemake.params.hmm  # noqa: F821
    in_file = snakemake.input[0]  # noqa: F821
    out_file = snakemake.output[0]  # noqa: F821
    write_alignments(hmmfile, in_file, out_file)
