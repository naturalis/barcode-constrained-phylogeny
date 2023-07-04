import logging
import tempfile

from Bio import SeqIO, AlignIO
from Bio.AlignIO import read as read_alignment
from subprocess import run
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
    with open(outfile, "w") as output:
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
            # Check if the first count is higher
            if count_1 >= count_2:
                # Add alignment to list
                alignments.append(alignment_original)
            else:
                logger.debug('Record %s is in reverse complement.' % record.id)
                # Add alignment to list
                alignments.append(alignment_reverse)
                # Log counter
                log_reverse += 1
        logger.info("There were %i sequences that were reverse complemented before alignment." % log_reverse)
        # Rewrite stockholm file to a fasta alignment file
        AlignIO.write(alignments, output, "fasta")


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
        run(['hmmalign', '-o', temp_stockholm.name, hmmfile, temp_fasta.name])
        # Read the stockholm alignment
        alignment = read_alignment(temp_stockholm.name, "stockholm")
    # Return number of "*" from stockholm alignment and the alignment itself
    return alignment.column_annotations['posterior_probability'].count('*'), alignment


if __name__ == '__main__':
    hmmfile = snakemake.params.hmm  # noqa: F821
    in_file = snakemake.input[0]  # noqa: F821
    out_file = snakemake.output[0]  # noqa: F821
    write_alignments(hmmfile, in_file, out_file)