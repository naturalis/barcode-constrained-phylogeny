import logging
import tempfile

from Bio import SeqIO, AlignIO
from Bio.AlignIO import read as read_alignment
from subprocess import run
logging.basicConfig(level=snakemake.params.log_level)  # noqa: F821
logger = logging.getLogger(__name__)

# TODO add documentation and logging

def write_alignments(hmmfile, seqfile, outfile):
    """
    Read a fasta file per record, align record using hmmer with and without reverse complementing the sequence,
    check which one to keep(highest "8" count returned from with align_score()),
    and write all sequences together to a FASTA file.
    :param hmmfile:
    :param seqfile:
    :param outfile:
    :return:
    """
    with open(outfile, "w") as output:
        alignments = []
        for r in SeqIO.parse(seqfile, "fasta"):
            count_1, alignment_original = align_score(r, hmmfile)
            r.seq = r.seq.reverse_complement()
            count_2, alignment_reverse = align_score(r, hmmfile)
            if count_1 >= count_2:
                alignments.append(alignment_original)
            else:
                logger.info('Record %s is in reverse complement.' % r.id)
                alignments.append(alignment_reverse)
        AlignIO.write(alignments, output, "fasta")

def align_score(r, hmmfile):
    """
    Uses a Hidden Markov Model to align a sequence using hmmalign. This gives an stockholm file as output.
    Returns the number of "*" in the stockholm file. Reads the output file and saves its contents in the
    variable 'alignment'.
    :param r:
    :param hmmfile:
    :return: int number of "*" in a stockholm file, content of the hmmalignment output file
    """
    with tempfile.NamedTemporaryFile(mode='w+') as tempf, tempfile.NamedTemporaryFile(mode='w+') as tempf_rc:
        SeqIO.write(r, tempf.name, 'fasta')
        run(['hmmalign', '-o', tempf_rc.name, hmmfile, tempf.name])
        alignment = read_alignment(tempf_rc.name, "stockholm")
    return alignment.column_annotations['posterior_probability'].count('*'), alignment


if __name__ == '__main__':
    # TODO change to snakemake inputs
    hmmfile = snakemake.params.hmm  # noqa: F821
    in_file = snakemake.input[0]  # noqa: F821
    out_file = snakemake.output[0]  # noqa: F821
    write_alignments(hmmfile, in_file, out_file)