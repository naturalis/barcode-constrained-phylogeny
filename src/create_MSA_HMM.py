import logging
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.AlignIO import read as read_alignment, convert
from subprocess import run
import tempfile

logging.basicConfig(level=snakemake.params.log_level)  # noqa: F821
logger = logging.getLogger(__name__)
fasta_dir = snakemake.params.fasta_dir  # noqa: F821


def reverse_complement(record):
    """Return the reverse complement of a sequence record."""
    return SeqRecord(seq=record.seq.reverse_complement(), id=record.id)


def align_score(hmmfile, seqfile, outfile):
    """Run hmmalign on a sequence file and return the score."""
    run(['hmmalign', '-o', outfile, hmmfile, seqfile])
    alignment = read_alignment(outfile, "stockholm")
    return alignment.column_annotations['posterior_probability']


def main(args):
    """
    Main function that aligns the sequences. At this stage of the pipeline, the input file consists of sequences that:
        - are all the same marker (for which an HMM is provided)
        - all exceed the same minlength
        - could have <3 sequences in them (not usable for phylogeny)
        - might still be reverse-complemented
    In this function, we attempt to get all sequences to the same orientation and align them in the process by using
    a Hidden Markov Model (HMM) as produced by `hmmbuild` from the `hmmr` package. We do this by aligning the input
    sequence and its reverse complement, then count the number of columns that have the highest ranked bin of
    posterior probabilities, which are annotated with `*` in the Stockholm temp files produced by `hmmalign`.
    :param args: command line arguments dictionary as parsed by ArgParse
    """

    # Read the sequences from the input file
    sequences = list(SeqIO.parse(args.input, 'fasta'))
    best_alignments = []
    total_sequences = len(sequences)

    # Process each sequence
    for i, record in enumerate(sequences, start=1):

        with tempfile.NamedTemporaryFile(mode='w+') as tempf, tempfile.NamedTemporaryFile(mode='w+') as tempf_rc:

            # Write the sequence to a temporary file and align it
            SeqIO.write([record], tempf.name, 'fasta')
            score = align_score(args.hmm, tempf.name, 'temp.sto').count('*')

            # Write the reverse complement to a temporary file and align it
            SeqIO.write([reverse_complement(record)], tempf_rc.name, 'fasta')
            score_rc = align_score(args.hmm, tempf_rc.name, 'temp_rc.sto').count('*')

            # Choose the alignment with the higher score
            if score_rc > score:
                best_alignments.append(read_alignment('temp_rc.sto', 'stockholm'))
                logging.info(f"Using reverse complement for {i} ({score_rc})")
            else:
                best_alignments.append(read_alignment('temp.sto', 'stockholm'))
                logging.info(f"Using input direction for {i} ({score})")

        logging.info(f"Processed {i}/{total_sequences} sequences")

    # Convert the alignments to FASTA format and write them to the output file
    with tempfile.NamedTemporaryFile(mode='w+') as tempf:
        convert(best_alignments, 'stockholm', tempf.name, 'fasta')
        with open(args.output, 'w') as out_f:
            out_f.write(tempf.read())


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Align sequences with HMMER')
    parser.add_argument('--input', type=str, help='Input fasta file')
    parser.add_argument('--output', type=str, help='Output fasta file')
    parser.add_argument('--hmm', type=str, help='HMM file')
    parser.add_argument('--loglevel', type=str, default='INFO', help='Log level (default: INFO)')
    args = parser.parse_args()

    main(args)
