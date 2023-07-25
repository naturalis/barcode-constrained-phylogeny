import logging
from subprocess import run
import os
logging.basicConfig(level=snakemake.params.log_level)  # noqa: F821
logger = logging.getLogger(__name__)

def create_large_file(sequences, directory):
    """Create a file containing several sequences.
    A cut off is used, because only a few hundred sequences are necessary to create the file.
    :param sequences: The file to be created.
    """
    # iterate over some files in directory
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)   # All the files in the family directory
        try:
            with open(f'{sequences}', "a+") as output, open(f, "r") as file:
                for line in file:
                    output.write(line)
            if len(os.listdir(directory)) > 30 and filename == os.listdir(directory)[30]:
                break
        except FileNotFoundError:
            print("File coult not be found.")
        file.close()
    logger.info(f"File {sequences} is made.")
    output.close()


def create_alignment_file(inputfile, msa_alignment):
    """Create an alignment file which can be used to create a HMM.
    Clustalw is used to create an alignment.
    An inputfile with different sequences is used to create an alignment file.
    :return: Nothing, create a file instead.
    """
    try:
        # Create file msa_alignment based on the inputfile
        run(['clustalw', f"-infile={inputfile}", f"-outfile={msa_alignment}"])
    except FileNotFoundError:
        print("File cannot be found.")
    except:
        print("An unknown error occurred.")
    logger.info(f"An alignment file {msa_alignment} was created.")


def perform_hmmbuild(hmmfile, alignment_file):
    """Create an alignment file using HMMER tool.
    The hmmfile is used to produce the alignment file.
    :param hmmfile: file containing the HMM.
    :param alignment_file: Multiple sequence alignment file to be produced.
    :return: Nothing, but create a file instead.
    """
    try:
        run(['hmmbuild', hmmfile, alignment_file])  # Create
    except FileNotFoundError:
        print("File cannot be found.")
    except:
        print("An unknown error occurred.")
    logger.info(f"Creation of hmm file {hmmfile} was successfull.")


if __name__ == '__main__':
    dir = '../data/fasta/family/'
    #dir = snakemake.input.[0] # noqa F821  # Not sure why it does not work
    sequence_file = snakemake.output[0]  # noqa:  F821
    alignment_file = snakemake.output[1]  # noqa: F821
    hmmfile = snakemake.params.hmm  # noqa: F821

    create_large_file(sequence_file, dir)
    create_alignment_file(sequence_file, alignment_file)
    perform_hmmbuild(hmmfile, alignment_file)