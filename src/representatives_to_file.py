from Bio import SeqIO
import os


def get_representatives_as_list(representatives):
    """Open file with the representatives.
    Read the file. 
    Add the representatives to a list.
    :return: list with the representatives"""
    list = []
    with open(representatives, "r") as input:
        lines = input.readlines()
        for line in lines:
            line = line.rstrip("\n")
            list.append(line)
    return list


def representatives_to_file(representatives, inputfile, output):
    """The sequences with the highest distance (which are also in de opentree of life db) are written to a fasta file.
    Read the inputfile.
    Compare the ott_id with the ott in the representatives list.
    If they are the same write to file.
    :return: nothing. Write to a file instead.
    """

    with open(inputfile, "r") as input:
        with open(output, "a+") as outputs:
            for record in SeqIO.parse(input, "fasta"):
                # Iterate over the representatives to find their corresponding sequences
                for header in representatives:
                    print(record.id)
                    id = record.id.split("|")
                    print(id[1])
                    if id[1] == header:
                        # Write to fasta file
                        SeqIO.write(record, outputs, 'fasta')


if __name__ == '__main__':
    family_fasta = snakemake.input[0]   # noqa: F821
    repr = snakemake.input[1]  # noqa: F821
    outputfile = snakemake.output[0] # noqa: F821
    representatives = get_representatives_as_list(repr)
    representatives_to_file(family_fasta, representatives, outputfile)

