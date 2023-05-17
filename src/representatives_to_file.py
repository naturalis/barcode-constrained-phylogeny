from Bio import SeqIO
def representatives_to_file(representatives, inputfile, output):
    """The sequences with the highest distance (which are also in de opentree of life db) are written to a fasta file.
    """

    with open(inputfile, "r") as input:
        with open(output, "a") as outputs:
            for record in SeqIO.parse(input, "fasta"):
                # Iterate over the representatives to find their corresponding sequences
                for header in representatives:
                    if record.id == header:
                        # Write to fasta file
                        SeqIO.write(record, outputs, 'fasta')


if __name__ == '__main__':
    inputfile = snakemake.input[0]  # noqa: F821
    representatives = snakemake.input[1] # noqa: F821
    outputfile = snakemak.output[2] # noqa: F821
    representatives_to_file(representatives, inputfile, outputfile)

