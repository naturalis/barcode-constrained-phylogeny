from Bio import SeqIO
import os


def get_representatives_as_list(representatives):
    list = []
    with open(representatives, "r") as input:
        lines = input.readlines()
        for line in lines:
            line = line.rstrip("\n")
            list.append(line)
    return list
def loop_over_fam(path,rep):
    for family in path:
        print(family)
        if family == path[0]:
            mode = "w"
            print("write")
        else:
            mode = "a"
        if family == path[0]:
            print("same")
            print(rep)
            representatives_to_file(rep, "data/fasta/{}/{}.fasta".format(family,family), outputfile, mode)


def representatives_to_file(representatives, inputfile, output, mode):
    """The sequences with the highest distance (which are also in de opentree of life db) are written to a fasta file.
    """

    with open(inputfile, "r") as input:
        with open(output, mode) as outputs:
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
    path2 = os.getcwd()  # Get current working directory
    path = snakemake.input[0]  # noqa: F821
    total_path = os.path.join(path2, path)
    family_list = os.listdir(os.path.join(path2, path))
    repr = snakemake.input[1]  # noqa: F821
    representatives = get_representatives_as_list(repr)
    outputfile = "data/fasta/backbone/rep.fasta"
    loop_over_fam(family_list, representatives)
    #outputfile = snakemake.output[0]     # noqa: F821

    #representatives_to_file(representatives, inputfile, outputfile)

