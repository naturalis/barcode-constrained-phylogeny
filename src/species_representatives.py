""""Getting ott_ids from a list of representatives.

Representatives are in ott_barcode format.
Getting the first part (thus only getting the ott) and writing it to an output file.
The outputfile is species.txt
This file is required for using megatree pruner
"""


def get_species(representatives):
    with open(outputfile, "a") as output:
        for header in representatives:
            header = header.split("_")
            output.write(header[0]+"\n")


if __name__ == '__main__':
    # For personal use, not for snakemake
    outputfile = "representatives.txt"
    representatives = []
    get_species(representatives)