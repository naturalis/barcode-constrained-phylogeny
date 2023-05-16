from Bio import SeqIO


def change_newick_backbone(file):
    """Obtain the newick created by megatree pruner.
    :param file: Megatree pruner newick file.
    :return: The string containing newick format.
    """
    with open(file, "r") as input:
        backbone = input.readline()
        return backbone


def make_dict(file):
    """Open the fasta file containing the representatives.
    Use the file to make a dict containing the representatives per ott.
    :param file: fasta file containing the representatives.
    :return: dictionary containing the representatives and the otts.
    """
    rep_dict = {}
    with open(file, "r") as rep:
        for record in SeqIO.parse(rep, "fasta"):
            # Make dict
            val = record.id
            ott = record.id.split("_")
            rep_dict[ott[0]] = val
    return rep_dict


def alter_newick_string(otol_newick, dict, file):
    """Use dictionary to alter the current newick string.
    Replace the otts with ott_barcode_id, which is the same format as in the fasta files.
    :param dict: To retrieve the values from the dict (consisting of ott and barcode id)
    :return: nothing, create a file with the altered newick string instead.
    """
    with open(file, "w+") as output:
        otol_newick = otol_newick.replace(":0", "")
        for key in dict.keys():
            otol_newick = otol_newick.replace(key, str(dict[key]))
        output.write(otol_newick)


if __name__ == '__main__':
    # Mock up for expected snakemake
    otol_newick = snakemake.input[0]
    representatives = snakemake.input[1] # in fasta format
    altered_newick = snakemak.output[0]
    newick = change_newick_backbone(otol_newick)
    dict = make_dict(representatives)
    alter_newick_string(newick, dict, altered_newick)