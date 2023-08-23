from Bio import SeqIO


def replace_newick(alignment_input, newick_input, newick_output):
    """Iterate over the entries in the file {alignemnt_input}.
    From the file change the header. 
    Open the newick file replace the :0 and ; with "". 
    If there is more than one entry with the same ott use ( ) to group them together.
    :param dict: To retrieve the values from the dict (consisting of ott and barcode id)
    :return: nothing. Write output to file instead.
    """
    headers_dict = {}
    for record in SeqIO.parse(alignment_input, "fasta"):
        opentol_id = record.id.split('_')[0]
        headers_dict.setdefault(opentol_id, []).append(record.id)
    with open(newick_input, "r") as input:
        with open(newick_output, "w+") as output:
            ott = input.readline()
            ott = ott.replace(":0", "")
            ott = ott.replace(";", "")
            for key in headers_dict.keys():
                if len(headers_dict[key]) > 1:
                    ott = ott.replace(key, "(" + str(headers_dict[key]) + ")")
                else:
                    ott = ott.replace(key, str(headers_dict[key]))
                ott = ott.replace("[", "")
                ott = ott.replace("]", "")
            output.write(ott)


if __name__ == '__main__':
    alignment_input = snakemake.input[0] # noqa: F821
    newick_input = snakemake.input[1] # noqa: F821
    newick_output = snakemake.output[0] # noqa: F821

    replace_newick(alignment_input, newick_input, newick_output)
