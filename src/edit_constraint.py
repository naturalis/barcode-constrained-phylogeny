from Bio import SeqIO


def replace_newick(alignment_input, newick_input, newick_output):
    """
    :param dict: To retrieve the values from the dict (consisting of ott and barcode id)
    :return:
    """
    headers_dict = {}
    for record in SeqIO.parse(alignment_input, "fasta"):
        opentol_id = record.id.split('_')[0]
        headers_dict.setdefault(opentol_id, []).append(record.id)
    with open(newick_input, "r") as input:
        with open(newick_output, "w+") as output:
            ott = input.readline()
            ott = ott.replace(":0", "")
            for key in headers_dict.keys():
                ott = ott.replace(key, "(" + str(headers_dict[key]) + ")")
                ott = ott.replace("[", "")
                ott = ott.replace("]", "")
            output.write(ott)



if __name__ == '__main__':
    alignment_input = snakemake.input[0]
    newick_input = snakemake.input[1]
    newick_output = snakemake.output[0]

    replace_newick(alignment_input, newick_input, newick_output)