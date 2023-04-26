import argparse
import os
import sqlite3
from Bio import SeqIO
import re

### WITH FASTA

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()
parser.add_argument('-db', default="../../../data/databases/BOLD_COI-5P_barcodes.db",
                    help="Name of the the database file: {file_name}.db")
args = parser.parse_args()
print(args.db)


conn = sqlite3.connect(args.db)
# Create a cursor
cursor = conn.cursor()
file_list = os.listdir(par_path + "../../fasta/alignment/")


def replace_with_ott():
    # Get dict with different ott's
    """Replace the header in the phylip file.
    Uses the phylip file with the barcode id as input and creates a new file with an ott id as header.
    """
    dir_path= "../../fasta/alignment/"
    fasta_file = "Abacionidae_NT.fasta"
    new_file = "species.txt"
    ott = []
    dict = {}
    count = 0
    species = []
    with open("{}/{}".format(dir_path,fasta_file), "r+") as input: # File with barcode id as header
        with open("{}".format(new_file), "w+") as output:  # File with ott id as header
            for line in input:
                count += 1
                if line.__contains__(">"):
                    id = line.strip(">")
                    result = cursor.execute(
                        "SELECT taxon.opentol_id, barcode.barcode_id FROM taxon LEFT JOIN barcode ON barcode.taxon_id = "
                        "taxon.taxon_id WHERE barcode.barcode_id = {}".format(id))
                    results = result.fetchall()
                    if results[0][0] not in ott:
                        ott.append(results[0][0])
                    if results[0][0] in dict:
                        dict[results[0][0]].append(results[0][0]+"_"+ str(results[0][1]))
                    else:
                        dict[results[0][0]] = []
                    if results[0][0] not in species:
                        species.append(results[0][0])
            # Write list to string
            species = "\n".join(species)
            # Writing string to the outputfile
            output.write(species)
    output.close()
    return dict



def replace_newick(dict):
    """
    :param dict: To retrieve the values from the dict (consisting of ott and barcode id)
    :return:
    """
    with open("tree.nwk", "r") as input:
        with open("final_temp_subtree.nwk", "w+") as output:
            ott = input.readline()
            ott = ott.replace(":0", "")
            for key in dict.keys():
                ott = ott.replace(key, "(" + str(dict[key]) + ")")
                ott = ott.replace("[", "")
                ott = ott.replace("]", "")
            output.write(ott)


if __name__ == '__main__':
    dict = replace_with_ott()
    replace_newick(dict)



