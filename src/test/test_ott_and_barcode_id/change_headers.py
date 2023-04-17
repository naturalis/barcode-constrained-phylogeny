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
#print(par_path)
file_list = os.listdir(par_path + "../../fasta/alignment/")
#print(file_list)

def replace_with_ott():
    """Replace the header in the phylip file.
    Uses the phylip file with the barcode id as input and creates a new file with an ott id as header.
    """
    dir_path= "../../fasta/alignment/"
    test_path=""
    fasta_file = "Abacionidae_NT.fasta"
    new_file = "Aba.fasta"
    ott = []
    dict = {}
    count = 0
    with open("{}/{}".format(dir_path,fasta_file), "r+") as input: # File with barcode id as header
        with open("{}".format(new_file), "w+") as output:  # File with ott id as header
            for line in input:
                count += 1
                if line.__contains__(">"):
                    id = line.strip(">")
                    result = cursor.execute(
                        "SELECT taxon.opentol_id, taxon FROM taxon LEFT JOIN barcode ON barcode.taxon_id = "
                        "taxon.taxon_id WHERE barcode.barcode_id = {}".format(id))
                    results = result.fetchall()
                    output.write(">" +results[0][0] + "_" + str(count) + "\n")
                    if results[0][0] not in ott:
                        ott.append(results[0][0])
                    if results[0][0] in dict:
                        dict[results[0][0]].append(results[0][0]+"_"+str(count))
                    else:
                        dict[results[0][0]] = []
                else:
                    output.write(line)

    output.close()
    print(dict)
    return dict

def replace_newick(dict):

    with open("temp.nwk", "r") as input:
        with open("temp_subtree", "w") as output:
            ott = input.readline()
            for key in dict.keys():
                ott = ott.replace(key,str(dict[key]))
                ott = ott.replace("[","")
                ott = ott.replace("]", "")
            output.write(ott)
            print(ott)


if __name__ == '__main__':
    #get_length()
    dict = replace_with_ott()
    replace_newick(dict)
