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
                        "SELECT taxon.opentol_id, barcode.barcode_id FROM taxon LEFT JOIN barcode ON barcode.taxon_id = "
                        "taxon.taxon_id WHERE barcode.barcode_id = {}".format(id))
                    results = result.fetchall()
                    output.write(">" +results[0][0] + "_" + str(results[0][1]) + "\n")
                    if results[0][0] not in ott:
                        ott.append(results[0][0])
                    if results[0][0] in dict:
                        dict[results[0][0]].append(results[0][0]+"_"+ str(results[0][1]))
                    else:
                        dict[results[0][0]] = []
                else:
                    output.write(line)

    output.close()
    print(dict)
    return dict

def remove_newick(dict):
    """Remove ott which are not found in the alignment file
    :param dict: contains the ott in the alignment file
    :return: writes a newick in a file
    """
    not_in_alignment = []
    with open("temp.nwk", "r") as input:
        with open("test.nwk", "w+") as output:
            line = input.readline()
            line = line.split(",")
            for ott in line:
                ott = ott.replace("(","")
                ott = ott.replace(")"," ")
                ott = ott.split(" ")
                for i in range(len(ott)):
                    if ott[i] not in dict:
                        not_in_alignment.append(ott[i])
            line = ", ".join(line)
            print(line)
            for i in not_in_alignment:
                line = line.replace(i, "")
            print(line)
            output.write(line)



def replace_newick(dict):
    """
    :param dict:
    :return:
    """
    with open("test.nwk", "r") as input:
        with open("temp_subtree", "w") as output:
            ott = input.readline()
            for key in dict.keys():
                ott = ott.replace(key,str(dict[key]))
                ott = ott.replace("[","")
                ott = ott.replace("]", "")
            output.write("(")
            output.write(ott)
            output.write(");")
            print(ott)


if __name__ == '__main__':
    #get_length()
    dict = replace_with_ott()
    remove_newick(dict)
    replace_newick(dict)
