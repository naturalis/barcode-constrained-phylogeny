import argparse
import os
import sqlite3
from Bio import SeqIO

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
print(file_list)

def get_length():
    # Check if same length
    fasta_file = "Abacionidae_NT.fasta"
    with open("{}".format(fasta_file), "r+") as input:
        for line in input:
            line = line.rstrip("\n")
            line = line.split(" ")
            #print(line[1])
            print(len(line[1]))

def replace_with_ott():
    """Replace the header in the phylip file.
    Uses the phylip file with the barcode id as input and creates a new file with an ott id as header.
    """
    dir_path= "/"
    phylip_file = "Abacionidae_NT.fasta.raxml.reduced.phy"
    new_file = "Aba.phy"
    ott = []
    with open("{}/{}".format(dir_path,phylip_file), "r+") as input: # File with barcode id as header
        with open("{}/{}".format(dir_path, new_file), "w+") as output:  # File with ott id as header
            header = input.readline()   # To ignore the first line
            output.write(header)
            for line in input:
                line = line.rstrip("\n")
                line = line.split(" ")
                result = cursor.execute(
                    "SELECT taxon.opentol_id, taxon FROM taxon LEFT JOIN barcode ON barcode.taxon_id = "
                                        "taxon.taxon_id WHERE barcode.barcode_id = {}".format(line[0]))
                results = result.fetchall()
                output.write(results[0][0]+ "_" + results[0][1] + " " +line[1] + "\n")
                if results[0][0] not in ott:
                    ott.append(results[0][0])
    output.close()
    print(ott)



if __name__ == '__main__':
    #get_length()
    replace_with_ott()
