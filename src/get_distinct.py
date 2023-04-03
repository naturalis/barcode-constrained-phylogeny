import argparse
import os
import sqlite3
from Bio import SeqIO

# Testfile to see if the sequences from the same ott are the same or distinct
par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()
parser.add_argument('-db', default="../data/databases/BOLD_COI-5P_barcodes.db",
                    help="Name of the the database file: {file_name}.db")
args = parser.parse_args()

conn = sqlite3.connect(args.db)
# Create a cursor
cursor = conn.cursor()

file_list = os.listdir(par_path + "/src/fasta/alignment/")

def add_taxon():
    # file_list[1:2] is only for 'Abacionidae_NT.fasta' to test
    for fasta_file in file_list[1:2]:
        print(fasta_file)
        with open("fasta/alignment_d/{}".format(fasta_file), "w+") as outputs:
            for record in SeqIO.parse('fasta/alignment/' + fasta_file, "fasta"):
                result = cursor.execute("SELECT taxon.opentol_id, taxon FROM taxon LEFT JOIN barcode ON barcode.taxon_id = "
                                        "taxon.taxon_id WHERE barcode.barcode_id = {}".format(record.id))
                results = result.fetchall()
                record.description = ''
                record.id = '{}_{}_{}'.format(record.id, results[0][0], results[0][1])
                # print(record)
                #print(outputs)
                SeqIO.write(record, outputs, 'fasta')
            outputs.close()

def get_tot_ott():
    dict = {}; fasta = ""; all_ott = []; all_fasta=[]
    for fasta_file in file_list[1:2]:
        with open("fasta/alignment/{}".format(fasta_file), "r") as outputs:
            for line in outputs:
                if line.__contains__(">"):
                    # If contains header
                    # Split and add as a key
                    header = line.split("_")
                    print(header[1])
                    all_ott.append(header[1])   # To check the amount of original ott
    total = 0
    my_list = list(set(all_ott))
    print(my_list)

def add_all():
    dict = {};
    fasta = "";
    all_ott = [];
    all_fasta = []
    for fasta_file in file_list[1:2]:
        with open("fasta/alignment_d/{}".format(fasta_file + "_distinct.fasta"), "w+") as distinct:
            with open("fasta/alignment_d/{}".format(fasta_file), "r") as outputs:
                for line in outputs:
                    # If line contains > add to key
                    # Check if fasta empty
                    # If not empty add as value
                    if line.__contains__(">"):
                        if dict.get(line) is None:  # No key yet
                            dict[line]=fasta
                    else:
                        fasta += line
    print(dict.keys())
    return dict

def get_key(dict, val):
    for key, value in dict.items():
        if val == value:
            return key

    return "key doesn't exist"

if __name__ == '__main__':
  add_taxon()
  get_tot_ott()
  #dict = add_all()
