import argparse
import os
import sqlite3
from Bio import SeqIO
import argparse
import os
import sqlite3
import subprocess
from io import StringIO
from Bio.Phylo.Applications import RaxmlCommandline
import opentree
from Bio import SeqIO
from Bio import Phylo

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
    """Make dict. Per ott get all fastas.
    """
    dict = {}; fasta = ""; all_ott = []; all_fasta=[]; dict_name_syn={}
    for fasta_file in file_list[1:2]:
        with open("fasta/alignment_d/{}".format(fasta_file), "r") as outputs:
            for line in outputs:
                if line.__contains__(">"):
                    # If contains header
                    # Split and add as a key
                    header = line.split("_")
                    all_ott.append(header[1])   # To check the amount of original ott
                    # Add header to dict if not already in dict
                    if fasta != "":
                        if dict.get(header[1]) is None:  # No key yet
                            dict[header[1]] = []
                            dict[header[1]].append(fasta)
                            print(dict.get(header[1]))
                            fasta = ""
                        else:   # Key exists
                            #print(dict[header[1]])
                            CONTAINS=True
                            for value in dict.values():
                                for v in value:
                                    if v != fasta:
                                        CONTAINS=False
                            if CONTAINS:
                                dict[header[1]].append(fasta)
                                fasta = ""
                else:
                    fasta += line
    total = 0
    my_list = list(set(all_ott))
    #print(my_list)
    for v in dict.values():
        print(len(v))
        for e in v:
            print(len(e))
    #with open("try.fasta", "w") as o:
    #    print(len(dict.values()))
    return dict

def distinct_file(dict):
    for fasta_file in file_list[1:2]:
        with open("fasta/alignment_d/{}".format(fasta_file+ "_distinct.fasta"), "w+") as outputs:
            for k,v in dict.items():
                for val in v:

                    outputs.write(">" + k + "\n")
                    outputs.write(val)


def tree():
    fasta_file = " fasta/alignment_d/Abacionidae_NT.fasta_distinct.fasta"
    RAXML = "raxmlHPC.exe"
    subprocess.run(
        "{} -p 100 -m  GTRGAMMA -n {} -s {}  -w raxml/{}".format(RAXML, "run2", fasta_file,
                                                                                   "run2"))
    #with open("raxml/run2/RAxML_bestTree.run2") as f:
    #    tree = f.read()
    #tree = Phylo.read(StringIO(tree), "newick")
    # print(tree)
    #Phylo.draw(tree)




if __name__ == '__main__':
  add_taxon()
  dict = get_tot_ott()
  distinct_file(dict)
  tree()
  #dict = add_all()
