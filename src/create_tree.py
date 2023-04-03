import argparse
import os
import sqlite3
import subprocess
import opentree
from Bio import SeqIO

# TODO When rewriting for snakemake use file name as input without extensions
#File declarations

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()
parser.add_argument('-db', default="BOLD_COI-5P_barcodes.db",
                    help="Name of the the database file: {file_name}.db")
args = parser.parse_args()

# RAXML in var
RAXML = "raxmlHPC.exe"



def change_ids(fasta_file):
    """
    Change the headers of every sequence in a fasta file from >{barcode_id} to >{opentol_id}_{barcode_id}.
    Saves corrected fasta file in new directory fasta/alignments_c
    :param fasta_file: name of the fasta file
    """
    # Make directory for alignments with corrected headers
    os.makedirs("fasta/alignment_c/", exist_ok=True)
    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect(args.db)
    # Create a cursor
    cursor = conn.cursor()
    with open("fasta/alignment_c/{}".format(fasta_file),"w") as outputs:
        for record in SeqIO.parse('fasta/alignment/' + fasta_file, "fasta"):
            result = cursor.execute("SELECT taxon.opentol_id FROM taxon LEFT JOIN barcode ON barcode.taxon_id = "
                                  "taxon.taxon_id WHERE barcode.barcode_id = {}".format(record.id))
            id = result.fetchall()[0][0]
            record.description = ''
            record.id = '{}_{}'.format(id, record.id)

            SeqIO.write(record, outputs, 'fasta')


if __name__ == '__main__':
    os.makedirs('constraint_trees/', exist_ok=True)
    os.makedirs('raxml/', exist_ok=True)
    # Change when trying a new raxml run
    run_name = 'with_constraint_test'

    # Includes both NT and AA alignments in file list
    file_list = os.listdir(par_path + "/src/fasta/alignment/")

    # file_list[1:2] is only for 'Abacionidae_NT.fasta' to test
    for fasta_file in file_list[1:2]:
        os.makedirs('raxml/{}'.format(run_name), exist_ok=True)
        # Get ott_id in fasta header
        change_ids(fasta_file)

        # Get opentol id per name
        sep = '_'
        taxon = fasta_file
        taxon = taxon.split(sep, 1)[0]
        ott_id = opentree.OT.get_ottid_from_name(taxon)

        # Make constraint tree using family name (derived from file name)
        with open('constraint_trees/{}.nwk'.format(fasta_file.strip('.fasta')), 'w') as outfile:
            outfile.write(str(opentree.OT.synth_subtree(ott_id=ott_id, label_format="id").tree))

        # Run raxml commandline
        subprocess.run("{} -p 100 -m  GTRCAT -n {} -s fasta/alignment_c/{} -g temp_subtree.nwk -w raxml/{}".format(
            RAXML, run_name, fasta_file, run_name))
