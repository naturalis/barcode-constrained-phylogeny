import os
import sqlite3
import subprocess

from Bio import AlignIO
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

#File declarations
par_path = os.path.abspath(os.path.join(os.pardir))
os.makedirs('phylip/', exist_ok=True)
os.makedirs('trees/upgma/', exist_ok=True)
import opentree


def convert_file(fasta_file, phylip_file):
    print('converting files')
    # Convert fasta to a phylip type file
    #TODO Change this in family_fasta.py so that there are no single sequence fasta files.
    count = 0
    with open (fasta_file, "r") as f:
        for record in SeqIO.parse(f, 'fasta'):
            count += 1
            if count > 1:
                input_f = SeqIO.parse(fasta_file, "fasta")
                SeqIO.write(input_f, phylip_file, "phylip")
                return True
    return False


def create_distance_matrix(phylip_file):
    # Create distance matrix from alignment
    alignment = AlignIO.read(open(phylip_file), 'phylip')
    constr = DistanceTreeConstructor()
    calc = DistanceCalculator('identity')
    distance_matrix = calc.get_distance(alignment)
    return constr, distance_matrix


def create_tree(constr, distance_matrix, family):
    # Create UPGMA tree and display it
    # TODO fix what is on the tree (names etc.)
    tree = constr.upgma(distance_matrix)
    # Write the tree to a file in Newick format
    with open(f"trees/upgma/{family}_tree.txt", "w") as f:
        f.write(tree.format('newick'))
    return tree

def create_tree_raxml(file_name):
    print('creating trees')
    print(file_name)
    command = "raxmlHPC -s {} -n {} -f d -g {} -m GTRCAT -p 12345"
    file = 'phylip/alignment_Acanthosomatidae.phylip'
    outfile = '{}.tre'.format(file_name)


    # Set up the OpenTree API client
    tree_constraint = opentree.OT.synth_subtree(ott_id=215958)


    subprocess.run(command.format(file_name + '.fasta', 'tree', tree_constraint), shell=True)

def change_ids(fasta_file):
    os.makedirs('corrected/', exist_ok=True)
    from Bio import SeqIO
    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect('BOLD_COI_barcodes.db')
    print(fasta_file)
    # Create a cursor
    cursor = conn.cursor()
    with open("corrected/{}".format(fasta_file), "w") as outputs:
        for record in SeqIO.parse('alignment/' + fasta_file, "fasta"):
            result = cursor.execute("SELECT taxon.opentol_id FROM taxon LEFT JOIN barcode ON "
                                   "barcode.taxon_id = taxon.taxon_id WHERE barcode.barcode_id = {}".format(record.id))
            id = result.fetchall()[0][0]
            record.id = '{}_{}'.format(record.id, id)
            SeqIO.write(record, outputs, 'fasta')

if __name__ == '__main__':
    file_list = os.listdir(par_path + "/src/alignment/")
    # for fasta_file in file_list:
    #     phylip_file = f"phylip/{fasta_file.rstrip('.fasta')}.phylip"
    #     passed = convert_file(fasta_file, phylip_file)
    #     if passed:
    #         constructor, distance_matrix = create_distance_matrix(phylip_file)
    #         tree = create_tree(constructor, distance_matrix, fasta_file.rstrip('.fasta'))
    fasta_file = 'alignment_Acanthosomatidae.fasta'
    phylip_file = f"phylip/{fasta_file.rstrip('.fasta')}.phylip"
    change_ids(fasta_file)
    # convert_file('corrected/alignment_Acanthosomatidae.fasta', phylip_file)
    create_tree_raxml('corrected/'+ fasta_file.rstrip('.fasta'))
    # change_ids(fasta_file)
