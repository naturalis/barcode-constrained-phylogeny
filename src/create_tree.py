import os

from Bio import AlignIO
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

#File declarations
par_path = os.path.abspath(os.path.join(os.pardir))
os.makedirs('phylip/', exist_ok=True)
os.makedirs('trees/upgma/', exist_ok=True)


def convert_file(fasta_file, phylip_file):
    # Convert fasta to a phylip type file
    #TODO Change this in family_fasta.py so that there are no single sequence fasta files.
    count = 0
    with open ('fasta/alignment/' + fasta_file, "r") as f:
        for record in SeqIO.parse(f, 'fasta'):
            count += 1
            if count > 1:
                input_f = SeqIO.parse('fasta/alignment/' + fasta_file, "fasta")
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


if __name__ == '__main__':
    file_list = os.listdir(par_path + "/src/fasta/alignment/")
    for fasta_file in file_list:
        phylip_file = f"phylip/{fasta_file.rstrip('.fasta')}.phylip"
        passed = convert_file(fasta_file, phylip_file)
        if passed:
            constructor, distance_matrix = create_distance_matrix(phylip_file)
            tree = create_tree(constructor, distance_matrix, fasta_file.rstrip('.fasta'))
