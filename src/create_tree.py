from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import SeqIO

#File declarations
input_fasta = "../data/family_NT.fasta"
input_phylip = "../data/family_NT.phylip"


def convert_file():
    # Convert fasta to a phylip type file
    input_f = SeqIO.parse(input_fasta, "fasta")
    SeqIO.write(input_f, input_phylip, "phylip")


def create_distance_matrix():
    # Create distance matrix from alignment
    alignment = AlignIO.read(open(input_phylip), 'phylip')
    constr = DistanceTreeConstructor()
    calc = DistanceCalculator('identity')
    distance_matrix = calc.get_distance(alignment)
    return constr, distance_matrix


def create_tree(constr, distance_matrix):
    # Create UPGMA tree and display it
    # TODO fix what is on the tree (names etc.)
    tree = constr.upgma(distance_matrix)
    print(tree)
    return tree


if __name__ == '__main__':
    convert_file()
    constructor, distance_matrix = create_distance_matrix()
    tree = create_tree(constructor, distance_matrix)
    Phylo.draw(tree)
