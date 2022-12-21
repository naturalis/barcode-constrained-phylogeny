import os
import sys
import subprocess
import time
from io import StringIO

import Bio.Phylo
import pandas as pd
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO, Phylo
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib
import matplotlib.pyplot as plt
import Bio.Phylo
import matplotlib.pyplot as plt
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor


def align():
    os.mkdirs('alignment/', exists_ok=True)
    fasta = open('fasta/family/family.txt', 'r')
    # !BUG, writing alignmnet.fa takes longer, so it is not present yet when running
    # the rest of the code (works now because it exists in directory)
    # Make Multiple Sequence alignment using MUSCLE
    for i in fasta.readlines()[1:-1]:
        print("aligning family %s" %i.strip())
        muscle_cline = MuscleCommandline(input="fasta/family/%s.fasta" % i.strip(), out="alignment/alignment_%s.fasta" % i.strip())
        print(muscle_cline)
        child = subprocess.run(str(muscle_cline), shell=True)

        while not os.path.exists(os.path.join('alignment/', 'alignment_%s.fasta' %i.strip())):
            time.sleep(1)
        print('Done')

    fasta.close()



def tree_notworking():
    # # Read alignment file
    alignment = AlignIO.read("alignment/alignment_Acanthuriformes.fasta", "fasta")

    # Calculate distance with identity method
    calculator = DistanceCalculator('identity')

    # Make distance matrix
    distance_matrix = calculator.get_distance(alignment)

    #Construct phylogenetic tree
    constructor = DistanceTreeConstructor()
    upgmatree = constructor.upgma(distance_matrix)

    # Visualize phylogenetic tree
    matplotlib.rc('font', size=6)
    tree = Bio.Phylo.draw(upgmatree)
    plt.savefig("plot.png")

def tree():
    # Open the file containing FASTA-formatted sequences
    with open('fasta/family/family.txt', 'r') as fasta:
        # Iterate over the lines in the file (skipping the first and last lines)
        # First is .fasta (no familyname) and last family.txt with all the
        # family names
        for i in fasta.readlines()[1:-1]:
            # Strip leading and trailing whitespace from the line
            sequence_name = i.strip()
            print(sequence_name)

            # Read the alignment file for the current sequence
            alignment = AlignIO.read(
                f"alignment/alignment_{sequence_name}.fasta", "fasta")

            # If the alignment contains more than one sequence...
            if len(alignment) > 1:
                # Create a DistanceCalculator object with the identity method
                calculator = DistanceCalculator('identity')

                # Calculate a distance matrix using the alignment
                distance_matrix = calculator.get_distance(alignment)

                # Create a DistanceTreeConstructor object
                constructor = DistanceTreeConstructor()

                # Construct a phylogenetic tree using the UPGMA method
                upgmatree = constructor.upgma(distance_matrix)

                # Write the tree to a file in Newick format
                with open(f"upgma/{sequence_name}_tree.txt", "w") as f:
                    f.write(upgmatree.format('newick'))

                # Set the font size for the tree plot
                matplotlib.rc('font', size=6)

                # Draw the tree and save it to a file
                Bio.Phylo.draw(upgmatree, do_show=False)
                plt.savefig(f"png/{sequence_name}_plot.png")
                print(f"{sequence_name} tree saved.")

                # Close the plot
                plt.close()


# align()
tree()
