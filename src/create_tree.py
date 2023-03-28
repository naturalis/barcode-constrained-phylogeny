import os
import subprocess
from Bio import AlignIO
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.PAML import codeml

#File declarations
par_path = os.path.abspath(os.path.join(os.pardir))
os.makedirs('phylip/', exist_ok=True)
os.makedirs('trees/upgma/', exist_ok=True)
RAXML = "RAXMLHPC.exe"

# RAXML in var


if __name__ == '__main__':
    subprocess.run("{}  -p 100 -m  GTRCAT -n run1 -s fasta/alignment/Abacionidae_NT.fasta ".format(RAXML))
