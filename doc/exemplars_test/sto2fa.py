import sys
from Bio import AlignIO


input_file_path = sys.argv[1]
output_file_path = sys.argv[2]
alignment = AlignIO.read(input_file_path, "stockholm")
AlignIO.write(alignment, output_file_path, "fasta")
