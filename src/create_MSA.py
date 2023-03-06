import os
import subprocess

# Specify paths
par_path = os.path.abspath(os.path.join(os.pardir))
masce_path = "../masce/macse_v2.06.jar"
input_fasta = par_path + "\\data\\family.fasta"     # Input file used


def align_seq():
    # Uses input FASTA file en generates an alignment in AA and NT
    alignseq = "java -jar {} -prog alignSequences -seq {} ".format(masce_path, input_fasta)
    subprocess.run(alignseq, shell=True)


if __name__ == '__main__':
    align_seq()
