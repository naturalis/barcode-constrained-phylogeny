import os
import subprocess
import re
import shutil
import datetime

# Specify paths
par_path = os.path.abspath(os.path.join(os.pardir))
masce_path = "../data/macse_v2.06.jar"
file_list = os.listdir(par_path + "/src/fasta/family/")

def align_seq():
    os.makedirs('fasta/alignment', exist_ok=True)
    #TODO change from two family test data to everything?
    for fasta_file in file_list[0:2]:
        output = 'fasta/alignment/' + fasta_file.rstrip('.fasta')
        input = str('fasta/family/' + fasta_file)
        # Uses input FASTA file en generates an alignment in AA and NT
        alignseq = "java -jar {} -prog alignSequences -seq {} -out_NT {}_NT.fasta -out_AA {}_AA.fasta".format(masce_path, input, output, output)
        subprocess.run(alignseq, shell=True)
        remove_exclamation_mark(output)     # Call function to remove exclamation mark


def remove_exclamation_mark(output):
    """Remove eclamation mark from file.
    :param output: File to be changed.
    :return: None, but create a file without exclamation marks.
    """
    subprocess.run("sed s/!/-/g {}_NT.fasta > {}_temp_NT.fasta".format(output, output), shell=True)
    print("{}_NT.fasta".format(output))
    os.remove("{}_NT.fasta".format(output))
    os.rename("{}_temp_NT.fasta".format(output), "{}_NT.fasta".format(output))


if __name__ == '__main__':
    align_seq()