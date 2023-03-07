import os
import subprocess

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


if __name__ == '__main__':
    align_seq()
