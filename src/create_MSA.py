import os
import subprocess

from Bio.Align.Applications import MuscleCommandline

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
        alignseq = "java -jar {} -prog alignSequences -seq {} -out_NT 2 -out {}_NT.fasta -out_AA 2 -out {}_AA.fasta -gap_char -".format(masce_path, input, output, output)
        subprocess.run(alignseq, shell=True)

def align():
    os.makedirs('alignment', exist_ok=True)
    for fasta_file in file_list[0:2]:


        print("aligning family %s" %fasta_file.rstrip('.fasta'))
        muscle_cline=MuscleCommandline(input="fasta/family/%s"%fasta_file,out="alignment/alignment_%s"%fasta_file,fasta=True)
        print(muscle_cline)
        #muscle_cline()
        child = subprocess.run(str(muscle_cline),shell=True)

        #whilenotos.path.exists(os.path.join('alignment/','alignment_%s.fasta'%i.strip())):
        #time.sleep(1)
        #print('Done')

if __name__ == '__main__':
    # align_seq()

    align()
    
