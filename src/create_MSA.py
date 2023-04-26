import os
import subprocess

def align_seq(in_file, out_file_nt, out_file_aa):
    os.makedirs('data/fasta/alignment', exist_ok=True)
    #TODO change from two family test data to everything?
    # Uses input FASTA file en generates an alignment in AA and NT
    out_file_nt_temp = out_file_nt.split('_')[0] + "_temp_NT.fasta"
    out_file_aa_temp = out_file_nt.split('_')[0] + "_temp_AA.fasta"
    alignseq = "java -jar {} -prog alignSequences -seq {} -out_NT {} -out_AA {}".format(masce_path, in_file,
                                                                                        out_file_nt_temp,
                                                                                        out_file_aa_temp)
    subprocess.run(alignseq, shell=True)
    remove_exclamation_mark(out_file_nt, out_file_nt_temp) # Call function to remove exclamation mark
    remove_exclamation_mark(out_file_aa, out_file_aa_temp)


def remove_exclamation_mark(out_file, outfile_temp):
    """Remove eclamation mark from file.
    :param output: File to be changed.
    :return: None, but create a file without exclamation marks.
    """
    subprocess.run("sed s/!/-/g {} > {}".format(outfile_temp, out_file), shell=True)
    os.remove("{}".format(outfile_temp))


if __name__ == '__main__':
    masce_path = snakemake.input[0]
    in_file = snakemake.input[1]
    out_file_nt = snakemake.output[0]
    out_file_aa = snakemake.output[1]
    align_seq(in_file, out_file_nt, out_file_aa)
    
