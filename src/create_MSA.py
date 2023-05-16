import os
import subprocess


def align_seq(in_file, out_file):
    os.makedirs('../data/fasta/alignment', exist_ok=True)
    # Uses input FASTA file en generates an alignment in AA and NT
    out_file_temp = out_file.split('_')[0] + "_temp_" + out_file.split('_')[1]
    # check if NT or AA is in filename
    acid = out_file.split('_')[1][:2]
    alignseq = "java -jar {} -prog alignSequences -seq {} -out_{} {}".format(masce_path, in_file, acid,
                                                                                        out_file_temp)
    subprocess.run(alignseq, shell=True)

    # Call function to remove exclamation mark
    remove_exclamation_mark(out_file, out_file_temp)


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
    out_file = snakemake.output[0]
    align_seq(in_file, out_file)