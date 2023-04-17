import tarfile


def tar_unzip(tu_inputfile):
    # Open file at input directory
    file = tarfile.open(tu_inputfile)

    # Extract file to output directory
    file.extractall("data2/")

    file.close()

if __name__ == '__main__':
    inputfile = snakemake.input[0]

    tar_unzip(inputfile)