import argparse
import os
import tarfile

# Working directory
par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()
# Directory location of tarzip file
parser.add_argument('-indir', default=par_path+"/data/",
                    help="Input folder where tar.gz file is located.")

# Desired output directory
parser.add_argument('-outdir', default=par_path+"/data/",
                    help="Output folder for unzipped files.")

# Tarzip file name
parser.add_argument('-file', default="BOLD_Public.30-DEC-2022.tar.gz",
                    help="Name of the file to be unzipped.")

# Get all arguments in args variable
args = parser.parse_args()


def tar_unzip():
    # Open file at input directory
    file = tarfile.open(args.indir + args.file)

    # Extract file to output directory
    file.extractall(args.outdir)

    file.close()


if __name__ == '__main__':
    tar_unzip()