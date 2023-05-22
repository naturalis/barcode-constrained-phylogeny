import glob

configfile: "src/config.yaml"

fasta_files = glob.glob("../data/family_fastas/*.fasta")

ott_database_file = "src/data/databases/outfile.db"

# For this test snakemake. The files to create the backbone are tested.
# Data structure used for this testfile:
# Barcode-constrained-phylogeny
    # data
        # families
            # Family data (trees, fasta, alignment)cr

rule create_matrix:
    # Matrix is based on the output from the tree
    # Create csv files for family
    # Create the csv file for family in directory data/fasta/family/matrices/
    input: "data/fasta/"
    threads: 4
    script:
        "distance_matrix.py"
        # Loop over family directories
        # Example: barcode-constrained-phylogeny/data/families/familyA/familyA.bestTree
        # mock code
            # Loop over possible directories in barcode-constrained-phylogeny/data/families/

rule create_submatrices:
    # Based on input from create_matrix
    input:
        "data/fasta/",
        ott_database_file
    threads:
        config["cpu_cores"]
    script:
        "create_submatrices.py"


