import glob

configfile: "src/config.yaml"

fasta_files = glob.glob("../data/family_fastas/*.fasta")

ott_database_file = "data/databases/outfile.db"


rule create_matrix:
    # Matrix is based on the output from the tree
    # Create csv files for family
    # Create the csv file for family in directory data/fasta/family/{familyname}/
    input: "data/fasta/family/"
    threads: 4
    script:
        "distance_matrix.py"
        # Loop over family directories
        # Example: barcode-constrained-phylogeny/data/families/familyA/familyA.bestTree
        # mock code
            # Loop over possible directories in barcode-constrained-phylogeny/data/families/

rule create_submatrices:
    # Based on input from create_matrix n
    input:
        "data/fasta/family/",
        ott_database_file
    threads:
        config["cpu_cores"]
    script:
        "create_submatrices.py"

rule write_representatives_to_file:
    input:
        "data/fasta/family/"
    threads:
        config["cpu_cores"]
    script:
        "get_representatives.py"

