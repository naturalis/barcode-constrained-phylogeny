import glob

configfile: "src/config.yaml"

fasta_files = glob.glob("../data/family_fastas/*.fasta")

ott_database_file = "data/databases/outfile.db"


rule create_matrix:
    # Matrix is based on the output from the tree
    # Create csv files for family
    # Create the csv file for family in directory data/fasta/family/{familyname}/
    input: "data/raxml/{family}/{family}.fasta.raxml.bestTree"
    output: "data/fasta/distance_matrix/{family}.txt"
    threads: 4
    script:
        "distance_matrix.py"


rule create_submatrices:
    # Based on input from create_matrix n
    input:
        "data/fasta/distance_matrix/{family}.txt",
        ott_database_file
    output:
        "data/fasta/submatrix/{family}.txt"
    threads:
        config["cpu_cores"]
    script:
        "create_submatrices.py"

rule get_representatives:
    input:
        "data/fasta/distance_matrix/{family}.txt",
    output:
        "data/fasta/representatives/rep.txt"
    threads:
        config["cpu_cores"]
    script:
        "get_representatives.py"

rule write_representatives:
    input:
        "data/fasta/family/{family}.fasta",
        "data/fasta/representatives/rep.txt"
    output:
        "data/fasta/representatives/representatives.fasta"
    threads:
        config["cpu_cores"]
    script:
        "representatives_to_file.py"

rule insert_backbone:
    input:
        "data/raxml/Aethionem/{family}.fasta.raxml.bestTree",
        "data/raxml/representatives/representatives.fasta.raxml.bestTree"
    output:
        "data/representatives/altered_backbone.nwk"
    threads:
        config["cpu_cores"]
    script:
        "insert_backbone.py"


#rule test:
#    input: "data/fasta/submatrix/Aethionem.txt"
#    shell: "touch {input}"

