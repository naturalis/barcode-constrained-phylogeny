import glob

configfile: "src/config.yaml"

fasta_files = glob.glob("../data/family_fastas/*.fasta")

ott_database_file = "src/data/databases/outfile.db"

rule representatives_to_file:
    input:
        "data/fasta/",
        "data/fasta/backbone/representatives.txt"
    #output: "data/fasta/backbone/rep.fasta"
    threads: 4
    script:
        "representatives_to_file.py"

rule replace_newick:
    input:
        "data/fasta/backbone/rep.bestTree",
        "data/fasta/backbone/rep.fasta"
    output:
        "data/fasta/backbone/rep_altered.bestTree"
