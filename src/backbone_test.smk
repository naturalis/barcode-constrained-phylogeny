import glob

configfile: "src/config.yaml"

fasta_files = glob.glob("../data/family_fastas/*.fasta")

ott_database_file = "src/data/databases/outfile.db"

rule representatives_to_file:
    input:
        "data/fasta/family/",
        "data/fasta/family/backbone/representatives.txt"
    #output: "data/fasta/backbone/rep.fasta"
    threads: 4
    script:
        "representatives_to_file.py"

rule replace_newick:
    input:
        "data/fasta/family/backbone/rep.bestTree",
        "data/fasta/family/backbone/rep.fasta"
    output:
        "data/fasta/family/backbone/rep_altered.bestTree"
