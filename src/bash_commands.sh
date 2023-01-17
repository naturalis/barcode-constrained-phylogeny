# sudo apt install muscle
source environment/bin/activate
cd src
python3.9 retrieve_bold.py
cd fasta/family
ls | sed 's/.fasta//g' > family.txt
cd ..
cd ..
python3.9 msa_alignment.py