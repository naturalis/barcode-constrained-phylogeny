import os
import subprocess
import argparse

def process_taxon_files(fasta_dir, tmp_file):
    i = 0
    with open(tmp_file, 'w') as tmp:
        for taxon in os.listdir(fasta_dir):
            taxon_path = os.path.join(fasta_dir, taxon)
            if os.path.isdir(taxon_path):
                i += 1
                infile = os.path.join(taxon_path, 'unaligned.fa')

                with open(infile, 'r') as f:
                    lines = f.readlines()

                filtered_lines = []
                for j in range(0, len(lines), 2):
                    if 'ott' in lines[j]:
                        header = lines[j].split('|')[2]
                        sequence = lines[j + 1]
                        filtered_lines.append(f'>{header.strip()}\n{sequence}')

                if i == 1:
                    tmp.writelines(filtered_lines)
                else:
                    with open(tmp_file, 'a') as tmp_append:
                        tmp_append.writelines(filtered_lines)

def make_blast_db(tmp_file, database):
    subprocess.run(['makeblastdb', '-in', tmp_file, '-dbtype', 'nucl', '-out', database, '-parse_seqids'])

def main():
    parser = argparse.ArgumentParser(description='Create BLAST database from FASTA files.')
    parser.add_argument('-f', '--fasta_dir', required=True, help='Directory containing FASTA files')
    parser.add_argument('-t', '--tmp_file', required=True, help='Temporary file to store filtered sequences')
    parser.add_argument('-d', '--database', required=True, help='Output BLAST database file')
    args = parser.parse_args()

    process_taxon_files(args.fasta_dir, args.tmp_file)
    make_blast_db(args.tmp_file, args.database)
    os.remove(args.tmp_file)

if __name__ == '__main__':
    main()