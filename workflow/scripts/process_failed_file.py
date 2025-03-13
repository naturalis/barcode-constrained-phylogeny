import argparse
from Bio import SeqIO


def remove_outgroup_records(failed_file, outgroup_file, output_file):
    # Read outgroup sequences
    outgroup_records = set()
    if outgroup_file:
        with open(outgroup_file, "r") as outgroup:
            for record in SeqIO.parse(outgroup, "fasta"):
                outgroup_records.add(record.id)

    # Filter failed file records
    with open(failed_file, "r") as failed, open(output_file, "w") as output:
        for record in SeqIO.parse(failed, "fasta"):
            if record.id not in outgroup_records:
                # Ensure a sequence is written in one line
                record.seq = record.seq.replace("-", "")
                SeqIO.write(record, output, "fasta-2line")  # Use "fasta-2line" to ensure one-line sequences


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process failed file by removing outgroup records and ensuring one-line sequences.")
    parser.add_argument("--failed", required=True, help="Path to the failed file (FASTA format).")
    parser.add_argument("--outgroup", required=False, help="Path to the outgroup file (FASTA format).")
    parser.add_argument("--output", required=True, help="Path to the output file (FASTA format).")
    args = parser.parse_args()

    remove_outgroup_records(args.failed, args.outgroup, args.output)