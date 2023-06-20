from Bio import SeqIO
from Bio.AlignIO import read as read_alignment
from subprocess import run

def check_orientation(hmmfile, seqfile):
    with open("../data/hmm/updated_fasta.fa", "w") as outputs:
        seqs = []
        for r in SeqIO.parse(seqfile, "fasta"):
            count_1 = hmmalign_singular(r, hmmfile)
            r.seq = r.seq[::-1]
            count_2 = hmmalign_singular(r, hmmfile)
            if count_1 >= count_2:
                r.seq = r.seq[::-1]
            seqs.append(r)
        SeqIO.write(seqs, outputs, 'fasta')

def hmmalign_singular(r, hmmfile):
    SeqIO.write(r, '../data/hmm/temp.fasta', 'fasta')
    run(['hmmalign', '-o', '../data/hmm/temp.aln', hmmfile, '../data/hmm/temp.fasta'])
    alignment = read_alignment('../data/hmm/temp.aln', "stockholm")
    return alignment.column_annotations['posterior_probability'].count('*')

def hmmalign_file(outfile, hmmfile, seqfile):
    run(['hmmalign', '-o', outfile, hmmfile, seqfile])


if __name__ == '__main__':
    check_orientation('../data/hmm/COI-5P.hmm', '../data/fasta/family/Abronicidae.fasta')
    hmmalign_file('../data/hmm/Abbronicidae.aln', '../data/hmm/COI-5P.hmm', '../data/hmm/updated_fasta.fa')