DB=$1
HITS=$2
IN=$3
OUT=$4

# get the distinct process IDs from the FASTA header, third pipe-separated word
grep '>' ${IN} | cut -f3 -d'|' | sort | uniq > ${IN}.ids

# now alias the IDs, which is supposed to be for better run time performance
# (and it suppresses a warning)
blastdb_aliastool -seqid_file_in ${IN}.ids -seqid_file_out ${IN}.ids.aliased

# run blast on the infile, omit the infile IDs in the query, get top 10 hits for each query, write to tabular report
blastn \
  -query ${IN} \
  -db ${DB} \
  -out ${IN}.blast \
  -outfmt 6 \
  -max_target_seqs 10 \
  -negative_seqidlist ${IN}.ids.aliased

# parse the blast report:
# - IDs of hits are in col 2
# - sort, get distinct IDs with occurrence counts, sort ascending on the counts, get most common $HITS
# - remove the occurrence count and write to file, one ID per line
cut -f 2 ${IN}.blast \
  | sort | uniq -c | sort -k 1 | tail -${HITS} \
  | sed 's/^[[:space:]]*[0-9][0-9]*[[:space:]]*//' > ${IN}.hits

# extract the sequences. This is a multi-line FASTA file. To make the constraint tree we need to map the
# process IDs from this file back to ott IDs and get the ott IDs from the ingroup. Also, we'll use this
# file to get the IDs for the outgroup when running raxml.
blastdbcmd -db ${DB} -entry_batch ${IN}.hits -out ${OUT}
