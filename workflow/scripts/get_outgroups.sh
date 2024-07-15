DB=$1
HITS=$2
IN=$3
OUT=$4

# This shell script, `get_outgroups.sh`, is responsible for identifying potential outgroup sequences for a given set of
# sequences. The script performs the following steps:
# 1. Extracts the distinct process IDs from the FASTA header of the input file.
# 2. Aliases the IDs for better runtime performance and to suppress warnings.
# 3. Runs a BLAST search on the input file against a specified database, excluding the input IDs in the query, and
#    retrieves the top 10 hits for each query.
# 4. Parses the BLAST report to get the distinct IDs of the hits along with their occurrence counts, sorts them in
#    ascending order based on the counts, and retrieves the most common hits specified by the user.
# 5. Removes the occurrence count and writes the hit IDs to a file, one ID per line.
# 6. Extracts the sequences of the hits from the BLAST database and writes them to an output file.
#
# The script uses command line arguments for the database file to query, number of hits to retrieve, input FASTA file,
# and output file for the hit sequences. It is invoked by the Snakefile as a shell command in the rule `get_outgroups`.

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
