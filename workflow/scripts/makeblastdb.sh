DATABASE=$1
FASTADIR=$2
ITERATIONS=$3
TMP=$4

# doing this sequentially to avoid race conditions in BLAST indexing
for i in $(seq 1 $ITERATIONS); do

  # using the output from family_fasta
  INFILE=${FASTADIR}/${i}-of-${ITERATIONS}/unaligned.fa

  # only keep records with ott IDs, reformat the headers to retain the process ID, write to $TMP
  # start new file on first iteration, then append
  if [ "${i}" -eq 1 ]; then
    egrep -A 1 'ott[0-9]+' ${INFILE} | awk -F'|' '/^>/ {print ">"$3; next} {print}' > ${TMP}
  else
    egrep -A 1 'ott[0-9]+' ${INFILE} | awk -F'|' '/^>/ {print ">"$3; next} {print}' >> ${TMP}
  fi

done

# index $TMP with makeblastdb, growing the database
makeblastdb -in ${TMP} -dbtype nucl -out ${DATABASE} -parse_seqids

# remove the temporary file
rm ${TMP}



