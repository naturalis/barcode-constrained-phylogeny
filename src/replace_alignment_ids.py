import sqlite3

from Bio import SeqIO


def replace_with_ott(cursor, alignment_input, alignment_output):
    # Get dict with different ott's
    """Replace the header in the phylip file.
    Uses the phylip file with the barcode id as input and creates a new file with an ott id as header.
    """
    with open(alignment_output, "w") as output:
        for record in SeqIO.parse(alignment_input, "fasta"):
            result = cursor.execute("SELECT taxon.opentol_id, barcode.barcode_id FROM taxon LEFT JOIN barcode ON"
                                        " barcode.taxon_id = taxon.taxon_id WHERE barcode.barcode_id = {}".format(
                                                                                                            record.id))
            results = result.fetchall()
            opentol_id = results[0][0]
            record.id = '{}_{}'.format(opentol_id, record.id)
            record.description = ""
            SeqIO.write(record, output, 'fasta')


if __name__ == '__main__':
    database_file = snakemake.input[0]
    alignment_input = snakemake.input[1]



    alignment_output = snakemake.output[0]

    conn = sqlite3.connect(database_file)
    cursor = conn.cursor()

    replace_with_ott(cursor, alignment_input, alignment_output)