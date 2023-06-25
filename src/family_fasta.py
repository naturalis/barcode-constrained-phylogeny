import sqlite3
import os
import pandas as pd
import logging

logging.basicConfig(level=snakemake.params.log_level)  # noqa: F821
logger = logging.getLogger(__name__)
fasta_dir = snakemake.params.fasta_dir  # noqa: F821
maxseq = snakemake.params.maxseq  # noqa: F821
minseq = snakemake.params.minseq  # noqa: F821



def write_genera(family, fasta_dir, conn):
    """
    See write_families. Here the same is done for distinct genera within the provided family.
    :param family:
    :param conn:
    """

    # Iterate over distinct genera within family
    famname = (family,)
    gen = pd.read_sql_query("SELECT DISTINCT(genus) FROM taxon WHERE family = ?", conn, params=famname)
    for genus in set(gen['genus']):

        # Fetch processid, not null opentol_id, distinct nucraw within genus
        names = (family, genus,)
        genseq = pd.read_sql_query("""
            SELECT b.nucraw AS sequence, MIN(b.processid) AS processid, t.opentol_id
            FROM barcode b
            JOIN taxon t ON b.taxon_id = t.taxon_id
            WHERE t.family = ? AND t.genus = ? AND t.opentol_id IS NOT NULL
            GROUP BY sequence, t.opentol_id;""",
                                   conn, params=names)

        # Write data
        logger.info("Write to FASTA genus: %s", genus)
        file_name = f"{fasta_dir}/{family}-{genus}.fasta"
        with open(file_name, 'w') as f:
            for _, row in genseq.iterrows():
                line = f'>{row["processid"]}|{row["opentol_id"]}\n{row["sequence"]}\n'
                f.write(line)


def write_families(conn):
    """
    Takes the barcodes from the SQLite database and divides them into their
    taxonomic family names. For every family a fasta file is made named
    'fasta/family/{family name}.fasta'. All distinct barcodes in that family are
    put in the fasta file with the processid and opentol_id from the SQLite db and their
    nucleotide sequence in FASTA format.
    :param conn: Connection to SQLite database
    """
    # Make directory to put FASTA files in
    os.makedirs(fasta_dir, exist_ok=True)

    # Iterate over distinct families
    fam = pd.read_sql_query("SELECT DISTINCT(family) from taxon", conn)
    for family in set(fam['family']):

        # Fetch processid, not null opentol_id, distinct nucraw within family
        logger.info("Writing to FASTA family: %s", family)
        famname = (family,)
        famseq = pd.read_sql_query("""
            SELECT b.nucraw AS sequence, MIN(b.processid) AS processid, t.opentol_id
            FROM barcode b
            JOIN taxon t ON b.taxon_id = t.taxon_id
            WHERE t.family = ? AND t.opentol_id IS NOT NULL
            GROUP BY sequence, t.opentol_id;""",
                                   conn, params=famname)

        # Only write whole family if smaller than maxseq and higher than minseq
        if maxseq >= len(famseq) > minseq:
            file_name = f"{fasta_dir}/{family}.fasta"
            with open(file_name, 'w') as f:
                for _, row in famseq.iterrows():
                    line = f'>{row["processid"]}|{row["opentol_id"]}\n{row["sequence"]}\n'
                    f.write(line)
        elif len(famseq) <= minseq:
            file_name = f"{fasta_dir}/combined_families.fasta"
            with open(file_name, 'a+') as f:
                for _, row in famseq.iterrows():
                    line = f'>{row["processid"]}|{row["opentol_id"]}\n{row["sequence"]}\n'
                    f.write(line)
        elif len(famseq) > maxseq:
            logger.debug("Family %s has more than %s sequences", family, maxseq)
            write_genera(family, fasta_dir, conn)


if __name__ == '__main__':
    database_file = snakemake.input[0]  # noqa: F821

    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect(database_file)

    # Create a cursor
    cursor = conn.cursor()

    # Write barcodes to FASTA in family groups
    write_families(conn)

    # Close the connection
    conn.close()