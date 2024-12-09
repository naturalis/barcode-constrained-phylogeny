import errno
import sqlite3
import os
import pandas as pd
import argparse
import util
from pathlib import Path

levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'all']

"""
This script, `family_fasta.py`, is responsible for generating FASTA files for each family of a specified higher taxon 
from a SQLite database.

The script performs the following steps:
1. Connects to the SQLite database.
2. Retrieves distinct families and BINs (Barcode Index Numbers) for the higher taxon defined in the query restrictions.
3. Iterates over each unique family and creates a directory for each family.
4. For each BIN in a family, it fetches the longest sequence and writes it to the corresponding family's FASTA file.
5. The sequence written to the FASTA file is stripped of non-ACGT characters as HMMER (used in later steps of the 
   workflow) chokes on them.
6. The script continues this process until it has iterated over all families and their respective BINs.
7. Closes the connection to the database.

The script uses command line arguments for the database file to query, directory to write FASTA files to, taxonomic 
level to filter (e.g. order), taxon name to filter (e.g. Primates), number of chunks (families) to write to file, 
marker code (e.g. COI-5P), and log level. This script is invoked by the Snakefile as a shell command with the required
arguments in the rule `family_fasta`.
"""

def get_family_bins(q, conn):
    """
    Gets distinct families and bins for the higher taxon defined in the query restrictions
    :param q: Dictionary with query restrictions
    :param conn: Connection to SQLite database
    :return Pandas data frame
    """

    # Check if filter_level in config.yaml is usable
    if q['level'].lower() in levels:

        # Select all distinct family names that match config.yaml filters
        level = q['level']
        name = q['name']
        marker_code = q['marker_code']
        sql = f'''
            SELECT DISTINCT family, genus, species, bin_uri
            FROM barcode
            WHERE marker_code = '{marker_code}' 
              AND "{level}" = '{name}'
              AND family IS NOT NULL 
              AND family <> ''
            ORDER BY family ASC, bin_uri ASC;
        '''
        fam = pd.read_sql_query(sql, conn)

        # Check if this configuration contains any records at all
        if len(fam) == 0:
            raise Exception(f"No records found with {q}.")

    else:
        raise Exception(f"Filter level {q['level']} from config file does not exists as a column in the database")
    return fam


def write_bin(q, conn, outfile):
    """
    Writes the longest sequence for a BIN to file
    :param q: query object
    :param conn: DB connection
    :param fh: file handle
    :return:
    """

    # Fetch the longest sequence in the BIN
    logger.info(f"Writing longest sequence for BIN {q['bin_uri']} to FASTA")
    query = f'''
        SELECT b.processid, t.bin_uri, t.opentol_id, t.species, b.nuc, b.barcode_id
        FROM barcode b
        JOIN taxon t ON b.taxon_id = t.taxon_id
        WHERE
            t."{q["level"]}" = "{q["name"]}" AND
            t."{q["rank"]}" = "{q["taxon"]}" AND
            t.bin_uri = "{q["bin_uri"]}" AND
            b.marker_code = "{q["marker_code"]}" AND
            t.species IS NOT NULL AND
            t.species <> '' 
        ORDER BY
        b.nuc_basecount DESC LIMIT 1
        '''
    logger.debug(query)
    famseq = pd.read_sql_query(query, conn)

    # Append to file handle fh
    with open(outfile, "w") as fh:
        for _, row in famseq.iterrows():
            defline = f'>{row["barcode_id"]}|ott{row["opentol_id"]}|{row["processid"]}|{row["bin_uri"]}|{row["species"]}\n'
            fh.write(defline)

            # Strip non-ACGT characters (dashes, esp.) because hmmer chokes on them
            seq = row['nuc'].replace('-', '') + '\n'
            fh.write(seq)


if __name__ == '__main__':
    # Define and process command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-d', '--database', required=True, help='Database file to query')
    parser.add_argument('-f', '--fasta_dir', required=True, help='Directory to write FASTA files to')
    parser.add_argument('-l', '--level', required=True, help='Taxonomic level to filter (e.g. order)')
    parser.add_argument('-n', '--name', required=True, help='Taxon name to filter (e.g. Primates)')
    parser.add_argument('-L', '--limit', required=True, type=int, help='Fasta sequence limit, switch to lower rank if above (e.g. 200)')
    parser.add_argument('-m', '--marker', required=True, help='Marker code, e.g. COI-5P')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()
    database_file = args.database

    level = args.level.lower()
    if level not in levels:
        raise Exception(f"Filter level {level} from config file does not exists as a column in the database")

    if levels.index(level) > levels.index("family"):
        raise Exception("Level filter value must be 'family' or higher rank")

    # Configure logger
    logger = util.get_formatted_logger('family_fasta', args.verbosity)

    try:
        os.makedirs(args.fasta_dir, exist_ok=True)
    except OSError as error:
        logger.error(error)
        exit(1)

    # Connect to the database (creates a new file if it doesn't exist)
    logger.info(f"Going to connect to database {args.database}")
    connection = sqlite3.connect(args.database)
    cursor = connection.cursor()

    # Get families and bins for configured level and name
    query = {
        'level': args.level,
        'name': args.name,
        'marker_code': args.marker
    }
    df = get_family_bins(query, connection)

    def write_fasta(query, rank, taxon, family_bin_uris, higher_ranks=None):
        logger.info(f"Writing {taxon} ({rank})")

        # Make directory and open file handle
        align_file = os.path.join(args.fasta_dir, "taxon", taxon, "unaligned.fa")
        try:
            os.makedirs(Path(align_file).parent, exist_ok=True)
        except OSError as error:
            logger.error(error)
            exit(1)

        # Iterate over bins in family
        for bin_uri in family_bin_uris:
            logger.debug(f"Writing {bin_uri}")
            query['bin_uri'] = bin_uri
            query["rank"] = rank
            query["taxon"] = taxon
            write_bin(query, connection, align_file)

    family_set = {}
    genus_set = {}
    # Iterate over distinct families
    with open(os.path.join(args.fasta_dir, "taxon_fasta.tsv"), "w") as fw:
        unique_families = df['family'].unique()
        split_families = []
        for family in unique_families:
            family_bin_uris = df[df['family'] == family]['bin_uri'].unique()
            if len(family_bin_uris) > args.limit:
                split_families.append(family)
                fw.write(f"{family}\tfamily\t{len(family_bin_uris)}\tTrue\t\n")
                continue

            write_fasta(query, "family", family, family_bin_uris)
            fw.write(f"{family}\tfamily\t{len(family_bin_uris)}\tFalse\t\n")

        split_genera = []
        if split_families:
            for family in split_families:
                unique_genera = df[df['family'] == family]['genus'].unique()
                for genus in unique_genera:
                    if not genus:
                        continue
                    genus_bin_uris = df[(df['family'] == family) & (df['genus'] == genus)]['bin_uri'].unique()
                    if len(genus_bin_uris) > args.limit:
                        split_genera.append(genus)
                        fw.write(f"{genus}\tgenus\t{len(genus_bin_uris)}\tTrue\t{family}\n")
                        continue

                    write_fasta(query, "genus", genus, genus_bin_uris)
                    fw.write(f"{genus}\tgenus\t{len(genus_bin_uris)}\tFalse\t{family}\n")

        if split_genera:
            for genus in split_genera:
                if not genus:
                    continue
                unique_species = df[df['genus'] == genus]['species'].unique()
                for species in unique_species:
                    species_bin_uris = df[(df['genus'] == genus) & (df['species'] == species)]['bin_uri'].unique()
                    if len(species_bin_uris) > args.limit:
                        raise NotImplementedError("This is impossible!")
                    write_fasta(query, "species", species, species_bin_uris)
                    fw.write(f"{species}\tspecies\t{len(species_bin_uris)}\tFalse\t{genus}\n")

    # Close the connection
    connection.close()
