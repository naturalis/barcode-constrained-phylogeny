import errno
import sqlite3
import os
import pandas as pd
import argparse
import util


def get_family_bins(q, conn):
    """
    Gets distinct families and bins for the higher taxon defined in the query restrictions
    :param q: Dictionary with query restrictions
    :param conn: Connection to SQLite database
    :return Pandas data frame
    """

    # Check if filter_level in config.yaml is usable
    if q['level'].lower() in ['kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'all']:

        # Select all distinct family names that match config.yaml filters
        level = q['level']
        name = q['name']
        marker_code = q['marker_code']
        sql = f'''
            SELECT DISTINCT family, bin_uri
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


def write_bin(q, conn, fh):
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
            t.family = "{q["family"]}" AND
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
    parser.add_argument('-c', '--chunks', required=True, help="Number of chunks (families) to write to file")
    parser.add_argument('-m', '--marker', required=True, help='Marker code, e.g. COI-5P')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()
    database_file = args.database

    # Configure logger
    logger = util.get_formatted_logger('family_fasta', args.verbosity)

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

    # Iterate over distinct families
    index = 1
    for family in df['family'].unique():
        logger.info(f"Writing {family}")

        # Make directory and open file handle
        subdir = os.path.join(args.fasta_dir, f"{index}-of-{args.chunks}")
        try:
            os.mkdir(subdir)
        except OSError as error:
            logger.warning(error)

        with open(os.path.join(subdir, 'unaligned.fa'), 'w') as handle:

            # Iterate over bins in family
            family_bin_uris = df[df['family'] == family]['bin_uri'].unique()
            for bin_uri in family_bin_uris:
                logger.debug(f"Writing {bin_uri}")
                query['bin_uri'] = bin_uri
                query['family'] = family
                write_bin(query, connection, handle)

        index += 1

    # Close the connection
    connection.close()
