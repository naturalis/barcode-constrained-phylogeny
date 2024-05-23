import sqlite3
import util
import argparse
import subprocess


def load_tsv(tsv_file, db_file, temp_table_name, table_name):
    """
    Loads the provided TSV file into the database under the provided table name
    :param tsv_file: tab-separated file
    :param db_file: SQLite database file
    :param temp_table_name: database table name
    :param table_name: database table name
    :return:
    """

    # Construct the SQLite import command
    logger.info(f"Going to import TSV file {tsv_file} into database {db_file} in the {temp_table_name} table")
    command = f".mode tabs\n.import {tsv_file} {temp_table_name}\n"

    # Execute the command using sqlite3 CLI via subprocess
    proc = subprocess.run(['sqlite3', db_file], input=command, text=True, capture_output=True)

    # Check if the operation was successful
    if proc.returncode == 0:
        logger.info("TSV file has been successfully imported into the database")
    else:
        logger.error(f"An error occurred: {proc.stderr}")

    # Get the headers
    with open(tsv_file, 'r') as file:
        first_line = file.readline().strip()
        column_headers = [f'"{item}"' for item in first_line.split('\t')]

    # Copy over
    logger.info("Going to copy into permanent table")
    joined = ', '.join(column_headers)
    statement = f'INSERT INTO {table_name}({joined}) SELECT {joined} FROM {temp_table_name};'
    logger.debug(statement)
    database_cursor.execute(statement)


def create_barcode_table(tsv_file, table_name, add_keys):
    """
    Creates the barcode table from the provided TSV file under the provided table name. This operation takes all the
    provided TSV column header names and turns these into database columns of type TEXT. Optionally, if the table is not
    temporary, a barcode_id column is created as an autoincrementing integer primary key, and a taxon_id column that is
    a foreign key to the taxon table.
    :param tsv_file: tab-separated file
    :param table_name: table name
    :param add_keys: add primary and foreign key if true
    :return:
    """

    # Read the first line of the TSV file to get the column headers
    logger.info(f"Going to create {table_name} table from {tsv_file}")
    with open(tsv_file, 'r') as file:
        first_line = file.readline().strip()
        column_headers = first_line.split('\t')  # Assuming a tab-separated values file

    # Construct the SQL statement
    # Start with 'barcode_id INTEGER PRIMARY KEY AUTOINCREMENT'
    create_table_statement = f'CREATE TABLE IF NOT EXISTS {table_name} ('

    # Add each column header as a TEXT type column
    column_headers[column_headers.index('COLLECTORS')] = 'COLLECTORS2' # double column names in the curated BOLD
    column_headers[column_headers.index('COLLECTION_DATE')] = 'COLLECTION_DATE2'
    column_headers[column_headers.index('COUNTRY')] = 'COUNTRY2'
    column_headers[column_headers.index('SITE')] = 'SITE2'
    column_headers[column_headers.index('COORD')] = 'COORD2'
    for i, header in enumerate(column_headers):
        if not add_keys and i == len(column_headers) - 1:
            create_table_statement += f'"{header}" TEXT);'
        else:
            create_table_statement += f'"{header}" TEXT, '

    # If the table is the definitive one, not the temporary one, add a primary and foreign key
    if add_keys:
        logger.info("Adding primary and foreign key")
        create_table_statement += '''
            barcode_id INTEGER PRIMARY KEY AUTOINCREMENT,
            taxon_id INTEGER, FOREIGN KEY(taxon_id) REFERENCES taxon(taxon_id));
        '''

    logger.debug(create_table_statement)
    database_cursor.execute(create_table_statement)


def create_taxon_table(table_name):
    """
    Creates the table that holds normalized taxonomy information under the provided table name. This needs to be
    executed before creating the barcode table because the latter has a foreign key constraint to it.
    :param table_name: table name to create.
    :return:
    """

    # Execute SQL statement to create the taxon table
    logger.info(f"Going to create {table_name} table")
    database_cursor.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name}(
            taxon_id INTEGER PRIMARY KEY AUTOINCREMENT,
            kingdom TEXT NOT NULL,
            phylum TEXT NOT NULL,
            class TEXT NOT NULL,
            "order" TEXT NOT NULL,
            family TEXT NOT NULL,
            subfamily TEXT NOT NULL,
            genus TEXT NOT NULL,
            species TEXT NOT NULL,            
            bin_uri TEXT NOT NULL,
            opentol_id INTEGER,
            UNIQUE(kingdom, phylum, class, "order", family, subfamily, genus, species, bin_uri));
        """)


def index_taxonomy(table):
    """
    Applies indexes on the taxonomy columns in the barcode table
    :return:
    """

    # Iterate over taxonomy columns and compute index
    for column in ['kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'species', 'bin_uri']:
        logger.info(f'Going to index {column} column of {table} table')
        database_cursor.execute(f'CREATE INDEX IF NOT EXISTS {table}_{column}_idx ON {table} ("{column}");')


def normalize_taxonomy():
    """
    Copies the distinct (normalized) taxonomy tuples from barcode table into taxonomy table
    :return:
    """

    logger.info("Going to copy normalized taxonomy tuples from barcode table into taxonomy table")
    database_cursor.execute("""
        INSERT INTO taxon(kingdom, phylum, class, "order", family, subfamily, genus, species, bin_uri)
        SELECT DISTINCT kingdom, phylum, class, "order", family, subfamily, genus, species, bin_uri FROM barcode;
    """)


def update_fk():
    """
    Updates the taxon_id foreign key column in the barcode table by joining all denormalized records in it
    with the normalized tuple in the taxon table. Applies indexes to the PK and FK.
    :return:
    """

    logger.info("Going to update foreign key taxon_id in barcode table")
    database_cursor.execute("""
        UPDATE barcode
        SET taxon_id = (
            SELECT taxon_id
            FROM taxon
            WHERE
                barcode.kingdom = taxon.kingdom AND
                barcode.phylum = taxon.phylum AND
                barcode.class = taxon.class AND
                barcode."order" = taxon."order" AND
                barcode.family = taxon.family AND
                barcode.subfamily = taxon.subfamily AND
                barcode.genus = taxon.genus AND
                barcode.species = taxon.species AND
                barcode.bin_uri = taxon.bin_uri
        );        
    """)
    logger.info('Going to index PK/FK taxon_id on taxon and barcode tables')
    database_cursor.execute('CREATE INDEX "taxon_taxon_id_idx" ON "taxon"("taxon_id")')
    database_cursor.execute('CREATE INDEX "barcode_taxon_id_idx" ON "barcode"("taxon_id")')


if __name__ == '__main__':

    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-i', '--intsv', required=True, help='Input BOLD tsv file')
    parser.add_argument('-o', '--outdb', required=True, help="Output SQLite database file")
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Instantiate logger
    logger = util.get_formatted_logger('create_database', args.verbosity)

    # Connect to the database
    logger.info('Going to connect to database')
    connection = sqlite3.connect(args.outdb)
    database_cursor = connection.cursor()

    # Create database tables
    create_taxon_table('taxon')
    create_barcode_table(args.intsv, 'barcode_temp', False)
    create_barcode_table(args.intsv, 'barcode', True)

    # Load TSV into temporary barcode table, unindexed, then copy over
    load_tsv(args.intsv, args.outdb, 'barcode_temp', 'barcode')
    database_cursor.execute('CREATE INDEX IF NOT EXISTS barcode_nuc_idx ON barcode(nuc)')
    database_cursor.execute('CREATE INDEX IF NOT EXISTS barcode_processid_idx ON barcode(processid)')
    database_cursor.execute('CREATE INDEX IF NOT EXISTS barcode_marker_code_idx ON barcode(marker_code)')

    database_cursor.execute('CREATE INDEX IF NOT EXISTS barcode_processid_idx ON barcode(processid)')
    database_cursor.execute('CREATE INDEX IF NOT EXISTS barcode_marker_code_idx ON barcode(marker_code)')

    # Index the barcode table's taxonomy, copy its distinct tuples into the taxon table, then index the latter
    index_taxonomy('barcode')
    normalize_taxonomy()
    index_taxonomy('taxon')

    # Update the barcode table's foreign key
    update_fk()

    # Commit everything
    logger.info('Cleaning up temporary table and header records')
    database_cursor.execute('DROP TABLE barcode_temp;')
    database_cursor.execute("DELETE FROM barcode where kingdom='kingdom';")
    database_cursor.execute("DELETE FROM taxon where kingdom='kingdom';")
    connection.commit()


