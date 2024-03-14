import sqlite3
import util
import argparse
import subprocess


def load_tsv(tsv_file, db_file, table_name):
    """
    Loads the provided TSV file into the database under the provided table name
    :param tsv_file: tab-separated file
    :param db_file: SQLite database file
    :param table_name: database table name
    :return:
    """

    # Construct the SQLite import command
    logger.info(f"Going to import TSV file {tsv_file} into database {db_file} in the {table_name} table")
    command = f".mode tabs\n.import {tsv_file} {table_name}\n"

    # Execute the command using sqlite3 CLI via subprocess
    proc = subprocess.run(['sqlite3', db_file], input=command, text=True, capture_output=True)

    # Check if the operation was successful
    if proc.returncode == 0:
        logger.info("TSV file has been successfully imported into the SQLite database.")
    else:
        logger.error(f"An error occurred: {proc.stderr}")


def create_barcode_table(tsv_file, table_name):
    """
    Creates the barcode table from the provided TSV file under the provided table name. This operation takes all the
    provided TSV column header names and turns these into database columns of type TEXT. In addition, a barcode_id
    column is created as an autoincrementing integer primary key, and a taxon_id column that is a foreign key to the
    taxon table.
    :param tsv_file: tab-separated file
    :param table_name: table name
    :return:
    """

    # Read the first line of the TSV file to get the column headers
    logger.info(f"Going to create {table_name} table from {tsv_file} + extra columns")
    with open(tsv_file, 'r') as file:
        first_line = file.readline().strip()
        column_headers = first_line.split('\t')  # Assuming a tab-separated values file

    # Construct the SQL statement
    # Start with 'barcode_id INTEGER PRIMARY KEY AUTOINCREMENT'
    create_table_statement = f'CREATE TABLE IF NOT EXISTS {table_name} (barcode_id INTEGER PRIMARY KEY AUTOINCREMENT, '

    # Add each column header as a TEXT type column
    for header in column_headers:
        create_table_statement += f'{header} TEXT, '

    # Postfix with 'taxon_id INTEGER'
    create_table_statement += 'taxon_id INTEGER, FOREIGN KEY(taxon_id) REFERENCES taxon(taxon_id));'

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
        CREATE TABLE {table_name}(
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
            UNIQUE(kingdom, phylum, class, "order", family, subfamily, genus, species))
        """)


def index_taxonomy(table):
    """
    Applies indexes on the taxonomy columns in the barcode table
    :return:
    """

    # Iterate over taxonomy columns and compute index
    for column in ['kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'species']:
        logger.info(f'Going to index {column} column of {table} table')
        database_cursor.execute(f'CREATE INDEX {column}_idx ON {table} ("{column}")')


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
    with the normalized tuple in the taxon table
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
    connection = sqlite3.connect(args.database)
    database_cursor = connection.cursor()

    # Create database tables
    create_taxon_table('taxon')
    create_barcode_table(args.intsv, 'barcode')

    # Load TSV into barcode table, unindexed
    load_tsv(args.intsv, args.outdb, 'barcode')

    # Index the barcode table's taxonomy, copy its distinct tuples into the taxon table, then index the latter
    index_taxonomy('barcode')
    normalize_taxonomy()
    index_taxonomy('taxon')

    # Update the barcode table's foreign key
    update_fk()


