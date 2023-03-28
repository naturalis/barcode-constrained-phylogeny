import argparse
import os
import sqlite3

par_path = os.path.abspath(os.path.join(os.pardir))

# User arguments
parser = argparse.ArgumentParser()

parser.add_argument('-db', default="/data/databases/BOLD_COI-5P_barcodes.db",
                    help="Name of the the database file: {file_name}.db")

args = parser.parse_args()


def make_tables(conn, cursor):
    """Create taxon and barcode tables in the database.
    :param conn: Connection object to the database.
    :param cursor: Cursor object to execute SQL commands.
    """
    # Create taxon table
    cursor.execute("""CREATE TABLE IF NOT EXISTS taxon (
        taxon_id INTEGER PRIMARY KEY,
        taxon TEXT,
        kingdom TEXT NOT NULL,
        family TEXT NOT NULL
        )
    """)
    # Create barcode table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS barcode (
        barcode_id INTEGER PRIMARY KEY,
        processid TEXT,
        marker_code TEXT,
        nucraw TEXT,
        country TEXT,
        taxon_id TEXT,
        FOREIGN KEY (taxon_id) REFERENCES taxon(taxon_id)
    )""")
    # Commit the changes
    conn.commit()


def make_distinct(conn, cursor):
    """Inserts (the distinct) data from temporary tables into the main tables. Opentol_id column is added to taxon
     table and drops taxon_id to the barcode table as foreignkeys. The temporary tables are dropped.
    :param conn: Connection object to the database.
    :param cursor: Cursor object to execute SQL commands.
    """
    # Select only the distinct taxon entries from taxon_temp, insert into taxon
    cursor.execute("""INSERT INTO taxon (taxon, kingdom, family)
     SELECT DISTINCT * FROM taxon_temp""")

    # Get taxon_id from taxon table as foreign key insert
    cursor.execute("""
     INSERT INTO barcode (processid, marker_code, nucraw, country, taxon_id) 
     SELECT DISTINCT barcode_temp.processid, barcode_temp.marker_code,
     barcode_temp.nucraw, barcode_temp.country, taxon.taxon_id
     FROM barcode_temp INNER JOIN taxon ON barcode_temp.taxon = taxon.taxon""")

    cursor.execute("""ALTER TABLE taxon ADD opentol_id varchar(10)""")


    # Drop old tables
    cursor.execute("""DROP TABLE taxon_temp""")
    cursor.execute("""DROP TABLE barcode_temp""")

    # Commit the changes
    conn.commit()


if __name__ == '__main__':
    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect(args.db)

    # Create a cursor
    cursor = conn.cursor()

    # Make SQLite tables
    make_tables(conn, cursor)

    # Keep only distinct barcodes and taxon
    make_distinct(conn, cursor)

    # Close the connection
    conn.close()
