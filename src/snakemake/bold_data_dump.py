import csv
import sqlite3
import pandas as pd



def extract_bold(conn, bold_tsv, marker):
    """
     It reads a TSV file with a snapshot of the BOLD database in chunks, selects rows that match the user's
    arguments for the desired marker and kingdom, and writes them to two tables in the
    SQLite database (taxon_temp and barcode_temp) with a subset of their columns.
    :param conn: Connection object to the database.
    :param args: Command line arguments
    :param desired_marker: String representing the desired marker
    """
    for chunk in pd.read_csv(bold_tsv, quoting=csv.QUOTE_NONE,
                             low_memory=False, sep="\t", chunksize=10000):
        # Keep rows that match user arguments
        if marker == "COI-5P":
            df = chunk.loc[
                (chunk['marker_code'] == marker) & (chunk["kingdom"] == "Animalia")]

        else:
            # variable marker matk_rcbl is split into two seperate markers
            marker_1 = marker.split('_')[0]
            marker_2 = marker.split('_')[1]
            df = chunk.loc[
                ((chunk['marker_code'] == marker_1) & (chunk["kingdom"] == "Plantae")) |
                (chunk['marker_code'] == marker_2) & (chunk["kingdom"] == "Plantae")]

        # Keep stated columns, do not keep rows where NAs are present
        df_temp = df[['taxon', 'kingdom', 'family']].dropna()

        # Add rows to SQLite table (makes table if not exitst yet)
        df_temp.to_sql('taxon_temp', conn, if_exists='append',
                       index=False)

        # Keep stated columns
        df_temp = df[['processid', 'marker_code', 'nucraw', 'country', 'taxon']]

        # Add rows to SQLite table (makes table if not exitst yet)
        df_temp.to_sql('barcode_temp', conn, if_exists='append', index=False)

        conn.commit()

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
    database = snakemake.output[0]
    bold_tsv = snakemake.input[0]
    marker = snakemake.params.marker

    # Make connection to the database
    conn = sqlite3.connect(database)

    # Create a cursor
    cursor = conn.cursor()

    # Dump BOLD data into DB in temporary tables
    print("BOLD dump into database")
    extract_bold(conn, bold_tsv, marker)

    # Make new tables with different names
    print("Make new tables")
    make_tables(conn, cursor)

    # Drop duplicates
    print("Make distinct")
    make_distinct(conn, cursor)

    # Close the connection
    conn.close()
