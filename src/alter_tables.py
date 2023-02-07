import sqlite3


def make_db(conn, cursor):
    # Create a table
    cursor.execute("""CREATE TABLE IF NOT EXISTS taxon (
        taxon TEXT,
        kingdom TEXT NOT NULL,
        family TEXT NOT NULL
        )
    """)
    # Create a table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS barcode (
        processid TEXT,
        marker_code TEXT,
        nucraw TEXT,
        country TEXT,
        taxon TEXT
    )""")
    # Commit the changes
    conn.commit()

def make_distinct(conn, cursor):
    cursor.execute("""INSERT INTO taxon SELECT DISTINCT * FROM taxon_temp""")
    cursor.execute("""INSERT INTO barcode SELECT DISTINCT * FROM barcode_temp""")
    cursor.execute("""DROP TABLE taxon_temp""")
    cursor.execute("""DROP TABLE barcode_temp""")
    # Commit the changes
    conn.commit()

if __name__ == '__main__':
    # Connect to the database (creates a new file if it doesn't exist)
    conn = sqlite3.connect("bold_COI_barcodes.db")

    # Create a cursor
    cursor = conn.cursor()

    # Make SQLite tables
    make_db(conn, cursor)
    make_distinct(conn, cursor)

    # Close the connection
    conn.close()