import sqlite3
import csv
from pathlib import Path
from sqlalchemy import Column, Integer, String, Float, Date, DateTime, Text, create_engine, ForeignKey, TypeDecorator
from sqlalchemy.orm import sessionmaker, scoped_session, joinedload, class_mapper
from sqlalchemy.ext.declarative import declarative_base
from contextlib import contextmanager
from sqlalchemy import and_
from datetime import datetime
import json

from sqlalchemy.orm import relationship

Base = declarative_base()


def _esc(col):
    return f'"{col}"'


class CustomDate(TypeDecorator):
    impl = String

    def process_bind_param(self, value, dialect):
        if value is not None:
            if isinstance(value, str):
                # If it's already a string, return it as is
                return value
            return value.strftime('%Y-%m-%d')  # Store dates in ISO format
        return value

    def process_result_value(self, value, dialect):
        if value is not None:
            try:
                # Try parsing as '%d-%b-%Y' first
                return datetime.strptime(value, '%d-%b-%Y').date()
            except ValueError:
                try:
                    # If that fails, try parsing as 'YYYY-MM-DD'
                    return datetime.strptime(value, '%Y-%m-%d').date()
                except ValueError:
                    # If both fail, return the original string
                    return value
        return value


class Database:
    """
    This class manages the creation of the SQLite database, creation of object-relational classes from the JSON
    schema provided by BOLD's BCDM packages, creation of database tables based on these ORM classes, and the
    ingestion of data. Usage:

    # First time use:
    db = Database('db.sqlite', '/path/to/BOLD_Public.18-Dec-2023')
    db.make_schema() # instantiates the database tables
    db.load_bcdm() # loads TSV data into the database

    # To merge the megatree:
    db = Database('db.sqlite', '/path/to/BOLD_Public.18-Dec-2023')
    db.merge_megatree('node.sqlite')

    # All subsequent times:
    db = Database('db.sqlite', '/path/to/BOLD_Public.18-Dec-2023')
    node = db.Node('ott34728')
    barcode = db.Barcode(...)
    taxon = db.Taxon(...)
    """

    def __init__(self, db_file: Path, data_package: Path):
        self.db_file = db_file
        self.connection = sqlite3.connect(self.db_file)
        self.cursor = self.connection.cursor()
        self.data_package = data_package
        self.initialized = False
        self.engine = create_engine(f'sqlite:///{self.db_file}')
        self.current_session = None
        self.Session = scoped_session(sessionmaker(bind=self.engine))
        self.Barcode = None
        self.Taxon = None
        self.Node = Node
        self.make_orm()

    def start_session(self):
        """Start a new session."""
        if self.current_session is not None:
            raise RuntimeError("A session is already in progress")
        self.current_session = self.Session()
        return self.current_session

    def end_session(self):
        """End the current session."""
        if self.current_session is None:
            raise RuntimeError("No session is currently in progress")
        self.current_session.close()
        self.current_session = None

    @contextmanager
    def session_scope(self):
        """Provide a transactional scope around a series of operations."""
        session = self.Session()
        try:
            yield session
            session.commit()
        except:
            session.rollback()
            raise
        finally:
            session.close()

    def get_thing(self, model, search_dict=None):
        """
        Retrieves instance(s) of the specified model from the database based on the search criteria.

        :param model: The SQLAlchemy model class to query
        :param search_dict: A dictionary of search criteria
        :return: A list of model instances that match the criteria, or an empty list if none found
        """
        if self.current_session is None:
            raise RuntimeError("No session is currently in progress")

        query = self.current_session.query(model)

        # Eager load all relationships using class-bound attributes
        for relationship in class_mapper(model).relationships:
            query = query.options(joinedload(getattr(model, relationship.key)))

        if search_dict:
            filters = [getattr(model, k) == v for k, v in search_dict.items()]
            query = query.filter(and_(*filters))

        return query.all()

    def get_barcode(self, search_dict=None):
        """
        Retrieves Barcode instance(s) from the database based on the search criteria.

        :param search_dict: A dictionary of search criteria
        :return: A list of Barcode instances that match the criteria, or an empty list if none found
        """
        return self.get_thing(self.Barcode, search_dict)

    def get_taxon(self, search_dict=None):
        """
        Retrieves Taxon instance(s) from the database based on the search criteria.

        :param search_dict: A dictionary of search criteria
        :return: A list of Taxon instances that match the criteria, or an empty list if none found
        """
        return self.get_thing(self.Taxon, search_dict)

    def get_node(self, search_dict=None):
        """
        Retrieves Node instance(s) from the database based on the search criteria.

        :param search_dict: A dictionary of search criteria
        :return: A list of Node instances that match the criteria, or an empty list if none found
        """
        return self.get_thing(self.Node, search_dict)

    def make_orm(self) -> None:
        """
        Generates an ORM schema from a BCDM data package's JSON file
        :return: None
        """
        # Construct the expected path for the datapackage.json file
        folder_name = self.data_package.name
        json_file = self.data_package / f"{folder_name}.datapackage.json"

        if not json_file.exists():
            raise FileNotFoundError(f"Could not find {json_file}")

        with open(json_file, 'r') as f:
            json_schema = json.load(f)
            self.Taxon = self.create_dynamic_taxon_class(json_schema)
            self.Barcode = self.create_dynamic_barcode_class(json_schema)

    def make_schema(self) -> None:
        """
        Applies the ORM schema to the database
        :return: None
        """
        if self.Barcode is None or self.Taxon is None:
            raise ValueError("ORM classes have not been created. Call make_orm() first.")

        # Create tables for DynamicBarcode, DynamicTaxon, and Node
        Base.metadata.create_all(self.engine, tables=[
            self.Barcode.__table__,
            self.Taxon.__table__,
            self.Node.__table__
        ])

        self.initialized = True
        print("Schema created successfully.")

    def load_bcdm(self, bcdm_file: Path) -> None:
        """
        Populates the database with data from a BCDM data package's TSV file
        :param bcdm_file: a TSV file from a BCDM data package
        :return: None
        """
        self._create_temp_tables()
        self._populate_temp_tables(bcdm_file)
        self._finalize_permanent_tables()

    def _create_temp_tables(self):
        """
        Creates temporary tables by cloning the structure of permanent tables,
        excluding primary and foreign keys.
        """
        # For barcode table
        self.cursor.execute("PRAGMA table_info(barcode)")
        barcode_columns = [
            f'{_esc(row[1])} {row[2]}'
            for row in self.cursor.fetchall()
            if row[1] not in ('barcode_id', 'taxon_id')  # Exclude primary and foreign keys
        ]
        self.cursor.execute(f"""
            CREATE TEMPORARY TABLE temp_barcode (
                {', '.join(barcode_columns)}
            )
        """)

        # For taxon table
        self.cursor.execute("PRAGMA table_info(taxon)")
        taxon_columns = [
            f'{_esc(row[1])} {row[2]}'
            for row in self.cursor.fetchall()
            if row[1] != 'taxon_id'  # Exclude primary key
        ]
        self.cursor.execute(f"""
            CREATE TEMPORARY TABLE temp_taxon (
                {', '.join(taxon_columns)}
            )
        """)

    def _populate_temp_tables(self, bcdm_file):
        """
        Loads TSV data from BCDM file into temporary tables.
        :param bcdm_file: Path to the TSV file
        """
        # Get column names from the permanent tables
        self.cursor.execute("PRAGMA table_info(barcode)")
        barcode_columns = [row[1] for row in self.cursor.fetchall() if row[1] not in ('barcode_id', 'taxon_id')]

        self.cursor.execute("PRAGMA table_info(taxon)")
        taxon_columns = [row[1] for row in self.cursor.fetchall() if row[1] not in ('taxon_id', 'opentol_id')]

        # Read TSV and insert into temporary tables
        with open(bcdm_file, 'r') as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            batch_size = 100000
            barcode_batch = []
            taxon_batch = []

            for row in reader:
                barcode_row = [row[col] for col in barcode_columns]
                taxon_row = [row[col] for col in taxon_columns]

                barcode_batch.append(barcode_row)
                taxon_batch.append(taxon_row)

                if len(barcode_batch) >= batch_size:
                    self._insert_batch(barcode_batch, taxon_batch, barcode_columns, taxon_columns)
                    barcode_batch = []
                    taxon_batch = []

            # Insert any remaining rows
            if barcode_batch:
                self._insert_batch(barcode_batch, taxon_batch, barcode_columns, taxon_columns)

    def _insert_batch(self, barcode_batch, taxon_batch, barcode_columns, taxon_columns):
        """
        Insert a batch of records into temporary tables
        """
        # Escape column names with double quotes
        escaped_barcode_columns = [_esc(col) for col in barcode_columns]
        escaped_taxon_columns = [_esc(col) for col in taxon_columns]

        # Convert date strings to ISO format if needed
        for batch in [barcode_batch, taxon_batch]:
            for row in batch:
                for i, value in enumerate(row):
                    if isinstance(value, str):
                        try:
                            # Try parsing as '%d-%b-%Y'
                            date_obj = datetime.strptime(value, '%d-%b-%Y')
                            row[i] = date_obj.strftime('%Y-%m-%d')
                        except ValueError:
                            # If parsing fails, it might already be in 'YYYY-MM-DD' format or not a date
                            pass

        self.cursor.executemany(
            f"""INSERT INTO temp_barcode ({','.join(escaped_barcode_columns)}) 
               VALUES ({','.join(['?' for _ in barcode_columns])})""",
            barcode_batch
        )
        self.cursor.executemany(
            f"""INSERT INTO temp_taxon ({','.join(escaped_taxon_columns)}) 
               VALUES ({','.join(['?' for _ in taxon_columns])})""",
            taxon_batch
        )

    def _finalize_permanent_tables(self):
        """
        Transfers data from temporary tables to permanent tables,
        ensuring correct foreign key relationships and data normalization
        """

        # Get column names for taxon table. Exclude the opentol_id and the taxon_id.
        self.cursor.execute("PRAGMA table_info(taxon)")
        taxon_columns = [row[1] for row in self.cursor.fetchall() if row[1] not in ('taxon_id', 'opentol_id')]

        # Get column names for the barcode table. Exclude the barcode_id and taxon_id.
        self.cursor.execute("PRAGMA table_info(barcode)")
        barcode_columns = [row[1] for row in self.cursor.fetchall() if row[1] not in ('barcode_id', 'taxon_id')]

        # Insert unique taxa into the taxon table
        taxon_cols_str = ', '.join(map(_esc, taxon_columns))
        taxon_insert_query = f"""
            INSERT INTO taxon ({taxon_cols_str})
            SELECT DISTINCT {taxon_cols_str}
            FROM temp_taxon
        """
        self.cursor.execute(taxon_insert_query)

        # Prepare the JOIN conditions for the barcode insert query
        join_conditions = ' AND '.join(f'tx.{_esc(col)} = t.{_esc(col)}' for col in taxon_columns)

        # Prepare the INSERT statement for the barcode table
        barcode_cols_str = ', '.join(f't.{_esc(col)}' for col in barcode_columns)
        barcode_insert_query = f"""
            INSERT INTO barcode ({','.join(map(_esc, barcode_columns))}, taxon_id)
            SELECT {barcode_cols_str}, tx.taxon_id
            FROM temp_barcode t
            JOIN taxon tx ON {join_conditions}
        """
        self.cursor.execute(barcode_insert_query)

        # Commit changes and clean up
        self.connection.commit()
        self.cursor.execute("DROP TABLE temp_barcode")
        self.cursor.execute("DROP TABLE temp_taxon")

    def merge_megatree(self, megatree_file: Path) -> None:
        """
        Merges the megatree topology into the database
        :param megatree_file: A SQLite database file as produced by the megatree-loader process
        :return: None
        """

    @classmethod
    def create_dynamic_barcode_class(cls, json_schema):
        """
        Dynamically creates an ORM class for the barcode table by introspection on a provided JSON schema.
        :param json_schema: A parsed JSON schema object
        :return: An ORM class
        """

        # Initial attributes for a SQLAlchemy ORM class, i.e. the table name to connect to and a primary key
        class_attrs = {
            '__tablename__': 'barcode',
            'barcode_id': Column(Integer, primary_key=True),
            'taxon_id': Column(Integer, ForeignKey('taxon.taxon_id'), nullable=False, index=True),
            'taxon': relationship('Taxon')
        }

        # Navigate to the substructure containing the field definitions
        resources = json_schema['resources']
        schema_fields = next(
            resource['schema']['fields'] for resource in resources if resource['profile'] == 'tabular-data-resource')

        # Map between types in the JSON schema and Python boxed types
        type_mapping = {
            'string': String,
            'integer': Integer,
            'number': Float,
            'date': CustomDate,
            'datetime': DateTime,
            'geopoint': String,  # Assuming geopoint is stored as a string
        }

        # Iterate over fields, attach their names as class attributes and lookup their types
        for field in sorted(schema_fields, key=lambda x: x['index']):
            column_name = field['name']
            column_type = type_mapping.get(field['type'], Text)  # Default to Text for unknown types
            if field['name'] == 'processid':
                class_attrs[column_name] = Column(column_type, nullable=False, index=True)
            else:
                class_attrs[column_name] = Column(column_type)

        # Return the class, its superclass, and the attached class attributes
        return type('Barcode', (Base,), class_attrs)

    @classmethod
    def create_dynamic_taxon_class(cls, json_schema):
        """
        Dynamically creates an ORM class for the fields in the BCDM data corresponding with the Linnean lineage
        of a barcode record.
        :param json_schema:
        :return:
        """

        # Initial attributes for a SQLAlchemy ORM class, i.e. the table name to connect to and a primary key
        class_attrs = {
            '__tablename__': 'taxon',
            'taxon_id': Column(Integer, primary_key=True)
        }

        # Traverse the JSON to get to the field definitions
        resources = json_schema['resources']
        schema_fields = next(
            resource['schema']['fields'] for resource in resources if resource['profile'] == 'tabular-data-resource')

        # Only process those field definitions whose description includes 'Taxonomy: '
        for field in schema_fields:
            if "Taxonomy: " in field.get('description', ''):
                column_name = field['name']
                column_type = CustomDate if field['type'] == 'date' else Text
                class_attrs[column_name] = Column(column_type, nullable=False, index=True)

        # Add additional fields from the original Taxon class
        class_attrs['bin_uri'] = Column(Text, nullable=False, index=True)
        class_attrs['opentol_id'] = Column(Integer, index=True)

        return type('Taxon', (Base,), class_attrs)


class Node(Base):
    """ORM class for the colunms of database files created by DBTree"""
    __tablename__ = 'node'

    # Primary key
    id = Column(Integer, primary_key=True)

    # Foreign key reference to the parent node
    parent = Column(Integer)

    # Pre-order, autoincrementing integer index
    left = Column(Integer)

    # Post-order, autoincrementing integer index
    right = Column(Integer)

    # Tip name. This is going to be an OTT id, i.e. a string 'ott[0-9]+'
    name = Column(String(20))

    # Branch length
    length = Column(Float)

    # Patristic distance to the root
    height = Column(Float)
