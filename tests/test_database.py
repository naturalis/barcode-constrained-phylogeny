import pytest
from pathlib import Path
from bactria.database import Database
import shutil
import tempfile

# Path to the actual data package
DATA_PACKAGE_PATH = Path(__file__).parent / 'BOLD_Public.18-Dec-2023'


@pytest.fixture(scope="module")
def temp_data_package():
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as tmpdirname:
        temp_path = Path(tmpdirname)

        # Copy the entire data package to the temporary directory
        shutil.copytree(DATA_PACKAGE_PATH, temp_path / DATA_PACKAGE_PATH.name)

        yield temp_path / DATA_PACKAGE_PATH.name


@pytest.fixture(scope="module")
def db_instance(temp_data_package):
    # Create a temporary database file
    db_file = temp_data_package.parent / "test.db"

    # Initialize the Database instance
    db = Database(db_file, temp_data_package)
    db.make_schema()

    yield db

    # Clean up (close connections, etc.) if necessary
    db.connection.close()


def test_database_initialization(db_instance, temp_data_package):
    assert db_instance.db_file == temp_data_package.parent / "test.db"
    assert db_instance.data_package == temp_data_package
    assert db_instance.Barcode is not None
    assert db_instance.Taxon is not None


def test_make_orm(db_instance):
    # Check for some expected attributes in Barcode and Taxon classes
    # You may need to adjust these based on your actual schema
    assert hasattr(db_instance.Barcode, 'processid')
    assert hasattr(db_instance.Barcode, 'bin_uri')
    assert hasattr(db_instance.Taxon, 'kingdom')
    assert hasattr(db_instance.Taxon, 'phylum')
    assert hasattr(db_instance.Taxon, 'class')


def test_make_schema(db_instance):
    assert db_instance.initialized == True

    # Check if tables are created
    cursor = db_instance.connection.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = cursor.fetchall()
    assert ('barcode',) in tables
    assert ('taxon',) in tables
    assert ('node',) in tables


def test_load_bcdm(db_instance, temp_data_package):
    bcdm_file = temp_data_package / f"{temp_data_package.name}.tsv"
    db_instance.load_bcdm(bcdm_file)

    # Start a new session
    db_instance.start_session()

    try:
        # Check if data is loaded
        cursor = db_instance.connection.cursor()
        cursor.execute("SELECT COUNT(*) FROM barcode;")
        barcode_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM taxon;")
        taxon_count = cursor.fetchone()[0]

        bc = db_instance.get_barcode({ 'processid': 'AAASF001-17'})[0]

        assert barcode_count == 999  # 1000 lines in barcode TSV, including header
        assert taxon_count > 0  # There should be some distinct taxa
        assert bc.identification == 'Lutzomyia cruciata'

    finally:
        # End the session
        db_instance.end_session()



def test_invalid_data_package_path():
    with pytest.raises(FileNotFoundError):
        Database(Path("test.db"), Path("invalid_path"))


def test_make_schema_without_orm():
    db = Database.__new__(Database)  # Create instance without calling __init__
    db.Barcode = None
    db.Taxon = None
    with pytest.raises(ValueError):
        db.make_schema()
