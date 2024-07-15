import argparse
import csv
import gzip
import json
import logging
import requests
import sys
from typing import TextIO, Iterator, Tuple


def setup_logging(verbosity: int) -> None:
    """
    Set up logging configuration based on verbosity level.

    Args:
        verbosity (int): The verbosity level (0 for INFO, 1 for DEBUG, 2+ for more detailed DEBUG).
    """
    levels = [logging.INFO, logging.DEBUG, logging.DEBUG]
    level = levels[min(verbosity, len(levels) - 1)]

    logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')

    if verbosity >= 2:
        logging.getLogger().setLevel(logging.DEBUG)
        requests_log = logging.getLogger("requests.packages.urllib3")
        requests_log.setLevel(logging.DEBUG)
        requests_log.propagate = True


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: An object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Process CSV file and lookup taxon information.")
    parser.add_argument("file_path", help="Path to the input CSV file (can be .txt or .gz)")
    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="Increase output verbosity (e.g., -v for DEBUG, -vv for more detailed DEBUG)")
    return parser.parse_args()


def open_file(file_path: str) -> TextIO:
    """
    Open a file, handling both plain text and gzipped files.

    Args:
        file_path (str): The path to the file to be opened.

    Returns:
        TextIO: A file-like object for reading the content of the file.
    """
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    return open(file_path, 'r')


def process_csv(file: TextIO) -> Iterator[Tuple[str, str, str]]:
    """
    Process the CSV file and yield relevant information for each row.

    Args:
        file (TextIO): A file-like object containing the CSV data.

    Yields:
        Tuple[str, str, str]: A tuple containing taxonID, genus, and specificEpithet for each row.
    """
    reader = csv.DictReader(file, quotechar='"')
    for row in reader:
        yield row['taxonID'], row['genus'], row['specificEpithet']


def fetch_taxon_info(genus: str, specific_epithet: str) -> Tuple[str, str]:
    """
    Fetch taxon information from the BOLD Systems API.

    Args:
        genus (str): The genus name.
        specific_epithet (str): The specific epithet.

    Returns:
        Tuple[str, str]: A tuple containing the taxid and taxon name if found, or empty strings if not found.
    """
    url = f"https://boldsystems.org/index.php/API_Tax/TaxonSearch?taxName={genus}%20{specific_epithet}"
    logging.debug(f"Sending request to: {url}")
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        logging.debug(f"Received data: {json.dumps(data, indent=2)}")
        if data.get('total_matched_names', 0) > 0:
            match = data['top_matched_names'][0]
            return match.get('taxid', ''), match.get('taxon', '')
    else:
        logging.warning(f"Request failed with status code: {response.status_code}")
    return '', ''


def main() -> None:
    """
    Main function to orchestrate the script's operation.

    This function sets up logging, parses arguments, processes the input file,
    fetches taxon information for each record, and outputs the results as TSV.
    """
    args = parse_arguments()
    setup_logging(args.verbose)

    logging.info(f"Processing file: {args.file_path}")

    print("taxonID\tgenus\tspecificEpithet\ttaxid\ttaxon")

    try:
        with open_file(args.file_path) as file:
            for taxon_id, genus, specific_epithet in process_csv(file):
                logging.debug(f"Processing: taxonID={taxon_id}, genus={genus}, specificEpithet={specific_epithet}")
                taxid, taxon = fetch_taxon_info(genus, specific_epithet)
                print(f"{taxon_id}\t{genus}\t{specific_epithet}\t{taxid}\t{taxon}")
    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        sys.exit(1)

    logging.info("Processing complete")


if __name__ == "__main__":
    main()