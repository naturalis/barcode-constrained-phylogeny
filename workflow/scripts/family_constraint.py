import argparse
import util
import opentol
import sqlite3


"""
This script, `family_constraint.py`, is responsible for generating a constraint tree for a given family from a SQLite 
database and a set of FASTA files.

The script performs the following steps:
1. Connects to the SQLite database.
2. Extracts the Open Tree of Life (OpenTOL) IDs from the headers of the ingroup and outgroup FASTA files.
3. If there are no IDs, it creates a zero-byte file for the next step in the workflow.
4. If there are IDs, it retrieves a subtree from the OpenTOL Web Service API using the extracted IDs.
5. Writes the subtree to an output file in Newick format.

The script uses command line arguments for the ingroup FASTA file, outgroup FASTA file, output tree file, SQLite 
database file, and log level. The script is invoked by the Snakefile as a shell command with the required arguments in
the rule `family_constraint`.
"""


def extract_id_from_fasta(unaligned, outgroups):
    """
    Extracts the ott IDs from a FASTA file. This operates specifically on
    FASTA headers where the ott ID is the second element (i.e. at index 1)
    in the pipe-separated FASTA header.
    :param unaligned: the location of an ingroup FASTA file
    :param outgroups: the location of an outgroups FASTA file
    :return: a list of ott IDs
    """
    ids = []

    # Process the ingroup file
    with open(unaligned, 'r') as file:
        for line in file:
            if line.startswith('>'):
                parts = line.strip().split('|')
                if len(parts) > 1:
                    id_with_prefix = parts[1]
                    if id_with_prefix != "ottNone":
                        id_number = id_with_prefix.replace('ott', '')
                        ids.append(int(id_number))

    # Process the outgroup file
    with open(outgroups, 'r') as file:
        for line in file:
            if line.startswith('>'):
                pid = line.strip().removeprefix('>')
                sql = f"SELECT t.opentol_id FROM taxon t, barcode b WHERE t.taxon_id=b.taxon_id and b.processid='{pid}'"
                ott = conn.execute(sql).fetchone()
                if ott:
                    ids.append(ott[0])

    return ids


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-i', '--ingroup', required=True, help='FASTA file with the ingroup')
    parser.add_argument('-g', '--outgroups', required=True, help='FASTA file with outgroup taxa')
    parser.add_argument('-o', '--outtree', required=True, help="Output constraint tree")
    parser.add_argument('-d', '--database', required=True, help='SQLite database')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    logger = util.get_formatted_logger('family_constraint', args.verbosity)
    logger.info(f"Connecting to database {args.database}")
    conn = sqlite3.connect(args.database)

    ott_ids = extract_id_from_fasta(args.ingroup, args.outgroups)

    if len(ott_ids) == 0:
        logger.warning('No valid OTT IDs found in the input files.')
        with open(args.outtree, "a"):  # Create an empty output file
            pass

    else:
        try:
            # Try to fetch the subtree using the OpenToL API
            tree = opentol.get_subtree(ott_ids)
            if tree is None:
                logger.error('The API returned None, indicating no tree could be generated.')
                with open(args.outtree, "a"):  # Create an empty file
                    pass
            else:
                # Confirm tree has correct format before writing
                if hasattr(tree, 'as_string'):
                    logger.info(f'Writing tree to {args.outtree}')
                    with open(args.outtree, "w") as output_file:
                        output_file.write(tree.as_string(schema="newick"))
                else:
                    logger.error('The tree object does not support "as_string"; skipping file writing.')
                    with open(args.outtree, "a"):
                        pass

        except Exception as e:
            logger.error(f"Failed to fetch or write subtree: {e}")
            with open(args.outtree, "a"):  # Create an empty file on failure
                pass