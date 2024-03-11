import argparse
import util
import opentol
import sqlite3


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

    # process the ingroup file
    with open(unaligned, 'r') as file:
        for line in file:
            if line.startswith('>'):
                parts = line.strip().split('|')
                if len(parts) > 1:
                    # Extract the second element and remove 'ott' prefix
                    id_with_prefix = parts[1]
                    id_number = id_with_prefix.replace('ott', '')
                    if id_number != 'None':
                        ids.append(int(id_number))

    # process the outgroup file
    with open(outgroups, 'r') as file:
        for line in file:
            if line.startswith('>'):
                pid = line.strip().removeprefix('>')
                sql = f"SELECT t.opentol_id FROM taxon t, barcode b WHERE t.taxon_id=b.taxon_id and b.processid='{pid}'"
                ott = conn.execute(sql).fetchone()
                ids.append(ott[0])

    # Remove all 'None' entries
    cleaned_list = [item for item in ids if item != 'None']
    return cleaned_list


if __name__ == '__main__':
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-i', '--ingroup', required=True, help='FASTA file with the ingroup')
    parser.add_argument('-g', '--outgroups', required=True, help='FASTA file with outgroup taxa')
    parser.add_argument('-o', '--outtree', required=True, help="Output constraint tree")
    parser.add_argument('-d', '--database', required=True, help='SQLite database')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Configure logging
    logger = util.get_formatted_logger('family_constraint', args.verbosity)

    # Connect to the database (creates a new file if it doesn't exist)
    logger.info(f"Going to connect to database {args.database}")
    conn = sqlite3.connect(args.database)

    # Read input alignment, get ott IDs
    ott_ids = extract_id_from_fasta(args.ingroup, args.outgroups)

    # If we have no IDs at all, we write a zero byte file for run_raxml
    if len(ott_ids) == 0:
        logger.warning('There were zero OTT IDs in the input file')
        with open(args.outtree, "a"):
            pass

    else:
        # Get subtree from OpenToL WS API
        tree = opentol.get_subtree(ott_ids)

        # Write output
        logger.info(f'Going to write tree to {args.outtree}')
        with open(args.outtree, "w") as output_file:
            output_file.write(tree.as_string(schema="newick"))


