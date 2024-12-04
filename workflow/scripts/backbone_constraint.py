import argparse
import util
import os.path
import sqlite3
import opentol
import dendropy

from Bio import SeqIO


"""
This script, `backbone_constraint.py`, is responsible for generating a constraint tree for a given family from a SQLite 
database and a set of FASTA files.

The script performs the following steps:
1. Connects to the SQLite database.
2. Extracts the Open Tree of Life (OpenTOL) IDs from the headers of the FASTA files.
3. If there are no IDs, it creates a zero-byte file for the next step in the workflow.
4. If there are IDs, it verifies that all process IDs belong to the same family. If not, it raises an error.
5. If the process IDs belong to the same family, it retrieves a subtree from the OpenTOL Web Service API using the 
   extracted IDs.
6. Remaps the leaf labels of the input tree to the values provided in the pidmap. Possibly grafts additional child 
   nodes if one-to-many mapping occurs.
7. Writes the subtree to an output file in Newick format.

The script uses command line arguments for the input FASTA files, output tree file, SQLite database file, and log level.
The script is invoked by the Snakefile as a shell command with the required arguments in the rule `backbone_constraint`.
"""


def get_ott_ids(pids, conn):
    """
    Given a list of BOLD process IDs, checks to verify that these all belong to the same family
    and if so, maps them to OTT IDs. The OTT IDs come from the synthetic tree (i.e. not 
    'broken' taxa).
    :param pids: a list of BOLD process IDs
    :param conn: a sqlite3 database handle
    :return: a dictionary where key is OTT ID, value is list of process IDs
    """
    pids_for_ott = {}
    family = None

    # Step 1: Verify all pids are in the same family
    for pid in pids:
        query = 'SELECT t.family FROM barcode b JOIN taxon t ON b.taxon_id = t.taxon_id WHERE b.processid = ?'
        cursor = conn.execute(query, (pid,))
        result = cursor.fetchone()
        if result:
            if family is None:
                family = result[0]
            elif family != result[0]:
                raise ValueError(f"Process IDs {pids} belong to different families.")
    if family is None:
        raise ValueError(f"No valid family found for the given process IDs {pids}.")

    # Step 2 and 3: Find a representative node for the family
    node_query = """SELECT n.name FROM node n JOIN taxon t ON 'ott' || t.opentol_id = n.name WHERE t.family = ?
        AND t.opentol_id IS NOT NULL
        LIMIT 1
    """
    node_cursor = conn.execute(node_query, (family,))
    node_result = node_cursor.fetchone()

    if node_result:
        node_name = node_result[0]
        pids_for_ott[node_name] = pids

    return pids_for_ott


def process_exemplars(exemplar_files, conn):
    """
    Iterates over the provided FASTA files, attempting to map the sequence IDs (which are BOLD process IDs) to OpenToL
    IDs that can go into the backbone constraint.
    :param exemplar_files: list of FASTA files
    :param conn: sqlite3 database connection
    :return:
    """
    # Gather input files
    logger.info('Going to process input files')
    logger.debug(exemplar_files)
    pids_for_ott = {}
    extinct_pids = []
    for exemplar_file in exemplar_files:
        if not os.path.isfile(exemplar_file):
            logger.info(f'Location "{exemplar_file}" is not a file - skipping...')
            continue

        # Read process IDs from exemplar file into a list
        logger.info(f'Processing input file {exemplar_file}')
        process_ids = []
        for record in SeqIO.parse(exemplar_file, 'fasta'):
            process_ids.append(record.id)
        if not process_ids == []:
            # Get an OTT ID for the process IDs in the list
            dict_of_lists = get_ott_ids(process_ids, conn)
            if dict_of_lists:
                pids_for_ott.update(dict_of_lists)
            else:

                # The process IDs have no family that occurs in the OpenTree. In two cases now this
                # occurred because the focal exemplar file consisted entirely of extinct taxa 
                # (i.e. an entire extinct family, such as Megaladapidae). The least bad thing to do
                # with these is to remove them from the backbone.
                logger.warning(f'Input file with unconstrained sequences (maybe extinct?): {exemplar_file}')
                extinct_pids.extend(process_ids)
    return pids_for_ott, extinct_pids


def remap_tips(tree, pidmap):
    """
    Remaps the leaf labels of the input tree to the values provided in the pidmap. Possibly grafts additional child
    nodes if one-to-many mapping occurs.
    :param tree: a dendropy tree
    :param pidmap: a dictionary of lists
    :return:
    """
    logger.info('Going to remap backbone tree')
    logger.debug(pidmap)
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            name = node.taxon.label
            if str(name).startswith('ott'):
                if len(pidmap[name]) == 1:

                    # map ott to pid
                    node.taxon.label = pidmap[name][0]
                else:

                    # Iterate over all pids subtended by the ott
                    for process in pidmap[name]:

                        # Create a new dendropy taxon and node, associate them, append to parent
                        new_taxon = dendropy.Taxon(label=process)
                        new_node = dendropy.Node()
                        new_node.taxon = new_taxon
                        node.add_child(new_node)
                        logger.info(f'Added child {process} to {name}')


if __name__ == '__main__':
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-d', '--database', required=True, help='SQLite database file')
    parser.add_argument('-i', '--inaln', required=True, help='Input exemplar FASTA files')
    parser.add_argument('-o', '--outtree', required=True, help="Output constraint tree")
    parser.add_argument('-e', '--extinctpids', required=True, help='Putatively extinct PIDs')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Configure logging
    logger = util.get_formatted_logger('backbone_constraint', args.verbosity)

    # Configure database connection
    logger.info(f"Going to connect to database {args.database}")
    connection = sqlite3.connect(args.database)

    # Get one-to-many mapping from OTT IDs to process IDs and store extinct PIDs
    pidmap, extinctpids = process_exemplars(args.inaln.split(' '), connection)
    connection.close()
    if len(extinctpids) != 0:
        with open(args.extinctpids, 'w') as file:
            for pid in extinctpids:
                file.write(f"{pid}\n")

    # Fetch opentol subtree, remap tips
    ids = [int(str(item).removeprefix('ott')) for item in list(pidmap.keys())]
    logger.info(f'Getting subtree for tips {ids}')
    ott_tree = opentol.get_subtree(ids)
    remap_tips(ott_tree, pidmap)

    # Write output
    logger.info(f'Going to write tree to {args.outtree}')
    with open(args.outtree, "w") as output_file:
        output_file.write(ott_tree.as_string(schema="newick"))

