import argparse
import util
import os.path
import sqlite3
import opentol
import dendropy

from Bio import SeqIO


def get_ott_ids(pids, conn):
    """
    :param pids: a list of BOLD process IDs
    :param conn: a sqlite3 database handle
    :return:
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
    for exemplar_file in exemplar_files:
        if not os.path.isfile(exemplar_file):
            logger.info(f'Location "{exemplar_file}" is not a file - skipping...')
            continue

        # Ideally, these have OTT IDs that we can use
        logger.info(f'Processing input file {exemplar_file}')
        process_ids = []
        for record in SeqIO.parse(exemplar_file, 'fasta'):
            process_ids.append(record.id)

        # Check the database
        dict_of_lists = get_ott_ids(process_ids, conn)
        if dict_of_lists:
            pids_for_ott.update(dict_of_lists)
        else:
            logger.warning(f'Input file with unconstrained sequences (maybe extinct?): {exemplar_file}')
    return pids_for_ott


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
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Configure logging
    logger = util.get_formatted_logger('backbone_constraint', args.verbosity)

    # Configure database connection
    logger.info(f"Going to connect to database {args.database}")
    connection = sqlite3.connect(args.database)

    # Get one-to-many mapping from OTT IDs to process IDs
    pidmap = process_exemplars(args.inaln.split(' '), connection)
    connection.close()

    # Fetch opentol subtree, remap tips
    ids = [int(str(item).removeprefix('ott')) for item in list(pidmap.keys())]
    ott_tree = opentol.get_subtree(ids)
    remap_tips(ott_tree, pidmap)

    # Write output
    logger.info(f'Going to write tree to {args.outtree}')
    with open(args.outtree, "w") as output_file:
        output_file.write(ott_tree.as_string(schema="newick"))

