import logging
import argparse
import sqlite3

from Bio.AlignIO import read as read_alignment
from Bio.Phylo import read as read_newick, write as write_newick

logging.basicConfig()
logger = logging.getLogger('prep_raxml')


def preprocess_aln(inaln, connection):
    logger.info(f"Going to preprocess alignment {inaln}")
    mapping = {}
    alnmap = {}
    aln = read_alignment(inaln, 'phylip')
    for seq in aln:

        # Clean up ID, which starts with barcode_id
        barcode_id = seq.id.split('|')[0]

        # Map barcode ID to process_id and opentol_id
        tuple = connection.execute(f"""
            SELECT b.processid, t.opentol_id 
            FROM barcode b, taxon t 
            WHERE b.barcode_id={barcode_id} AND b.taxon_id=t.taxon_id
            """).fetchone()
        logger.info(tuple)

        # Skip if no OTT
        if tuple[1] is None:
            logger.info(f'Sequence {barcode_id} has no opentol_id, skipping')
            continue
        key = 'ott' + str(tuple[1])

        alnmap[tuple[0]] = seq.seq

        # Initialize dictionary for one-to-many mapping opentol_id => [ process_id, ... ]
        if key not in mapping:
            mapping[key] = []

        # Add seen process_id
        mapping[key].append(tuple[0])
    return mapping, alnmap


def write_data(intree, outtree, mapping, alnmap, outaln):
    logger.info(f"Going to reformat constraint tree {intree} to {outtree}")
    logger.info(mapping)
    tree = read_newick(intree, 'newick')

    # Map opentol_id to process_id, possibly adding tips if there are multiple process_ids
    # for this opentol_id (which means there are multiple BINs in this species)
    for tip in tree.get_terminals():
        if tip.name is not None:
            processes = mapping[tip.name]

            # Add and label tips if needed
            if len(processes) > 1:
                tip.split(n=len(processes))
                i = 0
                for child in tip.get_terminals():
                    child.name = processes[i]
                    logger.info(f'Added child {processes[i]} to {tip.name}')
                    i += 1
            else:
                tip.name = processes[0]
        else:
            logger.warning(f'Encountered None tip in {intree}, probably empty tree')

    # Remove branch lengths
    for clade in tree.find_clades():
        clade.branch_length = None
    write_newick(tree, outtree, 'newick')

    # Write alignment
    with open(outaln, mode='w+') as fastafh:
        for defline in alnmap:
            fastafh.write(f'>{defline}\n{alnmap[defline]}\n')

        # This keeps only the tips in the guide tree
        # for tip in tree.get_terminals():
        #     name = tip.name
        #     if name is not None and name in alnmap:
        #         seq = alnmap[name]
        #         fastafh.write(f'>{name}\n{seq}\n')


if __name__ == '__main__':
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    parser.add_argument('-t', '--intree', required=True, help='Input Newick tree')
    parser.add_argument('-a', '--inaln', required=True, help='Input aligned PHYLIP file')
    parser.add_argument('-o', '--outtree', required=True, help='Output Newick tree')
    parser.add_argument('-f', '--outaln', required=True, help="Output FASTA alignment")
    parser.add_argument('-d', '--db', required=True, help="SQLite database")
    args = parser.parse_args()

    # Configure logger
    logger.setLevel(args.verbosity)

    # Connect to the database (creates a new file if it doesn't exist)
    logger.info(f"Going to connect to database {args.db}")
    connection = sqlite3.connect(args.db)

    # Write alignment
    ott_process_map, alnmap = preprocess_aln(args.inaln, connection)

    # Write tree
    write_data(args.intree, args.outtree, ott_process_map, alnmap, args.outaln)
