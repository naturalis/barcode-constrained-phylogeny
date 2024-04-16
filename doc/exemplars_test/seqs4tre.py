import dendropy
import sqlite3
import argparse


def tiplabels(treefile):
    tree = dendropy.Tree.get(
        path=treefile,
        schema="newick",
        rooting='default-rooted'
    )
    return [tip.taxon.label for tip in tree.seed_node.leaf_nodes()]


def printseqs(pids, dbfile):
    connection = sqlite3.connect(dbfile)
    query = 'SELECT nuc FROM barcode WHERE processid = ?'
    for pid in pids:
        seq = connection.execute(query, (pid,)).fetchone()[0]
        seq = str(seq).replace('-', '')
        print(f'>{pid}')
        print(seq)


if __name__ == '__main__':
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-d', '--database', required=True, help='SQLite database file')
    parser.add_argument('-i', '--intree', required=True, help='Input tree file')
    args = parser.parse_args()

    # Read tree, get tip labels
    labels = tiplabels(args.intree)

    # Print sequences
    printseqs(labels, args.database)

