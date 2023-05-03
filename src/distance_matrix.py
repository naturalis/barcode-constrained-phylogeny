import numpy as np
from ete3 import Tree       # To use for distance matrix
import pandas as pd
import SeqIO


def create_matrix():
    tree = Tree("Aba.bestTree",
                format=1)  # Format indicates newick structure in file
    nodes = list(tree.get_leaves())  # Get all leaves
    names = []
    dmat = np.zeros((len(nodes), len(nodes)))

    for num in range(len(nodes)):
        names.append(nodes[num].name)
        for num2 in range(num, len(nodes)):
            d = tree.get_distance(nodes[num],
                                  nodes[num2],
                                  topology_only=False)
            dmat[num, num2] = dmat[num2, num] = round(d, 6)  # dmat makes matrix

    dist_df = pd.DataFrame(data=dmat, index=names, columns=names)  # Make into dataframe
    highest = dist_df.max()     # Get all highest from rows
    print(max(highest))     # Get highest values
    result_list = search_matrix(dist_df, {max(highest)})
    # Get row and column name
    get_corresponding_ott(dist_df,result_list[0][0], result_list[0][1])

    dist_df.to_csv("../data/tree_distancematrix.txt",
                   sep='\t',
                   header=True)  # Write distance matrix to txt file


def search_matrix(df_data: pd.DataFrame, search_set: set) -> list:
    nda_values = df_data.values
    tuple_index = np.where(np.isin(nda_values, [e for e in search_set]))
    return [(row, col, nda_values[row][col]) for row, col in zip(tuple_index[0], tuple_index[1])]

def get_corresponding_ott(df, col_pos, row_pos):
    representatives = []
    colname = df.columns[col_pos]
    rowname = df.rows[row_pos]
    representatives.append(colname)
    representatives.append(rowname)
    return representatives


def representatives_to_file(representatives):
    """Important to get the normal sequences in this file"""
    dir = "src/fasta/family/"
    filename = "Abacionidae.fasta"
    with open("{}/{}".format(dir, filename), "r") as input:
        with open("representatives.fasta", "w+") as output:
            for record in SeqIO.parse(input, "fasta"):
                print(record.id)
                for seq in representatives:
                    if record.id == seq:
                        SeqIO.write(record, output, 'fasta')


if __name__ == "__main__":
    create_matrix()


    """" 
    # Make connection to the database
    conn = sqlite3.connect(database)

    # Create a cursor
    cursor = conn.cursor()

    # Close the connection
    conn.close()
    create_matrix()
"""


