import numpy as np
from ete3 import Tree       # To use for distance matrix
import pandas as pd


def create_matrix():
    tree = Tree("test/test with ott and selfmade id/Aba.bestTree",
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
    result_list = search_coordinate(dist_df, {max(highest)})
    # Printing the rows and colomns that contain the highest value
    print(f"\n\n{'row':<4} {'col':<4} {'name':>10}")
    [print(f"{row:<4} {col:<4} {name:>10}") for row, col, name in result_list]

    dist_df.to_csv("../data/tree_distancematrix.txt",
                   sep='\t',
                   header=True)  # Write distance matrix to txt file

def search_coordinate(df_data: pd.DataFrame, search_set: set) -> list:
    nda_values = df_data.values
    tuple_index = np.where(np.isin(nda_values, [e for e in search_set]))
    return [(row, col, nda_values[row][col]) for row, col in zip(tuple_index[0], tuple_index[1])]


if __name__ == "__main__":
    create_matrix()



