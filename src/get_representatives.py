import numpy as np
import pandas as pd
import os


def loop_over_families(flist):
    for family in flist:
        if family != "backbone":
            print("not same",family)
            mode = "w" if family == flist[0] else "a"
            csv_file = "data/fasta/family/{}/altered_matrix_{}.txt".format(family, family)
            df = pd.read_csv(csv_file, delimiter="\t")
            get_highest(df, mode)

def get_highest(df, mode):
    highest = df.max()
    result_list = search_pos(df, {max(highest)})
    print(resultl_list)
    # Get row and column name
    representatives = get_corresponding_ott(df, result_list[0][0], result_list[0][1], mode)
    print(representatives)


def search_pos(df_data: pd.DataFrame, search: set) -> list:
    nda_values = df_data.values
    tuple_index = np.where(np.isin(nda_values, [e for e in search]))
    return [(row, col, nda_values[row][col]) for row, col in zip(tuple_index[0], tuple_index[1])]


def get_corresponding_ott(df, col_pos, row_pos, mode):
    representatives = ""
    colname = df.columns[col_pos]
    rowname = df.index[row_pos]
    representatives += colname + "\n"
    representatives += rowname + "\n"
    with open("data/fasta/family/backbone/submatrices_representatives.txt", mode) as output:
        output.write(str(representatives))
    return representatives

if __name__ == "__main__":
    path2 = os.getcwd()  # Get current working directory
    path = snakemake.input[0]  # noqa: F821
    total_path = os.path.join(path2, path)
    family_list = os.listdir(os.path.join(path2, path))
    print(family_list)


    #input = snakemake.input[0]  # noqa: F821
    loop_over_families(family_list)