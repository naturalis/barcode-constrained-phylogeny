import numpy as np
import pandas as pd
import os

def create_fam_list():
    par_path = os.path.abspath(os.path.join(os.pardir))
    family_list = os.listdir(par_path + "/src/test/matK/trees/")
    flist = []
    for family in family_list:
        family = family.split(".")
        flist.append(family[0])
    return flist

def read_xlsx(excel_file, family):
    df = pd.read_excel(excel_file,sheet_name= family,index_col=[0])
    return df

def loop_over_df(flist, df):
    for family in flist:
        df1 = df.get(family)
        representatives = get_highest(df1)
        print(representatives)

def get_highest(df):
    highest = df.max()
    result_list = search_pos(df, {max(highest)})
    # Get row and column name
    representatives = get_corresponding_ott(df, result_list[0][0], result_list[0][1])
    print(representatives)
    return representatives


def search_pos(df_data: pd.DataFrame, search: set) -> list:
    nda_values = df_data.values
    tuple_index = np.where(np.isin(nda_values, [e for e in search]))
    return [(row, col, nda_values[row][col]) for row, col in zip(tuple_index[0], tuple_index[1])]


def get_corresponding_ott(df, col_pos, row_pos):
    representatives = ""
    colname = df.columns[col_pos]
    rowname = df.index[row_pos]
    representatives += colname + "\n"
    representatives += rowname + "\n"
    with open("test/matK/matrix/submatrices_representatives.txt", "a+") as output:
        output.write(str(representatives))
    return representatives

if __name__ == "__main__":
    input = snakemake.input[0]
    flist = create_fam_list()
    df = read_xlsx(input, flist)
    loop_over_df(flist, df)