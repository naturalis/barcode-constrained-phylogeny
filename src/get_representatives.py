import numpy as np
import pandas as pd
import os
import logging


logging.basicConfig(snakemake.params.log_level) # noqa: F821
logger = logging.getLogger(__name__)


def loop_over_families(csv_file, representatives_file):
        logger.info("Getting the highest otts from distance matrix.")
        #mode = "w" if family == flist[0] else "a"
        #csv_file = "fasta/alignment/{}/altered_matrix_{}.txt".format(family, family)
        df = pd.read_csv(csv_file, delimiter="\t")
        get_highest(df, representatives_file)

def get_highest(df, representatives_file):
    """Declare the highest value in the dataframe.
    Search for the position of the highest value in the dataframe.
    Use that position to find the corresponding highest value.
    """
    highest = df.max()  # The highest value in the dataframe
    result_list = search_pos(df, {max(highest)})    # Search for position of the highest value
    print(result_list)
    logger.info(f"{highest} found at position {result_list}")
    # Get row and column name
    representatives = get_corresponding_ott(df, result_list[0][0], result_list[0][1], representatives_file)
    logger.info(f"Writing the highest otts {representatives} to representatives file")


def search_pos(df_data: pd.DataFrame, search: set) -> list:
    """Search in dataframe for the highest value and return the corresponding otts.
    return the values and correspongding otts if found, else return an empty list.
    """
    try:
        nda_values = df_data.values
        tuple_index = np.where(np.isin(nda_values, [e for e in search]))
        return [(row, col, nda_values[row][col]) for row, col in zip(tuple_index[0], tuple_index[1])]
    except:
        print("Unknown errror occurred.")
        logger.debug("Something went wrong.")
        return []



def get_corresponding_ott(df, col_pos, row_pos, representatives_file):
    """Get the corresponding ott from certain value.
    This is done using the position of the value.
    A file i
    """
    representatives = ""
    colname = df.columns[col_pos]
    rowname = df.index[row_pos]
    representatives += colname + "\n"
    representatives += rowname + "\n"
    with open(representatives_file, "a+") as output:
        output.write(str(representatives))
    return representatives

if __name__ == "__main__":
    #path2 = os.getcwd()  # Get current working directory
    #path = snakemake.input[0]  # noqa: F821
    #total_path = os.path.join(path2, path)
    #family_list = os.listdir(os.path.join(path2, path))
    #print(family_list)
    csv_file = snakemake.input[0]  # noqa: F821
    representatives_file = snakemake.input[1] # noqa: F821
    loop_over_families(csv_file, representatives_file)