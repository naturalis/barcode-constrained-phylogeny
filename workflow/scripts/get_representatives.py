import numpy as np
import pandas as pd
import os
import logging


logger = logging.getLogger(__name__)


def call_functions(csv_file, representatives_file):
        logger.info("Getting the highest otts from distance matrix.")
        df = pd.read_csv(csv_file, delimiter="\t")
        get_highest(df, representatives_file)

def get_highest(df, representatives_file):
    """Declare the highest value in the dataframe.
    Search for the position of the highest value in the dataframe.
    Use that position to find the corresponding highest value.
    """
    highest = df.max()  # The highest value in the dataframe
    row_index = search_pos(df, highest[0], df[highest[0]].max())    # Search for position of the highest value
    logger.info(f"{highest} found at position {row_index}")
    # Get row and column name
    representatives = get_corresponding_ott(highest, highest[0], row_index, representatives_file)
    logger.info(f"Writing the highest otts {representatives} to representatives file")


def search_pos(df_data: pd.DataFrame, col_name, value):
    """Search in dataframe for the highest value and return the corresponding otts.
    return the values and correspongding otts if found, else return an empty list.
    """
    try:
        row_index = df_data[df_data[col_name] == value].index[0]
        return row_index
    except:
        print("Unknown errror occurred.")
        logger.debug("Something went wrong.")
        return []



def get_corresponding_ott(df, colname, row_pos, representatives_file):
    """Get the corresponding ott from certain value.
    This is done using the position of the value.
    A file i
    """
    representatives = ""
    rowname = df.index[row_pos]
    representatives += colname + "\n"
    representatives += rowname + "\n"
    with open(representatives_file, "a+") as output:
        output.write(str(representatives))
    return representatives

if __name__ == "__main__":
    csv_file = snakemake.input[0]  # noqa: F821
    representatives_file = snakemake.output[0] # noqa: F821
    call_functions(csv_file, representatives_file)
