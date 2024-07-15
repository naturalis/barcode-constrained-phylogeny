import sqlite3
import pandas as pd
from Bio import Phylo
import io
import sys

# 136 different colours, for lepidoptera
colours = ["#ff7f0e", "#53777a", "#c02942", "#d6a0c8", "#51b467", "#ff69b4", "#e9d460", "#9b59b6", "#f0e68c", "#99d9ea",
"#ead1dc", "#42d695", "#f5b6c7", "#800000", "#a0c4ff", "#c7ce25", "#ffff00", "#9932cc", "#228b22", "#ff0000", "#556b2f",
"#8fbc8f", "#cd7f32", "#d8bfd8", "#2ca02c", "#ffe4e1", "#4169e1", "#e49b0f", "#915c83", "#ffdead", "#38793c", "#ffa500",
"#85e04b", "#f473b9", "#7cfc00", "#da70d6", "#300030", "#c6e2b0", "#c2dfff", "#00ced1", "#f08080", "#304ffe",
"#aeea00", "#ff007f", "#4e9256", "#c71585", "#d000b6", "#1e90ff", "#e6e6fa", "#333333", "#b2becc", "#90ee90",
"#6495ed", "#ffa07a", "#20b2aa", "#ff4000", "#5f9ea0", "#c3254e", "#d9ffb3", "#2980b9", "#ffc0cb", "#3e82f7",
"#b33939", "#ffebcd", "#808000", "#95odb5", "#c2c2f0", "#008080", "#f000ff", "#27308b", "#c0c0c0", "#a9a9a9",
"#7b68ee", "#f7f0e3", "#1c2833", "#9d2934", "#d2b48c", "#000000", "#8a2be2", "#ffff99", "#2e8b57", "#fff5ee",
"#4682b4", "#d2691e", "#fffacd", "#266168", "#ffc7ff", "#3c272b", "#b22222", "#ff9999", "#29abca", "#ffffcc",
"#436eee", "#d7ccc8", "#1f4b43", "#94474e", "#cab65b", "#536878", "#f4c430", "#242a2e", "#b0e0e6", "#e74c3c",  
"#89d38a", "#f5deb3", "#907bcb", "#ffd700", "#79a7ff", "#ffdd99", "#64b5f6", "#e066ff", "#a9e2b0", "#ff9e00",  
"#69d2e2", "#f9e0e3", "#c27ba0", "#f0f8ff", "#93775c", "#ffd1dc", "#428bca", "#e0c382", "#c850a6", "#f4eeac",
"#808b9e", "#ffe4b5", "#39acac", "#e6beff", "#b2d763", "#ffc6ff", "#38b3a5", "#e3cf57", "#99b3ff", "#f7cac9", 
"#636e72", "#ffcc99", "#29abe2", "#eeccff"]

def gather_fams(database, order):
    """This code takes the barcodes database, and returns all
    families in the order. 
    """
    connection = sqlite3.connect(database)
    cursor = connection.cursor()

    sql = f"SELECT record_id, family FROM barcode WHERE `order` = '{order}'"
    fam = pd.read_sql_query(sql, connection)
    fam[['record_id', 'marker']] = fam['record_id'].str.split('.', expand=True)
    print(fam)
    print(len(set(fam["family"])))

    return fam


def write_colours(fam):
    global colours
    colours = list(set(colours))

    # Create a dictionary mapping families to colors (assuming unique families)
    family_colors = dict(zip(fam['family'].unique(), colours[:len(fam['family'].unique())]))

    # Add the 'colour' column based on the family mapping
    fam['colour'] = fam['family'].apply(lambda x: family_colors.get(x))


    # #915|777 branch #00ff00 dashed 5

    clades = []
    f = open('results/pruned.tre','r')
    tree = Phylo.read(io.StringIO(f.read()),'newick')
    for clade in tree.find_clades():
        clades.append(str(clade))

    print(clades[0:5])
    documentation = open("resources/colors_styles.txt", "a")
    for index, row in fam.iterrows():
        if row['record_id'] in clades:
            documentation.write(f"{row['record_id']} range {row['colour']} {row['family']}\n")

if __name__ == '__main__':
    database = sys.argv[1] # bold database
    order = sys.argv[2] # the order

    fam = gather_fams(database, order)
    write_colours(fam)
