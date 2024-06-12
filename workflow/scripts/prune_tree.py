import sys
from Bio import Phylo
import io
import numpy as np
import matplotlib.pyplot as plt


def get_lengths(filename_tree):
    f = open(filename_tree,'r')
    tree = Phylo.read(io.StringIO(f.read()),'newick')
    lengths = []
    for clade in tree.find_clades():
        length = clade.branch_length
        if not length == None:
            lengths.append(float(length))
        else:
            print(None)
    return lengths


def boxplot_filter(lengths):
    # boxplot
    plt.boxplot(np.log10(lenghts))
    # IQR
    Q1 = np.percentile(lenghts, 25, method='midpoint')
    Q3 = np.percentile(lenghts, 75, method='midpoint')
    IQR = Q3 - Q1
    filter_at = Q3+1.5*IQR
    print("Boxplot:")
    print(f"The IQR is {IQR} and the max value is {filter_at}")
    print(f"Using this method would filter out {len([x for x in lenghts if x >= filter_at])} tips.")
    plt.savefig("doc/box_log10_branchlengths.png")
    plt.clf()
    plt.boxplot(lenghts)
    plt.savefig("doc/box_branchlengths.png")
    plt.clf()
    return filter_at


def above_threshold(lenghts, threshold=1):
    print(f"Threshold: set to {threshold}")
    print(f"Using this method would filter out {len([x for x in lenghts if x >= threshold])} tips.")
    histogram(lenghts,threshold)
    return threshold


def histogram(lenghts, threshold):
    # histogram
    data_0_1 = [x for x in lenghts if 0 <= x <= threshold]
    data_outliers = [threshold for x in lenghts if x > threshold] # > threshold is set to threshold value to compress plot

    # Create the histogram with logarithmic scale on x-axis
    plt.hist(data_0_1, bins=100, density=True, label=f"0-{threshold} Range")
    plt.hist(data_outliers, density=True, label=f"Outliers (>={threshold})")
    plt.title("Distribution of Data with Outliers (Logarithmic Scale)")
    plt.legend()
    plt.savefig("doc/hist_branchlengths.png")
    plt.clf()


def sd_filter(lenghts, sds=4):
    print(f"Standard deviation: set to {sds}")
    elements = np.array(lenghts)

    mean = np.mean(elements, axis=0)
    sd = np.std(elements, axis=0)

    filter_at = mean + sds * sd
    filtered = [x for x in lenghts if (x >= filter_at)]
    print(f"The filter threshold is set at {filter_at}")
    print(f"Using this method would filter out {len(filtered)} tips.")
    return filter_at


def prune_tree(filename_tree, filter_at):
    f = open(filename_tree,'r')
    tree = Phylo.read(io.StringIO(f.read()),'newick')
    pruned = []
    for clade in tree.find_clades():
        length = clade.branch_length
        if not length == None:
            if float(length) >= filter_at:
                try:
                    tree.prune(clade)
                    pruned.append(str(clade))
                except ValueError: # This is a clade, not a terminal
                    for terminal in clade.get_terminals():
                        if not str(terminal) in pruned: # to ensure it hasn't been pruned already
                            tree.prune(terminal)
                        pruned.append(str(terminal))
        else:
            print(None)
    return tree


if __name__ == '__main__':
    filename_tree = sys.argv[1]
    print(filename_tree)
    lenghts = get_lengths(filename_tree)
    filter_at_bx = boxplot_filter(lenghts)
    filter_at_thr = above_threshold(lenghts, threshold=1)
    filter_at_sd = sd_filter(lenghts, sds=5)
    pruned_tree = prune_tree(filename_tree, filter_at_sd)
    save_to = sys.argv[2]
    Phylo.write(pruned_tree, save_to, "newick")
