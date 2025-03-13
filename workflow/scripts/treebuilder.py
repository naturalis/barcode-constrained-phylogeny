import os
import sys
import logging
import argparse
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional
import pandas as pd
from Bio import AlignIO, Phylo
from Bio.Phylo.BaseTree import Tree, Clade

"""
BOLD Taxonomy Tree Builder

This script builds a taxonomic tree from BOLD process IDs listed in a FASTA alignment file.
It extracts the process IDs from the FASTA file, looks up their taxonomic information in a 
BOLD BCDM (Barcode Core Data Model) TSV file, and constructs a hierarchical tree with the
process IDs at the tips and taxonomic information at the internal nodes.

Inputs:
-------
1. FASTA Alignment File (-f, --fasta):
  - Standard FASTA format with BOLD process IDs as sequence identifiers
  - The process ID should be the first word in the FASTA defline
  - Example: >AAASF001-17

2. BOLD BCDM TSV File (-b, --bold):
  - Standard BOLD BCDM TSV file containing taxonomic information
  - Must include 'processid' column and taxonomic level columns (kingdom, phylum, etc.)
  - See BOLD documentation for complete format specifications:
    https://github.com/boldsystems-central/BCDM/blob/main/field_definitions.tsv

Output:
-------
- Newick format tree file with process IDs at the tips and optional taxonomic labels at
 internal nodes (-o, --output)

Output Customization:
--------------------
1. Internal Node Labels (-n, --nodelabels):
  - By default, internal node labels are removed from the output
  - Use --nodelabels to include taxonomic labels at internal nodes
  - Note: Many phylogenetic programs do NOT support internal node labels in Newick files
  - However, these labels can be useful when working with the Open Tree of Life web services
    to further resolve polytomies based on known phylogenetic relationships

2. Unbranched Internal Nodes (-c, --collapse):
  - Use --collapse to remove unbranching internal nodes (nodes with only one child)
  - Many phylogenetic analysis programs do NOT correctly process trees with unbranched
    internal nodes and may produce errors or unexpected results
  - Collapsing these nodes maintains the same tree topology from a phylogenetic perspective

Usage Examples:
--------------
Basic usage:
   python treebuilder.py -f sequences.fasta -b bold_data.tsv -o taxonomy_tree.nwk

With all options:
   python treebuilder.py -f sequences.fasta -b bold_data.tsv -o taxonomy_tree.nwk -c -n -v

Notes:
------
- Branch lengths are removed from the output tree as they are not meaningful in this context
- Missing taxonomic levels (marked as 'None' or empty in BOLD data) are skipped
- The script handles large BOLD data files fairly efficiently by processing in chunks
"""

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('bold-tree-builder')

def extract_process_ids(fasta_file: str) -> Set[str]:
    """
    Extract BOLD process IDs from a FASTA alignment file.

    :param fasta_file: Path to the FASTA alignment file
    :return: A set of BOLD process IDs
    """
    logger.info(f"Extracting process IDs from {fasta_file}")
    process_ids = set()

    try:
        alignment = AlignIO.read(fasta_file, "fasta")
        for record in alignment:
            process_id = record.id.strip()
            if process_id:
                process_ids.add(process_id)

        logger.info(f"Extracted {len(process_ids)} process IDs")
        logger.debug(f"Process IDs: {', '.join(list(process_ids)[:5])}...")
        return process_ids

    except Exception as e:
        logger.error(f"Error reading FASTA file: {e}")
        sys.exit(1)


def read_bold_taxonomy(bold_file: str, process_ids: Set[str]) -> pd.DataFrame:
    """
    Read BOLD taxonomy information for the specified process IDs.

    :param bold_file: Path to the BOLD BCDM TSV file
    :param process_ids: Set of BOLD process IDs to filter
    :return: DataFrame containing taxonomy information for the specified process IDs
    """
    logger.info(f"Reading taxonomy data from {bold_file}")

    try:
        # Define taxonomic levels in hierarchical order
        taxonomy_levels = [
            'kingdom', 'phylum', 'class', 'order',
            'family', 'subfamily', 'genus', 'species', 'subspecies'
        ]

        # Read the BOLD TSV file, focusing on taxonomic columns and processid
        cols_to_read = ['processid'] + taxonomy_levels

        # Read the file in chunks to handle large files efficiently
        df_chunks = pd.read_csv(
            bold_file,
            sep='\t',
            usecols=cols_to_read,
            dtype=str,
            chunksize=100000
        )

        # Combine chunks, filtering for our process IDs
        records = []
        for chunk in df_chunks:
            filtered_chunk = chunk[chunk['processid'].isin(process_ids)]
            records.append(filtered_chunk)

            # If we've found all our process IDs, we can stop reading
            if len(set.union(*[set(df['processid']) for df in records])) == len(process_ids):
                logger.debug("Found all requested process IDs, stopping further reading")
                break

        if not records:
            logger.error("No matching records found in BOLD data")
            sys.exit(1)

        taxonomy_df = pd.concat(records, ignore_index=True)
        found_ids = set(taxonomy_df['processid'])
        missing_ids = process_ids - found_ids

        if missing_ids:
            logger.warning(f"Could not find {len(missing_ids)} process IDs in BOLD data")
            logger.debug(f"Missing IDs: {', '.join(list(missing_ids)[:5])}...")

        logger.info(f"Found taxonomy data for {len(found_ids)} process IDs")
        return taxonomy_df

    except Exception as e:
        logger.error(f"Error reading BOLD taxonomy file: {e}")
        sys.exit(1)


def build_taxonomy_paths(taxonomy_df: pd.DataFrame) -> Dict[str, List[str]]:
    """
    Build taxonomy paths for each process ID.

    :param taxonomy_df: DataFrame containing taxonomy information
    :return: Dictionary mapping process IDs to their taxonomy paths
    """
    logger.info("Building taxonomy paths for each process ID")

    taxonomy_levels = [
        'kingdom', 'phylum', 'class', 'order',
        'family', 'subfamily', 'genus', 'species', 'subspecies'
    ]

    taxonomy_paths = {}

    for _, row in taxonomy_df.iterrows():
        process_id = row['processid']
        path = []

        # Build path, skipping 'None' or empty values
        for level in taxonomy_levels:
            taxon = row[level]
            if pd.notna(taxon) and taxon != 'None' and taxon.strip():
                path.append(taxon)

        # Add the process ID as the final element
        path.append(process_id)
        taxonomy_paths[process_id] = path

    logger.debug(f"Sample path: {list(taxonomy_paths.values())[0]}")
    return taxonomy_paths


def build_tree_from_paths(taxonomy_paths: Dict[str, List[str]]) -> Tree:
    """
    Build a BioPython tree from taxonomy paths.

    :param taxonomy_paths: Dictionary mapping process IDs to taxonomy paths
    :return: BioPython Tree object
    """
    logger.info("Building taxonomic tree")

    # Create the root node
    tree = Tree(rooted=True, root=Clade(name="root"))

    # Track nodes by path to handle shared ancestry
    path_to_clade = {"": tree.root}

    # Sort paths to ensure parent nodes are created before children
    process_ids = sorted(taxonomy_paths.keys())

    for process_id in process_ids:
        path = taxonomy_paths[process_id]
        current_path = ""
        parent_path = ""

        # Process each level in the path except the last (which is the process ID)
        for i, taxon in enumerate(path[:-1]):
            current_path = f"{current_path}/{taxon}" if current_path else taxon

            # Create new node if this path doesn't exist yet
            if current_path not in path_to_clade:
                new_clade = Clade(name=taxon)
                path_to_clade[parent_path].clades.append(new_clade)
                path_to_clade[current_path] = new_clade

            parent_path = current_path

        # Add the leaf node (process ID)
        leaf_clade = Clade(name=process_id)
        path_to_clade[current_path].clades.append(leaf_clade)

    logger.info(f"Tree built with {len(list(tree.find_clades()))} nodes")
    return tree


def fix_root_unifurcation(tree: Tree) -> None:
    """
    Fix the case where the root has only one child by making that child the new root.

    :param tree: BioPython Tree object
    """
    if len(tree.root.clades) == 1:
        logger.info("Root has only one child - fixing root unifurcation")
        # The child of the root becomes the new root
        old_root = tree.root
        new_root = old_root.clades[0]

        # Set the new root's properties
        tree.root = new_root

        # If we want to preserve the old root's name in some way
        if old_root.name and old_root.name != "root":
            # We could add the old root's name to the new root if needed
            if new_root.name:
                new_root.name = f"{old_root.name}_{new_root.name}"
            else:
                new_root.name = old_root.name

        logger.info("Root unifurcation fixed")


def collapse_unbranching_nodes(tree: Tree) -> None:
    """
    Remove internal nodes that have only one child using recursion.

    :param tree: BioPython Tree object
    """
    logger.info("Collapsing unbranching internal nodes")

    # Start with the root node
    original_node_count = len(list(tree.find_clades()))
    _collapse_clade_recursively(tree.root)

    new_node_count = len(list(tree.find_clades()))
    logger.info(f"Tree after collapsing: {new_node_count} nodes (removed {original_node_count - new_node_count} nodes)")


def _collapse_clade_recursively(clade: Clade) -> bool:
    """
    Recursively process a clade and its children to collapse unbranching nodes.
    Returns True if this clade should be removed.

    :param clade: BioPython Clade object
    :return: Whether this clade should be removed
    """
    # Process all children first (bottom-up approach)
    i = 0
    while i < len(clade.clades):
        child = clade.clades[i]
        should_remove = _collapse_clade_recursively(child)

        if should_remove:
            # Replace the child with its own children
            clade.clades.pop(i)
            clade.clades[i:i] = child.clades
            # Don't increment i since we need to process the new children
        else:
            i += 1

    # A clade should be collapsed if it has exactly one child
    # We don't collapse the root, even if it has one child
    return len(clade.clades) == 1 and clade.name != "root"


def write_tree_to_newick_old(tree: Tree, output_file: str) -> None:
    """
    Write the tree to a Newick file with interior labels.

    :param tree: BioPython Tree object
    :param output_file: Path to the output Newick file
    """
    logger.info(f"Writing tree to {output_file}")

    # Remove branch lengths from the tree: they're all 0.00000 and meaningless
    for clade in tree.find_clades():
        clade.branch_length = None

    try:
        Phylo.write(tree, output_file, "newick", branch_length_only=False)
        logger.info("Tree successfully written to Newick file")
    except Exception as e:
        logger.error(f"Error writing tree to file: {e}")
        sys.exit(1)


def write_tree_to_newick(tree: Tree, output_file: str, include_internal_labels: bool = False) -> None:
    """
    Write the tree to a Newick file with optional interior labels and no branch lengths.

    :param tree: BioPython Tree object
    :param output_file: Path to the output Newick file
    :param include_internal_labels: Whether to include labels for internal nodes
    """
    logger.info(f"Writing tree to {output_file}")

    # Create a copy of the tree to avoid modifying the original
    import copy
    tree_copy = copy.deepcopy(tree)

    # Set all branch lengths to None
    for clade in tree_copy.find_clades():
        clade.branch_length = None

        # Optionally remove internal node labels
        if not include_internal_labels and clade.clades:  # If it's an internal node
            clade.name = ""

    try:
        # First write to a temporary file
        temp_file = f"{output_file}.tmp"
        Phylo.write(tree_copy, temp_file, "newick")

        # Read the file and remove branch lengths
        with open(temp_file, 'r') as f:
            newick_str = f.read().strip()

        # Remove all branch lengths (patterns like :0.0 or :0.000)
        import re
        cleaned_newick = re.sub(r':[0-9.]+', '', newick_str)

        # Write the cleaned Newick string to the output file
        with open(output_file, 'w') as f:
            f.write(cleaned_newick)

        # Clean up temporary file
        os.remove(temp_file)

        logger.info("Tree successfully written to Newick file")
    except Exception as e:
        logger.error(f"Error writing tree to file: {e}")
        sys.exit(1)

def remove_internal_labels(tree: Tree) -> None:
    """
    Remove labels from all internal nodes, keeping only leaf labels.

    :param tree: BioPython Tree object
    """
    logger.info("Removing internal node labels")

    for clade in tree.find_clades():
        # If this is an internal node (has children), remove its name
        if clade.clades:
            clade.name = ""

    logger.info("Internal node labels removed")


if __name__ == "__main__":
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Build a taxonomic tree from BOLD process IDs in a Phylip file.")
    parser.add_argument("-f", "--fasta", required=True, help="Aligned FASTA file with process IDs")
    parser.add_argument("-b", "--bold", required=True, help="BOLD BCDM TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output Newick tree file")
    parser.add_argument("-c", "--collapse", action="store_true", help="Collapse unbranching internal nodes")
    parser.add_argument("-n", "--nodelabels", action="store_true", help="Provide internal node labels")
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    args = parser.parse_args()

    # Check if files exist
    if not os.path.isfile(args.fasta):
        parser.error(f"FASTA file not found: {args.fasta}")
    if not os.path.isfile(args.bold):
        parser.error(f"BOLD TSV file not found: {args.bold}")

    # Set logging level based on verbosity
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Extract process IDs from FASTA file
    process_ids = extract_process_ids(args.fasta)

    # Read BOLD taxonomy for these process IDs
    taxonomy_df = read_bold_taxonomy(args.bold, process_ids)

    # Build taxonomy paths
    taxonomy_paths = build_taxonomy_paths(taxonomy_df)

    # Build the tree
    tree = build_tree_from_paths(taxonomy_paths)

    # Optionally collapse unbranching internal nodes
    if args.collapse:
        collapse_unbranching_nodes(tree)
        # After collapsing, fix any root unifurcation
        fix_root_unifurcation(tree)

    # Write the tree to a Newick file
    write_tree_to_newick(tree, args.output, args.nodelabels)

    logger.info("Process completed successfully")
