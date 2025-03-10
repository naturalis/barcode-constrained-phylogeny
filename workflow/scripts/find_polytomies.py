#!/usr/bin/env python3
"""
Find Polytomies in Phylogenetic Trees

This script identifies polytomies (nodes with >2 children) in a phylogenetic tree,
often used as a preprocessing step for tree refinement. Polytomies are identified 
and filtered based on their depth in the tree and number of descendant tips.

Input:
    - Newick format phylogenetic tree file

Output:
    - Text file with a list of polytomies, including:
      * Tips contained in each polytomy
      * Depth of the polytomy in the tree
      * Size (number of tips)
      * Taxonomic name associated with the polytomy

Usage:
    python find_polytomies.py -t tree_file.tre [options] > polytomies.txt

Dependencies:
    - ETE3 (Python package for tree manipulation)
"""

import argparse
from ete3 import Tree
import re

def get_node_depth(node):
    """
    Get the depth of a node in the tree (distance from root)
    
    Args:
        node: An ETE3 TreeNode object
    
    Returns:
        int: Depth of the node (0 for root, increases with distance from root)
    """
    depth = 0
    current = node
    while current.up:
        depth += 1
        current = current.up
    return depth

def find_taxonomic_name(node):
    """
    Extract taxonomic name from node or its parent nodes
    
    This function tries to find a taxonomic name (not formatted as a barcode ID)
    associated with the node or its ancestors.
    
    Args:
        node: An ETE3 TreeNode object
    
    Returns:
        str: Taxonomic name if found, None otherwise
    """
    if node.name and node.name.strip() and not re.match(r'^[A-Z]+\d+-\d+$', node.name):
        return node.name
    
    # If no name at this node, try parent nodes until we find a taxonomic name
    parent = node.up
    while parent:
        if parent.name and parent.name.strip() and not re.match(r'^[A-Z]+\d+-\d+$', parent.name):
            return parent.name
        parent = parent.up
    
    return None

def find_polytomies(tree_file, min_depth=0, max_depth=None, min_size=3, max_size=None):
    """
    Identify and report polytomies in a phylogenetic tree
    
    Args:
        tree_file: Path to Newick format tree file
        min_depth: Minimum depth to consider (filter shallower polytomies)
        max_depth: Maximum depth to consider (filter deeper polytomies)
        min_size: Minimum number of tips in a polytomy
        max_size: Maximum number of tips in a polytomy
    
    Returns:
        None: Results are printed to stdout
    """
    # Load tree from the input file
    tree = Tree(tree_file, format=1)  # Adjust format if needed

    # Find polytomies (nodes with more than 2 children)
    all_polytomies = [node for node in tree.traverse() if len(node.children) > 2]
    
    # Filter polytomies by depth and size
    filtered_polytomies = []
    for node in all_polytomies:
        depth = get_node_depth(node)
        size = len(node.get_leaf_names())
        taxon = find_taxonomic_name(node)
        
        if (min_depth <= depth and (max_depth is None or depth <= max_depth) and
            size >= min_size and (max_size is None or size <= max_size)):
            filtered_polytomies.append((node, depth, size, taxon))
    
    print(f"Found {len(filtered_polytomies)} polytomies (filtered from {len(all_polytomies)} total)")

    # Extract tip labels for each polytomy
    for i, (node, depth, size, taxon) in enumerate(filtered_polytomies):
        tips = node.get_leaf_names()
        print(f"Polytomy {i+1}: {tips} (depth={depth}, size={size}, taxon={taxon})")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find polytomies in a constraint tree.")
    parser.add_argument("-t", "--tree", required=True, 
                        help="Path to the constraint tree file (Newick format)")
    parser.add_argument("--min-depth", type=int, default=0, 
                        help="Minimum depth for polytomies to include (default: 0)")
    parser.add_argument("--max-depth", type=int, default=None, 
                        help="Maximum depth for polytomies to include (default: None)")
    parser.add_argument("--min-size", type=int, default=3, 
                        help="Minimum number of tips in polytomies (default: 3)")
    parser.add_argument("--max-size", type=int, default=None, 
                        help="Maximum number of tips in polytomies (default: None)")
    
    args = parser.parse_args()
    
    find_polytomies(args.tree, args.min_depth, args.max_depth, args.min_size, args.max_size)

"""
Example usage:

# Find all polytomies with at least 3 tips
python find_polytomies.py -t constraint_tree.tre > all_polytomies.txt

# Find polytomies with 3-50 tips at taxonomic level depth 5 or deeper
python find_polytomies.py -t constraint_tree.tre --min-depth 5 --min-size 3 --max-size 50 > filtered_polytomies.txt

"""