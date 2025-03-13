#!/usr/bin/env python3
"""
Resolve Polytomies Using OpenTree of Life API

This script resolves polytomies (multifurcating nodes) in a phylogenetic tree by 
replacing them with better-resolved subtrees from the OpenTree of Life database.
The process involves:
1. Identifying polytomies in the input tree
2. Finding the taxonomic name associated with each polytomy
3. Querying the OpenToL API to get the taxonomic ID
4. Retrieving a resolved subtree for that taxon
5. Replacing the polytomy with the resolved subtree

Input:
    - Constraint tree (Newick format)
    - List of polytomies (.txt file, see output from find_polytomies.py)

Output:
    - Resolved constraint tree with improved resolution at polytomies

Usage:
    python resolve_polytomies.py -t tree.tre -p polytomies.txt -o resolved_tree.tre [options]

Dependencies:
    - ETE3 (Python tree manipulation)
    - Requests (API communication)
    - Internet connection (for OpenToL API access)
"""

import argparse
import re
import ast
from ete3 import Tree
import sys
import logging
import requests
import json
import tempfile
import os

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('resolve_polytomies')

def parse_polytomies_file(filename):
    """
    Parse the output file from find_polytomies.py
    
    Extracts lists of tip labels for each polytomy from the output format
    of find_polytomies.py.
    
    Args:
        filename: Path to polytomies file
        
    Returns:
        list: Lists of tip labels for each polytomy
    """
    polytomies = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    # Skip the first line (count summary)
    for i, line in enumerate(lines[1:], 1):
        # Extract the list part using regex with non-greedy matching to avoid capturing extra content
        match = re.search(r'Polytomy \d+: (\[.*?\]) \(depth=', line)
        if match:
            # Convert string representation of list to actual list
            try:
                tips = ast.literal_eval(match.group(1))
                polytomies.append(tips)
            except (SyntaxError, ValueError):
                logger.warning(f"Could not parse polytomy on line {i+1}")
    
    return polytomies

def find_node_with_tips(tree, tips):
    """
    Find the node in the tree that contains exactly these tips
    
    Args:
        tree: ETE3 Tree object
        tips: List of tip labels to find
        
    Returns:
        TreeNode: Node containing exactly these tips, or None if not found
    """
    for node in tree.traverse():
        if set(node.get_leaf_names()) == set(tips):
            return node
    return None

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

def get_ott_id(taxon_name):
    """
    Query OpenTOL TNRS API to get OTT ID for a taxonomic name
    
    Uses the Taxonomic Name Resolution Service to find the OpenTree Taxonomy ID
    for a given taxonomic name.
    
    Args:
        taxon_name: Taxonomic name to look up
        
    Returns:
        int: OpenTree Taxonomy ID if found, None otherwise
    """
    url = "https://api.opentreeoflife.org/v3/tnrs/match_names"
    payload = {
        "names": [taxon_name],
        "do_approximate_matching": True
    }
    
    try:
        logger.info(f"Looking up OTT ID for '{taxon_name}'")
        response = requests.post(url, json=payload)
        data = response.json()
        
        if 'results' in data and len(data['results']) > 0 and len(data['results'][0]['matches']) > 0:
            ott_id = data['results'][0]['matches'][0]['taxon']['ott_id']
            logger.info(f"Found OTT ID for '{taxon_name}': {ott_id}")
            return ott_id
        else:
            logger.warning(f"No OTT ID found for '{taxon_name}'")
            return None
    except Exception as e:
        logger.error(f"Error querying OTT ID: {e}")
        return None

def get_resolved_subtree(ott_id):
    """
    Get a resolved subtree from OpenTOL using an OTT ID
    
    Queries the OpenTree of Life API to get a resolved subtree for a given
    taxonomic ID.
    
    Args:
        ott_id: OpenTree Taxonomy ID
        
    Returns:
        str: Newick string of resolved subtree if found, None otherwise
    """
    url = "https://api.opentreeoflife.org/v3/tree_of_life/subtree"
    payload = {
        "ott_id": ott_id,
        "label_format": "name"
    }
    
    try:
        logger.info(f"Requesting subtree for OTT ID {ott_id}")
        response = requests.post(url, json=payload)
        data = response.json()
        
        if 'newick' in data:
            return data['newick']
        else:
            logger.warning(f"No subtree found for OTT ID {ott_id}")
            return None
    except Exception as e:
        logger.error(f"Error getting subtree: {e}")
        return None

def parse_newick_safely(newick_str):
    """
    Parse a Newick string safely with multiple methods
    
    Tries multiple approaches to parse potentially problematic Newick strings,
    including handling quoted names, special characters, and various formats.
    
    Args:
        newick_str: Newick format string to parse
        
    Returns:
        Tree: ETE3 Tree object if parsing succeeds, None otherwise
    """
    def sanitize_newick(text):
        """Clean up problematic characters in newick strings"""
        # Replace problematic patterns
        text = re.sub(r"'([^']*)'", r"\1", text)  # Remove single quotes
        text = re.sub(r'"([^"]*)"', r"\1", text)  # Remove double quotes
        text = re.sub(r'nr\.\s+', r"nr_", text)   # Fix "nr. " pattern
        text = re.sub(r'\s+sp\.\s+', r"_sp_", text)  # Fix " sp. " pattern
        text = re.sub(r'\s+', r"_", text)         # Replace spaces with underscores
        
        # Add quotes around all taxonomic names for safety
        text = re.sub(r'([A-Za-z][A-Za-z0-9_.:-]+)', r"'\1'", text)
        return text
    
    # Try multiple parsing approaches
    try:
        # First try simple parsing with quoted names
        return Tree(newick_str, format=1, quoted_node_names=True)
    except Exception:
        try:
            # Try format 0 (most flexible)
            return Tree(newick_str, format=0)
        except Exception:
            # Last resort: write to temp file with sanitized content
            sanitized = sanitize_newick(newick_str)
            
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp:
                temp.write(sanitized)
                temp_file = temp.name
                
            try:
                # Try to parse the sanitized newick from file
                tree = Tree(temp_file, format=1, quoted_node_names=True)
                os.remove(temp_file)
                return tree
            except Exception:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                
                # If we still can't parse it, log more info and return None
                logger.debug(f"Raw newick: {newick_str[:100]}...")
                logger.debug(f"Sanitized: {sanitized[:100]}...")
                return None

def resolve_polytomy(constraint_tree, polytomy_tips, min_depth=3):
    """
    Resolve a polytomy using the OpenTOL API
    
    Takes a set of tips defining a polytomy and attempts to replace the
    corresponding node with a better-resolved subtree from OpenTree of Life.
    
    Args:
        constraint_tree: ETE3 Tree object containing the polytomy
        polytomy_tips: List of tip labels in the polytomy
        min_depth: Minimum depth to consider (skip shallower nodes)
        
    Returns:
        tuple: (
            success (bool), 
            taxon_name (str or None), 
            ott_id (int or None), 
            newick_str (str or None)
        )
    """
    # Find the node containing these tips
    polytomy_node = find_node_with_tips(constraint_tree, polytomy_tips)
    if not polytomy_node:
        logger.warning(f"Could not find node for tips in tree")
        return False, None, None, None
    
    # Skip polytomies that are too shallow in the tree
    node_depth = get_node_depth(polytomy_node)
    if node_depth < min_depth:
        logger.info(f"Skipping shallow polytomy (depth {node_depth} < {min_depth})")
        return False, None, None, None
    
    # Get taxonomic name for this node
    taxon_name = find_taxonomic_name(polytomy_node)
    if not taxon_name:
        logger.warning(f"Could not find taxonomic name for node")
        return False, None, None, None
    
    # Get OTT ID for this taxon
    ott_id = get_ott_id(taxon_name)
    if not ott_id:
        # Try parent taxon if this one isn't found
        if polytomy_node.up:
            parent_name = find_taxonomic_name(polytomy_node.up)
            if parent_name and parent_name != taxon_name:
                logger.info(f"Trying parent taxon '{parent_name}'")
                ott_id = get_ott_id(parent_name)
    
    if not ott_id:
        logger.warning(f"Could not find OTT ID for taxon '{taxon_name}'")
        return False, taxon_name, None, None
    
    # Get resolved subtree
    newick_str = get_resolved_subtree(ott_id)
    if not newick_str:
        logger.warning(f"Could not get subtree for OTT ID {ott_id}")
        return False, taxon_name, ott_id, None
    
    # Replace polytomy with resolved subtree
    try:
        # Parse the resolved subtree with proper handling of quoted names
        resolved_subtree = parse_newick_safely(newick_str)
        
        if not resolved_subtree:
            logger.warning(f"Failed to parse resolved subtree")
            return False, taxon_name, ott_id, newick_str
        
        # Check if the resolved subtree has a better structure
        if len(resolved_subtree.get_children()) <= len(polytomy_node.get_children()):
            logger.info(f"Resolved subtree doesn't improve structure, skipping")
            return False, taxon_name, ott_id, newick_str
        
        # Replace the polytomy with the resolved subtree
        parent = polytomy_node.up
        if parent:
            # Remove old node
            polytomy_node.detach()
            
            # Add new subtree (fix: removed 'pos' parameter)
            parent.add_child(resolved_subtree)
            logger.info(f"Successfully replaced polytomy with resolved subtree")
            return True, taxon_name, ott_id, newick_str
        else:
            # Special handling for root node
            logger.warning(f"Cannot replace root node")
            return False, taxon_name, ott_id, newick_str
    except Exception as e:
        logger.error(f"Error replacing polytomy: {e}")
        return False, taxon_name, ott_id, newick_str

def main():
    """
    Main function: Parse arguments and process polytomies
    """
    parser = argparse.ArgumentParser(description="Resolve polytomies using OpenTOL API")
    parser.add_argument("-t", "--tree", required=True,
                        help="Path to the original constraint tree")
    parser.add_argument("-p", "--polytomies", required=True, 
                        help="Path to polytomies file (output from find_polytomies.py)")
    parser.add_argument("-o", "--output", required=True,
                        help="Path for the resolved constraint tree output")
    parser.add_argument("-l", "--limit", type=int, default=None,
                        help="Limit the number of polytomies to process (for testing)")
    parser.add_argument("--min-size", type=int, default=3,
                        help="Minimum size of polytomy to process (default: 3)")
    parser.add_argument("--max-size", type=int, default=50,
                        help="Maximum size of polytomy to process (default: 50)")
    parser.add_argument("--min-depth", type=int, default=3,
                        help="Minimum depth in tree for polytomy to process (default: 3)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose output")
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Load the constraint tree
    try:
        logger.info(f"Loading constraint tree from {args.tree}")
        constraint_tree = Tree(args.tree, format=1)
    except Exception as e:
        logger.error(f"Failed to load tree: {e}")
        sys.exit(1)
    
    # Parse the polytomies file
    logger.info(f"Parsing polytomies from {args.polytomies}")
    polytomies = parse_polytomies_file(args.polytomies)
    logger.info(f"Found {len(polytomies)} polytomies")
    
    # Process each polytomy
    processed = 0
    resolved = 0
    no_ott_id = []
    failed_parsing = []
    other_fails = []
    
    for i, tips in enumerate(polytomies):
        # Check if we've reached the limit
        if args.limit and i >= args.limit:
            break
        
        # Skip polytomies that are too large or too small
        if len(tips) < args.min_size or len(tips) > args.max_size:
            logger.info(f"Skipping polytomy {i+1} with {len(tips)} tips (size out of range)")
            continue
        
        logger.info(f"Processing polytomy {i+1}/{len(polytomies)} with {len(tips)} tips")
        
        # Try to resolve this polytomy
        success, taxon_name, ott_id, newick_str = resolve_polytomy(constraint_tree, tips, min_depth=args.min_depth)
        processed += 1
        
        if success:
            resolved += 1
        elif taxon_name and not ott_id:
            no_ott_id.append(taxon_name)
        elif ott_id and newick_str and not success:
            failed_parsing.append(taxon_name)
        elif taxon_name:
            other_fails.append(taxon_name)
    
    # Save the resolved tree
    logger.info(f"Saving resolved tree to {args.output}")
    constraint_tree.write(format=1, outfile=args.output)
    
    # Generate a summary report
    logger.info(f"\n===== SUMMARY =====")
    logger.info(f"Total polytomies processed: {processed}")
    logger.info(f"Successfully resolved: {resolved} ({(resolved/processed*100) if processed > 0 else 0:.1f}%)")
    logger.info(f"Taxa with no OTT ID: {len(no_ott_id)}")
    logger.info(f"Taxa with parsing errors: {len(failed_parsing)}")
    logger.info(f"Taxa that failed for other reasons: {len(other_fails)}")
    logger.info(f"===================\n")

if __name__ == "__main__":
    main()

"""
Example usage:

# Basic usage
python resolve_polytomies.py -t constraint_tree.tre -p polytomies.txt -o resolved_tree.tre

# Focus on smaller polytomies at deeper taxonomic levels
python resolve_polytomies.py -t constraint_tree.tre -p polytomies.txt -o resolved_tree.tre --min-depth 5 --min-size 3 --max-size 25

# Process only a limited number of polytomies (for testing)
python resolve_polytomies.py -t constraint_tree.tre -p polytomies.txt -o resolved_tree.tre -l 10 -v

Workflow:
1. Generate polytomies file: python find_polytomies.py -t tree.tre --min-depth 4 > polytomies.txt
2. Resolve polytomies: python resolve_polytomies.py -t tree.tre -p polytomies.txt -o resolved_tree.tre
3. Use resolved tree as a constraint tree for phylogenetic analyses (RAxML, etc.)
"""