#!/usr/bin/env python3
"""
Format Constraint Tree with OTT IDs and Bifurcation

Converts a taxonomic-labeled constraint tree to a fully bifurcating tree
with OTT IDs at internal nodes and barcode IDs at tips.

Input:
    - Constraint tree with taxonomic names at internal nodes
    
Output:
    - Fully bifurcating constraint tree with OTT IDs at internal nodes and no branch lengths
"""

import argparse
import re
from ete3 import Tree
import logging
import requests
import os
import sys

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('format_constraint_tree')

# Cache for OTT IDs to avoid redundant API calls
ott_cache = {}

def get_ott_id(taxon_name):
    """Get OTT ID for a taxonomic name using OpenTree API"""
    if taxon_name in ott_cache:
        return ott_cache[taxon_name]
        
    # Clean up taxon name - replace underscores with spaces
    clean_name = taxon_name.replace("_", " ")
    
    url = "https://api.opentreeoflife.org/v3/tnrs/match_names"
    payload = {
        "names": [clean_name],
        "do_approximate_matching": True
    }
    
    try:
        logger.debug(f"Looking up OTT ID for '{clean_name}'")
        response = requests.post(url, json=payload)
        data = response.json()
        
        if 'results' in data and len(data['results']) > 0 and len(data['results'][0]['matches']) > 0:
            ott_id = data['results'][0]['matches'][0]['taxon']['ott_id']
            ott_cache[taxon_name] = f"ott{ott_id}"
            return f"ott{ott_id}"
        else:
            logger.warning(f"No OTT ID found for '{clean_name}'")
            ott_cache[taxon_name] = None
            return None
    except Exception as e:
        logger.error(f"Error querying OTT ID: {e}")
        return None

def is_barcode_id(name):
    """Check if a name looks like a barcode ID"""
    # Most barcode IDs follow patterns like ABCDE123-45
    return bool(re.match(r'^[A-Z]+\d+-\d+$', name))

def format_tree(input_tree, output_tree, force_bifurcate=True):
    """
    Format tree by replacing taxonomic names with OTT IDs and ensuring bifurcation
    
    Args:
        input_tree: Path to input tree file
        output_tree: Path to output tree file
        force_bifurcate: Whether to force the tree to be fully bifurcating
    """
    # Load tree
    tree = Tree(input_tree, format=1)
    logger.info(f"Loaded tree with {len(tree)} tips")
    
    # Process internal nodes
    internal_nodes = 0
    resolved_nodes = 0
    unnamed_count = 0
    
    # First pass: get OTT IDs for all taxonomic names
    for node in tree.traverse():
        # Skip tips - keep their barcode IDs
        if node.is_leaf():
            continue
            
        if node.name and node.name.strip():
            internal_nodes += 1
            
            # Skip nodes that already have OTT format
            if node.name.startswith('ott'):
                resolved_nodes += 1
                continue
                
            # Get OTT ID for this taxonomic name
            ott_id = get_ott_id(node.name)
            if ott_id:
                node.name = ott_id
                resolved_nodes += 1
            else:
                # If no OTT ID found, create a placeholder
                unnamed_count += 1
                node.name = f"unnamed{unnamed_count}"
    
    # Second pass: force bifurcation if requested
    if force_bifurcate:
        polytomies = 0
        resolved_polytomies = 0
        
        # Find all polytomies
        for node in tree.traverse():
            if len(node.children) > 2:
                polytomies += 1
                
                # Sort children by the number of descendants (smallest first)
                node.children.sort(key=lambda n: len(n.get_leaves()))
                
                # Create a ladder-like structure with the sorted children
                while len(node.children) > 2:
                    # Take the two smallest children
                    child1 = node.children[0]
                    child2 = node.children[1]
                    
                    # Remove them from the node
                    node.remove_child(child1)
                    node.remove_child(child2)
                    
                    # Create a new internal node to hold them
                    new_node = Tree()
                    new_node.name = f"unnamed{unnamed_count}"
                    unnamed_count += 1
                    
                    # Add the two children to the new node
                    new_node.add_child(child1)
                    new_node.add_child(child2)
                    
                    # Add the new node back to the original node
                    node.add_child(new_node)
                
                resolved_polytomies += 1
        
        logger.info(f"Resolved {resolved_polytomies} of {polytomies} polytomies")
    
    # Remove branch lengths
    for node in tree.traverse():
        node.dist = 0
    
    # Write the tree in Newick format
    newick = tree.write(format=9)  # format 9 is Newick with internal node names
    
    # Write to file
    with open(output_tree, 'w') as f:
        f.write(newick)
    
    # Verify the tree structure
    bifurcating = all(len(node.children) <= 2 for node in tree.traverse() if not node.is_leaf())
    
    logger.info(f"Processed {internal_nodes} internal nodes")
    logger.info(f"Found OTT IDs for {resolved_nodes} nodes ({resolved_nodes/internal_nodes*100:.1f}% if internal_nodes else 0)")
    logger.info(f"Created {unnamed_count} unnamed internal nodes")
    logger.info(f"Tree is{'fully' if bifurcating else 'NOT'} bifurcating")
    logger.info(f"Wrote formatted tree to {output_tree}")

def main():
    parser = argparse.ArgumentParser(description="Format constraint tree with OTT IDs and ensure bifurcation")
    parser.add_argument("-i", "--input", required=True,
                        help="Path to input tree file")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to output tree file")
    parser.add_argument("--no-bifurcate", action="store_true",
                        help="Don't force bifurcation (keep polytomies)")
    parser.add_argument("--batch-size", type=int, default=100,
                        help="Batch size for OTT ID lookups (default: 100)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose output")
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    format_tree(args.input, args.output, not args.no_bifurcate)

if __name__ == "__main__":
    main()