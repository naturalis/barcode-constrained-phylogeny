import requests
import logging
import argparse
import dendropy

# Function to modify tip labels and remove interior node labels
def modify_tree_labels(tree):
    """
    Processes the tip and node labels. If tip, remove everything but the
    ott ID. If node, remove the label. If an interior node has been
    elevated to a tip (which can be recognized by the mrca prefix in the
    label), prune it.
    :param tree: a dendropy tree
    :return:
    """

    # Iterate over all nodes in the tree
    for node in tree.preorder_node_iter():

        # Check if it is a terminal node (leaf)
        if node.is_leaf():

            # Modify tip labels to keep only the 'ott+number' part
            label_parts = node.taxon.label.split(' ')
            ott_part = label_parts[-1]

            # Check if MRCA and remove it.
            # TODO: this situation arises when a taxon is 'broken' according to OpenToL.
            # For example, requesting _Alouatta seniculus_ results in this because the
            # species subtends two ssp _A. seniculus_ as well as _A. sara_. The solution
            # appears to be to parse the mrca string, tease out the ott IDs, then re-request,
            # so that we can at least figure out the species names (the initial output
            # doesn't give them) and map that back to our input.
            if ott_part.startswith('mrca'):
                logger.warning(f'Tree has MRCA tip in it: {ott_part}')
                node.parent_node.remove_child(node)
            else:
                node.taxon.label = ott_part
        else:

            # Remove interior node labels
            node.label = None


def fetch_induced_subtree(ids):
    """
    Places a request to the OpenToL induced subtree web service endpoint.
    Parameterized by a list of ott IDs. This service call sometimes fails,
    when the parameter set includes IDs not in the subtree, which results
    in a warning being emitted. In that case, the return value includes
    a list of the unknown IDs. The caller can then remove these from the
    input list and try again.
    :param ids: a list of ott IDs
    :return: a JSON data structure
    """
    # The API endpoint URL
    url = "https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree"

    # The headers to indicate we are sending JSON data
    headers = {
        "Content-Type": "application/json",
    }

    # The data to be sent with the request, as a Python dictionary
    data = {
        "ott_ids": ids
    }
    logger.debug(data)

    # Make the POST request
    response = requests.post(url, json=data, headers=headers)

    # Check if the request was successful
    if response.status_code == 200:
        return response.json()
    else:
        logger.warning(f"Failed to retrieve data - will try again")
        return response.json()


def extract_id_from_fasta(file_path):
    """
    Extracts the ott IDs from a FASTA file. This operates specifically on
    FASTA headers where the ott ID is the second element (i.e. at index 1)
    in the pipe-separated FASTA header.
    :param file_path: the location of a FASTA file
    :return: a list of ott IDs
    """
    ids = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                parts = line.strip().split('|')
                if len(parts) > 1:
                    # Extract the second element and remove 'ott' prefix
                    id_with_prefix = parts[1]
                    id_number = id_with_prefix.replace('ott', '')
                    ids.append(id_number)

    # Remove all 'None' entries
    cleaned_list = [item for item in ids if item != 'None']
    return cleaned_list


def postprocess_tree(newick):
    """
    Given a newick string as input, parses the string into a dendropy tree,
    cleans up the tip labels (via modify_tree_labels) and removes unbranched
    interior nodes
    :param newick: a newick string
    :return: a dendropy tree object
    """
    logger.debug(f'Post-processing result tree {newick}')

    # Read tree from string with dendropy
    tree_obj = dendropy.Tree.get(
        data=newick,
        schema="newick"
    )

    # Clean the labels
    modify_tree_labels(tree_obj)

    # Remove unbranched internals
    tree_obj.suppress_unifurcations()

    return tree_obj


if __name__ == '__main__':
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Required command line arguments.')
    parser.add_argument('-i', '--inaln', required=True, help='Input exemplar FASTA files')
    parser.add_argument('-o', '--outtree', required=True, help="Output constraint tree")
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Configure logging
    logging.basicConfig()
    logger = logging.getLogger('backbone_constraint')
    logger.setLevel(args.verbosity)

    # Read input alignment, get ott IDs
    ott_ids = extract_id_from_fasta(args.inaln)

    # If we have no IDs at all, we write a zero byte file for run_raxml
    if len(ott_ids) == 0:
        logger.warning('There were zero OTT IDs in the input file')
        with open(args.outtree, "a"):
            pass

    else:
        # Get JSON and process it
        result = fetch_induced_subtree(ott_ids)

        # If it has unknowns in it, run it again
        if 'unknown' in result:
            unknown_ott_ids = [ item.removeprefix('ott') for item in result['unknown'].keys() ]
            logger.warning(f'Request included unkown ott IDs: {unknown_ott_ids}')
            cleaned_list = [item for item in ott_ids if item not in unknown_ott_ids]
            result = fetch_induced_subtree(cleaned_list)

        # Clean up the tip labels, remove internal node labels, remove unbranched internals
        tree = postprocess_tree(result['newick'])

        # Write output
        with open(args.outtree, "w") as output_file:
            output_file.write(tree.as_string(schema="newick"))

