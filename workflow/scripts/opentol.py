import util
import dendropy
import requests

logger = util.get_formatted_logger(__name__, 'INFO')


def _modify_tree_labels(tree):
    """
    Processes the tip and node labels. Removes everything but the
    ott or mrca ID.
    :param tree: a dendropy tree
    :return: None
    """

    # Iterate over all nodes in the tree
    for node in tree.preorder_node_iter():

        # Check if it is a terminal node (leaf)
        if node.is_leaf():

            # Modify tip labels via taxon to keep only the 'ott+number' or mrca part
            label_parts = node.taxon.label.split(' ')
            ott_part = label_parts[-1]
            node.taxon.label = ott_part

        else:
            # Modify tip labels to keep only the 'ott+number' or mrca part
            label_parts = node.label.split(' ')
            ott_part = label_parts[-1]
            node.label = ott_part


def _iterate_requests(ids):
    """
    Iterates through OpenToL API requests, whittling down the list of unkown OTT IDs
    :param ids: list of OTT IDs
    :return: dendropy tree object
    """

    # Get JSON and process it
    result = _opentol_request(ids)

    # If it has unknowns in it, run it again
    while 'unknown' in result and len(ids) != 0:
        unknown_ott_ids = [int(item.removeprefix('ott')) for item in result['unknown'].keys()]
        ids = [item for item in ids if item not in unknown_ott_ids]
        result = _opentol_request(ids)
    logger.debug(f'OpenToL returned a useable result: {result}')

    # Read tree from string with dendropy
    return result


def _opentol_request(ids):
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
    headers = {"Content-Type": "application/json"}

    # The data to be sent with the request, as a Python dictionary
    data = {"ott_ids": ids}
    logger.debug(f'Requested OpenToL IDs serialized in JSON as: {data}')

    # Make the POST request
    response = requests.post(url, json=data, headers=headers)

    # Check if the request was successful
    if response.status_code == 200:
        json_response = response.json()
        logger.debug(f'Successful request. Response serialized as: {json_response}')
        return json_response
    else:
        json_response = response.json()
        logger.debug(f'Request resulted in (perhaps recoverable) errors. Response serialized as: {json_response}')
        return json_response


def get_subtree(idmap):
    json_result = _iterate_requests(idmap)

    # Parse the newick string, which may still have mrca nodes
    tree_obj = dendropy.Tree.get(
        data=json_result['newick'],
        schema="newick"
    )

    # Clean up all the resolved tips, which have taxon names (and underscores) in them
    _modify_tree_labels(tree_obj)

    # Check if there are 'broken' tips
    for ott_id in dict(json_result['broken']).keys():

        # Lookup the mrca placeholder and find it in the tree
        mrca_node = dict(json_result['broken'])[ott_id]
        target_node = tree_obj.find_node_with_label(mrca_node)

        # Might have to look for taxon instead
        if target_node is None:
            target_node = tree_obj.find_node_with_taxon_label(mrca_node)

        # Create a new dendropy taxon and node, associate them,
        new_taxon = dendropy.Taxon(label=ott_id)
        new_node = dendropy.Node()
        new_node.taxon = new_taxon
        target_node.add_child(new_node)

        logger.debug(f'Grafted leaf {ott_id} as child of {mrca_node}')

    # Remove unbranched internal nodes
    tree_obj.suppress_unifurcations()

    # Remove interior node labels
    for node in tree_obj.preorder_node_iter():
        if not node.is_leaf():
            node.label = None

    return tree_obj
