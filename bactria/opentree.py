import requests
import bactria.logger
import bactria.config
from typing import List, Dict, Any
from bactria.tree import OpenTOLTree


class OpenTOLClient:
    """
    A client for interacting with the Open Tree of Life (OpenTOL) API.

    This class provides methods to interact with different endpoints of the OpenTOL API,
    including the TNRS (Taxonomic Name Resolution Service) and the induced subtree service.
    """

    def __init__(self, config: bactria.config.Config):
        """
        Initialize the OpenTOLClient.

        :param log_level: The logging level (default: 'INFO')
        """
        self.logger = bactria.logger.get_formatted_logger(__name__, config)
        self.base_url = "https://api.opentreeoflife.org/v3"
        self.endpoint = f"{self.base_url}/tnrs/match_names"
        self.headers = {"Content-Type": "application/json"}

    @property
    def endpoint(self) -> str:
        """Get the current API endpoint."""
        return self.endpoint

    @endpoint.setter
    def endpoint(self, new_endpoint: str):
        """
        Set a new API endpoint.

        :param new_endpoint: The new endpoint to use
        """
        self.endpoint = f"{self.base_url}/{new_endpoint}"
        self.logger.info(f"Endpoint changed to: {self.endpoint}")

    def handle_request(self, method: str, data: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Handle API requests and return parsed JSON response.

        :param method: HTTP method ('GET' or 'POST')
        :param data: Data to send with the request (for POST)
        :return: Parsed JSON response
        :raises requests.RequestException: If the request fails
        """
        try:
            if method.upper() == 'GET':
                response = requests.get(self.endpoint, headers=self.headers)
            elif method.upper() == 'POST':
                response = requests.post(self.endpoint, json=data, headers=self.headers)
            else:
                raise ValueError(f"Unsupported HTTP method: {method}")

            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            self.logger.error(f"Request failed: {str(e)}")
            raise

    def tnrs_match(self, names: List[str], kingdom: str, fuzzy: bool = False) -> Dict[str, Any]:
        """
        Perform a TNRS match for the given names.

        :param names: List of taxonomic names to match
        :param kingdom: Taxonomic context ('Animals' or 'Plants')
        :param fuzzy: Whether to use fuzzy matching (default: False)
        :return: TNRS match results
        """
        self.endpoint = "tnrs/match_names"
        data = {
            "names": names,
            "do_approximate_matching": fuzzy,
            "context": kingdom
        }
        self.logger.info(f"Performing TNRS match for {len(names)} names")
        return self.handle_request('POST', data)

    def get_induced_subtree(self, ott_ids: List[int]) -> OpenTOLTree:
        """
        Retrieve the induced subtree for the given OpenTOL IDs.

        :param ott_ids: List of OpenTOL IDs
        :return: OpenTOLTree object representing the induced subtree
        """
        self.endpoint = "tree_of_life/induced_subtree"
        data = {"ott_ids": ott_ids}
        self.logger.info(f"Retrieving induced subtree for {len(ott_ids)} OTT IDs")

        while True:
            result = self.handle_request('POST', data)
            if 'unknown' not in result or not ott_ids:
                break
            unknown_ott_ids = [int(item.removeprefix('ott')) for item in result['unknown'].keys()]
            ott_ids = [item for item in ott_ids if item not in unknown_ott_ids]
            data["ott_ids"] = ott_ids

        tree = OpenTOLTree.get(data=result['newick'], schema="newick")
        tree.process_tree(result.get('broken', {}))
        return tree
