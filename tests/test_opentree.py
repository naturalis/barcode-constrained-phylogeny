import pytest
import requests
import os
from unittest.mock import Mock, patch
from bactria.opentree import OpenTOLClient
from bactria.config import Config
from bactria.tree import OpenTOLTree


# Construct the path to the config file
CONFIG_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config', 'config.yaml')


@pytest.fixture
def config():
    cfg = Config()
    cfg.load_config(CONFIG_PATH)
    return cfg


@pytest.fixture
def client(config):
    return OpenTOLClient(config)


def test_initialization(client):
    assert isinstance(client, OpenTOLClient)
    assert client.base_url == "https://api.opentreeoflife.org/v3"
    assert client.headers == {"Content-Type": "application/json"}


def test_endpoint_property(client):
    assert client.endpoint == "https://api.opentreeoflife.org/v3/tnrs/match_names"
    client.endpoint = "new/endpoint"
    assert client.endpoint == "https://api.opentreeoflife.org/v3/new/endpoint"


@patch('requests.get')
@patch('requests.post')
def test_handle_request(mock_post, mock_get, client):
    mock_response = Mock()
    mock_response.json.return_value = {"test": "data"}
    mock_get.return_value = mock_response
    mock_post.return_value = mock_response

    # Test GET request
    result = client.handle_request('GET')
    assert result == {"test": "data"}
    mock_get.assert_called_once()

    # Test POST request
    data = {"key": "value"}
    result = client.handle_request('POST', data)
    assert result == {"test": "data"}
    mock_post.assert_called_once_with(client.endpoint, json=data, headers=client.headers)

    # Test unsupported method
    with pytest.raises(ValueError):
        client.handle_request('PUT')


@patch('requests.post')
def test_tnrs_match(mock_post, client):
    mock_response = Mock()
    mock_response.json.return_value = {"matches": [{"name": "Homo sapiens"}]}
    mock_post.return_value = mock_response

    result = client.tnrs_match(["Homo sapiens"], "Animals")
    assert result == {"matches": [{"name": "Homo sapiens"}]}
    mock_post.assert_called_once()


@patch('requests.post')
@patch('bactria.tree.OpenTOLTree.get')
def test_get_induced_subtree(mock_tree_get, mock_post, client):
    mock_response = Mock()
    mock_response.json.return_value = {
        "newick": "(A:1,B:1);",
        "broken": {"ott1": "mrca2"}
    }
    mock_post.return_value = mock_response

    mock_tree = Mock(spec=OpenTOLTree)
    mock_tree_get.return_value = mock_tree

    result = client.get_induced_subtree([1, 2, 3])
    assert result == mock_tree
    mock_post.assert_called_once()
    mock_tree_get.assert_called_once_with(data="(A:1,B:1);", schema="newick")
    mock_tree.process_tree.assert_called_once_with({"ott1": "mrca2"})


@patch('requests.post')
def test_handle_request_exception(mock_post, client):
    mock_post.side_effect = requests.RequestException("API Error")

    with pytest.raises(requests.RequestException):
        client.handle_request('POST', {"test": "data"})