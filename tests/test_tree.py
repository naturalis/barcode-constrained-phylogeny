import pytest
import dendropy
import os
from bactria.tree import OpenTOLTree
from bactria.config import Config

# Construct the path to the config file
CONFIG_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config', 'config.yaml')


@pytest.fixture
def config():
    cfg = Config()
    cfg.load_config(CONFIG_PATH)
    return cfg

@pytest.fixture
def simple_tree(config):
    newick = "(ott1:1,(ott2:1,(ott3:1,ott4:1):0.5):0.5);"
    return OpenTOLTree.get(string=newick, schema="newick", config=config)

def test_initialization(config):
    tree = OpenTOLTree(config)
    assert isinstance(tree, OpenTOLTree)
    assert isinstance(tree, dendropy.Tree)

def test_get_class_method(config):
    newick = "(A:1,B:1);"
    tree = OpenTOLTree.get(string=newick, schema="newick")
    assert isinstance(tree, OpenTOLTree)

def test_from_tree_class_method(config):
    newick = "(A:1,B:1);"
    dendro_tree = dendropy.Tree.get(string=newick, schema="newick")
    tree = OpenTOLTree.from_tree(dendro_tree)
    assert isinstance(tree, OpenTOLTree)

def test_process_tree(simple_tree):
    broken = {"ott1": "mrcaott2ott3"}
    simple_tree.process_tree(broken)
    assert all(not node.label for node in simple_tree.internal_nodes())

def test_simplify_tip_labels(simple_tree):
    simple_tree.simplify_tip_labels()
    assert all(len(leaf.taxon.label.split()) == 1 for leaf in simple_tree.leaf_nodes())

def test_alias_broken_tips(simple_tree):
    broken = {"ott1": "A"}
    simple_tree.alias_broken_tips(broken)
    assert any(leaf.taxon.label == "ott1" for leaf in simple_tree.leaf_nodes())

def test_remove_interior_labels(simple_tree):
    simple_tree.remove_interior_labels()
    assert all(not node.label for node in simple_tree.internal_nodes())

def test_remap_tips(simple_tree):
    pidmap = {"ott1": ["ott11"], "ott2": ["ott21", "ott22"], "ott3": ["ott31"], "ott4": ["ott41"]}
    simple_tree.remap_tips(pidmap)
    leaf_labels = [leaf.taxon.label for leaf in simple_tree.leaf_nodes()]
    assert set(leaf_labels) == {"ott11", "ott21", "ott22", "ott31", "ott41"}

def test_pick_exemplars_tallest(simple_tree):
    exemplars = simple_tree.pick_exemplars("tallest")
    assert len(exemplars) == 2
    assert set(exemplars) == {"ott1", "ott4"}

def test_pick_exemplars_shallowest(simple_tree):
    exemplars = simple_tree.pick_exemplars("shallowest")
    assert len(exemplars) == 2
    assert set(exemplars) == {"ott1", "ott2"}

def test_pick_exemplars_median(simple_tree):
    exemplars = simple_tree.pick_exemplars("median")
    assert len(exemplars) == 2
    assert set(exemplars) == {"ott1", "ott3"}

def test_get_pseudo_outgroup(simple_tree):
    outgroup = simple_tree.get_pseudo_outgroup()
    assert outgroup == {"ott1"}

def test_find_set_bipartition(simple_tree):
    bipartition = simple_tree.find_set_bipartition({"ott3", "ott4"})
    assert bipartition is not None
    assert set(taxon.label for taxon in bipartition.bipartition.leafset_taxa(simple_tree.taxon_namespace)) == {"ott3", "ott4"}

def test_non_bifurcating_root(config):
    newick = "(A:1,B:1,C:1);"
    tree = OpenTOLTree.get(string=newick, schema="newick", config=config)
    assert tree.pick_exemplars("tallest") is None
#    with pytest.raises(Exception):  # Adjust the exception type as needed
#        tree.get_pseudo_outgroup()