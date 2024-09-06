import dendropy
from typing import Dict
import bactria.logger
import bactria.config
import statistics


class OpenTOLTree(dendropy.Tree):
    """
    A subclass of dendropy.Tree specifically for handling Open Tree of Life (OpenTOL) trees.

    This class incorporates methods for processing OpenTOL trees, including modifying labels,
    handling broken tips, and removing interior labels.
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize the OpenTOLTree.
        :param config: bactria.config.Config object
        :param args: Positional arguments to pass to the superclass constructor
        :param kwargs: Keyword arguments to pass to the superclass constructor
        """
        super().__init__(*args, **kwargs)
        self.logger = bactria.logger.get_formatted_logger(__name__)

    @classmethod
    def get(cls, *args, **kwargs):
        """
        Override the get method to return an instance of OpenTOLTree instead of dendropy.Tree.

        :param args: Positional arguments to pass to dendropy.Tree.get()
        :param kwargs: Keyword arguments to pass to dendropy.Tree.get()
        :return: An instance of OpenTOLTree
        """
        config = kwargs.pop('config', None)
        tree = super().get(*args, **kwargs)
        return cls.from_tree(tree)

    @classmethod
    def from_tree(cls, tree: dendropy.Tree) -> 'OpenTOLTree':
        """
        Create an OpenTOLTree instance from an existing dendropy.Tree.

        :param tree: An existing dendropy.Tree instance
        :return: An OpenTOLTree instance with the same structure and data as the input tree
        """
        new_tree = cls(seed_node=tree.seed_node)
        new_tree._clone_from(tree, {})
        return new_tree

    def process_tree(self, broken: Dict[str, str]) -> None:
        """
        Process the tree by modifying labels and handling broken tips.

        :param broken: Dictionary of broken tips and their MRCA placeholders
        """
        self.logger.info("Processing OpenTOL tree")
        self.simplify_tip_labels()
        self.alias_broken_tips(broken)
        self.suppress_unifurcations()
        self.remove_interior_labels()

    def simplify_tip_labels(self) -> None:
        """
        Modify the labels of the tree to keep only the OTT or MRCA ID.
        """
        self.logger.debug('Modifying tree labels')
        for node in self.preorder_node_iter():
            if node.is_leaf():
                node.taxon.label = node.taxon.label.split()[-1]
            else:
                node.label = node.label.split()[-1] if node.label else None

    def alias_broken_tips(self, broken: Dict[str, str]) -> None:
        """
        Handle broken tips by grafting new leaves onto their MRCA placeholders.

        :param broken: Dictionary of broken tips and their MRCA placeholders
        """
        for ott_id, mrca_node in broken.items():
            target_node = self.find_node_with_label(mrca_node) or self.find_node_with_taxon_label(mrca_node)
            if target_node:
                new_taxon = dendropy.Taxon(label=ott_id)
                new_node = target_node.new_child(taxon=new_taxon)
                self.logger.debug(f'Grafted leaf {ott_id} as child of {mrca_node}')

    def remove_interior_labels(self) -> None:
        """
        Remove labels from interior nodes of the tree.
        """
        self.logger.debug('Removing interior node labels')
        for node in self.preorder_node_iter():
            if not node.is_leaf():
                node.label = None

    def remap_tips(self, pidmap):
        """
        Remaps the leaf labels of the input tree to the values provided in the pidmap. Possibly grafts additional child
        nodes if one-to-many mapping occurs.
        :param pidmap: a dictionary of lists
        :return:
        """
        self.logger.info('Going to remap backbone tree')
        self.logger.debug(pidmap)
        tips = self.leaf_nodes()
        for node in tips:
            name = node.taxon.label
            if str(name).startswith('ott'):
                if len(pidmap[name]) == 1:

                    # map ott to pid
                    node.taxon.label = pidmap[name][0]
                else:

                    # Iterate over all pids subtended by the ott
                    for process in pidmap[name]:
                        # Create a new dendropy taxon and node, associate them, append to parent
                        new_taxon = dendropy.Taxon(label=process)
                        new_node = dendropy.Node()
                        new_node.taxon = new_taxon
                        node.add_child(new_node)
                        self.logger.debug(f'Added child {process} to {name}')

    def pick_exemplars(self, strategy):
        """
        Picks two exemplar tips, either the two tallest, the two shallowest, or the two closest to the median
        :param strategy: How to pick exemplars
        :return:
        """
        self.logger.debug(f'Going to pick exemplars following {strategy} strategy')

        # List of exemplars to return
        representatives = []

        # Get root and its immediate children
        root = self.seed_node
        children = root.child_nodes()
        if len(children) == 2:
            for child in children:
                dists = []

                # Calculate distances to root from all leaves
                for tip in child.leaf_nodes():
                    dist = tip.distance_from_root()
                    dists.append({'dist': dist, 'id': tip.taxon.label})

                # Get shallowest tip
                dists = sorted(dists, key=lambda x: x['dist'])
                if str(strategy).capitalize().startswith('S'):
                    representatives.append(dists[0]['id'])

                # Get tallest tip
                elif str(strategy).capitalize().startswith('T'):
                    representatives.append(dists[-1]['id'])

                # Get median tip
                else:
                    median_dist = statistics.median([d['dist'] for d in dists])
                    closest = min(dists, key=lambda x: abs(x['dist'] - median_dist))
                    representatives.append(closest['id'])

        else:
            self.logger.warning(f'Ingroup root is not bifurcating, will approximate rooting')
            return None
        self.logger.debug(representatives)
        return representatives

    def get_pseudo_outgroup(self):
        """
        Given the provided rooted input tree, finds the smallest clade subtended by the root and returns for it
        the set of leaf labels of its descendants.
        :return: a set of leaf labels
        """
        self.logger.debug(f"Going to find pseudo outgroup in input tree")
        root = self.seed_node

        # Check if tree is bifurcating at the root
        if len(root.child_nodes()) != 2:
            self.logger.error("The tree is not bifurcating at the root.")

        # Find the smallest clade immediately subtended by the root
        smallest_clade = min(root.child_nodes(), key=lambda x: len(x.leaf_nodes()))

        # Return the leaves subtended by this clade
        return set([leaf.taxon.label for leaf in smallest_clade.leaf_nodes()])

    def find_set_bipartition(self, query_set):
        """
        Given an input set of leaf labels, finds the smallest edge bipartition of which the input
        is a subset.
        :param query_set: a set of leaf labels
        :return:
        """
        self.logger.debug(f"Going to find bipartition for provided input set")
        self.logger.debug(query_set)
        self.update_bipartitions()

        # Identify the smallest split that monophylizes tip_set
        smallest_split = None
        for edge in self.postorder_edge_iter():
            if edge.bipartition.leafset_bitmask is not None:
                leaves = set([taxon.label for taxon in edge.bipartition.leafset_taxa(self.taxon_namespace)])
                self.logger.debug(f'Processing bipartition {leaves}')
                if query_set.issubset(leaves):
                    self.logger.debug(f'Found smallest bipartition that has monophyletic {query_set}')
                    smallest_split = edge
                    break
        return smallest_split
