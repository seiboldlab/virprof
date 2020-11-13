"""This module contains classes for handling operations on the NCBI taxonomy"""

import os
import csv
import logging
from abc import ABC, abstractmethod
from collections import OrderedDict
from typing import Set, Sequence, Any, Iterator, Dict, Tuple, Callable, Optional

import graph_tool as gt  # type: ignore
import graph_tool.topology as gt_topology  # type: ignore
import graph_tool.search as gt_search  # type: ignore
import graph_tool.util as gt_util  # type: ignore
import networkx as nx  # type: ignore

LOG = logging.getLogger(__name__)


class Taxonomy(ABC):
    """Base class holding NCBI taxonomy

    The tree is implemented in subclasses TaxonomyGT and
    TaxonomyNX, using graph-tool and networkx libraries,
    respectively.

    Args:
      path: Either path to directory containing names.dmp and nodes.dmp or
            path to the binary dump of the tree.
    """
    #: Name of file containing nodes
    NODES_FN = "nodes.dmp"
    #: Name of file containing names
    NAMES_FN = "names.dmp"

    def __init__(self, path: Optional[str]) -> None:
        self.path = path

        if path is not None:
            if os.path.isdir(path):
                self.tree = self.load_tree()
            elif os.path.exists(path):
                self.tree = self.load_tree_binary()
            else:
                self.tree = self.load_tree()

    def _find_file(self, filename: str) -> Optional[str]:
        if self.path is None:
            return None
        path = os.path.join(self.path, filename)
        if os.path.exists(path):
            return path
        path = "".join((self.path, filename))
        if os.path.exists(path):
            return path
        raise RuntimeError("Unable to find {filename} with prefix {self.path}")


    @property
    def names_fn(self) -> Optional[str]:
        """Path to the names.dmp file"""
        return self._find_file(self.NAMES_FN)

    @property
    def nodes_fn(self) -> Optional[str]:
        """Path to the nodes.dmp file"""
        return self._find_file(self.NODES_FN)

    @abstractmethod
    def load_tree_binary(self) -> Any:
        """Load the tree from implementation specific binary"""
        raise NotImplementedError()

    @abstractmethod
    def save_tree_binary(self, path: str) -> None:
        """Save the tree to implementation specific binary"""
        raise NotImplementedError()

    @abstractmethod
    def load_tree(self) -> Any:
        """Load the tree from NCBI dmp files"""
        raise NotImplementedError()

    @property
    def root(self) -> int:
        """The ``tax_id`` of the root node

        Returns:
          Always 1
        """
        return 1

    def node_reader(self) -> Iterator[Tuple[int, Dict[str, str]]]:
        """Reader for NCBI names.dmp listing nodes and properties

        Returns:
           Generator of tax_id and property dictionary. The
           latter has ``name`` set to the scientific name from
           the NCBI file.
        """
        if self.names_fn is None:
            return
        LOG.info("Reading nodes from '%s'", self.names_fn)
        with open(self.names_fn, "r") as names_fd:
            reader = csv.DictReader(
                names_fd, delimiter='|',
                fieldnames=[
                    'tax_id', 'name', 'unique_name', 'name_class'
                ])
            for row in reader:
                if row['name_class'].strip() == 'scientific name':
                    tax_id = int(row['tax_id'])
                    data = dict((('name', row['name'].strip()),))
                    yield tax_id, data

    def edge_reader(self) -> Iterator[Tuple[int, int, str]]:
        """Reader for NCBI nodes.dmp, listing node and parent ids

        Returns:
           Generator over taxids ``(source, target)``
        """
        if self.nodes_fn is None:
            return
        LOG.info("Reading edges from '%s'", self.nodes_fn)
        with open(self.nodes_fn, "r") as nodes_fd:
            reader = csv.DictReader(
                nodes_fd, delimiter='|',
                fieldnames=[
                    'tax_id', 'parent_id', 'rank'
                ]
            )
            for row in reader:
                source_node = int(row['parent_id'])
                target_node = int(row['tax_id'])
                yield source_node, target_node, row['rank']

    @abstractmethod
    def get_name(self, tax_id: int) -> str:
        """Get the scientific name of a node

        Args:
          tax_id: ID of taxonomy node
        Returns:
          The NCBI taxonomy scientific name, or "Unknown"
          for ``tax_id``s not found or named in the tree
        """
        raise NotImplementedError()

    @abstractmethod
    def get_lineage(self, tax_id: int) -> str:
        """Get the lineage for a node

        Args:
          tax_id: ID of taxonomy node

        Returns:
          The names of all nodes from (excluding) the root node until
          (including) node given in ``tax_id``. Names are separated by
          semicolons. If the node ``tax_id`` is not found in the tree,
          "Unknown" is returned. Nodes not described are not mentioned
          in the linage string.
        """
        raise NotImplementedError()

    @abstractmethod
    def get_subtree_ids(self, name: str) -> Set[int]:
        """Get the ``tax_id``s for node ``name`` and all subnodes

        Args:
          name: Scientific name of node

        Returns:
          Set of ``tax_id``s of the named node and all subnodes.
        """
        raise NotImplementedError()

    def make_filter(self,
                    include: Optional[Sequence[str]] = None,
                    exclude: Optional[Sequence[str]] = None) -> Callable[[int], bool]:
        """Create a filtering function

        Args:
          include: List of strings matching NCBI taxonomy scientific
            names of nodes. Only tax_ids of nodes below any of the
            listed "clades" will be accepted by the generated filter
            function.

          exclude: List of strings matching NCBI taxonomy scientific
            names of nodes. IDs of nodes below any of the listed
            "clades" will be rejected by the generated filter
            function.

        Returns:
          A function that takes an integer ``tax_id`` and returns
          a boolean indicating whether the tax_id should be included.
        """
        if include is None:
            include_ids = set()
        else:
            include_ids = set(
                tid
                for name in include
                for tid in self.get_subtree_ids(name)
            )
        if exclude is None:
            exclude_ids = set()
        else:
            exclude_ids = set(
                tid
                for name in exclude
                for tid in self.get_subtree_ids(name)
            )
        if not include_ids:
            if not exclude_ids:
                def null_filter(_tax_id: int) -> bool:
                    """Always true"""
                    return True
                return null_filter

            def exclude_filter(taxid: int) -> bool:
                """True if taxid not excluded"""
                return taxid not in exclude_ids
            LOG.info("Made filter excluding %i tax_ids",
                     len(exclude_ids))
            return exclude_filter

        include_ids = include_ids.intersection(exclude_ids)

        def include_filter(taxid: int) -> bool:
            """True if ``taxid`` included"""
            return taxid in include_ids
        LOG.info("Made filter including %i tax_ids",
                 len(include_ids))
        return include_filter

    @staticmethod
    def is_null() -> bool:
        """Tests if this is an instance of TaxonmyNull (no taxonomy)"""
        return False


class TaxonomyGT(Taxonomy):
    """NCBI Taxonomy class using ``graph-tool`` to manage tree"""
    def load_tree(self) -> gt.Graph:
        LOG.info("Loading graph from '%s'", self.nodes_fn)
        tree = gt.load_graph_from_csv(
            self.nodes_fn,
            directed=True, eprop_types=[], string_vals=False,
            ecols=(1, 0),
            csv_options={'delimiter': '|'})
        tree.vp.name = tree.new_vertex_property("string")
        for tax_id, data in self.node_reader():
            tree.vp.name[tree.vertex(tax_id)] = data['name']
        tree.vp.rank = tree.new_vertex_property("string")
        for _, tax_id, rank in self.edge_reader():
            tree.vp.rank[tree.vertex(tax_id)] = rank.strip()
        return tree

    def load_tree_binary(self) -> gt.Graph:
        return gt.load_graph(self.path)

    def save_tree_binary(self, path: str) -> None:
        self.tree.save(path)

    def get_name(self, tax_id: int) -> str:
        try:
            return self.tree.vp.name[self.tree.vertex(tax_id)] or 'Unknown'
        except ValueError:
            return 'Unknown'

    def _get_path_to_root(self, target):
        res = [target]
        while int(target) != self.root:
            target = next(target.in_edges()).source()
            res.append(target)
        res.reverse()
        return res

    def get_lineage(self, tax_id: int) -> OrderedDict:
        # Get node and check that it's "real":
        try:
            target = self.tree.vertex(tax_id)
            # Check target has parent by tiggering StopIteration if it does not
            next(target.in_edges())
        except (ValueError, StopIteration):
            return OrderedDict((('no rank', 'Unknown'),))
        nodes = self._get_path_to_root(target)
        return OrderedDict((self.tree.vp.rank[node],  self.tree.vp.name[node])
                           for node in nodes[1:])


    def get_subtree_ids(self, name: str) -> Set[int]:
        vertices = gt_util.find_vertex(self.tree, self.tree.vp.name, name)
        if len(vertices) != 1:
            return set()
        ids = set(
            int(tgt) for _, tgt in
            gt_search.dfs_iterator(self.tree, vertices[0])
        )
        ids.add(int(vertices[0]))
        return ids


class TaxonomyNX(Taxonomy):
    """Implementation of Taxonomy class using NetworkX

    No longer actually used.
    """
    def load_tree_binary(self) -> nx.DiGraph:
        return nx.read_gpickle(self.path)

    def save_tree_binary(self, path: str) -> None:
        nx.write_gpickle(self.tree, path)

    def load_tree(self) -> nx.DiGraph:
        tree = nx.DiGraph()
        tree.add_nodes_from(self.node_reader())
        tree.add_edges_from(self.edge_reader())
        return tree

    def get_subtree_ids(self, name: str) -> Set[int]:
        for node, data in self.tree.nodes(data=True):
            if data.get('name') == name:
                taxid: int = node
                break
        else:
            return set()
        ids = set(*nx.descendants(self.tree, taxid))
        ids.add(taxid)
        return ids

    def get_lineage(self, tax_id: int) -> str:
        if tax_id not in self.tree:
            return "Unknown"
        path = nx.shortest_path(self.tree, self.root, tax_id)
        return '; '.join(self.tree.nodes[node].get('name')
                         for node in path[1:])

    def get_name(self, tax_id: int) -> str:
        if tax_id not in self.tree:
            return "Unknown"
        name: str = self.tree.nodes[tax_id].get('name')
        return name


class TaxonomyNull(Taxonomy):
    """Implementation of Taxonomy Class not doing anything

    Used when no taxonomy data is provided.
    """
    def load_tree_binary(self) -> None:
        pass

    def save_tree_binary(self, _path: str) -> None:
        pass

    def load_tree(self) -> None:
        pass

    def get_name(self, tax_id: int) -> str:
        return "Unknown"

    def get_lineage(self, tax_id: int) -> str:
        return "Unknown"

    def get_subtree_ids(self, name: str) -> Set[int]:
        return set()

    def make_filter(self,
                    include: Optional[Sequence[str]] = None,
                    exclude: Optional[Sequence[str]] = None) -> Callable[[int], bool]:
        def null_filter(_tax_id: int) -> bool:
            """Always return true"""
            return True
        return null_filter

    @staticmethod
    def is_null() -> bool:
        return True


def load_taxonomy(path: Optional[str], library: Optional[str] = None) -> Taxonomy:
    """Makes object of class Taxonomy

    If ``path`` is ``None``, returns an object of
    class `TaxonomyNull`, otherwise instanciates
    `TaxonomyGT`
    """
    if library is None:
        if path is None:
            return TaxonomyNull(None)
        elif path.endswith(".nx"):
            library = "networkx"
        else:
            library = "graph_tool"
    if path is None:
        raise RuntimeError(
            "If taxonomy library is specified, path cannot be None"
        )
    if library == "networkx":
        return TaxonomyNX(path)
    return TaxonomyGT(path)
