"""This module contains classes for handling operations on the NCBI taxonomy"""

import os
import csv
import logging
from abc import ABC, abstractmethod
from typing import (
    Set,
    Sequence,
    Any,
    Iterator,
    Iterable,
    Dict,
    Tuple,
    Callable,
    Optional,
    Union,
    List,
    cast,
)
from collections import Counter

import graph_tool as gt  # type: ignore
import graph_tool.search as gt_search  # type: ignore
import graph_tool.util as gt_util  # type: ignore

LOG = logging.getLogger(__name__)


class Taxonomy(ABC):
    """Base class holding NCBI taxonomy

    The tree is implemented in subclass TaxonomyGT using graph-tool library.

    Args:
      path: Either path to directory containing names.dmp and nodes.dmp or
            path to the binary dump of the tree.
    """

    #: Name of file containing nodes
    NODES_FN = "nodes.dmp"
    #: Name of file containing names
    NAMES_FN = "names.dmp"
    #: Name of file containing merged nodes
    MERGED_FN = "merged.dmp"

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

    @property
    def merged_fn(self) -> Optional[str]:
        """Path to the merged.dmp file"""
        return self._find_file(self.MERGED_FN)

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
        with open(self.names_fn, "r", encoding="UTF-8") as names_fd:
            reader = csv.DictReader(
                names_fd,
                delimiter="|",
                fieldnames=["tax_id", "name", "unique_name", "name_class"],
            )
            for row in reader:
                if row["name_class"].strip() == "scientific name":
                    tax_id = int(row["tax_id"])
                    data = dict((("name", row["name"].strip()),))
                    yield tax_id, data

    def edge_reader(self) -> Iterator[Tuple[int, int, str]]:
        """Reader for NCBI nodes.dmp, listing node and parent ids

        Returns:
           Generator over taxids ``(source, target, rank)``
        """
        if self.nodes_fn is None:
            return
        LOG.info("Reading edges from '%s'", self.nodes_fn)
        with open(self.nodes_fn, "r", encoding="UTF-8") as nodes_fd:
            reader = csv.DictReader(
                nodes_fd, delimiter="|", fieldnames=["tax_id", "parent_id", "rank"]
            )
            for row in reader:
                source_node = int(row["parent_id"])
                target_node = int(row["tax_id"])
                yield source_node, target_node, row["rank"]

    def merged_reader(self) -> Iterator[Tuple[int, int]]:
        """Reader for NCBI merged.dmp, listing taxids merged into another

        Returns:
           Iterator over taxids ``(old, new)``
        """
        if self.merged_fn is None:
            return
        LOG.info("Reading merged nodes from '%s'", self.merged_fn)
        with open(self.merged_fn, "r", encoding="UTF-8") as merged_fd:
            reader = csv.DictReader(merged_fd, delimiter="|", fieldnames=["old", "new"])
            for row in reader:
                yield int(row["old"]), int(row["new"])

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
    def get_lineage(self, tax_id: int) -> Sequence[str]:
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
    def get_lineage_ranks(self, tax_id: int) -> Sequence[str]:
        """Get list of ranks for all names in lineage of ``taxid``"""
        raise NotImplementedError()

    @abstractmethod
    def get_rank(self, tax_id: int, rank: str) -> str:
        """Get specific rank ``rank`` from lineage of ``taxid``"""
        raise NotImplementedError()

    @abstractmethod
    def get_lca(
        self, tax_ids: Union[Iterable[int], Iterable[Tuple[int, int]]]
    ) -> "LCA":
        """Create an LCA object"""
        raise NotImplementedError()

    @abstractmethod
    def get_taxid(self, name: str) -> Optional[int]:
        """Fetch NCBI taxid for given name

        Returns None if no match or multiple matches were found
        """
        raise NotImplementedError()

    @abstractmethod
    def get_siblings(self, tax_id: int) -> Set[int]:
        """Fetch sibling taxa at same rank"""
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

    def make_filter(
        self,
        include: Optional[Sequence[str]] = None,
        exclude: Optional[Sequence[str]] = None,
    ) -> Callable[[int], bool]:
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
                tid for name in include for tid in self.get_subtree_ids(name)
            )
        if exclude is None:
            exclude_ids = set()
        else:
            exclude_ids = set(
                tid for name in exclude for tid in self.get_subtree_ids(name)
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

            LOG.info("Made filter excluding %i tax_ids", len(exclude_ids))
            return exclude_filter

        include_ids = include_ids.intersection(exclude_ids)

        def include_filter(taxid: int) -> bool:
            """True if ``taxid`` included"""
            return taxid in include_ids

        LOG.info("Made filter including %i tax_ids", len(include_ids))
        return include_filter

    @staticmethod
    def is_null() -> bool:
        """Tests if this is an instance of TaxonmyNull (no taxonomy)"""
        return False


class LCA(ABC):
    """Lowest common ancestor

    This is implemented as a class to allow merging and querying of
    the set of sequences under comparison.
    """

    def __init__(self, taxonomy: "Taxonomy") -> None:
        self.taxonomy = taxonomy

    @abstractmethod
    def add(self, taxids: Iterable[int]) -> "LCA":
        """Add taxids to LCA object"""
        raise NotImplementedError()

    @abstractmethod
    def classify(self, quorum: float, majority: float) -> Tuple[int, float]:
        """Find lowest common ancestor (LCA)"""
        raise NotImplementedError()

    def __str__(self):
        return "no classification"

    def dump(self):
        """Dump LCA class content to STDOUT for debugging"""
        print("LCA dump: Not implemented")


class TaxonomyGT(Taxonomy):
    """NCBI Taxonomy class using ``graph-tool`` to manage tree"""

    def load_tree(self) -> gt.Graph:
        LOG.info("Loading graph from '%s'", self.nodes_fn)
        tree = gt.load_graph_from_csv(
            self.nodes_fn,
            directed=True,
            eprop_types=[],
            hashed=False,
            ecols=(1, 0),
            csv_options={"delimiter": "|"},
        )
        tree.vp.name = tree.new_vertex_property("string")
        for tax_id, data in self.node_reader():
            tree.vp.name[tree.vertex(tax_id)] = data["name"]
        tree.vp.rank = tree.new_vertex_property("string")
        for _, tax_id, rank in self.edge_reader():
            tree.vp.rank[tree.vertex(tax_id)] = rank.strip()
        for src, dst in self.merged_reader():
            tree.add_edge(dst, src)
        return tree

    def load_tree_binary(self) -> gt.Graph:
        LOG.info("Loading taxonomy from %s", self.path)
        return gt.load_graph(self.path)

    def save_tree_binary(self, path: str) -> None:
        self.tree.save(path)

    def get_vertex(self, tax_id: int) -> Optional[gt.Vertex]:
        """Retrieve graph-tools vertex for given ``tax-id``

        This finds the nearest named node.

        Returns ``None`` if the tax_id is not in the tree. This may
        ocurr for example for tax_id's that have been deleted.
        """
        # Get node and check that it's "real":
        try:
            vertex = self.tree.vertex(tax_id)
        except ValueError:
            # tax_id larger than largest in tree
            return None
        if not vertex.in_degree():
            # not connected at all, empty index in the middle
            return None
        while not self.tree.vp.name[vertex]:
            vertex = next(vertex.in_edges()).source()
            if not vertex.in_degree():
                return None
        return vertex

    def get_name(self, tax_id: int) -> str:
        try:
            return self.tree.vp.name[self.get_vertex(tax_id)]
        except ValueError:
            return "Unknown"

    def get_path_to_root(self, target):
        """Retrieve parent chain for vertex ``target``"""
        res = [target]
        while int(target) != self.root:
            target = next(target.in_edges()).source()
            res.append(target)
        res.reverse()
        return res

    def get_lineage(self, tax_id: int) -> Sequence[str]:
        target = self.get_vertex(tax_id)
        if target is None:
            return ["Unknown"]
        nodes = self.get_path_to_root(target)
        return [self.tree.vp.name[node] for node in nodes[1:]]

    def get_lineage_ranks(self, tax_id: int) -> Sequence[str]:
        """Get list of ranks for all names in lineage of ``taxid``"""
        target = self.get_vertex(tax_id)
        if target is None:
            return ["no rank"]
        nodes = self.get_path_to_root(target)
        return [self.tree.vp.rank[node] for node in nodes[1:]]

    def get_rank(self, tax_id: int, rank: str) -> str:
        """Get specific rank ``rank`` from lineage of ``taxid``"""
        target = self.get_vertex(tax_id)
        if target is not None:
            nodes = self.get_path_to_root(target)
            for node in reversed(nodes):
                if self.tree.vp.rank[node] == rank:
                    return self.tree.vp.name[node]
        return "Unknown"

    def get_lca(
        self, tax_ids: Union[Iterable[int], Iterable[Tuple[int, int]]]
    ) -> "LCA":
        """Create an LCA object"""
        return LcaGt(self).add(tax_ids)

    def get_taxid(self, name: str) -> Optional[int]:
        LOG.info("Searching for species '%s'...", name)
        vertices = gt_util.find_vertex(self.tree, self.tree.vp.name, name)
        if len(vertices) != 1:
            LOG.error("Species '%s' not found in taxonomy", name)
            return None
        LOG.info("Found taxid %i", vertices[0])
        return vertices[0]

    def get_siblings(self, tax_id: int) -> Set[int]:
        target = self.get_vertex(tax_id)
        if target is None:
            return set()
        parent = next(target.in_edges()).source()
        rank = self.tree.vp.rank[target]
        siblings = set()
        # We specifically do not go through the subnodes as the NCBI taxonomy
        # weirdly has entries below "unclassified speciesXYZ" ranked as species.
        for edge in parent.out_edges():
            dst = edge.target()
            if dst != target and self.tree.vp.rank[dst] == rank:
                siblings.add(int(dst))
        return siblings

    def get_subtree_ids(self, name: str) -> Set[int]:
        vertices = gt_util.find_vertex(self.tree, self.tree.vp.name, name)
        if len(vertices) != 1:
            return set()
        ids = set(int(tgt) for _, tgt in gt_search.dfs_iterator(self.tree, vertices[0]))
        ids.add(int(vertices[0]))
        return ids


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

    def get_lineage(self, tax_id: int) -> Sequence[str]:
        return ["Unknown"]

    def get_lineage_ranks(self, tax_id: int) -> Sequence[str]:
        return ["Unknown"]

    def get_rank(self, tax_id: int, rank: str) -> str:
        return "Unknown"

    def get_lca(self, tax_ids: Union[Iterable[int], Iterable[Tuple[int, int]]]):
        raise NotImplementedError()

    def get_taxid(self, name: str) -> Optional[int]:
        return None

    def get_siblings(self, tax_id: int) -> Set[int]:
        return set()

    def get_subtree_ids(self, name: str) -> Set[int]:
        return set()

    def make_filter(
        self,
        include: Optional[Sequence[str]] = None,
        exclude: Optional[Sequence[str]] = None,
    ) -> Callable[[int], bool]:
        def null_filter(_tax_id: int) -> bool:
            """Always return true"""
            return True

        return null_filter

    @staticmethod
    def is_null() -> bool:
        return True


class LcaGt(LCA):
    """Notes:

    - Depth is not the same across the tree, euks lack the classic
      6-depth classification scheme. The issue here is that the Euks
      may be "deeper" than an alternative classification from a
      different domain of life. May not be relevant because of lack
      of overlap.

    """

    def __init__(self, taxonomy: TaxonomyGT) -> None:
        super().__init__(taxonomy)
        self.taxonomy: TaxonomyGT = self.taxonomy
        self.cnt: List[Counter] = [Counter()]
        self.total = 0

    def add(self, taxids: Union[Iterable[int], Iterable[Tuple[int, int]]]) -> "LCA":
        """Add taxids to LCA object

        ``taxids`` must be an iterable either of integer taxids or of
        tuples of integer taxids and weights to be used for the LCA.
        """
        first, *rest = taxids
        if isinstance(first, tuple):
            self.add_one(first[0], first[1])
            for taxid, score in cast(Iterable[Tuple[int, int]], rest):
                self.add_one(taxid, score)
        else:
            self.add_one(first)
            for taxid in cast(Iterable[int], rest):
                self.add_one(taxid)
        return self

    def add_one(self, taxid: int, score: int = 1) -> "LCA":
        """Add a single taxid to LCA object"""
        leaf_vertex = self.taxonomy.get_vertex(taxid)
        if leaf_vertex is None:
            LOG.error("Taxonomy node %i not found. Maybe NT out of date?", taxid)
            return self
        self.total += score
        path = self.taxonomy.get_path_to_root(leaf_vertex)
        for depth, vertex in enumerate(path):
            while depth + 1 > len(self.cnt):
                self.cnt.append(Counter())
            self.cnt[depth][vertex] += score
        return self

    def classify(self, quorum=0.9, majority=0.9):
        """Find lowest common ancestor (LCA) vertex

        The ``quorum`` and ``majority`` parameters can be used to
        relax the algorithm from actual strict LCA to accomodate real
        life NCBI taxonomy entries.

        - The ``majority`` parameter reduces the fraction of input
          taxids that need to share the ancestor from 100% to default
          90%. This prevents underclassification in the presence of
          spurious input taxa, such as can be expected from a BLAST
          search.

        - The ``quorum`` parameter allows for incomplete input
          classifications, requiring that only ``majority`` fraction
          of all input taxa be classified to the depth resulting.

        Issues:

        - The quorum does not currently differentiate between
          conflicting and non-conflicting incomplete classifications
          but completely ignores those of insufficient depth.

        - The NCBI taxonomy has varying ranks per depth. Some parts of
          the taxonomy are split into more ranks by adding "clade"
          ranks and/or using ranks such as infraorder not used for the
          prokaryotic domains. This is currently ignored. A
          "infraorder" in the Eukaryotic domain may therefore be
          reported over a theoretical better Bacterial
          classification. However, this should only be possible for
          quorum*majority < 0.5 (naively, didn't think it through).

        Returns a tuple of taxonomy node and fraction of agreeing
        classifications at the reported depth.

        """
        # Step through from the deepest classifications to the most
        # shallow and return once we find a node that meets quorum and
        # majority criteria.
        for counter in reversed(self.cnt):
            if not counter:
                # Nothing at this level, skip
                continue
            if counter.total() < self.total * quorum:
                # We do not have quorum (not enough "votes" given for
                # this depth into tree), skip
                continue
            top, top_count = counter.most_common(1)[0]
            if top_count < counter.total() * majority:
                # We do not have majority (not enough votes within
                # this depth given to most common node), skip
                continue
            return top, top_count / counter.total()
        # Fallback that should really never run:
        return self.taxonomy.root, 0

    def __str__(self):
        return "no classification"

    def dump(self):
        for depth, counter in enumerate(self.cnt):
            if not counter:
                print(f"{depth} - not present")
                continue
            (top, top_count), *other = counter.most_common(2)
            name = self.taxonomy.tree.vp.name[top]
            rank = self.taxonomy.tree.vp.rank[top]
            if other:
                sec = other[0][0]
                sec_count = other[0][1]
                sec_name = self.taxonomy.tree.vp.name[sec]
                sec_rank = self.taxonomy.tree.vp.rank[sec]
                sec_print = f"({sec_rank} {sec_name} {sec_count}/{counter.total()})"
            else:
                sec_print = ""
            print(f"{depth} - {rank} {name} {top_count}/{counter.total()} {sec_print}")


def load_taxonomy(path: Optional[str], library: Optional[str] = None) -> Taxonomy:
    """Makes object of class Taxonomy

    If ``path`` is ``None``, returns an object of
    class `TaxonomyNull`, otherwise instanciates
    `TaxonomyGT`
    """
    if library is None:
        if path is None:
            return TaxonomyNull(None)
    if path is None:
        raise RuntimeError("If taxonomy library is specified, path cannot be None")
    return TaxonomyGT(path)
