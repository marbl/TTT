#!/usr/bin/env python3
"""
Tangle class to encapsulate all tangle-related data and reduce function parameters.
"""

from typing import Optional, Dict, Set, List, Tuple
import networkx as nx
from .node_id_mapper import NodeIdMapper


class Tangle:
    """
    Encapsulates all information about a tangle region in an assembly graph.
    
    A tangle represents a complex region in the assembly graph that requires
    specialized traversal to resolve structural variants or repetitive sequences.
    
    Core structural attributes (initialized in __init__):
        nodes (set): Oriented node IDs within the tangle (includes RC nodes as negatives)
        nor_nodes (set): Non-oriented/absolute node IDs (nodes without orientation)
        boundary_nodes (dict): Maps incoming boundary nodes → matching outgoing boundary nodes
        original_graph (nx.DiGraph): Original GFA graph with node attributes
        dual_graph (nx.DiGraph): Dual graph where nodes=junctions, edges=original nodes
        node_id_mapper (NodeIdMapper): Bidirectional mapping between node IDs and string names
    
    Solution/derived attributes (added progressively during processing):
        cleaned_tips (list): Tip nodes removed during pre-processing
        coverage_dict (dict): Node ID → coverage value mapping
        median_unique_coverage (float): Estimated coverage for unique nodes
        coverage_range (tuple/list): [low, high] coverage bounds for unique nodes
        multiplicities (dict): Node ID → multiplicity value (from MIP solution)
        multi_graph (nx.MultiDiGraph): Multiplied dual graph after MIP solution
        detected_coverage (float): Final calculated unique coverage from MIP
    """
    
    def __init__(self, nodes, nor_nodes, boundary_nodes, original_graph, dual_graph, node_id_mapper):
        """
        Initialize a Tangle with core structural information.
        
        Args:
            nodes (set): Oriented node IDs within the tangle
            nor_nodes (set): Non-oriented/absolute node IDs
            boundary_nodes (dict): Maps incoming → outgoing boundary nodes
            original_graph (nx.DiGraph): Original GFA graph
            dual_graph (nx.DiGraph): Dual graph representation
            node_id_mapper (NodeIdMapper): Node ID to name mapper
        """
        # Core structural data
        self.nodes: Set[int] = nodes
        self.nor_nodes: Set[int] = nor_nodes
        self.boundary_nodes: Dict[int, int] = boundary_nodes
        self.original_graph: nx.DiGraph = original_graph
        self.dual_graph: nx.DiGraph = dual_graph
        self.node_id_mapper: NodeIdMapper = node_id_mapper
        
        # Solution/derived attributes (to be set later)
        # These will be populated progressively as the pipeline executes
        self.cleaned_tips: List[int] = []
        self.coverage_dict: Dict[int, float] = {}
        self.median_unique_coverage: Optional[float] = None
        self.coverage_range: Optional[List[float]] = None
        self.multiplicities: Dict[int, int] = {}
        self.multi_graph: Optional[nx.MultiDiGraph] = None
        self.detected_coverage: Optional[float] = None
