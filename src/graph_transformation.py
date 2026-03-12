#!/usr/bin/env python3
import logging
import random
import networkx as nx
from .logging_utils import log_assert
from .node_id_mapper import NodeIdMapper

def rc_node(node):
    """RC nodes stored as negative"""
    return -node

def get_canonical_nodepair(oriented_node1, oriented_node2, original_graph, node_mapper):
    """Returns a canonical lexicographically smallest link for a junction"""
    to_node = oriented_node2
    from_node = oriented_node1
    #sink and source should not be in original graph, special case
    if to_node == 0 or from_node == 0:
        return (from_node, to_node)
    for node in original_graph.successors(oriented_node1):
        if node < to_node:
            to_node = node
    
    for node in original_graph.successors(rc_node(to_node)):
        if rc_node(node) < from_node:
            from_node = rc_node(node)
    if not (to_node in original_graph.successors(from_node)):
        logging.warning(f"Irregular junction {node_mapper.node_id_to_name_safe(from_node)} {node_mapper.node_id_to_name_safe(to_node)} {node_mapper.node_id_to_name_safe(oriented_node1)} {node_mapper.node_id_to_name_safe(oriented_node2)}")
        logging.warning(f"{node_mapper.node_id_to_name_safe(from_node)} -> {[node_mapper.node_id_to_name_safe(n) for n in original_graph.successors(from_node)]}")
        logging.warning(f"{node_mapper.node_id_to_name_safe(rc_node(to_node))} -> {[node_mapper.node_id_to_name_safe(n) for n in original_graph.successors(rc_node(to_node))]}")        
    return (from_node, to_node)

def get_canonical_rc_vertex(v, original_graph, node_mapper):
    """Returns a canonical reverse connection vertex for a given vertex"""
    return get_canonical_nodepair(-v[1], -v[0], original_graph, node_mapper)

def create_dual_graph(original_graph: nx.MultiDiGraph, node_mapper):
    """
    Transform to dual graph, vertices = junctions, edges = old nodes
    TODO: vertices should be not exactly junctions but junction+nodes since z-connections are legit
    """
    dual_graph = nx.MultiDiGraph()
    canonical_edges_set = set() # To store unique canonical edges (dual nodes)
    logging.info("Creating dual graph representation (nodes = canonical connections) using NetworkX...")

    # Iterate through each oriented node 'v' which acts as the junction
    for v in sorted(original_graph.nodes):
        # Determine the canonical edge C1 leading INTO node v
        predecessors = list(original_graph.successors(rc_node(v)))
        if not predecessors:
            # v is a source tip (no incoming edges for this orientation)
            C1 = ("TIP", v)
        else:
            # Pick any predecessor 'u' to define the canonical incoming edge
            # Note: get_canonical_edge expects (from, to), so use (pred, v)
            u = rc_node(predecessors[0])
            C1 = get_canonical_nodepair(u, v, original_graph, node_mapper)
        canonical_edges_set.add(C1)

        # Determine the canonical edge C2 leading OUT FROM node v
        successors = list(original_graph.successors(v))
        if not successors:
            # v is a sink tip (no outgoing edges for this orientation)
            C2 = (v, "TIP")
        else:
            # Pick any successor 'w' to define the canonical outgoing edge
            # Note: get_canonical_edge expects (from, to), so use (v, succ)
            w = successors[0]
            C2 = get_canonical_nodepair(v, w, original_graph, node_mapper)
        canonical_edges_set.add(C2)

        # Add the edge in the dual graph representing the transition through v
        # Ensure nodes C1 and C2 exist before adding the edge (add_nodes_from handles this later too)
        dual_graph.add_node(C1)
        dual_graph.add_node(C2)
        dual_graph.add_edge(C1, C2, original_node=v)
        logging.debug(f"Added edge from {C1} to {C2} in dual graph, original node: {node_mapper.node_id_to_name_safe(v)}")

    logging.info(f"Total oriented nodes in original GFA structure: {len(original_graph.nodes)}")
    logging.info(f"Total unique canonical edges (dual nodes): {len(canonical_edges_set)}")
    logging.info(f"Total nodes in final dual graph: {dual_graph.number_of_nodes()}") # Should match set size
    logging.info(f"Total edges (connections through original nodes) in dual graph: {dual_graph.number_of_edges()}")
    #logging.debug(f"Edges of dual graph: {list(dual_graph.edges())}")
    return dual_graph

def is_tangle_vertex(C, tangle_set):
    return C[0] in tangle_set or C[1] in tangle_set

def create_multi_dual_graph(tangle):
    """
    edge of multiplicity X -> X multiedges multiplicity 1
    """
    multi_dual_graph = nx.MultiDiGraph()
    logging.info("Creating multi-dual graph from dual graph and multiplicities...")

    # Ensure all nodes from the original dual graph are present
    for node in tangle.dual_graph.nodes():
        if is_tangle_vertex(node, tangle.nodes):
            multi_dual_graph.add_node(node)

    #adding start and sink nodes
    incomings = set(tangle.boundary_nodes.keys())
    outgoings = set(tangle.boundary_nodes.values())
    for b in tangle.boundary_nodes:
        incomings.add(-tangle.boundary_nodes[b])
        outgoings.add(-b)
    edges_added = 0

    for b in tangle.boundary_nodes.keys():
        for next_node in tangle.original_graph.successors(b):
            fw_edge = get_canonical_nodepair(b, next_node, tangle.original_graph, tangle.node_id_mapper)                                
            logging.debug(f"start {tangle.node_id_mapper.node_id_to_name_safe(b)} {tangle.node_id_mapper.node_id_to_name_safe(next_node)}")
            start = (0, b)
            multi_dual_graph.add_node(start)
            multi_dual_graph.add_edge(start, fw_edge, original_node=b, key = f"{b}_{edges_added}")
            edges_added += 1
            break 
    for b in tangle.boundary_nodes.values():
        for next_node in tangle.original_graph.predecessors(b):
            bw_edge = get_canonical_nodepair(next_node, b, tangle.original_graph, tangle.node_id_mapper)                                
            logging.debug(f"start {tangle.node_id_mapper.node_id_to_name_safe(next_node)} {tangle.node_id_mapper.node_id_to_name_safe(b)}")
            start = (b, 0)
            multi_dual_graph.add_node(start)
            multi_dual_graph.add_edge(bw_edge, start, original_node=b, key = f"{b}_{edges_added}")
            edges_added += 1
            break
    
    # Iterate through edges of the dual graph and add them to the multi-graph
    # based on the multiplicity of the original node they represent.
    for u, v, data in tangle.dual_graph.edges(data=True):
        if (not u in multi_dual_graph.nodes()) or (not v in multi_dual_graph.nodes()):
            continue
        original_node_oriented = data.get('original_node')
        if not original_node_oriented:
            logging.warning(f"Edge ({u}, {v}) in dual graph is missing 'original_node' attribute. Skipping.")
            continue
        
        original_node_base = original_node_oriented #abs(original_node_oriented)

        multiplicity = tangle.multiplicities[original_node_base]

        if multiplicity < 0:
            logging.error(f"Negative multiplicity {multiplicity} for node {tangle.node_id_mapper.node_id_to_name_safe(original_node_base)}. Treating as 0 for edge ({u} -> {v}).")
            exit(0)
        # Add the edge 'multiplicity' times to the multi-dual graph
        logging.debug(f"Adding {multiplicity} multiedges for {tangle.node_id_mapper.node_id_to_name_safe(original_node_oriented)}")
        for _ in range(multiplicity):
            # Add edge with the original node attribute
            multi_dual_graph.add_edge(u, v, original_node=original_node_oriented, key = str(original_node_oriented) + "_"+str(edges_added))
            edges_added += 1

    logging.info(f"Created multi-dual graph with {multi_dual_graph.number_of_nodes()} nodes.")
    logging.info(f"Added {edges_added} edges to multi-dual graph based on multiplicities (original dual graph had {tangle.dual_graph.number_of_edges()} unique edges).")
    for v in multi_dual_graph.nodes():
        logging.debug(f"{v}: in {multi_dual_graph.in_degree[v] } out {multi_dual_graph.out_degree[v] }")    

    return multi_dual_graph

def get_traversable_subgraph(tangle):
    """
    Supplementary for Euler path search - search itself moved to path_optimizer.py
    """
    multi_dual_graph = tangle.multi_graph
    #we may have some disconnected cycles in reverse-complement tangle, reversing them.
    start_vertices = []
    end_vertices = []
    for n in tangle.boundary_nodes:
        #border tips encoding
        start_vertices.append((0, n))
        end_vertices.append((n, 0))

    logging.info(f"Start and end vertices in the graph {start_vertices}, {end_vertices}")  
    s = ""
    for b in tangle.boundary_nodes:
        s += f"{tangle.node_id_mapper.node_id_to_name_safe(b)}->{tangle.node_id_mapper.node_id_to_name_safe(tangle.boundary_nodes[b])} "    
    logging.info(f"Boundary nodes and their connections: {s}")
    border_nodes_count = len(tangle.boundary_nodes)
    # Only 1-1 or 2-2 tangles for now
    log_assert(border_nodes_count == 1 or border_nodes_count == 2, f"Only 1-1 or 2-2 tangles are supported")
    log_assert(len(start_vertices) == border_nodes_count and len(end_vertices) == border_nodes_count, f"Start and end vertices count mismatch: {len(start_vertices)} vs {border_nodes_count} or {len(end_vertices)} vs {border_nodes_count}")
    start_vertices.sort()
    start_vertex = start_vertices[0]
    matching_end_vertices = []
    for i in range (len(start_vertices)):
        matching_end_vertex_node = tangle.boundary_nodes[start_vertices[i][1]]
        matching_end_vertices.append((matching_end_vertex_node, 0))

    for e in multi_dual_graph.edges(keys=True):
        logging.debug(f"Edge {e}")

    reachable_vertices = set()
    unreachable_edges = set()
    previous_reachable_count = -1
    MAX_ITERATIONS = 2000
    logging.info (f"Total vertices in multi-dual graph before transformation: {len(multi_dual_graph.nodes())}")
    for _ in range (MAX_ITERATIONS):
        reachable_vertices.clear()
        for ind in range (len(start_vertices)):
            reachable_vertices.update(nx.descendants(multi_dual_graph, start_vertices[ind]))
            reachable_vertices.add(start_vertices[ind])
        unreachable_vertices = set(multi_dual_graph.nodes()) - reachable_vertices
        unreachable_edges.clear()
        unreachable_original_nodes = set()
        for e in multi_dual_graph.edges(keys=True):
            if e[0] in unreachable_vertices:
                unreachable_edges.add(e)
                data = multi_dual_graph.get_edge_data(e[0], e[1], key=e[2])
                unreachable_original_nodes.add(data['original_node'])
        logging.info(f"Clearing unreachable edges: iteration {_} {len(unreachable_edges)} present")
        if len(reachable_vertices) == previous_reachable_count or len (unreachable_edges) == 0:
            logging.info(f"Reachable vertices count stabilized at {len(reachable_vertices)} after iteration {_}")
            break
        previous_reachable_count = len(reachable_vertices)
#all unreachable edges are combination of cycles, so we can reverse them
        for e in unreachable_edges:
            logging.debug(f"Reversing edge: {e}")
            data = multi_dual_graph.get_edge_data(e[0], e[1], key=e[2])
            logging.debug(f"Data {data}")
            multi_dual_graph.remove_edge(e[0], e[1], key = e[2])
            if e[2][0] == '-':
                new_key = e[2][1:]
            else:
                new_key = '-' + e[2]
            multi_dual_graph.add_edge(get_canonical_rc_vertex(e[1], tangle.original_graph, tangle.node_id_mapper), get_canonical_rc_vertex(e[0], tangle.original_graph, tangle.node_id_mapper), original_node=-int(data['original_node']), key = new_key)
            logging.debug(f"Added reversed edge from {get_canonical_rc_vertex(e[1], tangle.original_graph, tangle.node_id_mapper)} to {get_canonical_rc_vertex(e[0], tangle.original_graph, tangle.node_id_mapper)} with original node {-int(data['original_node'])} and key {new_key}")
    if _ == MAX_ITERATIONS - 1:
        logging.warning(f"Maximum iterations reached while trying to stabilize reachable vertices. Some traversal will still be found but it is not an expected situation. Please create a github issue to investigate this case.")
    
    mults = {}
    for e in multi_dual_graph.edges(keys=True):
        data = multi_dual_graph.get_edge_data(e[0], e[1], key=e[2])
        original_node = data['original_node']
        if not original_node in mults:
            mults[original_node] = 0
        mults[original_node] += 1
    logging.debug ("Final multiplicities after reassigning rc isolated cycles")
    for n in sorted(mults.keys()):
        logging.debug(f"Mult for {tangle.node_id_mapper.node_id_to_name_safe(n)} is {mults[n]} ")                   

    #TODO: possibly be fixed via considering rc loop switch in addition to regular loop rc switch.
    for i, v in enumerate(matching_end_vertices):
        descendants = nx.descendants(multi_dual_graph, start_vertices[i])
        if not v in descendants:            
            logging.error(f"One of the exit nodes {v} is not reachable from corresponding start vertex {start_vertices[i]}")
            logging.error(f"It is weird but possible situation for 2-2 tangles that mostly consists of a large inverted repeat, TTT currently cannot handle it properly")
            logging.error(f"Suggested workaround is to reconsider tangle borders to increase the tangle a bit (add additional diploid bulge to the tangle")
            exit(1)

    if (len(unreachable_edges) != 0):
        logging.warning(f"{len(unreachable_edges)} unreachable edges with nonzero multiplicities remains. They will be ignored in path finding.")
        logging.warning(f"This happens when after removal of edges with estimated multiplicity 0 some cycles of positive multiplicity are disconnected from the main graph component. Usually this indicates problems with graph structure or really uneven coverage")        
        out_str = "Unreachable edges: "
        for e in unreachable_edges:                    
            data = multi_dual_graph.get_edge_data(e[0], e[1], key=e[2])
            out_str += f"{tangle.node_id_mapper.node_id_to_name_safe(data['original_node'])} "
            multi_dual_graph.remove_edge(e[0], e[1], key = e[2])
        logging.warning(out_str)

    if border_nodes_count == 2:
        aux_node_str = "AUX"
        aux_int_id = tangle.node_id_mapper.parse_node_id(aux_node_str)
        multi_dual_graph.add_edge(matching_end_vertices[0], start_vertices[1], original_node=aux_int_id, key=f"{aux_int_id}_0")
        logging.info(f"Added auxiliary edge from {tangle.node_id_mapper.node_id_to_name_safe(matching_end_vertex_node)} to {tangle.node_id_mapper.node_id_to_name_safe(start_vertices[1][1])}")       

    reachable_subgraph = multi_dual_graph.subgraph(reachable_vertices)
    return reachable_subgraph, start_vertex
