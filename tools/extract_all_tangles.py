#!/usr/bin/env python3
import sys
from pathlib import Path

# Add ../tools to sys.path
sys.path.insert(0, str(Path(__file__).parent.parent ))

import src.input_parsing
import logging
from src.node_id_mapper import NodeIdMapper
import src.logging_utils
import networkx as nx
import os

cov_variation = 1.5
#bubble lens should not be significantly different
len_variation = 1.5
component_length_threshold = 1000000
path_node_count_threshold = 20
definitely_not_tangle = 1000000
max_tip_length = 30000




def max_nonbranching_from_node (graph, node):
    path = [node]    
    next_node  = node    
    while graph.out_degree(next_node) == 1:        
        next_node = list(graph.successors(next_node))[0]                    
        if graph.in_degree(next_node) != 1:
            break
        path.append(next_node)

    return path

def parse_alignment_file(alignment_file, node_id_mapper):
    alignments = {}
    with open(alignment_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            query_name = parts[0]
            target_name = parts[5]
            if int(parts[11]) < 20:
                continue
            nodes = src.input_parsing.parse_gaf_string(target_name, node_id_mapper)
            print (nodes)
            for node in nodes:
                abs_node = abs(node)
                if abs_node not in alignments:
                    alignments[abs_node] = []
                alignments[abs_node].append(nodes)
    return alignments


def get_tangle_components(gfa_file, coverage_file, alignment_file):
    node_id_mapper = NodeIdMapper()
    hifi_graph = src.input_parsing.parse_gfa(gfa_file, node_id_mapper)
    coverage_data = src.input_parsing.read_coverage_file(coverage_file, node_id_mapper)

    total_len = 0
    covs = []
    for node in hifi_graph.nodes:
        if node in coverage_data:
            cur_len =  hifi_graph.nodes[node]['length']
            total_len += cur_len
            #print (f"Node: {node}, Length: {hifi_graph.nodes[node]['length']}, Coverage: {coverage_data[node]}")
            covs.append([cur_len, coverage_data[node]])
    sorted_covs = sorted(covs, key=lambda x: x[1])
    cur_len = 0
    for i in range(len(sorted_covs)):
        cur_len += sorted_covs[i][0]
        if cur_len >= total_len // 2:
            median_cov = sorted_covs[i][1]
            break

    print(f"Median coverage: {median_cov}")
    graph_copy = hifi_graph.copy()
    deleted = True
    #trivial tip clipping
    while deleted:
        deleted = False
        for node in hifi_graph.nodes():
            if not (node in graph_copy.nodes()):
                continue
            # Tip clipping: remove tips shorter than max_tip_length
            if graph_copy.in_degree(node) == 0 and graph_copy.out_degree(node) == 1:
                path = max_nonbranching_from_node(graph_copy, node)
                path_length = sum(graph_copy.nodes[n]['length'] for n in path)
                if path_length < max_tip_length:
                    print(f"Removing tip: {', '.join([node_id_mapper.node_id_to_name_safe(n) for n in path])}")
                    for n in path:
                        graph_copy.remove_node(n)
                        graph_copy.remove_node(-n)
                    deleted = True
    deleted = True
    #trivial bulge removal
    while deleted:    
        deleted = False
        for node in hifi_graph.nodes:
            if node in graph_copy.nodes():
                next_nodes =  list (graph_copy.successors(node))
                if len(next_nodes) == 2:
                    next_next = [max_nonbranching_from_node(graph_copy, next_nodes[0]), max_nonbranching_from_node(graph_copy, next_nodes[1])]
                    end_nodes = []
                    plens = [sum(graph_copy.nodes[n]['length'] for n in path) for path in next_next]
                    for i in range(2):
                        end_nodes.append(list(graph_copy.successors(next_next[i][-1])))
                    #length/coverage check?
                    if plens[0] > plens[1] * len_variation or plens[1] > plens[0] * len_variation:
                        continue
                    if len (end_nodes[0]) == 1 and len (end_nodes[1]) == 1 and end_nodes[0][0] == end_nodes[1][0]:
                        deleted = True
                        print (f"Removing bulge: {', '.join([node_id_mapper.node_id_to_name_safe(n) for n in next_next[0] ])}")
                        for node in next_next[0]:
                            graph_copy.remove_node(node)
                            graph_copy.remove_node(-node)

    #remove long 1-1 chains
    to_delete = set()
    bubble_cleared_copy = graph_copy.copy()
    for node in graph_copy.nodes():          
        #skip checking midpath  
        if graph_copy.in_degree(node) == 1:
            prev_node = list(graph_copy.predecessors(node))[0]
            if graph_copy.out_degree(prev_node) == 1:
                continue
        cur_path = max_nonbranching_from_node(graph_copy, node)
        cur_len = 0
        for node in cur_path:
            cur_len += graph_copy.nodes[node]['length']
        if cur_len > component_length_threshold or len (cur_path) > path_node_count_threshold:
            print(f"Clearing bubble chain length {cur_len}")
            print (f"Bubble nodes: {', '.join([node_id_mapper.node_id_to_name_safe(n) for n in cur_path])}")            
            for node in cur_path:
                to_delete.add(node)
                to_delete.add(-node)
    for node_id in to_delete:
        graph_copy.remove_node(node_id)
    alignments = parse_alignment_file(alignment_file, node_id_mapper)
    tangle_components = list(nx.weakly_connected_components(graph_copy))
    for tangle in tangle_components:
        if len(tangle) == 1:
            continue
        neighbors = set()
        for node in tangle:
            for neighbor in bubble_cleared_copy.successors(node):
                if not (neighbor in tangle):
                    neighbors.add(neighbor)
            for neighbor in bubble_cleared_copy.predecessors(node):
                if not (neighbor in tangle):
                    neighbors.add(neighbor)
        
        filtered_neighbors = set()
        valid_tangle = True
        for n in neighbors:
            #TODO: verkko-based check instead of coverage
            abs_node = abs(n)
            if abs_node in alignments and len(alignments[abs_node]) == 1:
                filtered_neighbors.add(abs_node)
            else:
                if hifi_graph.out_degree(n) > 1 and not (list(hifi_graph.successors(n))[0] in tangle):
                    if hifi_graph.out_degree(n) == 2:
                        next_nodes = list(hifi_graph.successors(n))
                        print (f"fixing coverage, neighbor {node_id_mapper.node_id_to_name_safe(n)} coverage {coverage_data[abs(n)]} replacing with {', '.join([node_id_mapper.node_id_to_name_safe(x) for x in next_nodes])}")
                        for rn in next_nodes:
                            filtered_neighbors.add(abs(rn))
                    else:
                        print (f"Failed to fix coverage, neighbor {node_id_mapper.node_id_to_name_safe(n)} has out_degree {hifi_graph.out_degree(n)} and in_degree {hifi_graph.in_degree(n)}")
                        valid_tangle = False
                        break
                elif hifi_graph.in_degree(n) > 1 and not (list(hifi_graph.predecessors(n))[0] in tangle):
                    if hifi_graph.in_degree(n) == 2:
                        next_nodes = list(hifi_graph.predecessors(n))
                        print (f"fixing coverage, neighbor {node_id_mapper.node_id_to_name_safe(n)} coverage {coverage_data[abs(n)]} replacing with {', '.join([node_id_mapper.node_id_to_name_safe(x) for x in next_nodes])}")
                        for rn in next_nodes:
                            filtered_neighbors.add(abs(rn))
                    else:
                        print (f"Failed to fix coverage, neighbor {node_id_mapper.node_id_to_name_safe(n)} has in_degree {hifi_graph.in_degree(n)} and in_degree {hifi_graph.in_degree(n)}")
                        valid_tangle = False
                        break
                else:
                    print (f"Failed to fix coverage, neighbor {node_id_mapper.node_id_to_name_safe(n)} has in_degree {hifi_graph.in_degree(n)} and out_degree {hifi_graph.out_degree(n)}")
                    valid_tangle = False
                    break

        if valid_tangle:
            if len (filtered_neighbors) == 2 or len(filtered_neighbors) == 4:
                for first in filtered_neighbors:
                    if first in alignments and len(alignments[first]) == 1:
                        for second in filtered_neighbors:
                            if first < second and second in alignments and len(alignments[second]) == 1:
                                if alignments[first] == alignments[second]:
                                    print (f"Valid verkko traversal to compare between {first} and {second}")
                                else:
                                    print (f"Alignment splitted between {first} and {second}")
                    else:
                        print (f"Multiple alignments for {first}")
            else:
                print (f"Unexpected number of filtered neighbors: {len(filtered_neighbors)} {','.join([node_id_mapper.node_id_to_name_safe(node) for node in filtered_neighbors])}")
        print(f"Tangle component of length {sum(graph_copy.nodes[node]['length'] for node in tangle)}: total nodes {len(tangle)}, external connections: {len(filtered_neighbors)} {','.join([node_id_mapper.node_id_to_name_safe(node) for node in filtered_neighbors])}.")
        print (",".join([node_id_mapper.node_id_to_name_safe(node) for node in tangle]))
        print ("\n")


if __name__ == "__main__":
    #basedir = sys.argv[1]
    basedir = "/data/antipovd2/verkko2_paper/HG002/v2.2.1_hic/"
    graph = os.path.join(basedir, "2-processGraph", "unitig-unrolled-hifi-resolved.noseq.gfa")
    coverage = os.path.join(basedir, "2-processGraph", "unitig-unrolled-hifi-resolved.ont-coverage.csv")
    alignment_file = os.path.join(basedir, "node_realign", "scaffolds2utig1.gaf")
    get_tangle_components(graph, coverage, alignment_file)