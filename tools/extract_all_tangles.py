#!/usr/bin/env python3
import argparse
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

import edlib
from Bio import SeqIO

PATH_2_FASTA_SCRIPT = "/data/antipovd2/devel/utils/sequence_tools/path2fasta.py" 
REFERENCE = "/data/antipovd2/res/TTT_paper/HG002/hg002v1.1.hpc.fasta"
#qREFERENCE = "/data/antipovd2/data/refs/hg002v1.1.hpc.concat.fasta"
cov_variation = 1.5
#bubble lens should not be significantly different
len_variation = 1.5
component_length_threshold = 1000000
path_node_count_threshold = 20
definitely_not_tangle = 1000000
MAX_TIP_LENGTH = 30000

BASE_SUBDIR = "tangles"

def align_fasta_files(fasta1, fasta2, target_file):
    # Load queries and targets
    queries = list(SeqIO.parse(fasta1, "fasta"))
    targets = list(SeqIO.parse(fasta2, "fasta"))
    out_f = open(target_file, "w")

    for qrec in queries:
        best_target = None
        best_result = None

        for trec in targets:
            result = edlib.align(str(qrec.seq), str(trec.seq), mode="HW")  # semiglobal
            if best_result is None or result["editDistance"] < best_result["editDistance"]:
                best_result = result
                best_target = trec

        out_f.write(f"{qrec.id}\t{best_result['editDistance']}\n")


def max_nonbranching_from_node (graph, node):
    path = [node]    
    next_node  = node    
    while graph.out_degree(next_node) == 1:        
        next_node = list(graph.successors(next_node))[0]                    
        if graph.in_degree(next_node) != 1:
            break
        path.append(next_node)
    return path

def is_tip (graph, node):
    return graph.in_degree(node) == 0 or graph.out_degree(node) == 0

def parse_alignment_file(alignment_file, node_id_mapper):
    alignments = {}
    with open(alignment_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            query_name = parts[0]
            target_name = parts[1]
            
            nodes = src.input_parsing.parse_gaf_string(target_name, node_id_mapper)
            #logging.info (nodes)
            for node in nodes:
                abs_node = abs(node)
                if abs_node not in alignments:
                    alignments[abs_node] = []
                alignments[abs_node].append(nodes)
    return alignments


def get_tangle_components(gfa_file, coverage_file, alignments, node_id_mapper):
    hifi_graph = src.input_parsing.parse_gfa(gfa_file, node_id_mapper)
    coverage_data = src.input_parsing.read_coverage_file(coverage_file, node_id_mapper)

    total_len = 0
    covs = []
    for node in hifi_graph.nodes:
        if node in coverage_data:
            cur_len =  hifi_graph.nodes[node]['length']
            total_len += cur_len
            #logging.info (f"Node: {node}, Length: {hifi_graph.nodes[node]['length']}, Coverage: {coverage_data[node]}")
            covs.append([cur_len, coverage_data[node]])
    sorted_covs = sorted(covs, key=lambda x: x[1])
    cur_len = 0
    for i in range(len(sorted_covs)):
        cur_len += sorted_covs[i][0]
        if cur_len >= total_len // 2:
            median_cov = sorted_covs[i][1]
            break

    logging.info(f"Median coverage: {median_cov}")
    tangle_only_graph = hifi_graph.copy()
    deleted = True
    #trivial tip clipping
    while deleted:
        deleted = False
        for node in hifi_graph.nodes():
            if not (node in tangle_only_graph.nodes()):
                continue
            # Tip clipping: remove tips shorter than max_tip_length
            if tangle_only_graph.in_degree(node) == 0 and tangle_only_graph.out_degree(node) == 1:
                path = max_nonbranching_from_node(tangle_only_graph, node)
                path_length = sum(tangle_only_graph.nodes[n]['length'] for n in path)
                if path_length < MAX_TIP_LENGTH:
                    logging.debug(f"Removing tip: {', '.join([node_id_mapper.node_id_to_name_safe(n) for n in path])}")
                    for n in path:
                        tangle_only_graph.remove_node(n)
                        tangle_only_graph.remove_node(-n)
                    deleted = True
    deleted = True
    #trivial bulge removal
    while deleted:    
        deleted = False
        for node in hifi_graph.nodes:
            if node in tangle_only_graph.nodes():
                next_nodes =  list (tangle_only_graph.successors(node))
                if len(next_nodes) == 2:
                    next_next = [max_nonbranching_from_node(tangle_only_graph, next_nodes[0]), max_nonbranching_from_node(tangle_only_graph, next_nodes[1])]
                    end_nodes = []
                    plens = [sum(tangle_only_graph.nodes[n]['length'] for n in path) for path in next_next]
                    for i in range(2):
                        end_nodes.append(list(tangle_only_graph.successors(next_next[i][-1])))
                    #length/coverage check?
                    if plens[0] > plens[1] * len_variation or plens[1] > plens[0] * len_variation:
                        continue
                    if len (end_nodes[0]) == 1 and len (end_nodes[1]) == 1 and end_nodes[0][0] == end_nodes[1][0]:
                        deleted = True
                        logging.debug (f"Removing bulge: {', '.join([node_id_mapper.node_id_to_name_safe(n) for n in next_next[0] ])}")
                        for node in next_next[0]:
                            tangle_only_graph.remove_node(node)
                            tangle_only_graph.remove_node(-node)

    #remove long 1-1 chains
    nodes_between_tangles = set()
    simplified_graph = tangle_only_graph.copy()
    for node in tangle_only_graph.nodes():          
        #skip checking midpath  
        if tangle_only_graph.in_degree(node) == 1:
            prev_node = list(tangle_only_graph.predecessors(node))[0]
            if tangle_only_graph.out_degree(prev_node) == 1:
                continue
        cur_path = max_nonbranching_from_node(tangle_only_graph, node)
        cur_len = 0
        for node in cur_path:
            cur_len += tangle_only_graph.nodes[node]['length']
        if cur_len > component_length_threshold or len (cur_path) > path_node_count_threshold:
            logging.debug(f"Clearing bubble chain length {cur_len}")
            logging.debug (f"Bubble nodes: {', '.join([node_id_mapper.node_id_to_name_safe(n) for n in cur_path])}")            
            for node in cur_path:
                nodes_between_tangles.add(node)
                nodes_between_tangles.add(-node)
    for node_id in nodes_between_tangles:
        tangle_only_graph.remove_node(node_id)
    tangle_components = list(nx.weakly_connected_components(tangle_only_graph))

    neighbor_count = {}
    #should do not use simplified or tangle_only graph below!
    processed = set()
    valid_tangle_count = 0
    basedir = BASE_SUBDIR
    for tangle in tangle_components:
        #returning tip clipped/bubble removed nodes to tangle
        #to exclude double processing of rc component
        if list(tangle)[0] in processed:
            continue
        for node in tangle:
            processed.add(node)
            processed.add(-node)
        added = True
        while added:
            added = False
            for node in list(tangle):
                for neighbor in hifi_graph.successors(node):
                    if not neighbor in tangle and not neighbor in nodes_between_tangles:
                        tangle.add(neighbor)
                        added = True
                for neighbor in hifi_graph.predecessors(node):
                    if not neighbor in tangle and not neighbor in nodes_between_tangles:
                        tangle.add(neighbor)
                        added = True
        #Removing long tips from tangle
        tips = set()
        for node in tangle:
            if is_tip(hifi_graph, node) and hifi_graph.nodes[node]['length'] > MAX_TIP_LENGTH:
                tips.add(node)
        for tip in tips:
            tangle.remove(tip)

        if len(tangle) == 1:
            continue
        neighbors = set()
        for node in tangle:
            for neighbor in hifi_graph.successors(node):
                if not (neighbor in tangle):
                    neighbors.add(neighbor)
            for neighbor in hifi_graph.predecessors(node):
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
                logging.info(f"Multiple alignments ( {len(alignments.get(abs_node, []))}) for {node_id_mapper.node_id_to_name_safe(n)}")

                if hifi_graph.out_degree(n) > 1 and not (list(hifi_graph.successors(n))[0] in tangle):
                    if hifi_graph.out_degree(n) == 2:
                        next_nodes = list(hifi_graph.successors(n))
                        logging.info (f"fixing coverage, neighbor {node_id_mapper.node_id_to_name_safe(n)} coverage {coverage_data[abs(n)]} replacing with {', '.join([node_id_mapper.node_id_to_name_safe(x) for x in next_nodes])}")
                        for rn in next_nodes:
                            filtered_neighbors.add(abs(rn))
                    else:
                        logging.info (f"Failed to fix coverage, neighbor {node_id_mapper.node_id_to_name_safe(n)} has out_degree {hifi_graph.out_degree(n)} and in_degree {hifi_graph.in_degree(n)}")
                        valid_tangle = False
                        break
                elif hifi_graph.in_degree(n) > 1 and not (list(hifi_graph.predecessors(n))[0] in tangle):
                    if hifi_graph.in_degree(n) == 2:
                        next_nodes = list(hifi_graph.predecessors(n))
                        logging.info (f"fixing coverage, neighbor {node_id_mapper.node_id_to_name_safe(n)} coverage {coverage_data[abs(n)]} replacing with {', '.join([node_id_mapper.node_id_to_name_safe(x) for x in next_nodes])}")
                        for rn in next_nodes:
                            filtered_neighbors.add(abs(rn))
                    else:
                        logging.info (f"Failed to fix coverage, neighbor {node_id_mapper.node_id_to_name_safe(n)} has in_degree {hifi_graph.in_degree(n)} and in_degree {hifi_graph.in_degree(n)}")
                        valid_tangle = False
                        break
                else:
                    
                    logging.info (f"Failed to fix coverage, neighbor {node_id_mapper.node_id_to_name_safe(n)} has in_degree {hifi_graph.in_degree(n)} and out_degree {hifi_graph.out_degree(n)}")
                    valid_tangle = False
                    break
        local_count = len(filtered_neighbors)
        neighbor_count[len(filtered_neighbors)] = neighbor_count.get(len(filtered_neighbors), 0) + 1
        if valid_tangle:
            cur_dir = os.path.join(basedir, f"tangle_{valid_tangle_count}")
            os.makedirs(cur_dir, exist_ok=True)
            res_file = os.path.join(cur_dir, f"boundary.txt")
            if len(filtered_neighbors) == 2:
                with open(res_file, "w") as f:
                    f.write("\t".join([node_id_mapper.node_id_to_name_safe(x)[1:] for x in filtered_neighbors]) + "\n")
                f.close()
                valid_tangle_count += 1
            elif len(filtered_neighbors) == 4:                  
                nlist = list(filtered_neighbors)
                first = nlist[0]
                found = False
                for i in range (1, 4):
                    if first in alignments and nlist[i] in alignments:
                        if alignments[first] == alignments[nlist[i]] and len(alignments[first]) == 1:
                            logging.info (f"Found first verkko traversal to compare between {node_id_mapper.node_id_to_name_safe(first)} and {node_id_mapper.node_id_to_name_safe(nlist[i])}")
                            first_string = f"{node_id_mapper.node_id_to_name_safe(nlist[0])[1:]}\t{node_id_mapper.node_id_to_name_safe(nlist[i])[1:] }\n"
                            alts = []
                            for j in range (1,4):
                                if j != i:
                                    alts.append(nlist[j])
                            if alts[0] in alignments and alts[1] in alignments:
                                if alignments[alts[0]] == alignments[alts[1]] and len(alignments[alts[0]]) == 1:
                                    logging.info (f"Found second verkko traversal to compare between {node_id_mapper.node_id_to_name_safe(alts[0])} and {node_id_mapper.node_id_to_name_safe(alts[1])}")
                                    second_string = f"{node_id_mapper.node_id_to_name_safe(alts[0])[1:]}\t{node_id_mapper.node_id_to_name_safe(alts[1])[1:]}\n"
                                    with open(res_file, "w") as f:
                                        f.write(first_string)
                                        f.write(second_string)
                                    f.close()
                                    valid_tangle_count += 1
                                    found = True
                                    break                            
                if not found:
                    logging.info (f"Verkko didn't help to scaffold 2 haplo tangle, neighbors: {','.join([node_id_mapper.node_id_to_name_safe(node) for node in filtered_neighbors])}")

            elif len(filtered_neighbors) == 0:
                logging.info (f"Isolated tangle {valid_tangle_count} no incoming/outgoing, doing nothing")
            else:
                logging.info (f"Unexpected number of filtered neighbors: {len(filtered_neighbors)} {','.join([node_id_mapper.node_id_to_name_safe(node) for node in filtered_neighbors])}")
        else:
            logging.info (f"Invalid tangle failed to fix multicopy neighbours {local_count} {','.join([node_id_mapper.node_id_to_name_safe(node) for node in neighbors])}")

        logging.info(f"Tangle component of length {sum(hifi_graph.nodes[node]['length'] for node in tangle)}: total nodes {len(tangle)}, external connections: {len(filtered_neighbors)} {','.join([node_id_mapper.node_id_to_name_safe(node) for node in filtered_neighbors])}.")
        logging.info (",".join([node_id_mapper.node_id_to_name_safe(node) for node in tangle]))        
        logging.info ("\n")
    ssum = 0
    for i in range (0,20):
        logging.info (f"Tangles with {i} connections: {neighbor_count.get(i,0)}")
        ssum += neighbor_count.get(i,0)        
    logging.info(f"Total tangles detected: {ssum}")


def verify_alignment(alignments, node_id_mapper, outdir):
    res_gaf_file = os.path.join(outdir, "traversal.gaf")
    if not os.path.exists(res_gaf_file):
        logging.info(f"No traversal.gaf file in {outdir}, skipping verification")
        return
    logging.info(f"Verifying alignments in {res_gaf_file}")
    alignment_id = 0
    strings_to_output = []
    for_reeval = False
    for line in open(res_gaf_file):
        gaf_str = line.strip().split('\t')[1]
        alignment = src.input_parsing.parse_gaf_string(gaf_str, node_id_mapper)
        borders = [alignment[0], alignment[-1]]
        abs_borders = [abs(b) for b in borders]
        verkko_paths = [alignments.get(abs_borders[0], []), alignments.get(abs_borders[1], [])]
        if not verkko_paths[0] or not verkko_paths[1]:
            logging.info(f"Some of border nodes missing in verkko paths: {','.join([node_id_mapper.node_id_to_name_safe(b) for b in borders])}")
        elif verkko_paths[0] != verkko_paths[1]:
            logging.info(f"Verkko path do not match for border nodes: {','.join([node_id_mapper.node_id_to_name_safe(b) for b in borders])}")
        elif len (verkko_paths[0]) != 1:
            logging.info(f"Border nodes are not unique in Verkko paths: {','.join([node_id_mapper.node_id_to_name_safe(b) for b in borders])}")
        else:
            verkko_res_path = verkko_paths[0][0]
            start_id = 0
            end_id = len(verkko_res_path) - 1
            while start_id < end_id and not (abs(verkko_res_path[start_id]) in abs_borders):
                start_id += 1
            while end_id > start_id and not (abs(verkko_res_path[end_id]) in abs_borders):
                end_id -= 1

            if abs(verkko_res_path[start_id]) in abs_borders and abs(verkko_res_path[end_id]) in abs_borders:
                sub_path = verkko_res_path[start_id:end_id+1]
                if abs(sub_path[0]) == abs(borders[1]):
                    sub_path = [-n for n in reversed(sub_path)]
                has_gaps = False
                for n in sub_path:
                    nodename = node_id_mapper.node_id_to_name_safe(n)
                    if nodename.find("gap") != -1:
                        has_gaps = True
                        break
                verkko_subpath_string = "".join([node_id_mapper.node_id_to_name_safe(n) for n in sub_path])

                if has_gaps:
                    logging.info(f"Verkko's path has gaps for borders {','.join([node_id_mapper.node_id_to_name_safe(b) for b in borders])}, not evaluating")
                    logging.info (f"verkko_path {verkko_subpath_string}")
                    continue
                strings_to_output.append(f"verkko_{alignment_id}\t{verkko_subpath_string}\n")
                if sub_path[0] != borders[0] or sub_path[-1] != borders[1]:
                    logging.info (f"Testing error, orientation do not match borders {','.join([node_id_mapper.node_id_to_name_safe(b) for b in borders])}, verkko_path {''.join([node_id_mapper.node_id_to_name_safe(n) for n in verkko_res_path])}")
                else:                    
                    match = True
                    if len (alignment) != len(sub_path):
                        match = False
                        logging.info(f"Paths do not match in length {len(alignment)} vs {len(sub_path)}, borders {','.join([node_id_mapper.node_id_to_name_safe(b) for b in borders])}, verkko_path {verkko_subpath_string}")                            
                    for i in range(min (len(alignment), len(sub_path))):
                        if alignment[i] != sub_path[i]:
                            match = False
                            logging.info (f"Paths do not match for borders {','.join([node_id_mapper.node_id_to_name_safe(b) for b in borders])}, verkko_path {verkko_subpath_string}")
                            logging.info(f"First divergence at {node_id_mapper.node_id_to_name_safe(alignment[i])} vs {node_id_mapper.node_id_to_name_safe(sub_path[i])}, borders {''.join([node_id_mapper.node_id_to_name_safe(b) for b in borders])}, verkko_path {verkko_subpath_string}")
                            break
                    if match:
                        logging.info (f"Paths match for borders {','.join([node_id_mapper.node_id_to_name_safe(b) for b in borders])}")
                    else:
                        for_reeval = True
            else:
                logging.info (f"Testing error, failed to find borders {','.join([node_id_mapper.node_id_to_name_safe(b) for b in borders])}, in verkko_path correctly {''.join([node_id_mapper.node_id_to_name_safe(n) for n in verkko_res_path])}")    
        alignment_id += 1

    if for_reeval:
        verkko_subpath_gaf = os.path.join(outdir, "verkko_subpath.gaf")
        with open(verkko_subpath_gaf, "w") as verkko_subpath_gaf:
            for line in strings_to_output:
                verkko_subpath_gaf.write(line)


def realign_to_ref(graph, reference_path, outdir):
    verkko_gaf = os.path.join(outdir, "verkko_subpath.gaf")
    verkko_fasta = os.path.join(outdir, "verkko_subpath.fasta")
    if not os.path.exists(verkko_gaf):
        return
    TTT_gaf = os.path.join(outdir, "traversal.gaf")
    TTT_fasta = os.path.join(outdir, "traversal.full.fasta")
    verkko_res = os.path.join(outdir, "verkko.res")
    TTT_res = os.path.join(outdir, "TTT.res")
    # Convert GAF to FASTA
    convert_str = f"python {PATH_2_FASTA_SCRIPT} {graph} {verkko_gaf} > {verkko_fasta}"
    logging.info(f"Converting gaf->fasta {convert_str}")
    os.system(convert_str)


    convert_str = f"python {PATH_2_FASTA_SCRIPT} {graph} {TTT_gaf} > {TTT_fasta}"
    logging.info(f"Converting gaf->fasta {convert_str}")
    os.system(convert_str)

    # Align FASTA files
    minimap_str = f"minimap2 -x asm5 {reference_path} {verkko_fasta} -t 20 --secondary=no -c > {verkko_res}"
    logging.info(f"Aligning with minimap2: {minimap_str}")
    os.system(minimap_str)
    minimap_str = f"minimap2 -x asm5 {reference_path} {TTT_fasta} -t 20 --secondary=no -c > {TTT_res}"
    logging.info(f"Aligning with minimap2: {minimap_str}")
    os.system(minimap_str)
    
    #align_fasta_files(verkko_fasta, reference_path, verkko_res)
    #align_fasta_files(TTT_fasta, reference_path, TTT_res)
    

def verify_alignments(alignments, node_id_mapper, out_name):
    for dir in os.listdir(BASE_SUBDIR):
        outdir = os.path.join(os.getcwd(), BASE_SUBDIR, dir, out_name)
        verify_alignment(alignments, node_id_mapper, outdir)

def realign_all_to_ref(graph, reference_path, out_name):
    for dir in os.listdir(BASE_SUBDIR):
        logging.info(f"Processing tangle: {dir}")
        outdir = os.path.join(os.getcwd(), BASE_SUBDIR, dir, out_name)
        realign_to_ref(graph, reference_path, outdir)

def parse_aligns(align_file):
    aligns = {}
    for line in open(align_file):
        fields = line.strip().split("\t")
        query_name = fields[0].split('_')[1]
        if query_name in aligns.keys():
            aligns[query_name].append("broken")
        else:
            aligns[query_name] = fields
    return aligns

def extract_all_res(out_name):
    print (f"tangle\tverkko_chr\tverkko_start\tverkko_end\tscore_diff\tbroken alignment")
    for dir in os.listdir(BASE_SUBDIR):
        outdir = os.path.join(os.getcwd(), BASE_SUBDIR, dir, out_name)
        verkko_res = os.path.join(outdir, "verkko.res")
        TTT_res = os.path.join(outdir, "TTT.res")
        if not os.path.exists(verkko_res) or not os.path.exists(TTT_res):
            continue
        TTT_aligns = parse_aligns(TTT_res)
        verkko_aligns = parse_aligns(verkko_res)
        for al_id in TTT_aligns:
            if al_id not in verkko_aligns:
                continue
            # Compare alignments
            TTT_al = TTT_aligns[al_id]
            verkko_al = verkko_aligns[al_id]
            score_diff = int(TTT_al[14].split(":")[2]) - int(verkko_al[14].split(":")[2])
            broken = TTT_al[-1] == "broken" or verkko_al[-1] == "broken"
            out_line = f"{dir}_{al_id}\t{verkko_al[5]}\t{verkko_al[7]}\t{verkko_al[8]}\t{score_diff}\t{broken}"
            print(out_line)


def parse_arguments():
    #currently just for logging compatibilituy  
    parser = argparse.ArgumentParser(description="Test comparison against verkko")    
    args = parser.parse_args()
    args.log_level = "INFO"
    args.outdir = "."
    args.basename = "extract_tangle"
    return args


if __name__ == "__main__":
    #basedir = sys.argv[1]
    args = parse_arguments()
    src.logging_utils.setup_logging(args)
    basedir = "/data/antipovd2/verkko2_paper/HG002/v2.2.1_hic/"
    graph = os.path.join(basedir, "2-processGraph", "unitig-unrolled-hifi-resolved.noseq.gfa")
    graph_seq = os.path.join(basedir, "2-processGraph", "unitig-unrolled-hifi-resolved.gfa")
    coverage = os.path.join(basedir, "2-processGraph", "unitig-unrolled-hifi-resolved.ont-coverage.csv")
    alignment_file = os.path.join(basedir, "node_realign", "scaffolds2utig1.gaf")
    node_id_mapper = NodeIdMapper()
    alignments = parse_alignment_file(alignment_file, node_id_mapper)
    
    run_id = 4
    
    #get_tangle_components(graph, coverage, alignments, node_id_mapper)

    if False:
        for dir in os.listdir(BASE_SUBDIR):                
            boundary_file = os.path.join(os.getcwd(), BASE_SUBDIR, dir, "boundary.txt")
            if os.path.exists(boundary_file):
                outdir = os.path.join(os.getcwd(), BASE_SUBDIR, dir, f"TTT_{run_id}")
                script_file = os.path.join(BASE_SUBDIR, dir, "run_TTT.sh")                
                TTT_runstr = f"/data/antipovd2/devel/TTT/TTT.py --verkko-output /data/antipovd2/verkko2_paper/HG002/v2.2.1_hic/ --boundary-nodes {boundary_file} --outdir {outdir}"
                with open(script_file, "w") as f:
                    f.write("#!/bin/bash\n")                
                    f.write(f"{TTT_runstr}\n")
                os.system(f"chmod +x {script_file}")
                os.system(f"sbatch --time=3:00:00 --mem=40g --cpus-per-task=28 --partition=norm,quick {script_file}")
    
    verify_alignments(alignments, node_id_mapper, f"TTT_{run_id}")
    realign_all_to_ref(graph_seq, REFERENCE, f"TTT_{run_id}")
    extract_all_res(f"TTT_{run_id}")