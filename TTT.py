#!/usr/bin/env python3
import sys
import argparse
import logging
import math
import os
from src.path_optimizer import PathOptimizer
from src.path_supplementary import get_gaf_string
from src.alignment_scorer import AlignmentScorer
from src.logging_utils import setup_logging, print_warnings_summary
from src.node_id_mapper import NodeIdMapper
from src.MIP_optimizer import MIPOptimizer
from src.tangle import Tangle
from src.input_parsing import (
    parse_gfa, parse_gaf, read_tangle_nodes, get_oriented_boundaries, identify_tangle_nodes, new_identify_tangle_nodes, read_coverage_file,
    coverage_from_graph, verify_coverage, calculate_median_coverage, clean_tips, DETECTED_LOW_MEDIAN_COVERAGE_VARIATION, DETECTED_HIGH_MEDIAN_COVERAGE_VARIATION
)
from src.graph_transformation import (
    get_canonical_rc_vertex,
    create_dual_graph, create_multi_dual_graph,
    get_traversable_subgraph
)

# Create global instances
node_id_mapper = NodeIdMapper()

def optimize_paths(tangle, alignment_scorer: AlignmentScorer, args):
    """Optimize Eulerian paths."""
    best_path = None
    best_score = -1
    logging.info(f"Starting optimization with {args.num_initial_paths} initial paths, max {args.max_iterations} iterations per path.")
    rc_vertex_map = {}
    for vertex in tangle.multi_graph.nodes():
        rc_vertex_map[vertex] = get_canonical_rc_vertex(vertex, tangle.original_graph, tangle.node_id_mapper)
    subgraph_to_traverse, start_vertex = get_traversable_subgraph(tangle)    

    if not subgraph_to_traverse:
        logging.error(f"No Eulerian path found.")        
        sys.exit(1)
    pathOptimizer = PathOptimizer(subgraph_to_traverse, start_vertex, tangle.node_id_mapper, rc_vertex_map)

    for seed in range(args.num_initial_paths):
        logging.info(f"Generating initial path with seed {seed}.")
        pathOptimizer.generate_random_eulerian_path(seed)
        #TODO: most of this logic should go to PathOptimizer
        current_path = pathOptimizer.get_path()
        current_score = alignment_scorer.score_corasick(current_path)
        logging.info(f"Initial path score for seed {seed}: {current_score}.")

        iterations_since_improvement = 0
        for i in range(args.max_iterations):
            if iterations_since_improvement >= args.early_stopping_limit:
                logging.info(f"Early stopping for seed {seed} on iteration {i} after {iterations_since_improvement} iterations without improvement.")
                break
            
            new_path = pathOptimizer.get_random_change(seed * args.max_iterations + i)
            if new_path is None:
                logging.info(f"Not able to change a path at all, stopping")
                break
            new_score = alignment_scorer.score_corasick(new_path)
            if new_score > current_score:
                logging.info(f"Improved score for seed {seed} at iteration {i}: {current_score} -> {new_score}.")
                current_path = new_path
                pathOptimizer.set_path(new_path)
                current_score = new_score
                iterations_since_improvement = 0
            else:
                iterations_since_improvement += 1
        
        logging.info(f"Final score for seed {seed}: {current_score}.")
        if current_score > best_score:
            logging.info(f"New best path found for seed {seed} with score {current_score}.")
            logging.info(f"Path: {get_gaf_string(current_path, node_id_mapper)}")
            best_path = current_path
            best_score = current_score
    
    logging.info(f"Optimization completed. Best score: {best_score}.")        
    if best_path:
        logging.info("Path optimizing finished.")
        return best_path, best_score, pathOptimizer        
    else:
        logging.error("No valid path found during optimization.")
        sys.exit(1)

def print_final_path_info(best_path, pathOptimizer, tangle, alignment_scorer: AlignmentScorer, args):
    output_fasta = os.path.join(args.outdir, args.basename + ".hpc.fasta")
    output_gaf = os.path.join(args.outdir, args.basename + ".gaf")
    
    pathOptimizer.set_path(best_path)
    pathOptimizer.get_synonymous_changes(alignment_scorer)
    best_path_str = get_gaf_string(best_path, tangle.node_id_mapper)
    logging.info(f"Found traversal\t{best_path_str}")   
    # Output FASTA file
    logging.info(f"Writing best path to {output_fasta} and gaf to {output_gaf}")
    pathOptimizer.output_path(tangle.original_graph, output_fasta, output_gaf)
    alignment_scorer.not_satisfied_fraction(best_path)

    coverage_fraction = MIPOptimizer.ALWAYS_INCLUDE_COVERAGE_FRACTION

    cutoff = coverage_fraction * tangle.detected_coverage
    path_abs_ids = {abs(edge.original_node) for edge in best_path}
    tips_abs_ids = {abs(node) for node in tangle.cleaned_tips}
    non_tips_abs_ids = {abs(node) for node in tangle.nodes} - tips_abs_ids

    missing_tip_ids = {
        node_id for node_id in tips_abs_ids
        if tangle.coverage_dict.get(node_id, 0) >= cutoff and node_id not in path_abs_ids
    }
    missing_non_tip_ids = {
        node_id for node_id in non_tips_abs_ids
        if tangle.coverage_dict.get(node_id, 0) >= cutoff and node_id not in path_abs_ids
    }

    missing_tips_count = len(missing_tip_ids)
    missing_non_tips_count = len(missing_non_tip_ids)
    missing_abs_ids = sorted(missing_tip_ids | missing_non_tip_ids)

    missing_total = missing_tips_count + missing_non_tips_count
    if missing_total > 0:
        missing_names = ", ".join([tangle.node_id_mapper.node_id_to_unoriented_name(node_id) for node_id in missing_abs_ids])
        logging.warning(
            f"Missing high-coverage edges in final path: total={missing_total} (tips={missing_tips_count}, non_tips={missing_non_tips_count}), "
            f"threshold={coverage_fraction * tangle.detected_coverage:.2f}, edges={missing_names}"            
        )
        logging.warning(f"TTT was not able to include high-covered tips or some tip-like structures in the path because of the graph structure.")
    else:
        logging.info(
            f"No missing high-coverage edges in the final path, "
            f"threshold={coverage_fraction * tangle.detected_coverage:.2f}"
        )
        
    logging.info(f"Detected coverage for the path: {tangle.detected_coverage}")
    logging.info(f"Final path: {get_gaf_string(best_path, tangle.node_id_mapper)}")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Solve for integer multiplicities in a GFA tangle graph based on coverage.")
    parser.add_argument("--graph", required=False, help="Path to the GFA graph.")
    parser.add_argument("--alignment", required=False, help="Path to a file with graphaligner alignment")
    parser.add_argument("--outdir", required=True, type=str, help="Output directory for all result files (will be created if it doesn't exist)")
    parser.add_argument("--verkko-output", required=False, type=str, help="Path to dir with verkko results. Ovewrites --gfa, --alignment, --coverage-file with standart paths to hifi(utig1-) verkko's graph")
    parser.add_argument("--boundary-nodes", required=True, type=str, help="Path to a file listing boundary node pairs, tab-separated (required for 2-2 tangles).")

    parser.add_argument("--coverage", help="Path to a file with node coverages (verkko's format; newline separated pairs node_id coverage). If not provided, coverage will be filled from the GFA file.")
    parser.add_argument("--median-unique", type=float, help="Median coverage for reliable unique nodes in tangle. If not provided, autodetected.")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Set the logging level for log file (default: INFO).")

    parser.add_argument("--num-initial-paths", type=int, default=10, help="Number of initial paths to generate (default: 10).")
    parser.add_argument("--max-iterations", type=int, default=100000, help="Maximum iterations for path optimization (default: 100000).")
    #TODO: for quality 0 use random of alignments with highest score?
    parser.add_argument("--early-stopping-limit", type=int, default=15000, help="Early stopping limit for optimization (default: 15000).")
    parser.add_argument("--quality-threshold", type=int, default=20, help="Alignments with quality less than this will be filtered out, default 20")
    parser.add_argument("--basename", required=False, default="traversal", type=str, help="Basename for most of the output files, default `traversal`")
    parser.add_argument("--milp-time-limit", type=int, default=7200, help="Time limit for MILP solver in seconds (default: 7200 seconds = 2 hour).")
    args = parser.parse_args()
    if not args.verkko_output  and (not args.graph  or not args.alignment ):
        #logging not initialized yet
        sys.stderr.write("Either --verkko-output OR both --graph and --alignment are required options\n")
        return
    
    if args.verkko_output: 
        args.coverage = os.path.join(args.verkko_output, "2-processGraph", "unitig-unrolled-hifi-resolved.ont-coverage.csv")
        args.alt_coverage = os.path.join(args.verkko_output, "2-processGraph", "unitig-unrolled-hifi-resolved.hifi-coverage.csv")
        args.graph = os.path.join(args.verkko_output, "2-processGraph", "unitig-unrolled-hifi-resolved.gfa")
        args.alignment = os.path.join(args.verkko_output, "3-align", "alns-ont.gaf")
    else: 
        args.alt_coverage = None
    return args

def main():    
    args = parse_arguments()
    os.makedirs(args.outdir, exist_ok=True)
    setup_logging(args)
    logging.debug(f"args: {args}")
    logging.info("Reading files...")

    #TODO: separate function to clean Z connection, save them somewhere
    original_graph = parse_gfa(args.graph, node_id_mapper)
    if args.coverage:
        cov = read_coverage_file(args.coverage, node_id_mapper)
    else:
        cov = coverage_from_graph(original_graph)

    if args.alt_coverage:
        alt_cov = read_coverage_file(args.alt_coverage, node_id_mapper)
        verify_coverage(alt_cov, original_graph, node_id_mapper)
    # Verifying that coverage matches the graph
    # For rare cases coverage may be missing for some nodes, will update with median then 
    verify_coverage(cov, original_graph, node_id_mapper)
    #boundary nodes: map from incoming to outgoing

    # Here all sequence is stored in edges, junctions are new vertices
    dual_graph = create_dual_graph(original_graph, node_id_mapper)

    tangle_nodes, boundary_nodes = new_identify_tangle_nodes(args, original_graph, dual_graph, node_id_mapper)
    
    nor_nodes = {abs(node) for node in tangle_nodes}
    
    # Create Tangle object with core structural information
    tangle = Tangle(
        nodes=tangle_nodes,
        nor_nodes=nor_nodes,
        boundary_nodes=boundary_nodes,
        original_graph=original_graph,
        dual_graph=dual_graph,
        node_id_mapper=node_id_mapper
    )
    
    #TODO: do all operations on dual graph        
    tangle.cleaned_tips = clean_tips(tangle, node_id_mapper)
    tangle.dual_graph = create_dual_graph(original_graph, node_id_mapper)
    
    used_or_nodes = tangle.nodes.copy()
    for b in tangle.boundary_nodes:
        used_or_nodes.add(b)
        used_or_nodes.add(tangle.boundary_nodes[b])
        used_or_nodes.add(-b)
        used_or_nodes.add(-tangle.boundary_nodes[b])
    
    tangle.coverage_dict = cov
    tangle.coverage_range = calculate_median_coverage(tangle, args)
    median_unique = (tangle.coverage_range[0] * DETECTED_LOW_MEDIAN_COVERAGE_VARIATION)
    tangle.median_unique_coverage = median_unique
    
    filtered_alignment_file = os.path.join(args.outdir, f"{args.basename}.q{args.quality_threshold}.used_alignments.gaf")
    alignments = parse_gaf(args.alignment, used_or_nodes, filtered_alignment_file, args.quality_threshold, node_id_mapper)
    alignment_scorer = AlignmentScorer(alignments, original_graph, node_id_mapper)
    logging.info("Starting multiplicity counting...")
    # TODO: instead of a_values we just use coverage
    mip_optimizer = MIPOptimizer(node_id_mapper, time_limit=args.milp_time_limit)
    equations, nonzeros, a_values, boundary_values = mip_optimizer.generate_MIP_equations(tangle, directed=True)
    best_solution, score, detected_coverage = mip_optimizer.solve_MIP(equations, nonzeros, boundary_values, a_values, tangle.coverage_range, original_graph)
    if score > 0.2 and args.alt_coverage != None:
        logging.warning(f"High coverage inconsistency score {score}, trying to re-calculate coverage from alternative (hifi) coverage file {args.alt_coverage}")
        # Temporarily update tangle with alternative coverage
        tangle.coverage_dict = alt_cov
        tangle.coverage_range = calculate_median_coverage(tangle, args)
        median_unique = (tangle.coverage_range[0] * DETECTED_LOW_MEDIAN_COVERAGE_VARIATION)
        tangle.median_unique_coverage = median_unique
        equations, nonzeros, a_values, boundary_values = mip_optimizer.generate_MIP_equations(tangle, directed=True)
        best_solution_hifi, score_hifi, detected_coverage_hifi = mip_optimizer.solve_MIP(equations, nonzeros, boundary_values, a_values, tangle.coverage_range, original_graph)
        if score_hifi * 1.2 < score:
            logging.warning(f"Alternative coverage produced significantly better MIP score {score_hifi} vs {score}, using it")
            best_solution = best_solution_hifi.copy()
            score = score_hifi
            detected_coverage = detected_coverage_hifi
        else:
            logging.warning(f"Alternative coverage did not significantly improve MIP score {score_hifi} vs {score}, keeping original (ONT)")
            # Restore original coverage
            tangle.coverage_dict = cov
            tangle.coverage_range = calculate_median_coverage(tangle, args)
            tangle.median_unique_coverage = (tangle.coverage_range[0] * DETECTED_LOW_MEDIAN_COVERAGE_VARIATION)
    # Store MIP solution in tangle
    tangle.multiplicities = best_solution
    tangle.detected_coverage = detected_coverage
    
    # Define output filename for multiplicities CSV
    output_csv = os.path.join(args.outdir, args.basename + ".multiplicities.csv")
    mip_optimizer.write_multiplicities(output_csv, best_solution, cov)

    #Splitting edges according to multiplicities
    tangle.multi_graph = create_multi_dual_graph(tangle)

    best_path, best_score, pathOptimizer = optimize_paths(tangle, alignment_scorer, args)
    print_final_path_info(best_path, pathOptimizer, tangle, alignment_scorer, args)
    print_warnings_summary(args)

if __name__ == "__main__":
    main()
