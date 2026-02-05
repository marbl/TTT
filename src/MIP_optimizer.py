import logging
import pulp
import networkx as nx
from .logging_utils import log_assert
from .node_id_mapper import NodeIdMapper


class MIPOptimizer:

    ALWAYS_INCLUDE_COVERAGE_FRACTION = 0.7

    def __init__(self, node_mapper, time_limit):
        """
        Initialize the MIP optimizer with a node mapper.
        
        Args:
            node_mapper: NodeIdMapper instance for converting between node IDs and names
        """
        self.node_mapper = node_mapper
        self.time_limit = time_limit
    
    def normalized_score(self, solution, coverages, lengths, median_unique):
        total_deviation = 0
        total_length = 0
        for node_id in coverages:
            if node_id in solution and node_id in lengths:
                mult = solution[node_id] + solution.get(-node_id, 0)
                deviation = abs(mult * median_unique - coverages[node_id])
                length = lengths.get(node_id, 0)
                total_deviation += deviation * length
                total_length += length
                logging.debug(f"Node {self.node_mapper.node_id_to_name_safe(node_id)}: Deviation = {deviation}, Length = {length}")
        if total_length == 0:
            return float('inf')
        return total_deviation / (total_length * median_unique)

    def solve_MIP(self, equations, nonzeros, boundary_values, coverages, unique_coverage_range , original_graph):
        # last element in equation - amount of boundary nodes    
        #nonzeros = set()
        while True:
            prob = pulp.LpProblem("Minimize_Deviation", pulp.LpMinimize)    
            # Define multiplicity variables    

            keys = set()
            abs_keys = set()
            for eq in equations:
                keys.update(eq[0])
                keys.update(eq[1])
                for ind in (0, 1):
                    for j in eq[ind]:
                        abs_keys.add(abs(j))
            for key in abs_keys:
                log_assert(key in keys and -key in keys, 
                        f"One of the keys {self.node_mapper.node_id_to_name_safe(key)} or {self.node_mapper.node_id_to_name_safe(-key)} is not in equations {equations}")
                log_assert(key in coverages, 
                        f"Key {self.node_mapper.node_id_to_name_safe(key)} is not in coverages")

            # TODO: we can check for hairpin presence just here.
            x_vars = {i: pulp.LpVariable(f"x_{i}", cat="Integer") for i in keys}

            # Constraints: x_i = x_j + x_k + x_l
            # for lhs, rhs in equations.items():
            # prob += x_vars[lhs] == sum(x_vars[j] for j in rhs)
            for eq in equations:
                prob += sum(x_vars[j] for j in eq[0]) == sum(x_vars[j] for j in eq[1]) 

            for x_var in x_vars:
                prob += x_vars[x_var] >= 0
                if x_var > 0 and x_var in nonzeros:
                    prob += x_vars[x_var] + x_vars[-x_var] >= 1

            # boundary nodes set
            for x_var in boundary_values:
                prob += x_vars[x_var] == boundary_values[x_var]

            # Define continuous variables for absolute deviation |x_i - a_i|
            d_vars = {i: pulp.LpVariable(f"d_{i}", lowBound=0, cat="Continuous") for i in coverages}   
            lengths = {i: original_graph.nodes[i].get('length', 0) for i in coverages}
            #do not want boundary nodes to contribute to the objective, coverage is usually a bit lower within the tangle
            for i in boundary_values:
                if i in lengths:
                    lengths[i] = 0

            inv_unique_coverage = pulp.LpVariable("uniqueness_multiplier", 
                                                lowBound=1/unique_coverage_range[1], 
                                                upBound=1/unique_coverage_range[0], 
                                                cat="Continuous")
            for i in abs_keys:
                prob += d_vars[i] >= (x_vars[i] + x_vars[-i]) - coverages[i] * inv_unique_coverage
                prob += d_vars[i] >= coverages[i] * inv_unique_coverage - x_vars[i] - x_vars[-i]

            # Objective function: minimize sum of absolute deviations weighted by node lengths
            objective = pulp.lpSum(d_vars[i] * lengths[i] for i in abs_keys)
            #objective = pulp.lpSum(d_vars[i] for i in abs_keys)
            prob += objective

            logging.debug("Linear programming problem:")
            logging.debug(f"Objective: {prob.objective}")
            for constraint in prob.constraints.values():
                logging.debug(f"Constraint: {constraint}")
            for var in prob.variables():
                logging.debug(f"Variable: {var.name}, LowBound: {var.lowBound}, UpBound: {var.upBound}, Cat: {var.cat}")

            # MIP magic
            pulp.LpSolverDefault.msg = 1
            # No more than 1 hour
            prob.solve(pulp.GLPK(timeLimit=self.time_limit))
            
            result = {}    
            
            if pulp.LpStatus[prob.status] != "Optimal":
                logging.warning("MIP did not find an optimal solution.")
                if len (nonzeros) > 0:
                    #TODO: think about iterative removal strategy
                    #logging.warning (f"Removing the constraint that forced to include {self.node_mapper.node_id_to_name_safe(nonzeros[-1])} cov {coverages[nonzeros[-1]]} into traversing paths")
                    #nonzeros.pop()
                    logging.warning("Removing all forced inclusions and retrying.")
                    nonzeros = []
                    continue
                else:
                    break
            else:
                result = {i: int(pulp.value(x_vars[i])) for i in keys}
                
                # Convert result to use node names for logging
                result_with_names = {self.node_mapper.node_id_to_name_safe(i): result[i] for i in result}
                logging.info(f"Found multiplicities from MIP {result_with_names}")
                score = pulp.value(objective)  
                detected_coverage = 1/pulp.value(inv_unique_coverage)
                logging.info(f"Unique coverage for the best MIP solution {detected_coverage}, solution score {score}")
                if abs(detected_coverage - unique_coverage_range[0]) < 0.01 or abs(detected_coverage - unique_coverage_range[1]) < 0.01:
                    logging.info(f"Warning, detected best coverage is close to the allowed unique coverage borders{unique_coverage_range}")  
                normalized_score = self.normalized_score(result, coverages, lengths, detected_coverage)
                logging.info(f"Normalized score of the MIP solution: {normalized_score}")                                                                          
                return result, normalized_score
        logging.error("MIP did not find an optimal solution after all adjustment attempts.")
        exit(1)


    def generate_MIP_equations(self, tangle_nodes, nor_nodes, cov, median_unique, original_graph, boundary_nodes, directed=False):

        junction_equations = []    
        used = set()
        extended_tangle = tangle_nodes.copy()
        # format: used ins to corresponding used outs
        all_boundary_nodes = set()
        for b in boundary_nodes:
            all_boundary_nodes.add(abs(b))
            all_boundary_nodes.add(-abs(b))
            all_boundary_nodes.add(abs(boundary_nodes[b]))
            all_boundary_nodes.add(-abs(boundary_nodes[b]))
        extended_tangle.update(all_boundary_nodes)
        canonic_name = {}
        for node in extended_tangle:
            if directed:
                canonic_name[node] = node
            else:
                canonic_name[node] = abs(node)
        
        # Generate equations
        for from_node in extended_tangle:
            if from_node in used:
                continue
            arr = [[], []]
            back_node = ""
            bad_extension = 0
            for to_node in original_graph.successors(from_node):
                back_node = to_node
                if not directed:
                    used.add(-to_node)
                if to_node in extended_tangle:
                    arr[1].append(canonic_name[to_node])
                else:
                    bad_extension = to_node
                    break
               
            if bad_extension != 0:
                if not (from_node in all_boundary_nodes):
                    logging.error(f"Somehow jumped over boundary nodes {self.node_mapper.node_id_to_name_safe(from_node)} {self.node_mapper.node_id_to_name_safe(back_node)} {boundary_nodes}")
                    exit()
                else:
                    continue
            if back_node != "":
                for alt_start_node in original_graph.predecessors(back_node):
                    used.add(alt_start_node)
                    if alt_start_node in extended_tangle:
                        arr[0].append(canonic_name[alt_start_node])
                    else:
                        bad_extension = alt_start_node
                        break
            if bad_extension != 0:
                logging.error(f"Somehow jumped over boundary nodes (backwards) {self.node_mapper.node_id_to_name_safe(alt_start_node)} {self.node_mapper.node_id_to_name_safe(bad_extension)} {boundary_nodes}")
                exit()
            junction_equations.append(arr)

        
        must_use_nodes = []
        coverage = {}
        for node in nor_nodes:        
            coverage[node] = float(cov[node])  # / median_unique
            logging.debug(f"Coverage of {self.node_mapper.node_id_to_name_safe(node)} : {coverage[node]}")
            #if (coverage[node] / median_unique >= 0.3 and (coverage[node] >= 5)) or coverage[node]/median_unique >= 0.6:
            if (coverage[node] / median_unique >= self.ALWAYS_INCLUDE_COVERAGE_FRACTION):
                must_use_nodes.append((node, coverage[node]))
        must_use_nodes = sorted(must_use_nodes, key=lambda x: x[1], reverse=True)
        nonzeroes = [node for node, cov in must_use_nodes]
        boundary_values = {}
        for b in boundary_nodes:
            boundary_values[b] = 1
            boundary_values[boundary_nodes[b]] = 1
            boundary_values[-b] = 0
            boundary_values[-boundary_nodes[b]] = 0
            coverage[abs(b)] = float(cov[abs(b)]) 
            coverage[abs(boundary_nodes[b])] = float(cov[abs(boundary_nodes[b])])
        return junction_equations, nonzeroes, coverage, boundary_values
    
    def write_multiplicities(self, output_file, solution, cov):
        with open(output_file, 'w') as out_file:
            out_file.write("node\tcoverage\tmult\n")
            if len(solution) > 0:
                for node_id in cov.keys():                    
                    mult_value = solution.get(node_id, "X")
                    rev_mult_value = solution.get(-node_id, "X")
                    if mult_value != "X" and rev_mult_value != "X":
                        mult_value = int(mult_value) + int(rev_mult_value)
                    cov_value = cov.get(node_id, "N/A")

                    out_file.write(f"{self.node_mapper.node_id_to_unoriented_name(node_id)}\t{cov_value}\t{mult_value}\n")
        logging.info(f"Wrote multiplicity solutions to {output_file}")
