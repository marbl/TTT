#!/usr/bin/env python3
import logging
import ahocorasick


class AlignmentScorer:
    #TODO: deprioritize based on lengths and not just node counts
    #With othervise equivaltent solution we belive more in one with smaller inversion
    DEPRIORITIZE_RC_COEFFICIENT = 0.9

    #banning Z-connection-based links
    BANNED_Z_WEIGHT = -1000000

    #Small malus for diploid tangles and nodes present in only one haplotype with multiplicity > 1
    #R in hap1 and RR in hap2 is better than nothing in hap1 and RRR in hap2
    DEPRIORITIZE_ASSYMETRIC = -1
    def __init__(self, alignments, original_graph, node_id_mapper):
        
        self.automaton = ahocorasick.Automaton()
        self.pattern_counts = {}
        self.node_id_mapper = node_id_mapper
        used_nodes = set()
        #from lexicographical minimum of (pattern, rc_pattern) to max
        self.rc_patterns = {}
        for idx, alignment in enumerate(alignments):
            #logging.debug(f"Adding alignment {idx}: {alignment}")
            pattern_str = self.aln_to_string(alignment)
            for n in alignment: 
                used_nodes.add(n)
                used_nodes.add(-n)
            rc_nodes = [-n for n in alignment]
            rc_nodes.reverse()
            rc_pattern_str = self.aln_to_string(rc_nodes)
            
            #TOTHINK Possibly use lexicographical minimum of pattern and rc_pattern?

            if pattern_str not in self.pattern_counts:
                self.pattern_counts[pattern_str] = 0
                self.pattern_counts[rc_pattern_str] = 0
                self.automaton.add_word(pattern_str, pattern_str)
                self.automaton.add_word(rc_pattern_str, rc_pattern_str)     
                self.rc_patterns[pattern_str] = rc_pattern_str
                self.rc_patterns[rc_pattern_str] = pattern_str   
            self.pattern_counts[pattern_str] += 1
            self.pattern_counts[rc_pattern_str] += 1
        logging.info(f"{len(self.pattern_counts)} different alignment pattern used")
        total_banned_z = 0
        for u, v, data in original_graph.edges(data=True):
            #Z connection "overlaps" are negative
            if data['overlap'] < 0 and (u in used_nodes or v in used_nodes):
                logging.debug(f"banning connection between {u} and {v} because of Z vertice")
                total_banned_z += 1
                string_to_ban = self.aln_to_string([u,v])
                self.pattern_counts[string_to_ban] = self.BANNED_Z_WEIGHT
                self.automaton.add_word(string_to_ban, string_to_ban)
                #RC added from graph the same way
        logging.info (f"Banning {total_banned_z} Z connections")
        for pattern in self.pattern_counts:
            logging.debug(f"Pattern: {pattern} Count: {self.pattern_counts[pattern]}")
        self.automaton.make_automaton()
        logging.debug(f"automaton keys {list(self.automaton.keys())}")
        logging.info(f"Built automaton with {len(self.pattern_counts)} unique alignment patterns")

    def path_to_string(self, path):
        return "," + ",".join(str(edge.original_node) for edge in path) + ","
    
    def aln_to_string(self, aln):
        return "," + ",".join(str(node) for node in aln) + ","
    
    def score_corasick(self, path):
        path_str = self.path_to_string(path)
        found = set()
        for item in self.automaton.iter(path_str):
            pattern = item[1]
            found.add(pattern)
        score = 0
        for item in found:
            #only looking for one of (pattern, rc_pattern)
            if self.rc_patterns[item] in found:
                score += self.pattern_counts[item] * self.DEPRIORITIZE_RC_COEFFICIENT
                #TODO: possibly add paths length 1 for better deprioritization?
                logging.debug(f"deprioritizing {item} because of rc")
            else:
                score += self.pattern_counts[item] * 2
        #diploid tangles; RR + R is better than RRR + 0
        aux_pos = -1
        for idx in range (0, len(path)):
            if self.node_id_mapper.node_id_to_name_safe(path[idx].original_node) == "AUX":
                aux_pos = idx
                break
        if aux_pos != -1:
            edges = [{}, {}]
            for idx in range (len(path)):
                e = abs(path[idx].original_node)
                if idx < aux_pos:
                    edge_idx = 0
                else:
                    edge_idx = 1
                if not (e in edges[edge_idx]):
                    edges[edge_idx][e] = 0
                edges[edge_idx][e] += 1
            assymetric_nodes = 0
            for idx in range (2):
                for e in edges[idx]:
                    if edges[idx][e] > 1 and not (e in edges[idx - 1]):
                        assymetric_nodes += 1
            score += assymetric_nodes * self.DEPRIORITIZE_ASSYMETRIC
            logging.debug(f"Deprioritizing {assymetric_nodes} assymetric nodes")
        return score
    
    def not_satisfied_fraction(self, path):
        path_str = self.path_to_string(path)
        found = set()
        for item in self.automaton.iter(path_str):
            pattern = item[1]
            found.add(pattern)
        not_satisfied = 0
        total_reads = 0
        found_reads = 0
        used = set()
        for pattern in self.pattern_counts:
            if pattern in used:
                continue
            #not includiong banned Z connections
            if self.pattern_counts[pattern] < 0:
                continue
            rc_pattern = self.rc_patterns[pattern]
            if (pattern in found) or (rc_pattern in found):
                found_reads += self.pattern_counts[pattern]
            total_reads += self.pattern_counts[pattern]
            used.add(pattern)
            used.add(rc_pattern)
        logging.info(f"Not found reads fraction : {(total_reads - found_reads) / total_reads if total_reads > 0 else 0}")
