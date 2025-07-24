from typing import Any
import re,argparse, logging, sys, copy
from cocopaths import __version__
from .utils import (pairtablepath_to_dotbracket,path_to_pairtablepath, find_connected_modules,is_balanced_structure,acfp_terminal_input, structure_to_pairtable, is_legal_dotbracket)
import string

#______define_logger______#
logger = logging.getLogger('cocopaths')
console_handler = logging.StreamHandler()
formatter = logging.Formatter('# %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


class node():
    def __init__(self,name,prefix ="",suffix ="",complement= False,max_weight = 1,neighbors = None,connected = None):
        self.name = name
        self.prefix = prefix
        self.middle = ""
        self.suffix = suffix
        self.complement = complement
        self.max_weight = max_weight
        self.connected = connected
        self.neighbors = set() if neighbors is None else neighbors
        
    def __str__(self):
        neighbor_names = ', '.join(str(neighbor.name) for neighbor in self.neighbors)
        return f"Name: {self.name} Sequence: {self.prefix}{self.middle}{self.suffix} \nneighbors: {neighbor_names}, Complement: {self.complement}, Connected: {self.connected}, max_weight: {self.max_weight}\n"
    
    
class graph():
    def __init__(self):
        self.graph = {}
        self.edges = dict()
        self.edge_neighbors = dict()

    # -------- Graph creation methods --------
    def add_node(self, node):
         if node not in self.graph:
              self.graph[node] = []

    def add_edge(self, node1, node2):    
        if node2 not in node1.neighbors:
            node1.neighbors.add(node2)  
            node2.neighbors.add(node1)  

    def create_nodes_from_pairtable(self, pairtable):
        for x in range(1, pairtable[-1][0] + 1):
            new_node = node(x)
            self.add_node(new_node)

    def print_nodes(self):
        return "\n".join([str(node) for node in self.graph])

    def create_edges(self, acfp, nodes):
        for step in acfp[::-1]:
            for x in range(1, step[0]):
                if step[x] != 0:
                    node1 = nodes[x]
                    node2 = nodes[step[x]]
                    self.add_edge(node1, node2)
                else:
                    continue
    
    def get_edges(self):
        for node in self.graph:
            for neighbor in node.neighbors:
                self.edges[(node.name, neighbor.name)] = 1  # default weight of edge is 1

        # Remove duplicate edges.
        unique_edges = {}
        for edge in self.edges:
            if edge[1] >= edge[0]:
                unique_edges[edge] = self.edges[edge]
        self.edges = unique_edges

    def set_edge_neighbors(self):
        for edge in self.edges.keys():
            self.edge_neighbors[edge] = list()
            for neighbor in self.edges.keys():
                if neighbor[0] in edge or neighbor[1] in edge:
                    if neighbor != edge:
                        self.edge_neighbors[edge].append(neighbor)
        logger.debug(f"Edge Neighbors: {self.edge_neighbors}")

    # -------- Utility Methods --------
    def get_node_by_name(self, node_name: int) -> Any:
        """
        Returns the node object with the given name.
        """
        for node_obj in self.graph:
            if node_obj.name == node_name:
                return node_obj
        return None

    # -------- String-based Rule Helpers as Static Methods --------
    @staticmethod
    def apply_R0(curr_step, prev_step):
        """
        Rule 0: Extension.
        If the current structure is exactly the previous structure with an extra dot appended,
        then R0 applies.
        """
        return prev_step + "."
    
    @staticmethod
    def is_valid_dot_bracket(seq):
        """Checks if a dot‐bracket string is valid (balanced and properly nested)."""
        stack = []
        for char in seq:
            if char == '(':
                stack.append(char)
            elif char == ')':
                if not stack:
                    return False
                stack.pop()
        return len(stack) == 0

    @staticmethod
    def get_matching_index(structure, index):
        """Given a dot‐bracket structure, return the matching closing parenthesis index for the opening at index."""
        if structure[index] != '(':
            return None
        stack = []
        for i in range(index, len(structure)):
            if structure[i] == '(':
                stack.append(i)
            elif structure[i] == ')':
                if not stack:
                    return None
                start = stack.pop()
                if start == index:
                    return i
        return None

    def check_complementarity(structure):
        """
        Checks that every paired position in the structure corresponds to nodes with complementary values.
        Uses the pairtable obtained from the structure.
        Assumes that graph.node_list has been populated so that:
        node_list[i] corresponds to the node for structure position i (with node name = i+1).
        A valid pair has one node with complement value 0 and the other with value 1,
        and both nodes must be in the same connected module.
        """
        try:
            node_list = graph.node_list
        except AttributeError:
            raise RuntimeError("graph.node_list is not set. It must be assigned after node creation.")
        
        # Obtain the pairtable for the structure.
        # Here we assume that path_to_pairtablepath returns a list, and we take its first element.
        pairtable = structure_to_pairtable(structure)
        # For example, if structure = "(()).", then pairtable might look like:
        # [ [length], p1, p2, ..., pn ]
        #
        # We assume that pairtable[1]..pairtable[n] contain the pairing partner (or 0 if unpaired),
        # and that these indices are 1-indexed.
        
        for i in range(1, len(pairtable)):
            j = pairtable[i]
            # Only check each pair once (i.e. when i < j).
            if j != 0 and i < j:
                # Convert 1-indexed positions to 0-indexed list indices.
                node1 = node_list[i - 1]
                node2 = node_list[j - 1]
                logger.debug(f"for {structure = } Checking pair: position {i} with position {j}: "
                    f"node1.complement = {node1.complement}, node1.connected = {node1.connected}; "
                    f"node2.complement = {node2.complement}, node2.connected = {node2.connected}")
                # They are complementary if one has complement 0 and the other 1,
                # and they must be in the same connected module.
                if node1.complement == node2.complement or node1.connected != node2.connected:
                    return False
        return True


    @staticmethod
    def apply_R2_1(curr_step):
        """
        Try all possible ways to partition curr_step into:
        prefix + "(" + A + ")" + B + "." + suffix
        and produce new structures according to Rule R2.1.
        """
        outcomes = []
        n = len(curr_step)
        # We need at least 4 characters to match the pattern
        for i in range(n):
            if curr_step[i] != '(':
                continue
            for j in range(i+1, n):
                if curr_step[j] != ')':
                    continue
                # Now, consider all possible positions for the literal dot after the pair
                for k in range(j+1, n):
                    if curr_step[k] != '.':
                        continue
                    # Define the four segments:
                    prefix = curr_step[:i]
                    A = curr_step[i+1:j]
                    B = curr_step[j+1:k]
                    suffix = curr_step[k+1:]
                    # Check if A, B, prefix, suffix are valid dot-bracket structures
                    if not (graph.is_valid_dot_bracket(prefix) and
                            graph.is_valid_dot_bracket(A) and 
                            graph.is_valid_dot_bracket(B) and
                            graph.is_valid_dot_bracket(suffix)):
                        continue
                    new_structure = prefix + "." + A + "(" + B + ")" + suffix
                    if (graph.is_valid_dot_bracket(new_structure) and 
                        graph.check_complementarity(new_structure)):
                        outcomes.append(new_structure)
        return outcomes if outcomes else None

    @staticmethod
    def apply_R2_2(curr_step):
        """
        Rule 2.2:
        Transform a structure of the form:
            prefix + "." + A + "(" + B + ")" + suffix
        into:
            prefix + "(" + A + ")" + B + "." + suffix

        For example:
            curr_step = "..()()"
        can be interpreted as:
            prefix = "."
            "."      (the dot to be transformed, at index 1)
            A = ""  
            "(" + B + ")" = "()"  
            suffix = "()"
        yielding:
            new_structure = ".(" + "" + ")" + "" + "." + "()"  =>  ".().()"

        Returns a list of outcomes that also obey complementarity.
        """
        outcomes = []
        # Loop over every index in the string.
        for i in range(len(curr_step)):
            if curr_step[i] != '.':
                continue
            # Consider the substring starting at position i.
            substring = curr_step[i:]
            m = re.match(r'\.(?P<A>[^()]*)\((?P<B>[^()]*)\)(?P<suffix>.*)$', substring)
            if m:
                pre = curr_step[:i]
                A = m.group("A")
                B = m.group("B")
                suffix = m.group("suffix")
                new_structure = pre + "(" + A + ")" + B + "." + suffix
                # Optionally, check for validity and complementarity.
                if (graph.is_valid_dot_bracket(new_structure) and 
                    graph.check_complementarity(new_structure)):
                    outcomes.append(new_structure)
        return outcomes if outcomes else None

    @staticmethod
    def apply_R2_3(curr_step):
        """
        Rule 2.3:
        Transform a structure of the form:
            prefix + "(" + A + ")" + B + "." + suffix
        into:
            prefix + "(" + A + "." + B + ")" + suffix

        This implementation explores every possible partition of the content
        originally split as A (inside the parentheses) and B (between the closing
        parenthesis and the dot).

        For example, if curr_step = "()()." (with indices:
            0: '('
            1: ')'
            2: '('
            3: ')'
            4: '.'
        ),
        then one valid decomposition is:
            prefix = ""
            open_index = 0  (character '(' at pos0)
            close_index = 1 (character ')' at pos1)
            A = curr_step[1:1] = ""      (nothing between positions 0 and 1)
            B = curr_step[2:4] = "()"     (from pos2 to pos3)
            suffix = curr_step[5:] = ""    (after the dot)
        We then set T = A + B = "()", and consider every split:
            k = 0: A′ = "" and B′ = "()"  => new_structure = "(" + "" + "." + "()" + ")" = "(.())"
            k = 1: A′ = "(" and B′ = ")"  => new_structure = "(" + "(" + "." + ")" + ")" = "((.))"
            k = 2: A′ = "()" and B′ = ""  => new_structure = "(" + "()" + "." + "" + ")" = "(())."
        
        We also try other choices for the outer pair. For example, choosing
        open_index = 2 and close_index = 3 yields:
            prefix = "()"
            A = curr_step[3:3] = "" and B = ""  => new_structure = "()(" + "" + "." + "" + ")" = "()(. )",
        i.e. "()(.)"
        
        Only candidates that are valid dot‐bracket strings and obey complementarity
        (per graph.is_valid_dot_bracket and graph.check_complementarity) are returned.
        """
        outcomes = set()
        n = len(curr_step)
        # Find all positions of the literal dot that we want to shift inside a pair.
        dot_positions = [i for i, ch in enumerate(curr_step) if ch == '.']
        if not dot_positions:
            return None

        # Try each dot position (usually one, but we allow multiple)
        for dot_index in dot_positions:
            # Now, search for an outer pair (an opening "(" before the dot and its matching ")"
            # before the dot). We'll try every candidate.
            for open_index in range(dot_index):
                if curr_step[open_index] != '(':
                    continue
                for close_index in range(open_index + 1, dot_index):
                    if curr_step[close_index] != ')':
                        continue
                    # Decompose curr_step as:
                    #   prefix + "(" + A + ")" + B + "." + suffix
                    prefix = curr_step[:open_index]
                    A = curr_step[open_index + 1:close_index]
                    B = curr_step[close_index + 1:dot_index]
                    suffix = curr_step[dot_index + 1:]
                    # Merge A and B into one string T, then try every partition of T.
                    T = A + B
                    for k in range(len(T) + 1):
                        A_prime = T[:k]
                        B_prime = T[k:]
                        new_structure = prefix + "(" + A_prime + "." + B_prime + ")" + suffix
                        # Debug print:
                        logger.debug(f"Trying new_structure = {new_structure} "
                            f"(prefix='{prefix}', A='{A}', B='{B}', suffix='{suffix}', "
                            f"split at {k} yields A'='{A_prime}', B'='{B_prime}')")
                        if (graph.is_valid_dot_bracket(new_structure) and 
                            graph.check_complementarity(new_structure)):
                            outcomes.add(new_structure)
        outcomes = list(outcomes)
        logger.debug(f"Outcomes = {outcomes}")
        return outcomes if outcomes else None



    @staticmethod
    def apply_R2_4(curr_step):
        """
        Rule 2.4:
          Transform a structure of the form:
            prefix + "(" + A + "." + B + ")" + suffix
          into:
            prefix + "(" + A + ")" + B + "." + suffix
        Returns a list of outcomes (if any) that obey complementarity.
        """
        pattern = r'^(.*?)\((.*?)\.(.*?)\)(.*?)$'
        outcomes = []
        for m in re.finditer(pattern, curr_step):
            prefix = m.group(1)
            A      = m.group(2)
            B      = m.group(3)
            suffix = m.group(4)
            # (Optionally check validity of parts.)
            new_structure = prefix + "(" + A + ")" + B + "." + suffix
            if graph.is_valid_dot_bracket(new_structure):
                if graph.check_complementarity(new_structure): 
                    outcomes.append(new_structure)
        return outcomes if outcomes else None

    @staticmethod
    def apply_R2_5(curr_step):
        """
        Rule 2.5:
        Transform a structure of the form:
            prefix + "." + A + "(" + B + ")" + suffix
        into:
            prefix + "(" + "." + A + B + ")" + suffix

        For example, with curr_step = "..()()":
            If we choose the dot at index 1:
            - prefix = curr_step[:1] = "."
            - The substring from index 1 is ".()()"
                Applying the pattern:
                ^\.(?P<A>.*?)\((?P<B>.*?)\)(?P<suffix>.*)$
                yields:
                A = ""    (since immediately after the dot we see "(")
                B = ""    (from the inner "()")
                suffix = "()"  (the remainder)
            - The new structure is:
                prefix + "(" + "." + A + B + ")" + suffix
                which becomes:
                "." + "(" + "." + "" + "" + ")" + "()"  →  ".(.)()"
                
        Returns a list of outcomes that obey complementarity.
        """
        outcomes = []
        # Try each dot in the string as the candidate for transformation.
        for i in range(len(curr_step)):
            if curr_step[i] != '.':
                continue
            prefix = curr_step[:i]
            substring = curr_step[i:]
            m = re.match(r'^\.(?P<A>.*?)\((?P<B>.*?)\)(?P<suffix>.*)$', substring)
            if m:
                A = m.group("A")
                B = m.group("B")
                suffix = m.group("suffix")
                new_structure = prefix + "(" + "." + A + B + ")" + suffix
                # Check that the new structure is valid and complementary.
                if (graph.is_valid_dot_bracket(new_structure) and 
                    graph.check_complementarity(new_structure)):
                    outcomes.append(new_structure)
        return outcomes if outcomes else None


    @staticmethod
    def apply_R2_6(curr_step):
        """
        Rule 2.6:
          Transform a structure of the form:
            prefix + "(" + A + "." + B + ")" + suffix
          into:
            prefix + "." + A + "(" + B + ")" + suffix
        Returns a list of outcomes that obey complementarity.
        """
        pattern = r'^(.*?)\((.*?)\.(.*?)\)(.*?)$'
        outcomes = []
        for m in re.finditer(pattern, curr_step):
            prefix = m.group(1)
            A      = m.group(2)
            B      = m.group(3)
            suffix = m.group(4)
            new_structure = prefix + "." + A + "(" + B + ")" + suffix
            if (graph.is_valid_dot_bracket(new_structure) and 
                graph.check_complementarity(new_structure)):
                outcomes.append(new_structure)
        return outcomes if outcomes else None

    @staticmethod
    def apply_R2_7(curr_step):
        """
        Rule 2.7:
          Transform a structure of the form:
            prefix + "(" + A + ")" + B + "(" + C + ")" + suffix
          into:
            prefix + "(" + A + "(" + C + ")" + B + ")" + suffix
        Returns a list of outcomes that obey complementarity.
        """
        pattern = r'^(.*?)\((.*?)\)(.*?)\((.*?)\)(.*?)$'
        outcomes = []
        for m in re.finditer(pattern, curr_step):
            prefix = m.group(1)
            A      = m.group(2)
            B      = m.group(3)
            C      = m.group(4)
            suffix = m.group(5)
            new_structure = prefix + "(" + A + "(" + C + ")" + B + ")" + suffix
            if (graph.is_valid_dot_bracket(new_structure) and 
                graph.check_complementarity(new_structure)):
                outcomes.append(new_structure)
        return outcomes if outcomes else None

    @staticmethod
    def apply_R2_8(curr_step):
        """
        Rule 2.8 Transformation:
        
        Transform any occurrence of:
        
            pre  +  "("  +  A  +  "("  +  B  +  ")"  +  C  +  ")"  + suf
            
        into:
        
            pre  +  "("  +  A  +  ")"  +  B  +  "("  +  C  +  ")"  + suf
            
        The function scans the string 'curr_step' for every possible way
        that a substring of the form "( A ( B ) C )" appears and produces the
        corresponding transformed structure.
        
        For example:
            Input:  "()(())"  
            may be parsed as pre="()", transformation target="(())" where
            A = "", B = "", C = "".
            The new structure becomes: "()" + "(" + "" + ")" + "" + "(" + "" + ")"  → "()()()"
            
        Returns a list of outcomes (one for each occurrence) that result from applying
        the transformation.
        """
        outcomes = []
        n = len(curr_step)

        def find_matching_paren(s, i):
            """Given that s[i]=='(', return index of its matching ')'."""
            count = 0
            for j in range(i, len(s)):
                if s[j] == '(':
                    count += 1
                elif s[j] == ')':
                    count -= 1
                    if count == 0:
                        return j
            return -1  # no matching parenthesis found

        # Iterate over all possible starting positions for an outer '('.
        for i in range(n):
            if curr_step[i] != '(':
                continue
            outer_close = find_matching_paren(curr_step, i)
            if outer_close == -1:
                continue  # skip unmatched parentheses

            # Now the outer structure is curr_step[i : outer_close+1]
            # We want to see if it contains a substructure of the form: A ( B ) C
            # That is, inside the outer pair, we look for an inner '('.
            for k in range(i + 1, outer_close):
                if curr_step[k] != '(':
                    continue
                inner_close = find_matching_paren(curr_step, k)
                if inner_close == -1 or inner_close >= outer_close:
                    continue  # not a valid inner pair inside the outer pair

                # Extract A, B, and C.
                A = curr_step[i + 1 : k]
                B = curr_step[k + 1 : inner_close]
                C = curr_step[inner_close + 1 : outer_close]
                pre = curr_step[:i]
                suf = curr_step[outer_close + 1:]

                # Build the new structure according to the rule:
                # pre + "(" + A + ")" + B + "(" + C + ")" + suf
                new_structure = pre + "(" + A + ")" + B + "(" + C + ")" + suf

                # (Optional) Validate new_structure if needed:
                if graph.is_valid_dot_bracket(new_structure):
                    if  graph.check_complementarity(new_structure):
                        outcomes.append(new_structure)

        return outcomes
    
    @staticmethod
    def apply_R1(self, curr_step, prev_step):
        """
        Rule 1: Binding.
        Looks for the first candidate position (except the last) where the previous structure
        has a dot and the current structure has an opening parenthesis, and where the last
        character in prev_step is '.' and in curr_step is ')'.
        
        The new structure is obtained by taking prev_step and replacing the candidate
        position with '(' and the last character with ')'.
        
        Before finalizing, it checks:
        - That the candidate’s opening has a matching closing parenthesis.
        - That the resulting structure is a valid dot‐bracket string.
        - That the structure obeys complementarity.
        - And that the new binding does not “interrupt” any potential legal pairing of two
            unpaired, complementary domains (using check_unpaired_complements).
        
        Returns:
            The new structure (str) after binding, or None if no valid binding is detected.
        """
        # Ensure the last-character condition.
        if prev_step[-1] != '.' or curr_step[-1] != ')':
            return None

        candidate_index = None
        # Look for the first candidate position (except the last) where prev_step is '.' and curr_step is '('.
        for i in range(len(curr_step) - 1):
            if prev_step[i] == '.' and curr_step[i] == '(':
                candidate_index = i
                break

        if candidate_index is None:
            return None

        # Build the new structure by copying prev_step and "binding" candidate_index and the last position.
        new_structure_list = list(prev_step)
        new_structure_list[candidate_index] = '('  # Bind the candidate.
        new_structure_list[-1] = ')'               # Bind the last character.
        new_structure = ''.join(new_structure_list)

        # Ensure that the newly added opening at candidate_index has a matching closing parenthesis.
        matching_index = self.get_matching_index(new_structure, candidate_index)
        if matching_index is None:
            return None

        # Check that the new structure is a valid dot-bracket string and obeys complementarity.
        if not (graph.is_valid_dot_bracket(new_structure) and graph.check_complementarity(new_structure)):
            return None

        # Convert the new structure into a pairtable.
        new_pt = structure_to_pairtable(new_structure)
        # Use check_unpaired_complements to verify that the new binding does not interrupt
        # any potential legal pairing between unpaired, complementary domains.
        if not graph.check_unpaired_complements(self,[new_pt]):
            return None

        return new_structure

    @staticmethod
    def apply_R1_dfs(current):
        """
        DFS variant of Rule 1: Binding.
        For every pair of candidate positions (i, j) (with i < j) that are unbound ('.'),
        bind those positions by changing position i to '(' and position j to ')'.
        
        Only outcomes that yield a valid dot‐bracket structure and obey complementarity are returned.
        
        (Note: It is assumed that graph.check_complementarity(new_structure) verifies that the new
        pair is formed between complementary (and connected) nucleotides. If further checks are needed,
        for example using a function like graph.are_complementary(i, j), they can be inserted.)
        """
        if len(current) < 2:
            return None

        outcomes = []
        n = len(current)
        for i in range(n - 1):
            if current[i] != '.':
                continue
            for j in range(i + 1, n):
                if current[j] != '.':
                    continue
                new_structure_list = list(current)
                new_structure_list[i] = '('
                new_structure_list[j] = ')'
                new_structure = ''.join(new_structure_list)
                logger.debug(f"Trying pairing: i={i}, j={j} -> {new_structure}")
                # Validate the new structure:
                if (graph.is_valid_dot_bracket(new_structure) and 
                    graph.check_complementarity(new_structure)):
                    outcomes.append(new_structure)
        return outcomes if outcomes else None



    # -------- Main Rule Sequencing Method --------
    def dfs_R2(self, current, target, sequence, visited, memo, depth_limit):
        # Check memoized results.
        if current in memo:
            return memo[current]


        p_acfp = path_to_pairtablepath([current,target])

        if current == target:
            if self.check_unpaired_complements(p_acfp):
                memo[current] = sequence
                return sequence
            else:
                return None

        if depth_limit <= 0:
            memo[current] = None
            return None

        rules = [
            ("R1", self.apply_R1_dfs),
            ("R2.1", self.apply_R2_1),
            ("R2.2", self.apply_R2_2),
            ("R2.3", self.apply_R2_3),
            ("R2.4", self.apply_R2_4),
            ("R2.5", self.apply_R2_5),
            ("R2.6", self.apply_R2_6),
            ("R2.7", self.apply_R2_7),
            ("R2.8", self.apply_R2_8)
        ]
        logger.debug(f"\nCurrent structure: {current}")
        for rule_name, rule_func in rules:
            logger.debug("\nTrying rule:", rule_name)
            outcomes = rule_func(current)
            if outcomes is None:
                outcomes = []
            elif not isinstance(outcomes, list):
                outcomes = [outcomes]
            logger.debug(f"Outcomes = {outcomes}")

            for new_structure in outcomes:
                if not is_balanced_structure(new_structure):
                    continue

                    
                if new_structure is None or new_structure == current:
                    continue

                # For the first move, ensure the binding partner of the last character changes.
                new_pairtable = structure_to_pairtable(new_structure)
                current_pairtable = structure_to_pairtable(current)

                #if not sequence and new_pairtable[-1] == current_pairtable[-1]:
                #    continue

                if new_structure in visited:
                    continue

                visited.add(new_structure)
                new_sequence = sequence + [(rule_name, new_structure)]
                result = self.dfs_R2(new_structure, target, new_sequence, visited, memo, depth_limit - 1)
                if result is not None:
                    memo[current] = result
                    return result
                # Note: do NOT remove new_structure from visited here.
                # It has been fully explored and does not lead to a solution.
            # End for outcomes
        # End for rules
        memo[current] = None
        return None



    def find_rule_seq(self, curr_step, prev_step):
        # Apply R0 and R1 as before:
        modified_structure = prev_step
        rule_seq = []


        logger.debug(f"New step \n\n{curr_step = } {prev_step = }")
        r0 = self.apply_R0(curr_step, modified_structure)
        if r0 is not None:
            rule_seq.append(("R0",))
            if r0 != modified_structure:
                modified_structure = r0
            else:
                return rule_seq


        r1 = self.apply_R1(self,curr_step, modified_structure)
        if r1 is not None:
            rule_seq.append(("R1", r1))
            modified_structure = r1

      
         
        
                
        logger.debug(f"R1 post {modified_structure = }")

        # Now use DFS for R2.X rules:
        visited = set([modified_structure])
        memo = {}

        #First Rule must affect 
        if curr_step[-1] == "." and prev_step != curr_step[0:len(curr_step) -1]:
            return False
        dfs_result = self.dfs_R2(modified_structure, curr_step, [], visited, memo , depth_limit=20)
        logger.debug(f"{dfs_result = }")

        p_acfp = path_to_pairtablepath([prev_step,curr_step])
        if not self.check_unpaired_complements(p_acfp):
            return False
        
        if dfs_result is not None:
            rule_seq.extend(dfs_result)
            return rule_seq
        else:
            return False

    def is_favorable(self, acfp):
        """
        Checks whether an acfp is favorable by verifying that a valid sequence of rules exists
        to move from one folding step to the next.
        """

        if self.pseudoknot_attack_detection(path_to_pairtablepath(acfp)):
            return False
        
        rule_seq_overall = []
        for i in range(1, len(acfp)):
            curr_rule_seq = self.find_rule_seq(acfp[i], acfp[i-1])
            if curr_rule_seq:
                rule_seq_overall.append(curr_rule_seq)
            else:
                logger.debug(f"Path at step:{i} going from {acfp[i-1]} {acfp[i]} is not favorable\n{rule_seq_overall = }")
                return False
            if len(curr_rule_seq) == 0:
                return False
        logger.debug(f"{rule_seq_overall = }")
        return rule_seq_overall
    
    def get_current_edges(self,step_number):
        current_edges = []
        for edge in self.edges:
            if edge[1] <= step_number + 1:
                current_edges.append(edge)
        return current_edges
    """
    def get_weights(self, acfp):
        assigned_edges = {}
        for x, step in enumerate(acfp[1:], start=1):
            collected_edges = self.get_current_edges(x)
            logger.debug("--" * 50)
            logger.debug(f"\nAssigned Edges: {assigned_edges}\nCollected Edges: {collected_edges}\nCurrent Step: {step}")
            active_edges = []
            inactive_edges = []
            
            # Separate collected edges into active and inactive.
            for edge in collected_edges:
                if step[edge[0]] == edge[1] and step[int(edge[1])] == int(edge[0]):
                    active_edges.append(edge)
                else:
                    inactive_edges.append(edge)
            
            logger.debug(f"Inactive Edges: {inactive_edges}")
            inactive_edge_to_remove = []
            for inactive_edge in inactive_edges:
                logger.debug(f"Inactive Edge: {inactive_edge}")
                for neighbor in self.edge_neighbors[inactive_edge]:
                    logger.debug(f"Neighbor: {neighbor}")
                    if neighbor in assigned_edges and neighbor in active_edges and self.edges[inactive_edge] < self.edges[neighbor]:
                        logger.debug(f"Removed: {inactive_edge}")
                        inactive_edge_to_remove.append(inactive_edge)
            
            for edge in set(inactive_edge_to_remove):
                if edge in inactive_edges:
                    inactive_edges.remove(edge)
                #if edge not in assigned_edges:
                #    assigned_edges[edge] = 1
            
            # Sort active edges (e.g., by the second element)
            active_edges.sort(key=lambda x: x[1])
            logger.debug(f"Active Edges: {active_edges}")
            logger.debug(f"Inactive Edges: {inactive_edges}")
            
            # Process inactive edges using active edges.
            while len(inactive_edges) != 0:
                logger.debug(f"\nBegin weighting edges: {inactive_edges}")
                # If no active edges remain, assign weight 1 to all remaining inactive edges.
                if not active_edges:
                    for edge in inactive_edges:
                        if edge not in assigned_edges:
                            assigned_edges[edge] = 1
                            self.edges[edge] = 1
                    break
                
                current_edge = active_edges.pop()
                logger.debug(f"Current Edge: {current_edge}")
                
                inactive_flag = False
                for edge in inactive_edges:
                    # Check if the current active edge influences this inactive edge.
                    if edge[0] in current_edge or edge[1] in current_edge:
                        inactive_flag = True
                        logger.debug(f"Inactive edge influences edge {edge}")
                
                if current_edge not in assigned_edges and inactive_flag:
                    l_node = self.get_node_by_name(current_edge[0])
                    r_node = self.get_node_by_name(current_edge[1])
                    logger.debug(f"Left weight {l_node.max_weight, l_node}")
                    logger.debug(f"Right weight {r_node.max_weight, r_node}")
                    edge_weight = max(l_node.max_weight, r_node.max_weight) + 1
                    self.edges[current_edge] = edge_weight
                    logger.debug(f"Edge Weight: {edge_weight} for edge: {current_edge}")
                    
                    l_node.max_weight = edge_weight
                    r_node.max_weight = edge_weight
                    logger.debug(f"Updated Left weight {l_node.max_weight, l_node}")
                    logger.debug(f"Updated Right weight {r_node.max_weight, r_node}")
                    assigned_edges[current_edge] = self.edges[current_edge]
                
                edges_to_remove = []
                for edge in inactive_edges:
                    logger.debug(f"Current inactive edge {edge}")
                    if edge[0] in current_edge or edge[1] in current_edge:
                        edges_to_remove.append(edge)
                
                for edge in edges_to_remove:
                    inactive_edges.remove(edge)
                    logger.debug(f"Removed edge: {edge}")
                    if edge not in assigned_edges:
                        assigned_edges[current_edge] = self.edges[current_edge]
                        assigned_edges[edge] = self.edges[edge]
                        logger.debug(f"New edge assignment: {assigned_edges}")
                        
            logger.debug("--" * 50)
        logger.debug(self.edges)"""

    def get_weights(self, acfp):
        assigned_edges = {}
        for x, step in enumerate(acfp[1:], start=1):
            # Assume step is an ordered list of nodes;
            # the last element is the newly appended node.
            if len(step) < 2:
                # Not enough nodes to form an edge
                continue
            
            # Identify the edge that connects the new node with its predecessor.
            predecessor = step[-2]
            new_node = step[-1]
            edge = (predecessor, new_node)
            
            # Retrieve node objects.
            l_node = self.get_node_by_name(predecessor)
            r_node = self.get_node_by_name(new_node)
            
            # Determine the new edge weight.
            # If the edge already exists (i.e. being reapplied/replacing), increment its weight by 1.
            # Otherwise, set its weight based on the max_weight of the two nodes.
            if edge in self.edges:
                new_weight = self.edges[edge] + 1
            else:
                new_weight = max(l_node.max_weight, r_node.max_weight) + 1
            
            # Update the edge weight and update the nodes' max_weight.
            self.edges[edge] = new_weight
            l_node.max_weight = new_weight
            r_node.max_weight = new_weight
            
            assigned_edges[edge] = new_weight
            logger.debug(f"Step {x}: Edge {edge} assigned weight {new_weight}")
            
        logger.debug("Final edge weights: " + str(self.edges))
        
    def get_node_by_name(self, node_name):
        for node in self.graph:
            if node.name == node_name:
                return node
        return None

    def verify_weights(self,acfp):
        """
        Function to verify if the acfp can be formed based on the edge weights alone
        """
        logger.info("\n\nVerifying weights\n\n")
        assigned_edges = {}
        for x,step in enumerate(acfp[1:],start=1): 
            collected_edges = self.get_current_edges(x)
            logger.debug("--"*50)
            logger.debug(f"\nAssigned Edges: {assigned_edges}\n")
            logger.debug(f"\nCollected Edges: {collected_edges}")
            logger.debug(f"Current Step: {step}")
            active_edges = []
            inactive_edges = []
            
            # Get active and inactive edges
            for edge in collected_edges:
                if step[edge[0]] == edge[1] and step[int(edge[1])] == int(edge[0]):
            
                    active_edges.append(edge)
                else:
                    inactive_edges.append(edge)
            
            logger.debug(f"\n{active_edges = }")
            logger.debug(f"\n{inactive_edges = }")
            while len(inactive_edges) != 0:
                curr_edge = inactive_edges.pop()
                current_neighbors = [neighbor for neighbor in self.edge_neighbors[curr_edge] if neighbor in active_edges]
                max_weight = float("-inf")
                max_weight = max([self.edges[neighbor] for neighbor in current_neighbors], default=float("-inf"))

                if self.edges[curr_edge] > max_weight and self.edges[curr_edge] != 1:
                    logger.debug(f"{self.edges[curr_edge] = } , {max_weight = }")
                    logger.debug(f"\n{active_edges = }")
                    logger.debug(f"\n{inactive_edges = }")
                    logger.debug(f"Verify Weights: Folding path not possible in this step: {step}\n\n Goodbye")
                    return False 
        return True
          
    #_________Weights______#
    def get_domain_seq(self):
        """Returns the domain level sequence of the abstract folding path in form of a List.

        Returns:
            domain_seq(list): domain level sequence
        """
        self.create_domain_seq()
        domain_seq = []

        for x,node in enumerate(self.graph):
            domain_seq.append(str(" ".join(node.prefix) + " " + node.middle + " " + " ".join(node.suffix) + " " + "S" + str(x)))
        
        logger.debug(f"Resulting domain seq: {domain_seq}")
        logger.debug(f"{self.edges}")

        return domain_seq
    
    def create_domain_seq(self):
        logger.debug(f"************************Begin create_domain_seq****************************")
        domains = list(string.ascii_lowercase)
        domains += [x + z for x in string.ascii_lowercase for z in string.ascii_lowercase if 'l' not in (x, z) and 'm' not in (x, z)]
        
        visited_nodes = []
        for node in self.graph:
            current_node = node
            visited_nodes.append(current_node)
            logger.debug("____________________"*3)
            logger.debug(f"\n Current Node:  {current_node}\n")
            current_node.middle = "L" + str(node.connected)

            if current_node.max_weight > 1 and len(current_node.prefix) < current_node.max_weight - 1:
                prefix = domains[:current_node.max_weight - 1 - len(current_node.prefix)] + list(current_node.prefix)
                domains = domains[current_node.max_weight - 1 - len(current_node.prefix):]
                current_node.prefix = prefix

                suffix = list(current_node.suffix) + domains[:current_node.max_weight - 1 - len(current_node.suffix)] 
                domains = domains[current_node.max_weight - 1 - len(current_node.suffix):]
                current_node.suffix = suffix
            logger.debug(f"Current Node: {current_node}")
            
            for neighbor in current_node.neighbors:
                if neighbor not in visited_nodes and neighbor.connected == current_node.connected:
                    
                    logger.debug(f"graph edges :{self.edges}\n")
                    shared_weight = self.edges[(min(current_node.name, neighbor.name), max(current_node.name, neighbor.name))]
                    logger.debug(f"Shared Weight: {shared_weight}, Edges:{(current_node.name,neighbor.name)}")
                    

                    if shared_weight > 1:    
                        logger.debug(f"\n\nPrefix: {current_node.prefix}\n\n")
                        logger.debug(f"\n\nSuffix: {current_node.suffix}\n\n")
                        neighbor.prefix = current_node.suffix[:shared_weight - 1][::-1]
                        neighbor.suffix = current_node.prefix[- (shared_weight - 1):][::-1]
                    logger.debug(neighbor)

        for node in self.graph:
            logger.debug(f"Current Node: {node}")
            if node.complement:
                node.prefix = [str(domain).rstrip("*") + "*" for domain in node.prefix]
                node.middle += "*"
                node.suffix = [str(domain).rstrip("*") + "*" for domain in node.suffix]

        logger.debug("\nCreate Domain Seq Finished\n\n\n")
        logger.debug(self.print_nodes())

#_________Legal_folding_path_checks_______#
    def cycle_detection(self):
        visited = set()
        def dfs(node, parent):
            visited.add(node)
            for neighbor in node.neighbors:
                if neighbor != parent:
                    if neighbor in visited or dfs(neighbor, node):
                        return True
            return False

        for node in self.graph:
            logger.debug(node)
            if node not in visited:
                if dfs(node, None):
                    return True
        return False
    

    def pseudoknot_attack_detection(self, pairtable):
        for x in range(1, len(pairtable) - 1):
            secured_domains = []
            for i, j in enumerate(pairtable[x][1:], start=1):
                # Only consider if there's at least one index between i and j.
                if j > i + 1:
                    # Save the indices strictly between i and j.
                    domain_indices = list(range(i + 1, j))
                    secured_domains.append(domain_indices)
            
            logger.debug(f"secured domains {secured_domains}")
            

            for domain in secured_domains:
                if pairtable[x+1][-1] in domain:
                    return True
            
            """for indices in secured_domains:
                # Compare the elements at these indices in row x and row x+1.
                # Use a list comprehension to extract the elements.
                row_x_slice = [pairtable[x][k] for k in indices]
                row_x1_slice = [pairtable[x+1][k] for k in indices]
                if row_x_slice != row_x1_slice:
                    # Check that the boundaries in row x+1 (first and last indices of the domain)
                    # are nonzero.
                    #if pairtable[x+1][indices[0]] != 0 and pairtable[x+1][indices[-1]] != 0:
                        return True"""
        return False


    def check_unpaired_complements(self, acfp):
        """
        For each step in acfp (a pairtable where step[0] is the number of positions, positions are 1-indexed):
        - For every pair of unpaired positions x and y (i.e. step[x] == 0 and step[y] == 0, with x < y)
        - If the nodes at x and y are complementary and connected (i.e. they share the same 'connected'
            value but have opposite 'complement'),
        - Then check every index j between x and y.
            If no position j exists such that j is paired (step[j] != 0) and its partner lies outside the [x, y] range,
            then return False.
        If all complementary unpaired pairs are "interfered" with by an external pairing (or no such pair exists),
        return True.
        """
        for step in acfp:
            n = step[0]
            for x in range(1, n + 1):
                # Skip if x is paired.
                if step[x] != 0:
                    continue
                node_x = self.get_node_by_name(x)
                for y in range(x + 1, n + 1):
                    # Skip if y is paired.
                    if step[y] != 0:
                        continue
                    node_y = self.get_node_by_name(y)
                    # Check if x and y are complementary and in the same connected component.
                    if node_x.connected == node_y.connected and node_x.complement != node_y.complement:
                        interfering_found = False
                        # Check all indices between x and y.
                        for j in range(x + 1, y):
                            if step[j] != 0:
                                partner = step[j]
                                # If the partner is outside the [x, y] range, we mark that interference.
                                if partner < x or partner > y:
                                    interfering_found = True
                                    break
                        # If no interfering pairing is found, then this pair of unpaired, complementary domains is "isolated."
                        if not interfering_found:
                            return False
        return True

    def consistency_check(self,acfp):
        """Checks wether or not an aCFP is consistent. 

            1. Check wether or not a cycle is existant in the graph. 
            2. Check if the graph translation of the graph results in the targeted aCFP. 

        """

        logger.debug('\n Begin Consistency Check')
        #Cycle Check 
        if self.cycle_detection():
            logger.debug(f'{acfp = } is not consistent due to the occurence of a cycle in the path.')
            return False 

        pairtable_acfp = path_to_pairtablepath(acfp)
        self.get_weights(acfp=pairtable_acfp)

        if not self.verify_weights(acfp=pairtable_acfp):
            return False


        return True
    '''

    def consistency_check(self,acfp):
        """
        Checks whether an abstract cotranscriptional folding path (aCFP) is consistent.
        
        
        Args:
            acfp (list of str): The abstract cotranscriptional folding path, where each C_i is a string of length i.
            
        Returns:
            bool: True if the aCFP is consistent; False otherwise.
        """
        n = len(acfp)
        # Iterate over all possible ranges [i, j] (1-indexed).
        for i in range(1, n+1):
            for j in range(i+1, n+1):
                valid_subs = set()
                # Only constraints with length >= j are considered (they correspond to steps k from j to n).
                for k in range(j, n+1):
                    # Get the substring from C_k corresponding to indices [i, j].
                    # Python strings are 0-indexed so we take acfp[k-1][i-1:j]
                    sub = acfp[k-1][i-1:j]
                    if self.is_valid_dot_bracket(sub):
                        valid_subs.add(sub)
                if len(valid_subs) > 1:
                    print(f"Inconsistency found in domain [{i},{j}]: {valid_subs}")
                    return False
        return True'''
    
    def is_bipartite_dfs(self, node, visited_nodes, complement, connected, x):
        """
        Performs a depth-first search to check bipartiteness of a connected component.
        """
        visited_nodes[node] = True
        node.complement = complement

        for neighbor in node.neighbors:
            if neighbor.name not in connected:
                continue  # Only consider neighbors in the current connected component

            if neighbor not in visited_nodes:
                if not self.is_bipartite_dfs(neighbor, visited_nodes, not complement, connected, x):
                    return False
            else:
                if neighbor.complement == node.complement:
                    return False

        return True

    def bipartite_check(self, connected_components):
        """
        Checks whether the graph is bipartite given the connected components.
        Returns True if bipartite, False otherwise.
        """
        visited_nodes = {}
        logger.debug(f"Connected components: {connected_components}")

        for x, connected in enumerate(connected_components):
            for node in self.graph:
                if node.name in connected and node not in visited_nodes:
                    if not self.is_bipartite_dfs(node, visited_nodes, True, connected, x):
                        return False
        return True


def preprocessing(acfp):
    """Takes an acfp and creates a graph based on it. Using the graph a domain level sequence can be created. 

        Args:
            acfp(list): each entry corresponds to a structural constraint
        Returns:
            acfp_graph(obj.): graph object
    """
    if acfp[0][0] == ".":
        pairtable_acfp = path_to_pairtablepath(acfp)
    else:
        pairtable_acfp = acfp

    db_acfp = pairtablepath_to_dotbracket(pairtable_acfp)

    logger.debug(f"Pairtable Input: {pairtable_acfp}")
    acfp_graph = graph()
    acfp_graph.create_nodes_from_pairtable(pairtable_acfp)
    graph.node_list = sorted(list(acfp_graph.graph.keys()), key=lambda node: node.name)

    logger.info("Graph after adding all Nodes")
    logger.info(acfp_graph.print_nodes())

    nodes = list(acfp_graph.graph.keys())
    nodes.insert(0,0)

    acfp_graph.create_edges(pairtable_acfp,nodes)
    
    logger.info("Graph after adding edges:")
    logger.info(acfp_graph.print_nodes())
    acfp_graph.get_edges()

    acfp_graph.set_edge_neighbors()

    logger.info(f"Edges: {acfp_graph.edges}")

    connected_components = find_connected_modules(pairtable_acfp)
    logger.debug(f"Connected c: {connected_components}")
    logger.info("\nFollowing Nodes are connected:")
    for x,component in enumerate(connected_components):
        nodes_in_component = set()
        logger.debug(component)
        for node_name in component:
            acfp_graph.get_node_by_name(node_name).connected = x

        for node_index in component:
            nodes_in_component.add(nodes[node_index])

    #check if sls_translatable 

    checks = dict()

    if acfp_graph.bipartite_check(connected_components):
        logger.debug(f"{db_acfp} is assignable")
        checks["assignable"] = True
    else:
        checks["assignable"] = False
        #raise SystemExit(f"aCFP {acfp} is not assignable")

    if acfp_graph.is_favorable(db_acfp):
        logger.debug(f"{acfp} is favorable")
        checks["favorable"] = True
    else: 
        checks["favorable"] = False
        #raise SystemExit(f"Path {acfp} is not favorable")

    if acfp_graph.consistency_check(db_acfp):
        logger.debug(f"{acfp = } is consistent")
        checks["consistent"] = True
    else:
        checks["consistent"] = False
        #raise SystemExit("aCFP not consistent")

    # Now, check if any of the tests failed:
    failed_checks = [name for name, passed in checks.items() if not passed]

    if failed_checks:
        error_message = "The following checks failed: " + ", ".join(failed_checks)
        raise SystemExit(error_message)
    
    assignment = [[] for _ in acfp]
    i = 0
    for current_node in acfp_graph.graph:
        current_node.middle = "L" + str(current_node.connected)
        if current_node.complement:
            current_node.middle += "*"
        assignment[i] = str(current_node.middle)
        i += 1

    print(f"{assignment = }")
    return assignment


def translate_acfp(acfp):
    
    A = preprocessing(acfp)    
    pt_path = path_to_pairtablepath(acfp)
    S = [[] for _ in acfp]
    x = 0 

    for i in range(len(acfp)):
        S[i].append(A[i])
        pt = pt_path[i]
        if pt[-1] != -1 and pt[pt[-1]] != -1:
            S[i] = comp

def set_verbosity(console_handler,verbosity):
    if verbosity == 0:
        console_handler.setLevel(logging.CRITICAL)
    elif verbosity == 1:
        console_handler.setLevel(logging.ERROR)
    elif verbosity == 2:
        console_handler.setLevel(logging.INFO)
    elif verbosity >= 3:
        console_handler.setLevel(logging.DEBUG)

def input_parser(file_path):
    try:
        with open(file_path, 'r') as file:
            file_content = file.read()  # Add this line to print the file content
            lines = file_content.split('\n')  # Split lines based on newline characters
            lines = [line.strip() for line in lines]  # Remove leading and trailing whitespace
            # Filter lines based on the regex pattern
            lines = [line.strip() for line in lines if re.match(r'^[0-9,().\[\]]+$', line) and not line.startswith('#')]
            if not lines:
                raise SystemExit("File is empty. Here is your sequence: ")

            return lines

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        raise SystemExit("File not Found")
        return []
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return []

def main():
        #_________________Argparse_________________#
        parser = argparse.ArgumentParser(
            formatter_class = argparse.RawDescriptionHelpFormatter,
            description=f"""Copaths path_to_seq version {__version__}: generates a domain level sequence based on an abstract folding path 


The acfp must have following propperties to be translated to a domain level sequence:
            - Pseudoknot-Free
            - Structure must be fully saturated
            - If a structure was defined it cant be changed if no additional step influences the defined substructure"""
                        )

        parser.add_argument("-i","--input",help="Reads txt file as input if not specified user can input via console.")
        parser.add_argument("-v", "--verbose",action="count", default = 0,
            help = "Track process by writing verbose output to STDOUT during calculations.")
        
        args = parser.parse_args()

        set_verbosity(logger,args.verbose)
        

        #_________________Input_Parsing_________________#
        
        if args.input is None and not sys.stdin.isatty():
            # Handle input from stdin (piped data)
            acfp = sys.stdin.read().strip().split(',')
             
        elif args.input == None:

            acfp = acfp_terminal_input()
        
        else:
            acfp = input_parser(args.input)
            print(f"\n\nInput folding path:\n{acfp}\n\n")

        #_________________Translation_of_d_seq_________________#
        logger.debug(acfp)
        

        domain_seq = translate_acfp()



if __name__ == "__main__":

    main()
