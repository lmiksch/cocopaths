from typing import Any
import re,argparse, logging

try:
    from collections import deque
    from copaths.convert_functions import (path_to_pairtablepath, find_connected_modules,is_balanced_structure)
    import string
except:
    pass


#______define_logger______#

logger = logging.getLogger('copaths')
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

#_________Graph_creation________#

    def add_node(self,node):
         if node not in self.graph:
              self.graph[node] = []

    def add_edge(self,node1,node2):    
        if node2 not in node1.neighbors:
            node1.neighbors.add(node2)  
            node2.neighbors.add(node1)  

    def create_nodes_from_pairtable(self,pairtable):
        for x in range(1,pairtable[-1][0] + 1):
            new_node = node(x)
            self.add_node(new_node)

    def dfs(self, start_node):
        visited = set()

        def dfs_recursive(node):
            visited.add(node.name)
            logger.debug(node.name)

            for neighbor in node.neighbors:
                if neighbor.data not in visited:
                    dfs_recursive(neighbor)

        dfs_recursive(start_node)

    def print_nodes(self):
        return "\n".join([str(node) for node in self.graph])

    def create_edges(self,afp,nodes):
        for step in afp[::-1]:
            for x in range(1,step[0]):
                if step[x] != 0:
                    node1 = nodes[x]
                    node2 = nodes[step[x]]
                    
                    self.add_edge(node1, node2)
                else:
                    continue
    
    def get_edges(self):
        
        for node in self.graph:
            for neighbor in node.neighbors:
                self.edges[(node.name , neighbor.name)] = 1 #default weight of edge is 1 since they will always bind with middle domain


        #Remove duplicate edges
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

#_________Bipartite_check________#

    def is_bipartite_dfs(self, node, visited_nodes, complement, connected,x):
        visited_nodes[node] = True
        node.complement = complement
        
        
        for neighbor in node.neighbors:
            
            if neighbor.name not in connected:
                
                continue  # Skip neighbors not in the specified connected components

            if neighbor not in visited_nodes:
                node.connected = x
                
                if not self.is_bipartite_dfs(neighbor, visited_nodes, not complement, connected,x):
                    return False
            else:
                node.connected = x
                if neighbor.complement == node.complement:
                    return False

        return True

    def bipartite_check(self, connected_components):
        visited_nodes = {}
        for x, connected in enumerate(connected_components):
            for node in self.graph: 
                if node.name in connected and node not in visited_nodes:
                    node.connected = x
                    if not self.is_bipartite_dfs(node, visited_nodes, True, connected,x):
                        return False

        return True

    def get_current_edges(self,step_number):
        current_edges = []
        for edge in self.edges:
            if edge[1] <= step_number + 1:
                current_edges.append(edge)
        return current_edges
        
    def get_inequalities(self,afp):
        """Creates Inequalities based on the abstract folding path and the graph
        currently not used due to the get_weights function
        """
        
        inequalities = []

        for x,step in enumerate(afp): 
            collected_edges = self.get_current_edges(x)
            active_edges = []
            inactive_edges = []
            for edge in collected_edges:
                if step[edge[0]] == edge[1] and step[int(edge[1])] == int(edge[0]):
                    active_edges.append(edge)
                else:
                    inactive_edges.append(edge)

            active_edges.sort(key= lambda x: (x[1]))
            

            #Remove inactive edges which were allready defined
            if len(inequalities) > 0:    
                for ineq in inequalities:
                    result = [[act, inact] for inact in inactive_edges for act in active_edges if ineq[0] == act and ineq[1] == inact]
                    if result:
                        inactive_edges.remove(result[0][1])

            while len(inactive_edges) != 0: 
                
                current_edge = active_edges.pop()
                
                l_node = current_edge[0]
                r_node = current_edge[1]
    
                neighbor_edge = [edge for edge in inactive_edges if edge[0] == r_node or edge[1] == r_node]

                for edge in neighbor_edge:
                    inequalities.append([current_edge,edge])
                    inactive_edges.remove(edge)

                neighbor_edge = [edge for edge in inactive_edges if edge[0] == l_node or edge[1] == l_node]

                
                for edge in neighbor_edge:
                    inequalities.append([current_edge,edge])
                    inactive_edges.remove(edge)
            
        final_inequalities = []
        
        for ineq in inequalities:
            left = ineq[0]
            right = ineq[1]
            

            final_inequalities.append(f"{left} > {right}")

        return(final_inequalities)

#_________Weights_and_domainseq______#

    def get_weights(self,afp):
        
        assigned_edges = {}
        for x,step in enumerate(afp[1:],start=1): 
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
            """
            # Remove inactive edges which are allready taken care of due to the active edges
            for edge in assigned_edges:
                if edge in active_edges:
                    
                    for inactive_edge in inactive_edges:
                        logger.debug(f"inactive edge {inactive_edge} active edge {edge}")
                        if inactive_edge[0] in edge or inactive_edge[1] in edge:
                            if self.edges[inactive_edge] < self.edges[edge]:    
                                logger.debug(f"Removed: {inactive_edge}")
                                inactive_edges.remove(inactive_edge)
                                break"""

            
            for inactive_edge in inactive_edges: #doesnt go through all inactive edges
                logger.debug(f"Inactive Edge: {inactive_edge}")              
                for neighbor in self.edge_neighbors[inactive_edge]:
                    logger.debug(f"Neighbor: {neighbor}")
                    if neighbor in assigned_edges and neighbor in active_edges and self.edges[inactive_edge] < self.edges[neighbor]:
                            logger.debug(f"Removed: {inactive_edge}")
                            inactive_edges.remove(inactive_edge)
                            
            
            
            
            # Sort active edges 
            active_edges.sort(key= lambda x: (x[1]))
            logger.debug(f"Active Edges: {active_edges}")
            logger.debug(f"Inactive Edges: {inactive_edges}")  

            # Check for different substructures
            for edge in active_edges:
                if edge in assigned_edges:
                    logger.debug(f"Edge: {edge}, Edge weight: {self.edges[edge]}")
                    for inactive in inactive_edges:
                        logger.debug(f"Inactive Edge: {inactive}, Edge weight: {self.edges[inactive]}")
                        if self.edges[edge] < self.edges[inactive] and inactive_edge in assigned_edges:
                            no_assigned_neighbor_flag = False
                            for neighbor in self.edge_neighbors[inactive]:
                                logger.debug(f"Neighbor: {neighbor}Edge weight: {self.edges[neighbor]}")
                                if neighbor not in assigned_edges:
                                    no_assigned_neighbor_flag = True
                                    logger.debug(f"Neighbor: {neighbor}, assigend = {assigned_edges}")
                                

                            if not no_assigned_neighbor_flag:
                                raise SystemExit(f"\nFolding Path invalid. Following step: {step} not possible")

            while len(inactive_edges) != 0: 
                logger.debug(f"\nBegin weighting edges: {inactive_edges}")
                current_edge = active_edges.pop()
                logger.debug(f"Current Edge: {current_edge}")
                
                inactive_flag = False
                for edge in inactive_edges:
                    if edge[0] in current_edge or edge[1] in current_edge:
                        inactive_flag = True
                        logger.debug(f"Inactive edge influences edge {edge}")

                if current_edge not in assigned_edges and inactive_flag:

                    l_node = self.get_node_by_name(current_edge[0])
                    r_node = self.get_node_by_name(current_edge[1])
                    logger.debug(f"Left weight {l_node.max_weight,l_node}")
                    logger.debug(f"right weight {r_node.max_weight,r_node}")
                    edge_weight = max(l_node.max_weight,r_node.max_weight) + 1
                    self.edges[current_edge] = edge_weight
                    logger.debug(f"Edge Weight: {edge_weight} for edge: {current_edge}")

                    l_node.max_weight = edge_weight
                    r_node.max_weight = edge_weight
                    logger.debug(f"Left weight {l_node.max_weight,l_node}")
                    logger.debug(f"Right weight {r_node.max_weight,r_node}")
                    assigned_edges[current_edge] = self.edges[current_edge]

                    

                edges_to_remove = []

                for edge in inactive_edges:
                    logger.debug(f"current inactive edge {edge}")
                    if edge[0] in current_edge or edge[1] in current_edge:
                        edges_to_remove.append(edge)

                for edge in edges_to_remove:
                    inactive_edges.remove(edge)
                    logger.debug(f"Removed edge: {edge}")
                    if edge not in assigned_edges:
                        assigned_edges[current_edge] = self.edges[current_edge]
                        assigned_edges[edge] = self.edges[edge]
                        logger.debug(f"new edge: {assigned_edges}")


                
            logger.debug("--"*50)
        

        logger.debug(self.edges)
        
    def get_node_by_name(self, node_name):
        for node in self.graph:
            if node.name == node_name:
                return node
        return None

    def add_weights_to_graph(self,edge_weights):
        
        max_weights = {}
        for edge, weight in edge_weights.items():
            for node in edge:
                max_weights[node] = max(max_weights.get(node, 0), weight)
                self.edges[edge] = weight
        logger.debug(max_weights)

                
        for node in self.graph:
            if node.name in max_weights:
                node.max_weight = max_weights[node.name]
        logger.debug(self.print_nodes())
    
    def get_domain_seq(self):
        """Returns the domain level sequence of the abstract folding path in form of a List.

        Returns:
            domain_seq(list): domain level sequence
        """
        self.create_domain_seq()
        domain_seq = []

        for x,node in enumerate(self.graph):
            domain_seq.append(str(" ".join(node.prefix) + " " + node.middle + " " + " ".join(node.suffix) + " " + "l" + str(x)))
        
        logger.debug(f"Resulting domain seq: {domain_seq}")
        logger.debug(f"{self.edges}")

        return domain_seq
    
    def create_domain_seq(self):
        
        logger.debug(f"************************Begin create_domain_seq****************************")
        domains = list(string.ascii_lowercase)
        domains += [x + z for x in string.ascii_lowercase for z in string.ascii_lowercase if 'l' not in (x, z) and 'm' not in (x, z)]
        
        visited_nodes = {}
        for node in self.graph:
            
            current_node = node
            visited_nodes[current_node] = True
            logger.debug("____________________"*3)
            logger.debug(f"\n Current Node:  {current_node}\n")
            current_node.middle = "m" + str(node.connected)

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
                node.prefix = [str(domain) + "*" for domain in node.prefix]
                node.middle += "*"
                node.suffix = [str(domain) + "*" for domain in node.suffix]

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
            if node not in visited:
                if dfs(node, None):
                    return True

        return False
    
    def pseudoknot_attack_detection(self,pairtable):
            
            
            for x in range(1,len(pairtable)-1):
                secured_domains = []
                
                for i,j in enumerate(pairtable[x][1:],start=1):
                    if j > i+1:
                        secured_domains.append((i +1,j))

                

                for tuple in secured_domains:
                    if pairtable[x][tuple[0]:tuple[1]] != pairtable[x+1][tuple[0]:tuple[1]]:
                        return True
                    
            return False
    
    def different_substructures_detection(self,pairtable):
        """Detects if an occuring substructure was allready defined and can therfore not be used
        
        """


        for x,step in enumerate(pairtable[1:],start = 1):
            for i in range(1,len(step)):
                for struct in pairtable[1:x]: 
                    x = 2


      
"""Modules in Graph as nodes first initialize 


Args:
    afp(list): each entry in the list corresponds to one step in the abstract folding path the input can either be dot-bracket annotated or as a pairtable should be in list 

    
"""


def build_graph(afp):

    if afp[0][0] == ".":
        pairtable_afp = path_to_pairtablepath(afp)
    else:
        pairtable_afp = afp
    logger.debug(f"Pairtable Input: {pairtable_afp}")

    
    afp_graph = graph()

    if afp_graph.pseudoknot_attack_detection(pairtable_afp):
        raise SystemExit("Pseudoknot attack detected. Choose a path without Pseudoknot attacks")
    
    afp_graph.create_nodes_from_pairtable(pairtable_afp)

    
    logger.info("Graph after adding all Nodes")
    logger.info(afp_graph.print_nodes())

    nodes = list(afp_graph.graph.keys())
    nodes.insert(0,0)

    afp_graph.create_edges(pairtable_afp,nodes)
    
    
    logger.info("Graph after adding edges:")
    logger.info(afp_graph.print_nodes())
    afp_graph.get_edges()

    afp_graph.set_edge_neighbors()


    #Check if there are cycles present in the graph
    if afp_graph.cycle_detection():
        raise SystemExit("Cycle in graph detected. Currently no viable solutions for this case.")
    
    logger.info(f"Edges: {afp_graph.edges}")

    connected_components = find_connected_modules(pairtable_afp)
    logger.info("\nFollowing Nodes are connected:")
    for component in connected_components:
        nodes_in_component = set()
        for node_index in component:
            nodes_in_component.add(nodes[node_index])
        logger.debug(component)

    if not afp_graph.bipartite_check(connected_components):
        raise ImportError("Graph not Bipartite can't design domain level sequence. Check your input")

    #print("\nBipartite Check Complete:")
    #for x in afp_graph.graph:
        #print(x)
    
    afp_graph.get_weights(afp=pairtable_afp)


    return afp_graph

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
            lines = file.readlines()
            lines = [line.strip() for line in lines if re.match(r'^[0-9,().]*$', line) and not line.startswith('#')]
            
            if not lines:
                raise SystemExit("File is empty. Here is your sequence: ")

            return lines


    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return []

def main():
        #_________________Argparse_________________#
        parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            description="Copaths path_to_seq: generates a domain level sequence based on an abstract folding path"
        )

        parser.add_argument("-i","--input",help="Reads txt file as input if not specified user can input via console.")
        parser.add_argument("-v", "--verbose",type = int, default = 0,
            help = "Track process by writing verbose output to STDOUT during calculations.")
        
        args = parser.parse_args()

        #_________________Set_logger_verbosity________________#

        set_verbosity(logger,args.verbose)
        

        #_________________Input_Parsing_________________#
        
        
        if args.input == None:
            afp = []
            while True:
                print("\n")
                print(f"Current Input: {afp}")
                print("Please input a folding path in dot-bracket annotation or use '$' to finish and continue or use '@' to exit :")
                user_input = input()
                # Check for exit conditions
                if user_input == "@":
                    print("\nExiting copaths")
                    exit()
                elif user_input == "$":
                    print(f"\n\nFinal Input:\n{afp}\n\n")
                    break
                elif user_input == "r" or user_input == "R":
                    afp = []
                    print("Input cleared")
                    continue
                
                if is_balanced_structure(user_input):

                    # Check if the user input contains only ".", "(", and ")"
                    
                    if all(char == "." or char in ("(", ")") for char in user_input):
                        if len(user_input) == len(afp) + 1:
                            afp.append(user_input)
                        else:
                            print("Please add 1 character per step")
                    else:
                        print("Error: Invalid character in the folding path. Only '.', '(', and ')' are allowed.")
                else:
                    print("Structure is not balanced -> closing/opening brackets don't match")
                
        else:
            afp = input_parser(args.input)
            print(f"\n\nInput folding path:\n{afp}\n\n")


        #_________________Creation of Graph_________________#
        
        
        logger.debug(afp)

        afp_graph = build_graph(afp)
        

        print("Resulting Domain Level sequence: "," ".join(afp_graph.get_domain_seq()))

        



if __name__ == "__main__":

    main()