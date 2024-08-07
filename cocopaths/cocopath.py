from typing import Any
import re,argparse, logging, sys
from cocopaths import __version__
from .utils import (path_to_pairtablepath, find_connected_modules,is_balanced_structure,acfp_terminal_input)
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

    def print_nodes(self):
        return "\n".join([str(node) for node in self.graph])

    def create_edges(self,acfp,nodes):
        for step in acfp[::-1]:
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
                if not self.is_bipartite_dfs(neighbor, visited_nodes, not complement, connected,x):
                    return False
                
            else:
                if neighbor.complement == node.complement:
                    return False

        return True

    def bipartite_check(self, connected_components):
        visited_nodes = {}
        logger.debug(f"\n\n\n\n\n\n\nConnected_components {connected_components}\n\n\n\n\n\n\n\n")
        for x, connected in enumerate(connected_components):
            for node in self.graph: 
                if node.name in connected and node not in visited_nodes:
                    if not self.is_bipartite_dfs(node, visited_nodes, True, connected,x):
                        return False
        return True

    def get_current_edges(self,step_number):
        current_edges = []
        for edge in self.edges:
            if edge[1] <= step_number + 1:
                current_edges.append(edge)
        return current_edges
    
 
#_________Weights______#
    def get_weights(self,acfp):
        assigned_edges = {}
        for x,step in enumerate(acfp[1:],start=1): 
            collected_edges = self.get_current_edges(x)
            logger.debug("--"*50)
            logger.debug(f"\nAssigned Edges: {assigned_edges}\nCollected Edges: {collected_edges}\nCurrent Step: {step}")
            active_edges = []
            inactive_edges = []
            
            # Get active and inactive edges
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
                            
            for edge in inactive_edge_to_remove:
                inactive_edges.remove(edge)
                if edge not in assigned_edges:
                    assigned_edges[edge] = 1
            
            # Sort active edges 
            active_edges.sort(key= lambda x: (x[1]))
            logger.debug(f"Active Edges: {active_edges}")
            logger.debug(f"Inactive Edges: {inactive_edges}")  

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

            while len(inactive_edges) != 0:
                curr_edge = inactive_edges.pop()
                current_neighbors = [neighbor for neighbor in self.edge_neighbors[curr_edge] if neighbor in active_edges]
                max_weight = float("-inf")
                for neighbor in current_neighbors:
                    if self.edges[neighbor] > max_weight:
                        max_weight = self.edges[neighbor]

                if self.edges[curr_edge] > max_weight:
                    raise SystemExit(f"Folding path not possible in this step: {step}\n\n Goodbye")

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
            logger.debug(node)
            if node not in visited:
                if dfs(node, None):
                    return True
        return False
    
    def pseudoknot_attack_detection(self,pairtable):
            for x in range(1,len(pairtable)-1):
                secured_domains = []
                for i,j in enumerate(pairtable[x][1:],start=1):
                    if j > i+1 and pairtable[x][i + 1] != 0 and pairtable[x][j] != 0:
                        secured_domains.append((i +1,j-1))

                logger.debug(f"secured domains {secured_domains}")
                for tuple in secured_domains:
                    if pairtable[x][tuple[0]:tuple[1]] != pairtable[x+1][tuple[0]:tuple[1]]:
                        if pairtable[x + 1][tuple[0]] != 0 and pairtable[x + 1][tuple[1]] != 0:
                            return True
            return False

    def check_unpaired_complements(self,acfp):
        """
        Checks if there are two unpaired domains that are complementary in the given domain sequence.

        Parameters:
        - acfp (list): Abstract folding path in the form of a pairtable.

        Raises:
        - SystemExit: If two unpaired complementary domains are found.
        """
        for step in acfp: 
            for x in range(1,step[0]):
                if step[x] == 0: 
                    for y in range(x + 1,step[0]):
                        if self.get_node_by_name(x).connected == self.get_node_by_name(y).connected and step[y] == 0:
                            raise SystemExit(f"Two unpaired complements in the folding path: {acfp}. Please pair them.")

    def consistency_check(self,acfp):
        """Applies the consistency check on the input acfp. If it doesn't pass it raises a SystemExit. 

        Args:
            acfp(list): abstract folding path in the form of a list
        """
        for i in range(1,len(acfp)):
            if acfp[i][-1] == "." and acfp[i][0:i] != acfp[i-1]:
                raise SystemExit(f"CFolding path is not consitent in following step: {acfp[i]}-{acfp[i-1]}")

def build_graph(acfp):
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
    logger.debug(f"Pairtable Input: {pairtable_acfp}")
    acfp_graph = graph()
    acfp_graph.create_nodes_from_pairtable(pairtable_acfp)

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
        
    if not acfp_graph.bipartite_check(connected_components):
        raise ImportError("Graph not Bipartite can't design domain level sequence. Check your input\n\n Goodbye")

    if acfp_graph.pseudoknot_attack_detection(pairtable_acfp):
            raise SystemExit("SystemExit: Pseudoknot attack detected. Choose a path without Pseudoknot attacks")
    
    acfp_graph.check_unpaired_complements(pairtable_acfp)
    acfp_graph.consistency_check(acfp)
    graph.print_nodes
    acfp_graph.get_weights(acfp=pairtable_acfp)
    acfp_graph.verify_weights(acfp=pairtable_acfp)

    return acfp_graph

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

        #_________________Creation of Graph_________________#
        logger.debug(acfp)
        acfp_graph = build_graph(acfp)
        print(f"Resulting Domain Level sequence: {' '.join(acfp_graph.get_domain_seq())}")
        print('Length of domain seq = ', len(" ".join(acfp_graph.get_domain_seq()).split()))

if __name__ == "__main__":

    main()

    
