from typing import Any


try:
    from collections import deque
    from convert_functions import (path_to_pairtablepath, find_connected_modules)
    import string
    import inequality_solver
except:
    pass

class node():
    def __init__(self,name,prefix ="",suffix ="",complement= False,max_weight = 0,neighbors = None,connected = None):
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

    # dont know if necessary
    def set_middle_domain(self,middle_domain): 
         self.middle = middle_domain
    

class graph():


    def __init__(self):
        self.graph = {}
        self.edges = dict()
    
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
            print(node.name)

            for neighbor in node.neighbors:
                if neighbor.data not in visited:
                    dfs_recursive(neighbor)

        dfs_recursive(start_node)

    def print_nodes(self):
        for node in self.graph:
            print(node)

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
                self.edges[(node.name , neighbor.name)] = 0 #default weight of edge is 1 since they will always bind with middle domain


        #Remove duplicate edges
        unique_edges = {}
        for edge in self.edges:
            if edge[1] >= edge[0]:
                unique_edges[edge] = self.edges[edge]

        self.edges = unique_edges


    

    def is_bipartite_dfs(self, node, visited_nodes, complement, connected,x):
        visited_nodes[node] = True
        node.complement = complement
        
        #print("is bipartite\n")
        #print("\nVisited Nodes")
        #for x in visited_nodes.keys():
            #print("visitded nodes",x)
        #print("")
        for neighbor in node.neighbors:
            #print("node: ",node)
            #print("neighbor: ",neighbor)
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
           # print("\nFollowing component will be analyzed", connected)    
            for node in self.graph: 
                #print("nodename",node.name)
                #print("node name in connected",node.name in connected)
                if node.name in connected and node not in visited_nodes:
                    #print("here")
                    node.connected = x
                    if not self.is_bipartite_dfs(node, visited_nodes, True, connected,x):
                        #print("h2")
                        return False

        return True




    def get_current_edges(self,step_number):
        current_edges = []
        for edge in self.edges:
            if edge[1] <= step_number + 1:
                current_edges.append(edge)
        return current_edges


    def get_inequalities(self,afp):
        
        inequalities = []

        for x,step in enumerate(afp): 
            collected_edges = self.get_current_edges(x)
            print("-"*50)
            print("\nCollected Edges",collected_edges)
            print("Current Step:",step)
            active_edges = []
            inactive_edges = []
            
            for edge in collected_edges:
                if step[edge[0]] == edge[1] and step[int(edge[1])] == int(edge[0]):
                    active_edges.append(edge)
                else:
                    inactive_edges.append(edge)

            active_edges.sort(key= lambda x: (x[1]))
            print("Active Edges",active_edges)
            print("Inactive Edges",inactive_edges)  

            
            #Remove inactive edges which were allready defined
            if len(inequalities) > 0:    
                for ineq in inequalities:
                    result = [[act, inact] for inact in inactive_edges for act in active_edges if ineq[0] == act and ineq[1] == inact]
                    print(f"For ineq: {ineq} was following result found: {result}")
                    if result:
                        print(f"Remove following edge {result[0][1]} due to result being {result}")
                        inactive_edges.remove(result[0][1])

            while len(inactive_edges) != 0: 
                print("\nBegin ineq inactive edges:",inactive_edges)
                current_edge = active_edges.pop()
                print(f"Current Edge: {current_edge}")

        

                l_node = current_edge[0]
                r_node = current_edge[1]
                
                

                neighbor_edge = [edge for edge in inactive_edges if edge[0] == r_node or edge[1] == r_node]

                print(f"R neighbor edge {neighbor_edge}")
                for edge in neighbor_edge:
                    inequalities.append([current_edge,edge])
                    inactive_edges.remove(edge)

                neighbor_edge = [edge for edge in inactive_edges if edge[0] == l_node or edge[1] == l_node]

                print(f"L neighbor edge: {neighbor_edge}")
                for edge in neighbor_edge:
                    inequalities.append([current_edge,edge])
                    inactive_edges.remove(edge)
            print("-"*50)
                

            print(f"\nInequalities {inequalities}")

        final_inequalities = []
        print("\nResulting Inequalties: ",inequalities)
        for ineq in inequalities:
            left = ineq[0]
            right = ineq[1]
            

            final_inequalities.append(f"{left} > {right}")

        return(final_inequalities)
    
    def get_weights(self,afp):
        
        assigned_edges = {}
        for x,step in enumerate(afp): 
            collected_edges = self.get_current_edges(x)
            print("-"*50)
            print("\nCollected Edges",collected_edges)
            print("Current Step:",step)
            active_edges = []
            inactive_edges = []
            

            for edge in collected_edges:
                if step[edge[0]] == edge[1] and step[int(edge[1])] == int(edge[0]) and edge not in assigned_edges:
                    active_edges.append(edge)
                else:
                    inactive_edges.append(edge)

            active_edges.sort(key= lambda x: (x[1]))
            print("Active Edges",active_edges)
            print("Inactive Edges",inactive_edges)  

            

            while len(active_edges) != 0: 
                print("\nBegin weighting edges:",inactive_edges)
                current_edge = active_edges.pop()
                print(f"Current Edge: {current_edge}")

                l_node = self.get_node_by_name(current_edge[0])
                r_node = self.get_node_by_name(current_edge[1])
                print(f"Left weight {l_node.max_weight,l_node}")
                print(f"right weight {r_node.max_weight,r_node}")
                edge_weight = max(l_node.max_weight,r_node.max_weight) + 1
                self.edges[current_edge] = edge_weight
                print("Edge Weight: ", edge_weight)

                l_node.max_weight = edge_weight
                r_node.max_weight = edge_weight
                print(f"Left weight {l_node.max_weight,l_node}")
                print(f"right weight {r_node.max_weight,r_node}")
                assigned_edges[current_edge] = True
            


                
            print("-"*50)
        

        print(self.edges)
        

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
        print(max_weights)

                
        for node in self.graph:
            if node.name in max_weights:
                node.max_weight = max_weights[node.name]
        self.print_nodes()
    
    def get_domain_seq(self):
        """Returns the domain level sequence of the abstract folding path in form of a List.

        Returns:
            domain_seq(list): domain level sequence
        """

        domain_seq = []

        for node in self.graph:
            domain_seq.append(str(" ".join(node.prefix) + " " + node.middle + " " + " ".join(node.suffix)))
        
        return domain_seq

"""Modules in Graph as nodes first initialize 


Args:
    afp(list): each entry in the list corresponds to one step in the abstract folding path the input can either be dot-bracket annotated or as a pairtable should be in list 

    
"""


def build_graph(afp):

    if afp[0][0] == ".":
        pairtable_afp = path_to_pairtablepath(afp)
    else:
        pairtable_afp = afp
    print(pairtable_afp)


    afp_graph = graph()
    afp_graph.create_nodes_from_pairtable(pairtable_afp)

    #print(afp_graph)

    print("Graph after adding all Nodes")
    afp_graph.print_nodes()

    nodes = list(afp_graph.graph.keys())
    nodes.insert(0,0)

    afp_graph.create_edges(pairtable_afp,nodes)
    
    # Print the graph after adding edges
    print("Graph after adding edges:")
    afp_graph.print_nodes()
    afp_graph.get_edges()
    
    print("Edges: ", afp_graph.edges)

    connected_components = find_connected_modules(pairtable_afp)
    print("\nFollowing Nodes are connected:")
    for component in connected_components:
        nodes_in_component = set()
        for node_index in component:
            nodes_in_component.add(nodes[node_index])
        print(component)

    if not afp_graph.bipartite_check(connected_components):
        raise ImportError("Graph not Bipartite can't design domain level sequence. Check your input")

    print("\nBipartite Check Complete:")
    for x in afp_graph.graph:
        print(x)
    """
    inequalities = afp_graph.get_inequalities(afp=pairtable_afp)

    
    final_edges = afp_graph.get_current_edges(len(pairtable_afp))
    print("\n Final Inequalities",inequalities)
    inequalities_solutions = inequality_solver.inequality_solver(inequalities,final_edges)
    """
    afp_graph.get_weights(afp=pairtable_afp)
    print("\n")


    return afp_graph

def create_domain_seq(graph):
    """
    for connected in connected_components:
        print(connected_components)"""
    
    domains = list(string.ascii_lowercase)
    domains += [x + z for x in string.ascii_lowercase for z in string.ascii_lowercase if 'l' not in (x, z) and 'm' not in (x, z)]
    
    visited_nodes = {}
    for node in graph.graph:
        
        current_node = node
        visited_nodes[current_node] = True
        print("____________________"*3)
        print(f"\n Current Node:  {current_node}\n")
        current_node.middle = "m" + str(node.connected)

        if current_node.max_weight > 1 and len(current_node.prefix) < current_node.max_weight - 1:
            prefix = domains[:current_node.max_weight - 1 - len(current_node.prefix)] + list(current_node.prefix)
            domains = domains[current_node.max_weight - 1 - len(current_node.prefix):]
            current_node.prefix = prefix

            suffix = list(current_node.suffix) + domains[:current_node.max_weight - 1 - len(current_node.suffix)] 
            domains = domains[current_node.max_weight - 1 - len(current_node.suffix):]
            current_node.suffix = suffix
        print("current Node",current_node)
        
        for neighbor in current_node.neighbors:
            if neighbor not in visited_nodes and neighbor.connected == current_node.connected:
                
                print(f"graph edges :{graph.edges}\n")
                shared_weight = graph.edges[(min(current_node.name, neighbor.name), max(current_node.name, neighbor.name))]
                print(f"Shared Weight: {shared_weight}, Edges:{(current_node.name,neighbor.name)}")
                neighbor.prefix = current_node.suffix[:shared_weight][::-1]
                neighbor.suffix = current_node.prefix[:shared_weight][::-1]
                
                
                
                
                #break
            

                print(neighbor)
        
    

    #Need to implement complementary domains 

    for node in graph.graph:
        if node.complement:
            node.prefix = [str(domain) + "*" for domain in node.prefix]
            node.middle += "*"
            node.suffix = [str(domain) + "*" for domain in node.suffix]


    graph.print_nodes()
    print("done")


    def main():
        #write import
        print("main")



if __name__ == "__main__":
    #afp = [[1, 0], [2, 2, 1], [3, 0, 3, 2], [4, 4, 3, 2, 1]]

    afp = [[1,0],[2,2,1],[3,0,3,2],[4,4,3,2,1],[5,0,3,2,5,4],[6,6,3,2,5,4,1]]
    #afp = [[1,0],[2,2,1],[3,0,3,2],[4,4,3,2,1],[5,5,0,4,3,1]] currently 
    #afp = [".","()",".()","(())","(.())","(())()",".()(())",".(.(()))"]
    #afp = [".","()",".()","()()",".()()"]
    #afp = [".","()","().","()()","().()","()(())"]
    
    afp_graph = build_graph(afp)

    
    
    
    create_domain_seq(afp_graph)
    domain_seq = afp_graph.get_domain_seq()

    print("___________"*3)
    print("Resulting domain level sequence: \n",domain_seq)
    