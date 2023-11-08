# Cocopaths - A Compiler for contranscriptional folding pathways

[![codecov](https://codecov.io/gh/lmiksch/cocopaths/graph/badge.svg?token=6PVQSOEK8R)](https://codecov.io/gh/lmiksch/cocopaths)
## Getting Started
To install the package you need to clone the repository using git: 

```bash
# Clone this repository
$ git clone https://github.com/lmiksch/CoFPT

#Change your directory 
$ cd copaths
```

The package can then be installed via pip:

```bash
$ pip .
```

For development purpose:

```bash
$ pip -e .
```

Afterwards use pytest to verify that all tests are working:
```bash
$ pytest
```

After installation you can try if the script is installed by using:
```bash
$ copaths --help
```

The script gets called by using:
```bash
$ copaths 
```

the abstract folding path(AFP) can be either put in as a .txt file or directly in the terminal

The AFP must have following propperties to be translated to a domain level sequence:
  - Pseudoknot-Free
  - Structure must be fully saturated
  - If a structure was defined it can't be changed if no additional step influences the defined substructure



# Pseudocode for Copaths

## Initialize the graph

1. **Input**: Abstract Folding Path (AFP)
2. Create an empty graph G
3. Create an empty list of connected graphs connected_graphs
4. Create a node for each step in the AFP, naming them based on the step
   - Initialize each node with the following properties:
     - `node.prefix = ""`
     - `node.middle = ""`
     - `node.weight = 1`
     - `node.complement = False`
     - `node.connected = None`
     - `node.neighbors = []`

## Pseudoknot Detection

5. Check for Pseudoknot-attack in the AFP

### Node Creation and Edge Addition

6. For each step in reverse(AFP):
   - 7 If two nodes are paired in the AFP:
     - 8 Create an edge between the paired nodes and name it based on the nodes
     - 9 Add the edge to graph G

### Find Connected Graphs

10. Divide G into connected subgraphs and store them in connected_graphs

### Complementary Node Assignment

11. For each connected_graph in connected_graphs:
   - 12 Create a stack to keep track of nodes to visit
   - 13 Push the first node onto the stack
   - 14 Create an empty set to store visited nodes
   - 15 While the stack is not empty:
      - 16 Pop a node from the stack and mark it as visited
      - 17 Assign the node as complementary
      - 18 For each unvisited neighbor of the current node:
        - 19 Push the neighbor onto the stack
   - 20 If not all nodes are visited:
     - 21 Raise an error (complementary assignment not possible)

### Weight Assignment

22. Create a dictionary called assigned_edges
23. For each step in AFP:
   - 24 Create two empty lists: active_edges and inactive_edges
   - 25 For each edge in the graph:
     - 26 If the two nodes indicated by the edge are currently paired:
       - 27 Add the edge to active_edges
     - 28 Else:
       - 29 Add the edge to inactive_edges
   - 30 For each inactive_edge in inactive_edges (remove inactive edges that are already taken care of due to active edges):
     - 31 If the inactive_edge neighbors an active_edge in assigned_edges and the weight of the inactive_edge is less than active_edge:
       - 32 Remove inactive_edge from inactive_edges
   - 33 Sort active_edges based on occurrence in ascending order

34. For each active_edge in active_edges, that is also in assigned: # Different Substructure check:
    - 35 For each inactive_edge in inactive_edges, that is also in assigned:
      - 36 If the active_edge.weight < inactive_edge.weight:
          - 37 If every neighbor of the inactive_edge is assigned:
              - 38 Raise a SystemExit 

39. While inactive_edges are not empty:
   - 40 Pop an edge from active_edges (current_edge)
   - 41 For each edge in inactive_edges:
     - 42 If the edge neighbors the current_edge and the current_edge is not in assigned_edges:
       - 43 Get the weights of nodes from the current_edge
       - 44 Set the weight of the current_edge to the maximum weight of the nodes + 1
       - 45 Update the node weights to the maximum weight
       - 46 Add the current_edge and edge to assigned_edges
       - 47 Remove inactive_edge from inactive_edges if it neighbors the current_edge

### Domain Sequence Creation

48. Create a list called domains consisting of all possible 2-char combinations of the alphabet
49. Remove entries in domains that contain 'l' or 'm'
50. Create a list called visited_nodes
51. For each node in the graph G:
   - 52 Append the node to visited_nodes
   - 53 Set node.middle to 'm' + the index of connected_components
   - 54 If the weight of the node is greater than 1 and the length of node.prefix is less than node.weight - 1:
     - 55 Update node.prefix by adding domains[node.weight - 1 - len(node.prefix)] to the beginning
     - 56 Remove the used entries from domains
     - 57 Update node.suffix by adding domains[node.weight - 1 - len(node.suffix)] to the end
     - 58 Remove the used entries from domains
   - 59 For each neighbor in node.neighbors:
     - 60 If the neighbor is not in visited_nodes and the neighbor.connected is equal to node.connected:
       - 61 Calculate the shared_weight as the weight of the edge connecting both nodes
       - 62 If shared_weight is greater than 1:
         - 63 Update neighbor.prefix with the last (shared_weight - 1) characters of node.suffix in reverse order
         - 64 Update neighbor.suffix with the last (shared_weight - 1) characters of node.prefix in reverse order

65. For each node in G:
   - 66 If the node is complementary:
     - 67 Append a marker (e.g., '*') after every domain


74. Create an empty list called Domain_seq
75. For each node in the graph:
   - 76 Concatenate node.prefix, node.middle, and node.suffix into a single

