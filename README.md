# copaths

## Other name ideas
- CocoPaths - general package
- CocoFold - path to domain seq generator
- CocoTranslator
- copat
- CoFoldify
- CoTranscripter
- CoPathCompiler


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

afterwards use pytest to verify that all tests are working:
```bash
$ pytest
```


# Pseudocode for Copaths

## Initialize the graph

1. **Input**: Abstract Folding Path (AFP)
2. Create an empty graph G
3. Create an empty list of connected graphs connected_graphs
4. Create a node for each step in the AFP, naming them based on the step
5. Initialize each node with the following properties:
   - `node.prefix = ""`
   - `node.middle = ""`
   - `node.weight = 1`
   - `node.complement = False`
   - `node.connected = None`
   - `node.neighbors = []`

## Pseudoknot Detection

6. Check for Pseudoknot-attack in the AFP

## Node Creation and Edge Addition

7. For each step in reverse(AFP):

    8. If two nodes are paired in the AFP:
        
        9. Create an edge between the paired nodes and name it based on the nodes
        10. Add the edge to graph G

## Find Connected Graphs

11. Divide G into connected subgraphs and store them in connected_graphs

## Complementary Node Assignment

12. For each connected_graph in connected_graphs:

    13. Create a stack to keep track of nodes to visit
    14. Push the first node onto the stack
    15. Create an empty set to store visited nodes
    16. While the stack is not empty:

        17. Pop a node from the stack and mark it as visited
        18. Assign the node as complementary
        19. For each unvisited neighbor of the current node:

            20. Push the neighbor onto the stack

    21. If not all nodes are visited:

        22. Raise an error (complementary assignment not possible)

## Weight Assignment

23. Create a dictionary called assigned_edges
24. For each step in AFP:

    25. Create two empty lists: active_edges and inactive_edges
    26. For each edge in the graph:

        27. If the two nodes indicated by the edge are currently paired:

            28. Add the edge to active_edges
        29. Else:

            30. Add the edge to inactive_edges
    31. For each inactive_edge in inactive_edges (remove inactive edges that are already taken care of due to active edges):

        32. If the inactive_edge neighbors an active_edge in assigned_edges and the weight of the inactive_edge is less than active_edge:

            33. Remove inactive_edge from inactive_edges

    34. Sort active_edges based on occurrence in ascending order

35. For each active_edge in active_edges:

    36. If the edge is in assigned_edges:

        37. For each inactive_edge in inactive_edges:

            38. If the weight of active_edge is less than inactive_edge and inactive_edge is in assigned_edges:

                39. Set no_assigned_neighbor_flag to False
                40. For each neighbor of the inactive_edge:

                    41. If the neighbor is not in assigned_edges:

                        42. Set no_assigned_neighbor_flag to True

                43. If not no_assigned_neighbor_flag:

                    44. Raise a SystemExit

45. While inactive_edges are not empty:

    46. Pop an edge from active_edges (current_edge)
    47. For each edge in inactive_edges:

        48. If the edge neighbors the current_edge and the current_edge is not in assigned_edges:

            49. Get the weights of nodes from the current_edge
            50. Set the weight of the current_edge to the maximum weight of the nodes + 1
            51. Update the node weights to the maximum weight
            52. Add the current_edge and edge to assigned_edges
            53. Remove inactive_edge from inactive_edges if it neighbors the current_edge

## Domain Sequence Creation

54. Create a list called domains consisting of all possible 2-char combinations of the alphabet
55. Remove entries in domains that contain 'l' or 'm'
56. Create a list called visited_nodes
57. For each node in the graph G:

    58. Append the node to visited_nodes
    59. Set node.middle to 'm' + the index of connected_components
    60. If the weight of the node is greater than 1 and the length of node.prefix is less than node.weight - 1:

        61. Update node.prefix by adding domains[node.weight - 1 - len(node.prefix)] to the beginning
        62. Remove the used entries from domains
        63. Update node.suffix by adding domains[node.weight - 1 - len(node.suffix)] to the end
        64. Remove the used entries from domains

    65. For each neighbor in node.neighbors:

        66. If the neighbor is not in visited_nodes and the neighbor.connected is equal to node.connected:

            67. Calculate the shared_weight as the weight of the edge connecting both nodes
            68. If shared_weight is greater than 1:
            
                69. Update neighbor.prefix with the last (shared_weight - 1) characters of node.suffix in reverse order
                70. Update neighbor.suffix with the last (shared_weight - 1) characters of node.prefix in reverse order

71. For each node in G:

    72. If the node is complementary:

        73. Append a marker (e.g., '*') after every domain

74. Create an empty list called Domain_seq
75. For each node in the graph:

    76. Concatenate node.prefix, node.middle, and node.suffix into a single
