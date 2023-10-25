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


## Pseudocode for cocofold

#Initiation of Graph

- Take an Abstract Folding Path(AFP) as input
- Check for Pseudoknot attack in the AFP
- Create a node for each step in the AFP
- From the End to the Start of the AFP:
    - If two nodes are paired in the AFP:
        - Create an edge between them and name it based on the two nodes it connects
        - Add edge to graph
- Create a list where each sublist represents the connected graphs
- For each connected Graph:
    - start at the first node and then walk through the graph and assign every other node as complementary till all nodes are visited
    - If not possible raise error
 
#Weight assignment

- Create a dictionary called assigned edges
- For step in AFP 
    - collect active and inactive edges and put both in a list 
    - Remove inactive edges which neighbor active edges that are in assigned edges and 
    - sort active edges based on occurence ascending
    - While inactive edges != empty:
        - current edge = active edge.pop()
        - if an inactive edges neighbors a current edge:
            - get weight of nodes from current edges 
            - set current edgeweight to the max weight of the nodes + 1 
            - update the node weights to the max weight 
        - add current edge to assigned edges
        - remove inactive edges that neighbor current edge

#Domain seq creation

- 




