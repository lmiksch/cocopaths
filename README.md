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
- Make Graph etc. 
- Create a dictionary called assigned edges
- For step in AFP 
    - collect active and inactive edges and put both in a seperated list
    - Remove inactive edges which neighbor active edges that are in assigned edges
    - sort active edges 
    - While inactive edges != empty:
        - current edge = active edge.pop()
        - if an inactive edges neighbors a current edge:
            - get weight of nodes from current edges 
            - set current edgeweight to the max weight of the nodes + 1 
            - update the node weights to the max weight 
        - add current edge to assigned edges
        - remove inactive edges that neighbor current edge






