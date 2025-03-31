# Cocopaths - A Compiler for cotranscriptional folding pathways

[![codecov](https://codecov.io/gh/lmiksch/cocopaths/graph/badge.svg?token=6PVQSOEK8R)](https://codecov.io/gh/lmiksch/cocopaths)
## Getting Started
To install the package you need to clone the repository using git: 

```bash
# Clone this repository
$ git clone https://github.com/lmiksch/cocopaths

#Change your directory 
$ cd copaths
```

The package can then be installed via pip:

```bash
$ pip install .
```

For development purpose:

```bash
$ pip install -e .
```

Afterwards use pytest to verify that all tests are working:
```bash
$ pytest
```

After installation you can try if the script is installed by using:
```bash
$ cocopath --help
```

The script is dependent on following libraries: 

[Infrared Software](https://www.lix.polytechnique.fr/~will/Software/Infrared/Doc/index.html)

[ViennaRNA](https://www.tbi.univie.ac.at/RNA/)


The script gets called by using:
```bash
$ cocopath
```

the abstract folding path(AFP) can be either put in as a .txt file or directly in the terminal

The AFP must have following propperties to be translated to a domain level sequence:
  - Pseudoknot-Free
  - If a structure was defined it can't be changed if no additional step influences the defined substructure

# How to use cocopaths and cocosim 


Cocopaths and cocosim can be used in a pipeline to simulate the resulting domain level sequence

```bash
$ cocopath -i test.txt | cocosim 
```

In the examples directory, you can find sample inputs for both cocopaths and cocosim.

A typical command might look like this:

```bash
$ cocopath -i fp_1.txt | cocosim 
```

Using only the command without input will prompt you to enter the domain-level sequence in the terminal.

```bash
$ cocosim 
```

```bash
$ cocosim 
Please enter a domain level sequence:
```

cocopath also checks stdin and accepts if each step is seperated by a ','. 

## Cocosim 

Cocosim can be used separately from cocopaths. It accepts input either through the terminal or in .pil format. An example of the correct .pil format can be found in the examples folder. 

For more information about CocoSim, please refer to its dedicated GitHub repository: https://github.com/lmiksch/CocoSim

## CocoDesign

CocoDesign is currently designed to operate independently. It accepts either a .pil format file as input, where the user can specify the domain lengths, or an input via the terminal, with the default parameters for domain lengths.
CocoDesign initiates by validating the cotranscriptional folding path at the domain level. If the simulated path does not align with the abstract folding path, it will not continue unless the -force flag is activated. 


  
