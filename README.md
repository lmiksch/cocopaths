
# Cocopaths - A Compiler for Cotranscriptional Folding Pathways

CocoPaths is a toolkit for designing RNA sequences that follow predefined cotranscriptional folding pathways. It uses a domain-level abstraction to translate high-level folding paths into concrete sequences and includes three tools: CocoPath for compiling folding paths, CocoSim for simulating folding behavior, and CocoDesign for generating nucleotide-level sequences. 

For more information, see: https://utheses.univie.ac.at/detail/75721/

[![codecov](https://codecov.io/gh/lmiksch/cocopaths/graph/badge.svg?token=6PVQSOEK8R)](https://codecov.io/gh/lmiksch/cocopaths)

## Getting Started

To install the package, first clone the repository using Git:

```bash
# Clone this repository
git clone https://github.com/lmiksch/cocopaths

# Change into the project directory
cd cocopaths
```

Then install the package using pip:

```bash
pip install .
```

For development purposes:

```bash
pip install -e .
```

The package depends on the following libraries:

- [Infrared Software](https://www.lix.polytechnique.fr/~will/Software/Infrared/Doc/index.html)
- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/)

> ⚠️ Note: Infrared needs to be installed separately.

After installation, run pytest to verify that all tests are passing:

```bash
pytest
```

To verify that the script was installed correctly:

```bash
cocopath --help
```

You can invoke the script using:

```bash
cocopath
```

The abstract folding path (AFP) can be provided either as a `.txt` file or entered directly in the terminal.

## How to Use Cocopaths and Cocosim

Cocopaths and Cocosim can be used in a pipeline to simulate the resulting domain-level sequence:

```bash
cocopath -i test.txt | cocosim
```

In the `examples` directory, you can find sample inputs for both Cocopaths and Cocosim.

A typical command might look like this:

```bash
cocopath -i fp_1.txt | cocosim
```

Running the tool without input will prompt you to enter the domain-level sequence in the terminal:

```bash
cocosim
```

```
Please enter a domain level sequence:
```

Cocopath also checks `stdin` and accepts input where each step is separated by a comma.

## Cocosim

Cocosim can be used separately from Cocopaths. It accepts input either directly via the terminal or in `.pil` format. An example of the correct `.pil` format can be found in the `examples` folder.

For more information about CocoSim, please refer to its dedicated GitHub repository: https://github.com/lmiksch/CocoSim

## CocoDesign

CocoDesign is currently designed to operate independently. It accepts either a `.pil` format file as input—where the user can specify domain lengths—or input via the terminal with default domain lengths.

CocoDesign starts by validating the cotranscriptional folding path at the domain level. If the simulated path does not align with the abstract folding path, it will not continue unless the `--force` flag is used.

An example workflow might look like this:

```bash
cocodesign
```

You will be prompted to enter the domain-level sequence:

```
Please enter a domain level sequence:
```

After entering the domain sequence:

```
Please input a domain level sequence:
e* a* L0* b* f* Z0 c* b L0 a d* Z1 d a* L0* b* c Z2 f b L0 a e Z3
Please input the aCFP by which the domain level sequence was designed:

Please input the full folding path as space-separated dot-bracket structures.
Example: . () .() (())
. () .() (())
```

Afterward, a domain-level sequence will be generated.

Domain lengths can be adjusted using the `.pil` format to fine-tune each domain individually. Note that for the design logic to remain valid, domain lengths should not be changed drastically. It is recommended to first adjust the domain lengths using Cocosim to verify whether the design is viable at the domain level.