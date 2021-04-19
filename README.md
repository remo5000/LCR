## Synopsis

This is a fork of https://github.com/DeLaChance/LCR .
Please refer to `README.md.old` if you wish to cite the original author.

This project contains several methods for answering label-constraint reachability queries (LCR-queries).
The major contribution of this fork is the BackboneIndex.

## Code Example

The main part of the code is written in C++ and there are scripts in bash/sh and python3.

The C++-part has:
- a definition of a labelled graph (Graph)
- a set of indices to answer LCR-queries (Index/Unbounded)
- a number of tests (tests/Index/Unbounded) 
- scripts to run the experiments or generate queries for the experiments (experiments/Index/Unbounded).

The Python/sh parts consist of:
- a script to generate a synthetic graph under a given model with a specified number of vertices. (datagen)
- a script to download real graphs and label them randomly (experiments...install.sh)

## Installation

For C++ you need:
- g++ (>=9.3.0) or clang (>=12.0.0)
- OpenMP libraries

For python3 you need:
- snap-stanford

## Build

The code can be built by running:

```sh
cd LCRIndexing
sh ./build.sh
```

You can rebuild using `rebuild.sh`.

## Experiment Driver

To run experiments, we use `runExperiment.sh` as an entry point to `runExperiment.cc`.
