## Synopsis

This is the code + thesis for my graduation project at the Technical University of Eindhoven, July 2016.

This project contains several methods for answering label-constraint reachability queries (LCR-queries). Our major contribution is the LandmarkedIndex-method. 


“
To cite this work, please use the following reference.

Landmark indexing for evaluation of label-constrained reachability queries. Lucien Valstar, George Fletcher, Yuichi Yoshida. SIGMOD 2017. http://doi.acm.org/10.1145/3035918.3035955

For further details, please see the following work on which this is based.

Landmark indexing for scalable evaluation of label-constrained reachability queries. Lucien Valstar.  MSc thesis, Eindhoven University of Technology, 2016. https://pure.tue.nl/ws/files/46945317/855447-1.pdf
“



## Code Example

The main part of the code is written in C++ and there are scripts in bash/sh and Python (2.7).

The C++-part has a definition of a labelled graph (Graph), a set of indices to answer LCR-queries (Index/Unbounded), a number of tests (tests/Index/Unbounded) and to run the experiments or generate queries for the experiments (experiments/Index/Unbounded).

The Python-part consists of a script to generate a synthetic graph under a given model with a specified number of vertices.

For example

```python
python2.7 genGraph.py 1000 5 8 exp pa
```

generates a graph with 1000 vertices and roughly 5000 edges under the 'Preferential-Attachment' model. Any edge has one of 8 possible labels. The distribution of the edge labels is an exponential one.

## Installation

You need g++ (>=4.5) and python2.7.

The code can be built by running:

```sh
cd LCRIndexing
sh ./rebuild
```

on Mac or Linux. 

For Windows you should first remove the directories: build and .waf* (where * is an arbitary sequence). Then run:

```python
python2.7 waf configure
python2.7 waf build
```

## Contributors

Myself. Lucien Valstar: l + d + j + (last name in lower case letters) at gmail.com

My supervisors. George Fletcher: g dot h dot l dot (last name in lower case letters) at tue.nl

Yuichi Yoshida: y dot (last name in lower case letters) at gmail.com

