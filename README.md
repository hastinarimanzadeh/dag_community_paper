# Community detection in directed acyclic graphs

authors: Leo Speidel\*, Taro Takaguchi\*, Naoki Masuda (\*contributed equally)<br/>
paper: [European Physical Journal B 88, 203 (2015)](http://link.springer.com/article/10.1140%2Fepjb%2Fe2015-60226-y)

## Installation

1. Install igraph and glpk. On a Mac, you can use<br/> `brew install igraph` and `brew install glpk`
2. Open src/Makefile and adjust path to igraph
3. Change to directory src and type `make` 

Alternatively, if igraph and glpk are unavailable, remove -ligraph and -lglpk.<br/>
These are only required for optimal and optimal\_community, so you should be able to compile everything else.

## Usage

See ./src/readme.txt for details.

There are five binaries,

 - spectral: implements the spectral method for community detection
 - post\_processing: a heuristic post\_processing algorithm for improving modularity by reassigning communities of individual nodes
 - modularity: calculates the modularity of a community assignment
 - optimal: optimal assignment of nodes to communities (slow)
 - optimal\_community: same as optimal but different implementation (slow)

## Pre-processing

There is an R script in Rfunctions/pre\_proc\_graph.R which should prepare the required input files.
Output is edgelist.txt and layer.txt.

