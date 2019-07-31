#!/bin/bash

# Prepare input data
Rscript ./example_hepth.R

# Infer communities using spectral method
../src/spectral edgelist.txt 2 layer.txt 1 > hepth.membership

# Post-processing
../src/post_processing edgelist.txt 2 layer.txt hepth.membership 1 > hepth.membership.postprocessed

# Calculate modularity of partition
echo "t2"
../src/modularity edgelist.txt 0 layer.txt hepth.membership.postprocessed
../src/modularity edgelist.txt 1 layer.txt hepth.membership.postprocessed
../src/modularity edgelist.txt 2 layer.txt hepth.membership.postprocessed

