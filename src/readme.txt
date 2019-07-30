
This file lists/summarizes the contents and usage of the code in this folder

graph.h: header file for class graph. It includes all the essential functions

graph.cpp: The actual code for the class that is defined in graph.h

spectral.cpp: spectral method for Q^und, Q^dir and Q^dag
    usage: ./spectral [edgelist] [type] [layer] [resolution parameter] > [output file]
    where [type] = 0 for und, 1 for dir, and 2 for dag
    fine tuning is implemented, see the corresponding function in graph.cpp

post_processing.cpp: post processing algorithm for partitions
    usage: ./post_processing [edgelist] [type] [layer] [type of pp] > [output file]
    where [type of pp] = 2 for finding best community for node
          [type of pp] = 1 for merging communities

optimal.cpp: code for finding optimal partition using integer linear progr.
    usage: ./optimal [edgelist] [type] [layer] > [output file]

optimal_community.cpp: code for finding optimal partition using integer linear progr. as implemented in the igraph package for C++
    usage: same as optimal
      
Makefile: for compiling the code.

Some remarks:
In the paper we are describing the method in its most general form.
In particular, we allow for layers in which no link penetrates the layer.
This code cannot handle graphs with such layers.

Pre processing networks (DAGs) for this code:
There is a R function in the respective folder (pre_proc_graph.R), which checks whether the graph 
is a DAG and then creates layers for nodes using a leaf removal algorithm. It
also sorts node indices such that they satisfy the condition of layers.
they are enumerated also to respect layer  

Plotting community assignment of DAGs:
There is a R function in the respective folder (plot_ordered_graph), which
uses an edgelist, community assignment as input to create a figure. A special
layout as well as fruchtermann-reinfold is implemented - slight changes to
code are needed to switch between the two (just comment-out a line)


