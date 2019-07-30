#!/bin/sh

cp hepth_lcc_edgelist_reassign.txt edgelist.txt
cp hepth_lcc_date_reassign.txt layer.txt

echo "t0"
./modularity edgelist.txt 0 layer.txt hepth_lcc_community_t0.txt
./modularity edgelist.txt 1 layer.txt hepth_lcc_community_t0.txt
./modularity edgelist.txt 2 layer.txt hepth_lcc_community_t0.txt

echo "t1"
./modularity edgelist.txt 0 layer.txt hepth_lcc_community_t1.txt
./modularity edgelist.txt 1 layer.txt hepth_lcc_community_t1.txt
./modularity edgelist.txt 2 layer.txt hepth_lcc_community_t1.txt

echo "t2"
./modularity edgelist.txt 0 layer.txt hepth_lcc_community_t2.txt
./modularity edgelist.txt 1 layer.txt hepth_lcc_community_t2.txt
./modularity edgelist.txt 2 layer.txt hepth_lcc_community_t2.txt

echo "louvain"
./modularity edgelist.txt 0 layer.txt hepth_lcc_community_louvain.txt
./modularity edgelist.txt 1 layer.txt hepth_lcc_community_louvain.txt
./modularity edgelist.txt 2 layer.txt hepth_lcc_community_louvain.txt

