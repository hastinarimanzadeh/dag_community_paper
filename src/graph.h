#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cassert>
#include <unistd.h>

#define all(c) (c).begin(), (c).end()
#define iter(c) __typeof((c).begin())
#define rep(i, n) for (int i = 0; i < (int)n; i++)
#define tr(c, i) for (iter(c) i = (c).begin(); i != (c).end(); ++i)

#define EPS 1e-05 // threshold for termination of power method
#define EPS_Q 0.0 // threshold for termination of modularity maximization
#define ITER_MAX 1000 // the maximun number of iteration in a sigle execution of power method

class Graph;
typedef void(Graph::*fp)(std::vector<long double>&, std::vector<long double>&);

class Graph{

 private:
  int V, E;
  int number_community;
  int t_max;
  long double original_E;
  
  std::vector<int> nodes;
  std::vector< std::pair<int, int> > edges, new_edges;
  std::vector< std::vector<int> > adj;
  std::vector<int> col, row, ids, appear_time, community_id;
  std::vector<long double> k_in, k_out;

 public:
  static Graph readGraph(char filename[]);
  void show_graph(int type);
  
  int set_community_ids(char filename[]);
  int set_comm_id(int i, int base);
  long double modularity(int type, int max_comm_index);
  
  bool check_connectivity(int comm,std::vector< std::pair<int,int> > &parent,int &parent_lcc);
  
  void set_appeartime(char filename[]); 
  void set_null_model(int type, long double alpha);

  void product_undirected(std::vector<long double> &s, std::vector<long double> &new_s);
  void product_directed(std::vector<long double> &s, std::vector<long double> &new_s);
  void product_ordered(std::vector<long double> &s, std::vector<long double> &new_s);

  std::vector<long double> power_method(int type);
  
  Graph extract(int comm);
  long double recalculate_tunedQ(int type,std::vector<long double> &s,long double &sum_skout,long double &sum_skin,int index);
  void fine_tuning(int type,std::vector<long double> &s);
  
  void divide_community(int comm, std::vector<long double> &s);
  Graph divide_components(int comm, int type, long double &Q, Graph tmp_g);
  long double diff_modularity(int type, std::vector<long double> &s);

  long double innprod_k_out(std::vector<long double> &k_out, std::vector<long double> &x, int t);
  long double innprod_k_in(std::vector<long double> &k_in, std::vector<long double> &x, int t);
  std::vector<long double> multiply_k_out(std::vector<long double> &k_out, long double C, int t);
  std::vector<long double> multiply_k_in(std::vector<long double> &k_in, long double C, int t); 

  long double innprod(std::vector<long double> &x, std::vector<long double> &y);

  void community_assignment(std::vector< std::vector<int> > & community_membership,std::vector<int> & community_size, int &max_comm_size);
  void get_adjacent_communities(int k,std::vector<int> &community_member,std::vector<int> &community_size,std::vector<int> &adjacent_comms,int max_comm_index);
  long double isolate_node(int node_index,int type,int max_comm_index,int former_community,std::vector<int> &community_size,std::vector<int> &community_membership,int is_final);
  long double merge_communities(int comm_tobeincluded, int comm_toinclude, int type,std::vector<int> &community_size, std::vector< std::vector<int> > &community_membership,int is_final);
  void adjust_community_indices(int max_comm_index, std::vector<int> &community_size);
  
  int get_V();
  int get_original_E();
  int get_number_community();
  void show_community_id();
};

