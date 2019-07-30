#include "graph.h"

using namespace std;

int main(int argc, char* argv[]){

  /* Main function takes the command-line inputs    
   * as character strings separated by single spaces
   */
 
  if(argc != 5){
    fprintf(stderr, "usage: ./modularity [input network] [type] [input appear_time] [output community membership]\n");
    fprintf(stderr, "\t[type] 0:undirected, 1:directed, 2:ordered.\n");
    fprintf(stderr, "\tif type==2 and appear_time==NA, the indices of nodes are used as appear_time.\n");
    return 1;
  }

  /* CREATE A GRAPH OBJECT FROM THE EDGELIST WE JUST READ
   * read a link list from file
   * index: 0 ~ N-1
   * link_list -> sparse_matrix
   */
  Graph g = Graph::readGraph(argv[1]);

  int type = atoi(argv[2]); // convert char to int

  /*ordered: 
   *read appearance time of each node
   */
  if(type != 0){
    g.set_appeartime(argv[3]);
  }

  /* CALCULATE THE NULL MODEL SPECIFIC VECTORS:
   * link_list -> null model vectors
   * undirected: degree
   * directed: in-degree, out-degree
   * ordered: modified in-degree, modified out-degree 
   */
  int alpha = 1; //resolution parameter, not needed here
  g.set_null_model(type,alpha);
  
  
  //read community
  int max_comm_index = g.set_community_ids(argv[4]);
  
  if(0){
  //leo: create a pointer to file
  FILE *fp_r;
  //leo: buffer
  char buf[100];
  //leo: return error if path does not exist
  if((fp_r = fopen(argv[4], "r")) == NULL){
    fprintf(stderr, "file open error.-%s\n", argv[4]);
    exit(EXIT_FAILURE);
  }
  
  int u;
  vector<int> community_index;
  fgets(buf, 100, fp_r);
  
  while(fgets(buf, 100, fp_r) != NULL) {
    
    sscanf(buf, "%d", &u);
    community_index.push_back(u);
  }
  fclose(fp_r);
  }
  
  long double Q;
  Q = g.modularity(type,max_comm_index);
  
  //fprintf(stderr,"comm = %d \n",community_index[1]);
  fprintf(stderr,"Number of communities = %d \n", max_comm_index + 1);
  fprintf(stderr,"Q = %Lf \n", Q);
 
  
}