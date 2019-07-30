include "igraph.h"

using namespace std;


int main(){

  /*********************
   *COMMAND LINE INPUTS*
   *********************/ 

  /* Main function takes the command-line inputs    
   * as character strings separated by single spaces
   */

  if(argc != 4){
    fprintf(stderr, "usage: ./optimal [input network] [type] [input layer] > [output community membership]\n");
    fprintf(stderr, "\t[type] 0:undirected, 1:directed, 2:ordered.\n");
    fprintf(stderr, "\tif type==2 and layer==NA, the indices of nodes are used as layer.\n");
    return 1;
  }

 /**********************
  *READ GRAPH FROM FILE*
  **********************/


  /* Read graph from argv[1] using igraphs read_graph function
  */ 

  igraph_t graph;
  igraph_integer_t N=0;   //number of nodes to be read from egdelist, set to 0 if you want the whole edgelist
  igraph_bool_t directed;    //whether the graph is directed 1 or undirecte

  if(type == 0){
    directed = 0;
  }else{
    directed = 1;
  }

  //open edgelist file
  FILE *file_edgelist;
  file_edgelist=fopen(argv[1],"r");
  //read from edgelist
  igraph_read_graph_edgelist(&graph,file_edgelist,N,directed);
  //close egdelist file
  fclose(file_edgelist);
  //error handling; returns an error if graph initialization failed
  extern igraph_error_handler_t igraph_error_handler_abort;


  /*********************
   *MAXIMIZE MODULARITY*
   *********************/

  double modularity;
  igraph_vector_t membership;
  igraph_vector_init(&membership,no_of_nodes);

  igraph_community_optimal_modularity(&graph,&modularity,&membership,0);
  cout << "#: Q= "<< modularity << endl;

  for(int i=0; i<V; i++){
    cout << igraph_vector_e(&membership,i) << endl;
  }

  return 0;
}
