#include "graph.h"
#include <igraph.h>
#include <fstream> 
#include "igraph_interface.h" 
#include "igraph_structural.h"  
#include "igraph_community.h"  
#include "igraph_error.h"  
//#include "igraph_glpk_support.h"
//#include "igraph_interrupt.h"  
#include "igraph_centrality.h" 
//#include "config.h"  
#include <glpk.h> //solver

using namespace std;


int set_layer(char filename[],igraph_integer_t no_of_nodes,vector<int> &layer){

  FILE *fp_r;
  char buf[100];

  int t_max;

  if(strcmp(filename, "NA") == 0){
    
    layer.resize(no_of_nodes);
    rep(i, no_of_nodes){
      layer[i] = i;
    }
    t_max = no_of_nodes;
  
  }
  else{
    if((fp_r = fopen(filename, "r")) == NULL){
      //-disp
      //fprintf(stderr, "file open error-%s\n", filename);
      exit(EXIT_FAILURE);
    }

    layer.resize(no_of_nodes);

    int i = 0;
    int tmp_t;
    t_max = 0;
    while(fgets(buf, 100, fp_r) != NULL){
      sscanf(buf, "%d", &tmp_t);
      layer[i] = tmp_t;
      i++;

      if(t_max < tmp_t){
        t_max = tmp_t;
      }
    }
    fclose(fp_r);

    t_max++;
    assert(i == no_of_nodes);
  }
  return(t_max);
}



void calculate_kappa(int no_of_nodes, int L, igraph_vector_t indegree, igraph_vector_t outdegree, vector<int> layer, vector<double> &k_in,vector<double> &k_out){

 
    vector<long double> mu, lambda;
   
    k_in.resize(no_of_nodes);
    k_out.resize(no_of_nodes);
    mu.resize(L);
    lambda.resize(L);


    rep(i,no_of_nodes){
      k_in[i] = 0.0;
      k_out[i] = 0.0;
    }
    rep(i,L){
      mu[i] = 0.0;
      lambda[i] = 0.0;
    }


    rep(i,no_of_nodes){
      k_in[i] = (double) VECTOR(indegree)[i];
      k_out[i] = (double) VECTOR(outdegree)[i];
      //cerr << k_out[i] << endl;
    }

    // ordered: null model == modified in / out-degree
 
    int t;
     
    //calculate mu, lambda and sum_k
    rep(i, no_of_nodes){
      long double value = k_in[i] - k_out[i];
      for(t=layer[i]+1; t<L; t++){
        mu[t] += value;
        lambda[t] += value;
        //cerr << value << "," << mu[t] << endl;
      }

      if(layer[i] > 0){
        lambda[layer[i]] -= k_out[i];
      }

    }

   rep(t, L){
    
      if((t != 0)&&(t != L-1)&&(lambda[t] == 0)){
        fprintf(stderr, "there is a zero flux point at t = %d\n", t);
        exit(EXIT_FAILURE);
      }

    }

    vector<long double> a(L), b(L);    

    if(L == 1){
      fprintf(stderr, "t_max should be equal to or larger than 2.\n");
      exit(EXIT_FAILURE);
    }
    else if(L > 2){
      // if t_max == 2, the network is a directed bipartite graph
      // in which all the directed links go from nodes at t=1 to nodes at t=0.
      // the null model is equivalent to that in the directed case.

      // a is f(1,l_i) and b is f(l_i,L)
      a[0] = 1.0;
      a[1] = 1.0/mu[1];
      
      for(t=2; t<L; t++){
        a[t] = a[t-1]*(lambda[t-1]/mu[t]);
      }

      b[0] = 1.0;
   
      for(t=1; t<L-1; t++){
        b[t] = b[t-1]*(mu[t]/lambda[t]);
      }
      b[L-1] = b[L-2] * mu[L-1];
   
      rep(i, no_of_nodes){
        k_in[i] = b[layer[i]]*k_in[i];
        k_out[i] = a[layer[i]]*k_out[i];
      }
    }


       

    
}







int main(int argc, char* argv[]){

  /*Plan is to use igraph's optimal community and then to extend their algorithm with other null models
   *So we first implement optimal community
   */


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


  /*For now I leaf the code for making graph object with graph.h here
  */

  /* set seed    
   * In order to set different seed at every run, we use time and process id.
   */
  srandom(time(NULL)+(int)getpid());

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
  if(type == 2){
    g.set_appeartime(argv[3]);
  }

  /* CALCULATE THE NULL MODEL SPECIFIC VECTORS:
   * link_list -> null model vectors
   * undirected: degree
   * directed: in-degree, out-degree
   * ordered: modified in-degree, modified out-degree 
   */
  long double alpha = 1.0;
  g.set_null_model(type,alpha);



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



  /*****************************
   *BASIC QUANTITIES OF NETWORK*
   *****************************/
 
  igraph_integer_t no_of_nodes=(igraph_integer_t) igraph_vcount(&graph);
  igraph_integer_t no_of_edges=(igraph_integer_t) igraph_ecount(&graph);
  
  igraph_vector_t indegree;
  igraph_vector_t outdegree;
  //initialize vectors
  IGRAPH_VECTOR_INIT_FINALLY(&indegree, no_of_nodes); //what does finally mean?
  IGRAPH_VECTOR_INIT_FINALLY(&outdegree, no_of_nodes);

  IGRAPH_CHECK(igraph_strength(&graph, &indegree, igraph_vss_all(), 
        IGRAPH_IN, IGRAPH_LOOPS, 0));
  IGRAPH_CHECK(igraph_strength(&graph, &outdegree, igraph_vss_all(), 
        IGRAPH_OUT, IGRAPH_LOOPS, 0)); 


  //for Qdag we need lambda, mu and the strength kappa_in[i] and kappa_out[i]
  //then, if layer j is larger than layer i, we have P_ij = kappa_in[i]*kappa_out[j]
  vector<int> layer;
  vector<double> k_in, k_out;
  int L;
  L = set_layer(argv[3],no_of_nodes,layer);
  
  calculate_kappa(no_of_nodes, L, indegree, outdegree, layer, k_in, k_out);



 /***********************
  *DEFINE LINEAR PROGRAM*
  ***********************/
 
  
  int no_of_variables=no_of_nodes * (no_of_nodes+1)/2;  //N*(N-1)/2
  int i, j, k, l, st;
  int idx[] = {0, 0, 0, 0 };  //for some reason we need 4 entries here, not really sure why
  double coef[] = {0.0, 1.0, 1.0, -2.0 };
 
  double modularity;
  igraph_vector_t membership;
  igraph_vector_init(&membership,no_of_nodes);

  /*DEBUG 
    for(int i=0;i < no_of_nodes; i++){
    cout << igraph_vector_e(&indegree,i) << endl;
    }

    for(int i=0;i < no_of_nodes; i++){
    cout << igraph_vector_e(&outdegree,i) << endl;
    }
  */ 

  
  glp_prob *ip;  //pointer to lp
  glp_iocp parm;

  glp_term_out(GLP_ON);
  ip = glp_create_prob();
  IGRAPH_FINALLY(glp_delete_prob, ip);
  glp_set_obj_dir(ip, GLP_MAX); //problem is a maximization problem
  st=glp_add_cols(ip, no_of_variables); //not sure what this is doing, seems to add columns, but its value is always 1


  /* variables are binary */
  for (i=0; i<no_of_variables; i++) {
    glp_set_col_kind(ip, (1+i), GLP_BV);
  }

  #define IDX(a,b) ((b)*((b)+1)/2+(a)+1) //this maps the pair (a,b) to b(b+1)/2 + a
 
  /*DEBUG*/
  //for(int i = 0; i < no_of_nodes; i++){
  //  for(int j = 0; j < no_of_nodes; j++){
  //    cerr << "i = " << i << "j = " << j << "IDX(i,j) = " << IDX(i,j) << endl;
  //  }  
  //}
  //cerr << "no_of_variables = " << no_of_variables << " IDX(N,N) = " << IDX(no_of_nodes-1,no_of_nodes-1) << endl;


  /*****************
   *SET CONSTRAINTS*
   *****************/

  /* reflexivity */
  for(i=0; i<no_of_nodes; i++) {
    glp_set_col_bnds(ip, IDX(i,i), GLP_FX, 1.0, 1.0);
  }

  /* transitivity */
  for (i=0; i < no_of_nodes; i++) {
    for (j=i+1; j < no_of_nodes; j++) {
      for (k=j+1; k < no_of_nodes; k++) {
                 
        int newrow=glp_add_rows(ip, 3);
       
        glp_set_row_bnds(ip, newrow, GLP_UP, 0.0, 1.0);
        idx[1] = IDX(i,j); 
        idx[2] = IDX(j,k); 
        idx[3] = IDX(i,k);
        glp_set_mat_row(ip, newrow, 3, idx, coef);

        glp_set_row_bnds(ip, newrow+1, GLP_UP, 0.0, 1.0);
        idx[1] = IDX(i,j); 
        idx[2] = IDX(i,k); 
        idx[3] = IDX(j,k);
        glp_set_mat_row(ip, newrow+1, 3, idx, coef);

        glp_set_row_bnds(ip, newrow+2, GLP_UP, 0.0, 1.0);
        idx[1] = IDX(i,k); 
        idx[2] = IDX(j,k); 
        idx[3] = IDX(i,j);
        glp_set_mat_row(ip, newrow+2, 3, idx, coef);
 
      }
    }
  }

  /************************
   *SET OBJECTIVE FUNCTION*
   ************************/
  
   double c;

   if(type == 0){

   /* first part: -strength(i)*strength(j)/total_weight for every node pair */
    for (i=0; i<no_of_nodes; i++) {
      for (j=i+1; j<no_of_nodes; j++) {
        
        c = (double) -1.0*((VECTOR(indegree)[i]*VECTOR(outdegree)[j]) + (VECTOR(outdegree)[i]*VECTOR(indegree)[j])) / (2.0*no_of_edges);
        glp_set_obj_coef(ip, IDX(i,j), c);
            
      }
      /* special case for (i,i) */
      c = (double) -1.0*(VECTOR(indegree)[i]*VECTOR(outdegree)[i]) / (2.0*no_of_edges);
      glp_set_obj_coef(ip, IDX(i,i), c);
    }

   }else if(type == 1){
    
    /* first part: -strength(i)*strength(j)/total_weight for every node pair */
    for (i=0; i<no_of_nodes; i++) {
      for (j=i+1; j<no_of_nodes; j++) {
        
        c = (double) -1.0*((VECTOR(indegree)[i]*VECTOR(outdegree)[j]) + (VECTOR(outdegree)[i]*VECTOR(indegree)[j])) / (no_of_edges);
        glp_set_obj_coef(ip, IDX(i,j), c);
            
      }
      /* special case for (i,i) */
      c = (double) -1.0*(VECTOR(indegree)[i]*VECTOR(outdegree)[i]) / (no_of_edges);
      glp_set_obj_coef(ip, IDX(i,i), c);
    }
   }else if(type == 2){
    for (i=0; i<no_of_nodes; i++) {
      for (j=i+1; j<no_of_nodes; j++) {
        int i_tmp = i, j_tmp = j;
        if(layer[i] > layer[j]){
          j_tmp = i;
          i_tmp = j;
        }
        if(layer[i_tmp] < layer[j_tmp]){
          c = (double) -1.0 * k_in[i_tmp] * k_out[j_tmp];
          
          glp_set_obj_coef(ip, IDX(i,j), c);
         }
      }
    }
     //no need for IDX(i,i), because it is always zero
   }
 
   

    /* second part: add the adjacency matrix to the coefficient matrix */
    for (k=0; k<no_of_edges; k++) {

      i = IGRAPH_FROM(&graph, k);
      j = IGRAPH_TO(&graph, k);
      if(i > j){
        l = j;
        j = i;
        i = l;
      }
      if(type == 0){
        c = 2.0;
      }else{
        c = 1.0;
      }
      //cerr << c << "," << glp_get_obj_coef(ip, IDX(i,j)) << endl;
      glp_set_obj_coef(ip, IDX(i,j), c + glp_get_obj_coef(ip, IDX(i,j)));
      //cerr << glp_get_obj_coef(ip,IDX(i,j)) << endl;

    }
  

  /*DEBUG
  for(i = 0; i < no_of_nodes; i++){
    for(j = i+1;  j< no_of_nodes; j++){
      if(IDX(i,j) > no_of_variables) cerr << "IDX(i,j) = " << IDX(i,j) << endl;
        cerr << i << " " << j << " " <<  glp_get_obj_coef(ip,IDX(i,j)) << endl;
    }
  }
  */
  


  /**********
   *SOLVE LP* 
   **********/

  
  glp_init_iocp(&parm);
  parm.br_tech = GLP_BR_DTH;
  parm.bt_tech = GLP_BT_BLB;
  parm.presolve = GLP_ON;
  parm.binarize = GLP_ON;
  //parm.cb_func = igraph_i_glpk_interruption_hook;
  //IGRAPH_GLPK_CHECK(glp_intopt(ip, &parm), "Modularity optimization failed");
  if(glp_intopt(ip, &parm) == 1) cerr << "optimization failed" << endl; 
    

  /****************
   *RETURN RESULTS* 
   ****************/

  if(type == 0){
    modularity = glp_mip_obj_val(ip) / (2*no_of_edges);  
  }else{
    modularity = glp_mip_obj_val(ip) / no_of_edges;
  }

  long int comm=0;   /* id of the last community that was found */
  IGRAPH_CHECK(igraph_vector_resize(&membership, no_of_nodes));
  for (i=0; i<no_of_nodes; i++) {

    for (j=0; j<i; j++) {
      int val=(int) glp_mip_col_val(ip, IDX(j,i));
      if (val==1) {
        VECTOR(membership)[i]=VECTOR(membership)[j];
        break;
      }
    }
    if (j==i) {   /* new community */
      VECTOR(membership)[i]=comm++;
    }

  }

  cout << "#: Q= "<< modularity <<  endl;

  for(int i=0; i<no_of_nodes; i++){
    cout << igraph_vector_e(&membership,i) << endl;
  }

  /*DEBUG
  igraph_community_optimal_modularity(&graph,&modularity,&membership,0);
  cout << "#: Q= "<< modularity << endl;

  for(int i=0; i<V; i++){
    cout << igraph_vector_e(&membership,i) << endl;
  }
  */

  #undef IDX

  igraph_vector_destroy(&indegree);
  igraph_vector_destroy(&outdegree);
  glp_delete_prob(ip);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;


}




