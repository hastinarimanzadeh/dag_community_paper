#include "graph.h"

using namespace std;

int main(int argc, char* argv[]){

  /* Main function takes the command-line inputs    
   * as character strings separated by single spaces
   */
 
  if(argc != 5){
    fprintf(stderr, "usage: ./spectral [input network] [type] [input appear_time] [resolution limit parameter] > [output community membership]\n");
    fprintf(stderr, "\t[type] 0:undirected, 1:directed, 2:ordered.\n");
    fprintf(stderr, "\tif type==2 and appear_time==NA, the indices of nodes are used as appear_time.\n");
    return 1;
  }

  
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
  long double alpha = sqrt(atof(argv[4]));
  g.set_null_model(type,alpha);
  
  /******debug******/
  /*
  vector<double> v(g.get_V(), 1.0);
	fprintf(stderr, "Q before decomp.:%Lf\n", g.diff_modularity(type, v));
	return 1;
  */
  /******debug******/   

  g.show_graph(type); //should be muted after debugging.
  
  
  
  
  
  
  /**********************
   *FIRST BIPARTITIONING*
   **********************/
  
  //fprintf(stderr,"\n");
  
  /* POWER METHOD
   * power-method for finding the dominant eigenvector
   * g.get_V() returns the number of nodes of g, (note that V is private)  
   */
  vector<long double> s(g.get_V());
  long double Q; // modularity value
  s = g.power_method(type);
  
  
  /******debug******/
  /*
  fprintf(stderr, "s:");
	fp_r = fopen("index.txt", "r");
	target\index[i]=rep(i, g.get_V()){
    fprintf(stderr, "%Lf ", s[i]);
  }
  fprintf(stderr, "\n");
  */
  
  //fprintf(stderr, "g.diff: %Lf E: %d \n", g.diff_modularity(type, s), g.get_original_E());
  /******debug******/
  
  
  /* FINE TUNING
   * Fine tuning is implemented to tune the partition we obtain from the power method.
   * Basically we implemented the tuning technique of Newman, PNAS, 2006.
   */
  
  g.fine_tuning(type,s);
  
  
  /* CALCULATE Q 
   * (modularity of bipartition)
   */
  //Q = g.diff_modularity(type, s) / (4.0*g.get_original_E());
  Q = g.diff_modularity(type, s);
  

  //long double foo = Q; 
  //foo /= (4.0*g.get_original_E());
  
  //-disp
  //fprintf(stderr, "\nQ at 1st partition:%.10Lf, %.10Lf \n", Q/(4.0*g.get_original_E()),foo);
  
  /* BIPARTITIONING
   * If Q is positive, we accept the bipartition
   */
  if(Q > 0){
    g.divide_community(0, s);
  }
  else{
    //-disp
    //fprintf(stderr, "no bipartition is possible.\n");
    
    /*******************************************
     *Print result (all node in same community)*
     *******************************************/
    
    rep(i,g.get_V()){
      s[i] = 1;
    }
    g.divide_community(0,s);
    
    Q /= (4.0*g.get_original_E());
    printf("# Q:%Lf\n", Q);
    g.show_community_id();
    return 0;
  }
  
  
  /**************************
   *REAPEATED BIPARTITIONING*
   **************************/
  
  /* Set if(0) if you only want one bipartition
   */
  if(0==0){
 
  int comm = 0; // index of community that is about to be bipartitioned
  Graph tmp_g;  // subgraph of the community that is about to be bipartitioned
  long double diff_Q;  // change in modularity that a bipartition of community comm would cause

  while((comm < g.get_number_community())&&(g.get_number_community() < g.get_V())){
    //-disp
  //fprintf(stderr, "\ncomm:%d\n", comm);
    
    do{
      /* EXTRACT SUBGRAPH FOR COMMUNITY
       * g.extract extracts the subgraph for community comm from the original graph g
       */
      tmp_g = g.extract(comm);
      
      /* DIVIDE INTO CONNECTED COMPONENTS
       * The community may be disconnected, and divide_components divides tmp_g into its connected components. O(n log n), n = tmp_g.get_V()
       * Return value is a graph object corresponding to the giant component of tmp_g
       * NOTE: If communities are disconnected at this point, they will have to be divided into connected components anyway (in order to maximize Q).
       *       This step ensures that output consists of connected communities AND speeds up the algorithm, because it is faster that partitioning
       *       the graph via spectral partitioning.
       */
      tmp_g = g.divide_components(comm,type,Q,tmp_g);      
      
      /******debug******/
      //tmp_g.show_graph(type);
      /******debug******/
      
      /* s will be the dominant eigenvector
       */
      s.resize(tmp_g.get_V(), 0.0);
      
      
      /* POWER METHOD 
       * There is no need to apply the spectral method for communities with one element
       */
      if(tmp_g.get_V()==1){
        s[0] = 1.0;
      } else{
        s = tmp_g.power_method(type);
      }
      
      /******debug******/
      //rep(i,s.size()){
      //  s[i]=1;
      //}
	    /*
      fprintf(stderr, "s:");
      rep(i, tmp_g.get_V()){
        fprintf(stderr, "%Lf ", s[i]);
      }
      fprintf(stderr, "\n");
      */
      /******debug******/
 
      /* FINE TUNING
       */
      tmp_g.fine_tuning(type,s);
      /* DIVIDE INTO CONNECTED COMPONENTS
       */
      //tmp_g = g.divide_components(comm,type,Q,tmp_g);  
       
      
      /* CALCULATE CHANGE IN Q
       */
      //diff_Q = tmp_g.diff_modularity(type, s) / (4.0*g.get_original_E());
      diff_Q = tmp_g.diff_modularity(type, s);
      //fprintf(stderr,"diff_Q= %Lf \n",diff_Q);
      
      /******debug******/
      ///fprintf(stderr,"%Lf\n",diff_Q);
      //?? EPS is fine but how is EPS_Q also already declared..
      //Taro: EPS is the threshold for power method (which we actually don't use in the end)
      //Taro: and EPS_Q is the threshold for modularity maximization.
      /******debug******/
     
      /* BIPARTITIONING
       * If change in Q is positive, we accept the bipartitioning
       * I.e. EPS_Q = 0
       */
      if(diff_Q > EPS_Q){	
        g.divide_community(comm, s);
        Q += diff_Q;
        //-disp
        //fprintf(stderr,"diff_Q=%.10g\n",diff_Q);
        //-disp
        //fprintf(stderr, "\nQ at %dth partition:%.10g\n\n", g.get_number_community()-1, Q);
      }
      else{
        //-disp
        //fprintf(stderr, "no bipartition is possible.\n\n");
      }
    }while(diff_Q > EPS_Q);

    comm++;
  }
  }  
  
  Q /= (4.0*g.get_original_E());
 
  /*Taro: I want to separately get the results and error messages.
   *Taro: The results go to file via printf and redirection ">" in command line 
   *Taro: Error messages go to terminal via fprintf and stderr.
   */
   
   
  /* output results
   * print modularity value
   * print partition membership
   */
  printf("# Q:%Lf\n", Q);
  
  /* WRITE OUTPUT TO FILE
   */ 
  g.show_community_id();
  
 
  return 0;
}
