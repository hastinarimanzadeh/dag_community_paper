#include "graph.h"
//for random shuffle
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]){


  /* Main function takes the command-line inputs    
   * as character strings separated by single spaces
   */

  if(argc != 6){
    fprintf(stderr, "usage: ./check_singleton [input network] [type] [input appear_time] [input community membership file] [type of post processing] > [output community membership file]\n");
    fprintf(stderr, "\t[type] 0:undirected, 1:directed, 2:ordered.\n");
    fprintf(stderr, "\tif type==2 and appear_time==NA, the indices of nodes are used as appear_time.\n");
    return 1;
  }


  /****************************
   *CREATE GRAPH FROM EDGELIST*
   ****************************/ 

  Graph g = Graph::readGraph(argv[1]);



  /************************
   *SET TYPE OF MODULARITY*
   ************************/

  int type = atoi(argv[2]); // convert char to int

  /*ordered: 
   *read appearance time of each node
   */
  if(type != 0){
    g.set_appeartime(argv[3]);
  }



  /*****************
   *SET NULL MODELS*
   *****************/ 

  long double alpha = 1.0;
  g.set_null_model(type,alpha);



  /********************
   *READ COMMUNITY IDs*
   ********************/ 

  int max_comm_index = g.set_community_ids(argv[4]);

  long double Q = g.modularity(type,max_comm_index);

  //fprintf(stderr,"Q= %Lf\n",Q);

  /************
   *SINGLETONS*
   ************/ 

  long double increase_Q = 0.0;
  long double increase_Q_old = 0.0;

  /*PLAN*
   *make vector for size of communities
   *go through all communities and try to merge them with an adjacent community
   *i.e. we also need to make a vector of ajacent communities
   */

  vector< vector<int> > community_membership;
  vector<int> community_size(max_comm_index+1,0);
  community_membership.resize(max_comm_index+1);
  int max_comm_size;

  g.community_assignment(community_membership,community_size,max_comm_size);

  int current_size = 1;

  //tune code by constructing vector with size dist
  //if this code works, try to combine both post_processing steps

  int type_of_pp = atoi(argv[5]);

  if(type_of_pp == 1){
    while(current_size < max_comm_size + 1){
      rep(k,max_comm_index){
        if(community_size[k] == current_size){
          //fprintf(stderr,"debug\n");

          /***************************
           * get adjacent communities*
           ***************************/
          vector<int> adjacent_comms;
          g.get_adjacent_communities(k,community_membership[k],community_size,adjacent_comms,max_comm_index + 1);
          //int foo = (int) adjacent_comms.size();
          //fprintf(stderr,"number ad comm = %d \n",foo);
          
          /*********************************************
           * attempt to merge with adjacent communities*
           *********************************************/
          long double Q_temp = 0.0, Q_max = 0.0;
          int merge_comm_index = k;
          rep(i,adjacent_comms.size()){
            int is_final = 0;
            Q_temp = g.merge_communities(k,adjacent_comms[i],type,community_size,community_membership,is_final);
            if(Q_temp < Q_max){
              merge_comm_index = adjacent_comms[i];
              Q_max = Q_temp;
            }
          }
          
          /*********************************************************************
           * merge with community that gives the largest increase in modularity*
           *********************************************************************/ 
          if(merge_comm_index != k){
            int is_final = 1;
            Q_max = g.merge_communities(k,merge_comm_index,type,community_size,community_membership,is_final);
            increase_Q -= Q_max; 
            //fprintf(stderr,"Q_max = %Lf \n", Q_max);
            //fprintf(stderr,"toinclude = %d \n",merge_comm_index);
          }
        }
      }
      current_size++;
    }
  }

  //rep(i,max_comm_index+1){
  //  fprintf(stderr,"comm_size = %d\n",community_size[i]);
  //}

  if(type_of_pp == 2){

    max_comm_index++;
    //fprintf(stderr,"max_comm_index = %d\n",max_comm_index);
    community_membership.resize(max_comm_index+1);
    community_membership[max_comm_index].push_back(0);
    community_size.push_back(1);
    long double increase_Q_it = increase_Q;


    /*****************
     * random shuffle*
     *****************/
    std::random_device rd;
    std::mt19937 gen(rd());

    vector<int> perm(g.get_V(),0);
    rep(i,g.get_V()){
      perm[i] = i;
    }
    //shuffle vector perm
    shuffle(perm.begin(), perm.end(), gen);
    /*DEBUG*/
    //rep(i,perm.size()){
    //  fprintf(stderr,"perm[%d]=%d \n",i,perm[i]);
    //}

    //do{ 
      increase_Q_it = increase_Q;
      //fprintf(stderr,"check\n");
      rep(j,g.get_V()){
        //fprintf(stderr,"i = %d\n",i);
   

        /*************************************
         *set i equal to the permutated value*
         *************************************/ 
        int i = perm[j];
        //fprintf(stderr,"i=%d\n",i);


        /***************
         * isolate node*
         ***************/
        //set former_community = max_comm_index if you want to isolate node
        int is_final = 0;
        int former_community = max_comm_index;

        long double Q_tmp = g.isolate_node(i,type,max_comm_index,former_community,community_size,community_membership[max_comm_index],is_final);
        long double increase_Q_old = increase_Q;
        increase_Q += Q_tmp;

        //fprintf(stderr,"comm_member = %d \n",community_membership[max_comm_index][0]);
        /***********************************
         *get list of adjacent communities*
         ***********************************/        
        vector<int> adjacent_comms;
        g.get_adjacent_communities(max_comm_index,community_membership[max_comm_index],community_size,adjacent_comms,max_comm_index + 1);
        int foo = (int) adjacent_comms.size();
        //fprintf(stderr,"number ad comm = %d \n",foo);
        
        /*********************************************
         * attempt to merge with adjacent communities*
         *********************************************/ 
        long double Q_temp = 0.0, Q_max = 0.0;
        int merge_comm_index = max_comm_index;
        rep(k,adjacent_comms.size()){
          is_final = 0;
          //fprintf(stderr,"adjacent_comms[0] = %d\n",adjacent_comms[0]);
          //fprintf(stderr,"size = %d\n",community_size[adjacent_comms[0]]);
          Q_temp = g.isolate_node(i,type,max_comm_index,adjacent_comms[k],community_size,community_membership[max_comm_index],is_final);
          if(Q_temp < Q_max){
            merge_comm_index = adjacent_comms[k];
            Q_max = Q_temp;
          }
        }


        /********************************************************************
         *merge with community that gives the largest increase in modularity*
         ********************************************************************/ 
        assert(merge_comm_index != max_comm_index);
        is_final = 1;
        Q_max = g.isolate_node(i,type,max_comm_index,merge_comm_index,community_size,community_membership[max_comm_index],is_final);
        increase_Q -= Q_max; 
        community_size[max_comm_index] = 1;
        //if(increase_Q - increase_Q_old < 0) fprintf(stderr,"Q_diff=%Lf \n",increase_Q-increase_Q_old);
        assert(increase_Q - increase_Q_old >= -1.0*EPS); //we want to set 0 here but somtimes it is equal to -0
        //fprintf(stderr,"Q_max = %Lf \n", Q_max);
        //fprintf(stderr,"toinclude = %d \n",max_comm_index);
        //fprintf(stderr,"increase_Q = %Lf, mod = %Lf\n",Q + increase_Q/(4.0*g.get_original_E()),g.modularity(type,max_comm_index));
        assert(increase_Q >= 0);
      }
    //}while(increase_Q > increase_Q_it);
  }
  
  
  /*DEBUG*/
  //rep(i,max_comm_index+1){
  //  fprintf(stderr,"comm_size = %d\n",community_size[i]);
  //}
  //fprintf(stderr,"max_comm_index = %d\n",max_comm_index);
  
  /***************************
   * adjust community indices*
   ***************************/
  g.adjust_community_indices(max_comm_index,community_size);


  /* output results
   * print modularity value
   * print partition membership
   */
  Q += increase_Q/(4.0*g.get_original_E());

  printf("#Q: %Lf\n",Q); 

  /* WRITE OUTPUT TO FILE
  */ 
  g.show_community_id();

  //fprintf(stderr,"Q = %Lf\n",Q);    

  return 0;


}
