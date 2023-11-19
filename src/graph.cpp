#include "graph.h"

using namespace std;

Graph
Graph::readGraph(char filename[]) {

  //leo: create a pointer to file
  FILE *fp_r;
  //leo: buffer
  char buf[100];
  //leo: return error if path does not exist
  if((fp_r = fopen(filename, "r")) == NULL){
    fprintf(stderr, "file open error.-%s\n", filename);
    exit(EXIT_FAILURE);
  }

  //graph object
  Graph g;

  g.V = 0;
  g.E = 0;

  int u, v;

  //leo: read graph from file
  while(fgets(buf, 100, fp_r) != NULL) {

    sscanf(buf, "%d %d", &u, &v);

    //number of edges
    g.E++;

    //add node to V if node is new
    if(u > g.V){
      g.V = u;
    }
    if(v > g.V){
      g.V = v;
    }
    //add edge connecting u to v to graph
    g.edges.push_back(make_pair(u, v));
  }
  fclose(fp_r);
  //??why do we increase the number of nodes at this point and not in while..
  //Taro: Because I assume that the indices of nodes begin with 0.
  //Taro: Therefore the last node has index V-1, whereas the number of nodes is equal to V. 
  g.V++;
  //??what is the difference of E to original_E and why is it a long double? (weights?..)
  //Taro: original_E is the number of edges in the original (whole) network
  //Taro: g.E is the number of edges in the (extracted) subgraph.
  g.original_E = (long double)g.E;


  // adjacency list (undirectionized)
  g.adj.resize(g.V);

  //??what is tr?
  //Taro: see the preamble of graph.h.
  tr(g.edges, it){
    u = it->first, v = it->second;
    if (u != v) g.adj[u].push_back(v);             
    //?? no self loops?
    //Taro: remove self-loops if any in the data
    if (u != v) g.adj[v].push_back(u);
  }
  rep(i, g.V){
    sort(all(g.adj[i]));
  }


  // sparse matrix representation (undirectionized)
  //leo:　not checked in detail yet
  g.col.resize(2*g.E);
  g.row.resize(g.V+1);
  int count=0;
  g.row[0] = 0;
  rep(i, g.V){
    rep(j, (int)g.adj[i].size()){
      g.col[count] = g.adj[i][j];
      count++;
    }
    g.row[i+1] = count;
  }


  // set node ids
  g.ids.resize(g.V);
  rep(i,g.V){
    g.ids[i] = i;
  }

  // set initial community id
  g.community_id.resize(g.V, 0);
  g.number_community = 1;

  return g;
}


////////////////////////
void
Graph::show_graph(int type){
  //-disp
  //fprintf(stderr, "V:%d\n", V);
  //-disp
  //fprintf(stderr, "E:%d original_E:%Lf\n", E, original_E);
  /*
     fprintf(stderr, "Adjacency list:\n");
     rep(i, V){
     fprintf(stderr, "[%d] ", ids[i]);
     rep(j, adj[i].size()){
     fprintf(stderr, "%d ", ids[adj[i][j]]);
     }
     fprintf(stderr, "\n");
     }
     */
  /*
     rep(i, V){
     int count = 0;
     rep(j, V){
     if(count == adj[i].size()){
     fprintf(stderr, "0, ");
     }
     else{
     if(adj[i][count] == j){
     fprintf(stderr, "1, ");
     count++;
     }
     else{
     fprintf(stderr, "0, ");
     }
     }
     }
     }
     fprintf(stderr, "\n");
     */
  /*
     fprintf(stderr, "Sparse matrix representation:\n");
     fprintf(stderr, "col:");
     rep(j, col.size()){
     fprintf(stderr, "%d ", col[j]);
     }
     fprintf(stderr, "\n");
     fprintf(stderr, "row:");
     rep(j, row.size()){
     fprintf(stderr, "%d ", row[j]);
     }
     fprintf(stderr, "\n");
     */
  /*
     if(k_in.size() == V){
     fprintf(stderr, "degree sequence:\n");
     if(type == 0){
     rep(i, V){
  //fprintf(stderr, "[%d] %d\n", i, (int)k_in[i]);
  fprintf(stderr, "%d, ", (int)k_in[i]);
  }
  fprintf(stderr, "\n");
  }
  else{
  rep(i, V){
  fprintf(stderr, "[%d] %Lf %Lf\n", i, k_in[i], k_out[i]);
  }
  }
  }
  */
  /*
     if(appear_time.size() == V){
     fprintf(stderr, "Appear time:\n");
     rep(i, V){
     fprintf(stderr, "[%d] %d\n", ids[i], appear_time[i]);
     }
     }
     */
  /*
     rep(i, V){
     rep(j, adj[i].size()){
     printf("%d %d\n", i, adj[i][j]);
     }
     }
     printf("\n");
     */
}


////////////////////////
int
Graph::set_community_ids(char filename[]){


  FILE *fp_r;
  //leo: buffer
  char buf[100];
  //leo: return error if path does not exist
  if((fp_r = fopen(filename, "r")) == NULL){
    fprintf(stderr, "file open error.-%s\n", filename);
    exit(EXIT_FAILURE);
  }


  int comm_index;
  int max_comm_index = 0;
  int count = 0;

  fgets(buf, 100, fp_r); 

  //read community indices from file
  while(fgets(buf, 100, fp_r) != NULL) {

    sscanf(buf, "%d", &comm_index);
    community_id[count] = comm_index;
    //fprintf(stderr, "count= %d, comm_index=%d \n", count, comm_index);
    count++;  
    if(comm_index > max_comm_index) max_comm_index = comm_index;
  }
  fclose(fp_r);

  /******DEBUG*******/
  /*
     for(int i=0; i < 10; i++){ 
     fprintf(stderr,"%d,%d\n",i,community_id[i]);
     }
     */
  /******DEBUG*******/

  //fprintf(stderr,"V= %d, count = %d \n", V,count); 
  assert(count == V);


  return(max_comm_index);

}


///////////////////////
int
Graph::set_comm_id(int i,int base){

  //base can be thought of as max_comm_index+1

  int index = 0;
  int remain, num;
  vector<int> digits;
  //vector<int> check(base,0);
  //int sum = 0;
  int order = 0;

  while(i != 0){
    remain = i % base;
    i = i/base;
    digits.push_back(remain);
    if(remain > order) return 1;
    order++;
    //if(check[remain] == 0){
    //  check[remain] = 1;
    //  sum += 1;
    //}
  }

  //if(sum != base) return 1;

  rep(k,digits.size()){
    community_id[k] = digits[k];
    //fprintf(stderr,"digits[%d]=%d \n",k,digits[k]);
  }

  return 0;
}


long double
Graph::modularity(int type, int max_comm_index){

  vector< vector<int> > comm;

  comm.resize(max_comm_index + 1); //comm[i]: membership vector of community i

  rep(i, V){
    //fprintf(stderr,"%d \n", community_id[i]);
    comm[community_id[i]].push_back(i);
  }

  long double Q = 0;
  long double p = 0;

  //directed is not implemented correctly yet!!

  rep(k, max_comm_index + 1){
    //fprintf(stderr,"Community Size: %d \n", comm[k].size());
    rep(i, comm[k].size()){
      rep(j, comm[k].size()){
        //Q += A_ij - P_ij
        long double a = 0;
        rep(l, adj[comm[k][i]].size()){
          if(adj[comm[k][i]][l] == comm[k][j]) a = 1.0;
        }
        if(type == 0){
          p = k_in[comm[k][i]]*k_out[comm[k][j]]/(2.0*original_E);
          Q += a-p;
        }else if(type == 1){
          p = (k_in[comm[k][i]]*k_out[comm[k][j]] + k_out[comm[k][i]]*k_in[comm[k][j]])/original_E;
          //if(appear_time[comm[k][i]] < appear_time[comm[k][j]]){
          Q += a-p;
          //}else{
          //  Q -= p;
          //}
        }else if(type == 2){
          if(appear_time[comm[k][i]] < appear_time[comm[k][j]]) Q += a - (k_in[comm[k][i]]*k_out[comm[k][j]]);
          if(appear_time[comm[k][i]] > appear_time[comm[k][j]]) Q += a - k_out[comm[k][i]]*k_in[comm[k][j]]; 

        }
      }
    }
  }

  if(type != 2){
    Q /= 2.0*original_E;
  }else{
    Q /= 2.0*original_E;
  }


  return(Q);

}


////////////////////////
bool
Graph::check_connectivity(int comm,vector< pair<int,int> > &parent,int &parent_lcc){

  //check whether community comm is connected
  //if not, assign each connected component a new community index

  vector<int> size_of_component(V,1); //size_of_component storages the size of the component with parent i
  parent_lcc = 0;
  int number_components = V; //counts the number of disconnected components
  rep(i,V){
    parent.push_back(make_pair(i,i)); //first entry is parent, second entry is next node in same component
  }
  int u,v,foo; 

  tr(new_edges, it){
    u = it->first, v = it->second; 

    //if parents are different but nodes are connected by link, we have to merge the components
    if(parent[u].first != parent[v].first){
      //union
      number_components--; //merging two components reduces number of components by 1
      //adapt parent of larger component
      if(size_of_component[parent[u].first] > size_of_component[parent[v].first]){
        //change all parents for which parent is same as for u 
        //update component size
        size_of_component[parent[u].first] += size_of_component[parent[v].first];
        //update parent_lcc
        if(size_of_component[parent[u].first] > size_of_component[parent_lcc]) parent_lcc = parent[u].first;

        parent[v].first = parent[u].first;

        //We want to create a linked loop that connects all nodes in one component
        //link u to v and element before v to element next to u
        //i.e. foo will become parent[v].second
        foo = parent[u].second;
        //link u to v
        parent[u].second = v;

        //update parent index for all nodes in the smaller component, i.e. follow loop until we return to v
        int it = parent[v].second;
        int it_old = v;
        //one loop
        while(it != v){
          assert( parent[it].first != parent[u].first);
          parent[it].first = parent[u].first;
          it_old = it;
          it = parent[it].second;
        }
        //link element before v to element next to u
        parent[it_old].second = foo;
      }else{
        //analogue to above with u and v interchanged
        size_of_component[parent[v].first] += size_of_component[parent[u].first];
        if(size_of_component[parent[v].first] > size_of_component[parent_lcc]) parent_lcc = parent[v].first;

        parent[u].first = parent[v].first;
        foo = parent[v].second;
        parent[v].second = u;
        int it = parent[u].second;
        int it_old = u;

        while(it != u){
          assert( parent[it].first != parent[v].first);
          parent[it].first = parent[v].first;
          it_old = it;
          it = parent[it].second;
        }
        assert(parent[it_old].second == u);
        parent[it_old].second = foo;
      }
    }
  }



  //-disp
  //fprintf(stderr,"Number of components = %d\n",number_components);
  assert(number_components > 0);

  bool is_connected = 0;
  if(number_components > 1) is_connected = 1;

  return(is_connected);

}


/////////////////////////////////////
void
Graph::set_appeartime(char filename[]){

  FILE *fp_r;
  char buf[100];

  if(strcmp(filename, "NA") == 0){
    appear_time.resize(V);
    rep(i, V){
      appear_time[i] = i;
    }
    t_max = V;
  }
  else{
    if((fp_r = fopen(filename, "r")) == NULL){
      //-disp
      //fprintf(stderr, "file open error-%s\n", filename);
      exit(EXIT_FAILURE);
    }

    appear_time.resize(V);

    int i = 0;
    int tmp_t;
    t_max = 0;
    while(fgets(buf, 100, fp_r) != NULL){
      sscanf(buf, "%d", &tmp_t);
      appear_time[i] = tmp_t;
      i++;

      if(t_max < tmp_t){
        t_max = tmp_t;
      }
    }
    fclose(fp_r);

    //Taro: increment t_max with the samme reason we increment g.V in readGraph
    t_max++;
    assert(i == V);
  }
}


///////////////////////////////
void
Graph::set_null_model(int type,long double alpha){

  //calculates k(appa)^in and k(appa)^out for the different models
  assert((type==0)||(type==1)||(type==2));

  k_in.resize(V, 0.0);
  k_out.resize(V, 0.0);

  if(type == 0){
    // undirected
    rep(i, V){
      k_in[i] = (long double)adj[i].size(); // degree after discarding link directions
      k_out[i] = k_in[i];
      //fprintf(stderr, "[%d] k_in:%d k_out:%d\n", i, (int)k_in[i], (int)k_out[i]);
    }
  }
  else if(type == 1){
    // directed
    rep(e, (int)edges.size()){
      k_in[edges[e].second]++; // in-degree
      k_out[edges[e].first]++; // out-degree
    }
    /*
       rep(i, V){
       fprintf(stderr, "[%d] k_in:%d k_out:%d\n", i, (int)k_in[i], (int)k_out[i]);
       }
       */
  }
  else if(type == 2){
    // ordered
    //fprintf(stderr, "t_max:%d\n", t_max);

    // ordered: null model == modified in / out-degree
    rep(e, (int)edges.size()){
      k_in[edges[e].second]++; // in-degree
      k_out[edges[e].first]++; // out-degree
    }

    //debug
    /*
       int t_min = *min_element(appear_time.begin(), appear_time.end());
       fprintf(stderr, "t_min:%d\n", t_min);
       rep(i, V){
       if(appear_time[i] == t_min){
       fprintf(stderr, "[%d] k_in:%Lf k_out:%Lf k_in-k_out:%Lf\n", i, k_in[i], k_out[i], k_in[i]-k_out[i]);
       }
       }*/

    int t;

    //vector<long double> sum_k_in(t_max, 0), sum_k_out(t_max, 0);
    vector<long double> mu(t_max, 0.0), lambda(t_max, 0.0);

    //calculate mu, lambda and sum_k
    rep(i, V){
      long double value = k_in[i] - k_out[i];
      for(t=appear_time[i]+1; t<t_max; t++){
        mu[t] += value;
        lambda[t] += value;
      }

      if(appear_time[i] > 0){
        lambda[appear_time[i]] -= k_out[i];
      }

      //sum_k_in[appear_time[i]] += k_in[i];
      //sum_k_out[appear_time[i]] += k_out[i];
    }

    /*
    //debug
    rep(i,t_max){
    if(mu[i] <= 0.0){
    fprintf(stderr, "t_max:%d\n", t_max);
    fprintf(stderr, "mu[%d]:%Lf\n", i, mu[i]);
    }
    //assert(mu[i] > 0.0);
    }
    */

    //fprintf(stderr, "t mu lambda sum_k_in sum_k_out\n");
    rep(t, t_max){
      //fprintf(stderr, "t:%d mu:%Lf lambda:%Lf\n", t, mu[t], lambda[t]);
      //printf("%d %Lf %Lf %Lf %Lf\n", t, mu[t], lambda[t], sum_k_in[t], sum_k_out[t]);
      //fprintf(stderr, "t:%d sum_k_in:%Lf sum_k_out:%Lf\n", t, sum_k_in[t], sum_k_out[t]);

      if((t != 0)&&(t != t_max-1)&&(lambda[t] == 0)){
        fprintf(stderr, "there is a zero flux point at t = %d\n", t);
        exit(EXIT_FAILURE);
      }
    }

    vector<long double> a(t_max), b(t_max);    

    if(t_max == 1){
      fprintf(stderr, "t_max should be equal to or larger than 2.\n");
      exit(EXIT_FAILURE);
    }
    else if(t_max > 2){
      // if t_max == 2, the network is a directed bipartite graph
      // in which all the directed links go from nodes at t=1 to nodes at t=0.
      // the null model is equivalent to that in the directed case.

      ///??didn't quite understand why a and b can be calculated this way
      /*
         a[0] = 1.0;
         for(t=1; t<t_max-1; t++){
         a[t] = a[t-1]*(1.0 + sum_k_out[t]/lambda[t]);
         }
         a[t_max-1] = 0.0;

         b[t_max-1] = 1.0;
         for(t=t_max-2; t>=1; t--){
         b[t] = b[t+1]*(1.0 + sum_k_in[t]/lambda[t]);
         }
         b[0] = 0.0;
         */

      /*
         vector<double> cumprod_lambda(t_max,1.0), cumprod_mu(t_max,1.0);

         for(t = 1; t < t_max-1; t++){
         cumprod_lambda[t] = cumprod_lambda[t-1] * lambda[t];
         cumprod_mu[t] = cumprod_mu[t-1] * mu[t];
         fprintf(stderr,"%Lf, %Lf\n",cumprod_lambda[t], cumprod_mu[t]);
         }
         cumprod_mu[t_max-1] = cumprod_mu[t_max-2] * mu[t_max-1];

         a[0] = 1.0;
         for(t = 1; t < t_max; t++){
         a[t] = cumprod_lambda[t-1]/cumprod_mu[t];
         }

         cumprod_lambda[t_max-1] = 1.0;
         cumprod_mu[t_max-1] = mu[t_max-1];
         for(t = t_max-2 ; t >= 1; t--){
         cumprod_lambda[t] = cumprod_lambda[t+1] * lambda[t];
         cumprod_mu[t] = cumprod_mu[t+1] * mu[t];
         }
         cumprod_lambda[0] = 1.0;
         cumprod_mu[0] = 1.0;

         b[t_max-1] = 1.0;
         for(t = t_max-2;t >= 0; t--){
         b[t] = cumprod_lambda[t+1]/cumprod_mu[t+1];
         }
         */


      // a is f(1,l_i) and b is f(l_i,L)
      a[0] = 1.0;
      a[1] = 1.0/mu[1];
      //fprintf(stderr,"a[1] = %Lf \n",lambda[1]/(mu[1]*mu[2]));  
      for(t=2; t<t_max; t++){
        a[t] = a[t-1]*(lambda[t-1]/mu[t]);
        //fprintf(stderr,"a[%d] = %Lf \n",t,a[t]);  
      }

      //for(t=1; t<t_max; t++){
      //  a[t] /=mu[1];
      //}

      /*
         b[t_max-1] = 1.0;
         b[t_max-2] = 1.0/mu[t_max-1];
         for(t=t_max-3; t>=0; t--){
         b[t] = b[t+1]*(lambda[t+1]/mu[t+1]);
      //fprintf(stderr,"b[%d] = %Lf \n",t,b[t]);  
      }
      */


      b[0] = 1.0;
      //fprintf(stderr,"a[1] = %Lf \n",lambda[1]/(mu[1]*mu[2]));  
      for(t=1; t<t_max-1; t++){
        b[t] = b[t-1]*(mu[t]/lambda[t]);
        //fprintf(stderr,"a[%d] = %Lf \n",t,a[t]);  
      }
      b[t_max-1] = b[t_max-2] * mu[t_max-1];

      //for(t=1; t<t_max-1; t++){
      //  b[t] = b[t]*mu[t]/lambda[t];
      //  a[t] = a[t]/lambda[t];
      //}


      /*
         for(t=0;t<t_max;t++){
         fprintf(stderr, "t:%d a:%Lf b:%Lf\n", t, a[t], b[t]);
         }
         */

      //long double prod_lambda = 1.0;
      //long double prod_mu = mu[t_max-1];
      //long double f1n = 1.0/mu[t_max-1];
      ///
      //for(t=1; t<t_max-1; t++){
      //  f1n *= (lambda[t]/mu[t]);
      //}
      //f1n /= mu[t_max-1];

      ///f1n *= original_E; //this is canceled out below in product_order so essentially we don't need it here
      //printf(stderr, "original_E:%Lf prod_lambda:%Lf prod_mu:%Lf\n", original_E, prod_lambda, prod_mu);
      //fprintf(stderr, "f1n:%Lf\n", f1n);

      rep(i, V){
        ///k_in[i] = f1n*a[appear_time[i]]*k_in[i];
        k_in[i] = b[appear_time[i]]*k_in[i];
        k_out[i] = a[appear_time[i]]*k_out[i];
        //fprintf(stderr, "[%d] k_in:%.15g k_out:%.15g\n", i, k_in[i], k_out[i]);
      }
    }
    ///debug
    ///fprintf(stderr, "[%d] k_in:%Lf k_out:%Lf\n", 1, k_in[10322], k_out[11191]); 
  }

  rep(i,V){
    k_out[i] = alpha * k_out[i]; //for detectability
    k_in[i] = alpha * k_in[i]; //for detectability
  }

  // calculate the lower bound of eigenvalues of the modularity matrix
  /*
     if(type == 0){
     lambda_inf = -2.0*(*max_element(k_in.begin(), k_in.end());
     }
     else if(type == 1){
     lambda_inf = -2.0*(*max_element(k_in.begin(), k_in.end()));
     }
     else if(type == 2){
     vector<long double> s0(V, 1.0), s1(V);
     product_ordered(s0, s1);

     rep(i, V){
     s1[i] -= 2.0*k_in[i];
     }

     lambda_inf = *min_element(s1.begin(), s1.end());
     }
     */
  // for check
  //lambda_inf = 0.0;
  //fprintf(stderr, "lambda_inf:%Lf\n", lambda_inf);
}


////////////////////////////////
vector<long double>
Graph::power_method(int type){

  //calculate the dominant eigenvector of B (of the corresponding method) 

  vector<long double> new_s(V,0);
  vector<long double> old_s(V);
  vector<long double> old2_s(V);
  vector<long double> diff_s(V);
  vector<long double> one(V,1);
  vector<long double> new_one(V);


  //?? probably product is pointer to function but why does this declaration work?
  fp product;
  if(type == 0){
    product = &Graph::product_undirected;
  }
  else if(type == 1){
    product = &Graph::product_directed;
  } 
  else if(type == 2){
    product = &Graph::product_ordered;
  } 

  long double norm;
  long double EPS_power;
  long double lambda_inf, lambda_new,lambda_old;
  lambda_new = 0.0;
  lambda_inf = random()/(RAND_MAX+1.0); //to make sure that when dominant eigenvalue is not unique, the method still works.
  // of course, lambda_max = -lambda_min-0.01 is still possible..
  //lambda_inf = 0.0;

  do{
    EPS_power = EPS;
    //shift eigenvalues
    lambda_inf += fabs(lambda_new);


    //fprintf(stderr, "lambda_inf:%Lf\n", lambda_inf);

    rep(i, V){
      //old_s[i] = 2.0*(long double)random()/(RAND_MAX+1.0)-1.0;
      //old_s[i] = 2.0*(random() % 2) -1.0;
      old_s[i] = random()/(RAND_MAX+1.0);
      //fprintf(stderr,"s[%d]=%Lf\n",i,old_s[i]);
    }

    int iter = 0;

    lambda_new = 0.0;

    do{
      //fprintf(stderr, "iter:%d\n", iter);
      lambda_old = lambda_new;
      /*
         fprintf(stderr, "old_s:");
         rep(i, min(V,20)){
         fprintf(stderr, "%Lf ", old_s[i]);
         }     
         fprintf(stderr, "\n");
         */
      // new_s = (A - P)*old_s

      ///leo: this points to the objects whose member function is being excecuted, i.e. in this case the graph

      //calculate old2_s=[B^(g)+|lambda|I]*s
      //leo: 	make sure the method converges to the leading eigenvector 
      //			i.e. with positive eigenvalue
      // B^(g) + |lambda|I (shift eigenvalue)
      (this->*product)(old_s, old2_s);
      (this->*product)(one,new_one);
      rep(i,V){
        old2_s[i] -= (new_one[i]*old_s[i]);
        old2_s[i] += (lambda_inf*old_s[i]);
      }

      // normalize old2_s
      norm = sqrt(innprod(old2_s, old2_s));
      if(norm == 0.0){
        //-disp
        //fprintf(stderr,"break\n");
        break;
      }
      else{
        rep(i, V){
          old2_s[i] /= norm;         
        }
      }
      //fprintf(stderr,"norm %d, %Lf\n",iter,norm);

      //calculate [B^(g)+|lambda|I]*s again
      (this->*product)(old2_s, new_s);
      (this->*product)(one,new_one);
      rep(i,V){
        //assert(new_one[i]==0);
        new_s[i] -= (new_one[i]*old2_s[i]);
        new_s[i] += (lambda_inf*old2_s[i]);
      }


      //Graph::product_ordered(old_s, new_s);
      /*
         if(iter == 0){
         fprintf(stderr, "null model matrix:\n");
         rep(i, V){
         rep(j, V){
         if(i < j){
         fprintf(stderr, "%Lf, ", k_in[i]*k_out[j]/original_E);
         }
         else{
         fprintf(stderr, "0, ");	   
         }
         }
         }
         fprintf(stderr, "\n");
         }
         */
      // new_s = (A - lambda_inf*I - P)*old_s


      //fprintf(stderr,"lambda_old=%Lf\n",lambda_new);



      // approx eigenvalue
      //lambda_new = innprod(new_s, new_s) / innprod(old_s, new_s);
      // changed this but since we don't use the eigenvalue it shoudn't have had any impact
      lambda_new = innprod(new_s, old2_s) / innprod(old2_s, old2_s);
      //fprintf(stderr, "lambda_old:%Lf lambda_new:%Lf\n", lambda_old, lambda_new);
      //fprintf(stderr, "%Lf\n", lambda_new);

      // normalize new_s and copy new_s to s
      norm = sqrt(innprod(new_s, new_s));
      //fprintf(stderr, "norm:%Lf\n", new_s[0]);
      if(norm == 0.0){
        fprintf(stderr,"break\n");
        break;
      }
      else{
        //fprintf(stderr, "new_s:");
        rep(i, V){
          new_s[i] /= norm;
          diff_s[i] = new_s[i] - old_s[i];
          old_s[i] = new_s[i];
        }
        //fprintf(stderr, "\n");
      }

      //convergence criterion
      //norm = fabs(lambda_new - lambda_old);

      norm = sqrt(innprod(diff_s, diff_s));
      //fprintf(stderr,"norm=%Lf\n",norm);
      //let it run at least ITER_MAX steps
      iter++;
      //if(iter < ITER_MAX) norm = 0.1;

      //to make sure that when the smallest element is 
      //eventually close to zero, we increase the accuracy a little bit.
      //if(minimum < EPS_power*100 && EPS_power > EPS*0.0001 && norm < 0.01){
      //  EPS_power = minimum*0.01;
      //}

    }while( norm > EPS_power);
    //}while(iter < 10*ITER_MAX);
    //}while(norm > EPS);

    //fprintf(stderr,"Power method is repeated until norm < EPS_power, where\n EPS_power * 10^10 = %Lf\n",EPS_power*1e10);
    //fprintf(stderr,"The element of new_s with minimal magnitude:\n s[%d] * 10^10 = %Lf\n",m,new_s[m]*1e10);
    //-disp
    //fprintf(stderr, "Number of iterations in power method: %d\n", iter);
}while(lambda_new < 0.0);
///debug
///fprintf(stderr, "%Lf\n", lambda_new);

//debug
/*
   vector<long double> check(V);
   vector<long double> lambda_check(V);
   vector<long double> one_check(V);
   long double mu_check=0;
   long double sigma_check=0;
   (this->*product)(new_s,check);
   (this->*product)(one,one_check);
   long double sign=0;
   rep(i,V){
   check[i] -= one_check[i];
//if(sign==0){
// sign = check[i]/fabs(check[i]);
//}
//check[i] *= sign;
if(new_s[i]!=0 && check[i]!=0){
lambda_check[i]=check[i]/new_s[i];
}
sigma_check += (lambda_check[i]*lambda_check[i]);
mu_check += lambda_check[i];
}
sigma_check /= V;
mu_check /= V;
fprintf(stderr,"Variance of lambda: %Lf\n",sigma_check-(mu_check*mu_check));
fprintf(stderr,"lambda=%Lf\n", lambda_check[0]);
*/
//fprintf(stderr,"lambda=%Lf\n", lambda_new-lambda_inf);
//rep(i,V){
//  fprintf(stderr,"s[%d] = %Lf , ",i,new_s[i]);
//}
//fprintf(stderr,"\n");

return new_s;
}


/////////////////////////////////////////////////////
void
Graph::product_undirected(vector<long double> &s, vector<long double> &new_s){

  //leo: calculate new_s = Bs (= As-k(k^Ts)/2M )	
  //leo: if false then error is returned
  assert((int)s.size() == V);
  long double kx = 0.0;
  rep(i, V){
    new_s[i] = 0.0;

    // Ax
    for(int j=row[i]; j<row[i+1]; j++){
      new_s[i] += s[col[j]];
    }

    //(k^T s)
    kx += (k_out[i]*s[i]);
  }

  // k*(k^T s)/(2*E)
  long double C = kx/(2.0*original_E);
  rep(i, V){
    new_s[i] -= (C*k_out[i]);
  }

}


//////////////////////////////////////////////
void
Graph::product_directed(vector<long double> &s, vector<long double> &new_s){
  //leo: s is old_s, new_s is new_s
  //leo: this calculates  (A+A^T)*s+k_in*(k_out^T s)/M + k_out*(k_in^T s)/M

  assert((int)s.size() == V);

  new_s.resize(V, 0.0);

  long double k_out_x, k_in_x;
  k_out_x = k_in_x = 0.0;
  rep(i, V){
    // Ax
    //leo: new_s=A*s (possible to write in this form because A has only entries in {0,1})
    new_s[i] = 0.0;
    for(int j=row[i]; j< row[i+1]; j++){
      new_s[i] += s[col[j]];
    }

    // k_in*(k_out^T s) + k_out*(k_in^T s)
    k_out_x += k_out[i]*s[i];
    k_in_x += k_in[i]*s[i]; 
  }

  long double C_out = k_out_x/original_E;
  long double C_in = k_in_x/original_E;
  rep(i, V){
    new_s[i] -= (k_in[i]*C_out + k_out[i]*C_in);
  }
}


/////////////////////////////////////////////
void
Graph::product_ordered(vector<long double> &s, vector<long double> &new_s){
  //leo: this calculates Ax+sum(kappa^in(l)(kappa^out(l)^T x)) (l layer)


  assert((int)s.size() == V);

  // sum_{t=1}^{t/max} (k_in(t)*(k_out(t)^T s) + k_out(t)*(k_in(t)^T s)
  long double k_out_x, k_in_x;
  vector<long double> p(V, 0.0);
  //vector<long double> v1(V), v2(V);
  int t;
  int t_min = *min_element(appear_time.begin(), appear_time.end());

  for(t=t_min; t<t_max-1; t++){

    // inner product (k_out(t)^T s) and (k_in(t)^T s)
    k_out_x = k_in_x = 0.0;
    rep(i, V){
      if(appear_time[i] > t){
        k_out_x += k_out[i]*s[i];
      }
      else if(appear_time[i] == t){
        k_in_x += k_in[i]*s[i];
      }
    }

    // vector multiplied by a scalar
    // (k_in(t)*(k_out(t)^T s) and k_out(t)*(k_in(t)^T s)
    rep(i, V){
      if(appear_time[i] > t){
        p[i] += k_out[i]*k_in_x;
      }
      else if(appear_time[i] == t){
        p[i] += k_in[i]*k_out_x;
      }
    }
  }

  rep(i, V){
    new_s[i] = 0.0;
    // Ax
    for(int j=row[i]; j< row[i+1]; j++){
      new_s[i] += s[col[j]];
    }

    //leo devided by original_E because k_in was multiplied by original_E in null model
    // sum_{t=0}^{t_max-2} (k_in(t)*(k_out(t)^T s) + k_out(t)*(k_in(t)^T s)
    //new_s[i] -= p[i]/original_E;
    new_s[i] -= p[i];
  }
}


/////////////////////////
Graph
Graph::extract(int comm){

  Graph tmp_g;

  //fprintf(stderr,"edge %d %d\n",this->edges[0].first,this->edges[0].second);

  // Node i in the original network has index j in the extracted network.  
  map<int, int> new_index;
  vector<int> old_index(this->V);
  tmp_g.ids.resize(this->V);
  tmp_g.V = 0;
  rep(i, this->V){
    if(this->community_id[i] == comm){
      tmp_g.ids[tmp_g.V] = this->ids[i];
      new_index[i] = tmp_g.V;
      old_index[tmp_g.V] = i;
      tmp_g.V++;
    }
  }
  old_index.resize(tmp_g.V);
  tmp_g.ids.resize(tmp_g.V);


  // copy a part of adjacency list
  // set E
  tmp_g.adj.resize(tmp_g.V);
  tmp_g.original_E = this->original_E;
  tmp_g.E = 0;
  vector<int> res;
  rep(i, tmp_g.V){
    res.clear();

    set_intersection(this->adj[old_index[i]].begin() , this->adj[old_index[i]].end(), old_index.begin(), old_index.end(), back_inserter(res));

    tmp_g.adj[i].resize((int)res.size());
    rep(j, res.size()){
      tmp_g.adj[i][j] = new_index[res[j]];
      assert((0 <= new_index[res[j]])&&( new_index[res[j]] < tmp_g.V));
    }

    tmp_g.E += tmp_g.adj[i].size();
  }
  assert(tmp_g.E % 2 == 0);
  tmp_g.E /= 2;

  //fprintf(stderr,"edge %d %d\n",this->edges[0].first,this->edges[0].second);
  //create edgelist for tmp_g (this is implemented for use in check_connectivity)
  tmp_g.new_edges.resize(tmp_g.E);
  int u,v,link_count = 0;
  tr(this -> edges,it){
    u = it -> first, v = it -> second;
    if(u==v) fprintf(stderr,"u = %d, v= %d \n",u,v);
    assert(u != v);
    if((this -> community_id[u] == comm) && (this -> community_id[v] == comm)){
      /*
         vector< int >::iterator cIter1 = find(old_index.begin(),old_index.end(),u);
         vector< int >::iterator cIter2 = find(old_index.begin(),old_index.end(),v);
         if( (cIter1 == old_index.end()) || (cIter2 == old_index.end())){
         fprintf(stderr,"u or v is not in old_index\n");
         }
         */
      tmp_g.new_edges[link_count].first = new_index[u];
      tmp_g.new_edges[link_count].second = new_index[v];
      link_count++;
    }
  }
  assert(link_count == tmp_g.E);
  //fprintf(stderr,"link_count= %d , tmp_g.E= %d \n",link_count,tmp_g.E);


  // convert the adjacency list to a sparse matrix representation
  // copy community_id, k_out, k_in
  tmp_g.col.resize(2*tmp_g.E);
  tmp_g.row.resize(tmp_g.V+1);
  tmp_g.community_id.resize(tmp_g.V);
  tmp_g.k_in.resize(tmp_g.V);
  tmp_g.k_out.resize(tmp_g.V);
  int count=0;
  tmp_g.row[0] = 0;
  rep(i, tmp_g.V){
    rep(j, (int)tmp_g.adj[i].size()){
      tmp_g.col[count] = tmp_g.adj[i][j];
      count++;
    }
    tmp_g.row[i+1] = count;

    tmp_g.community_id[i] = comm;

    tmp_g.k_in[i] = this->k_in[old_index[i]];
    tmp_g.k_out[i] = this->k_out[old_index[i]];
  }  

  // copy a part of appear_time
  if((int)this->appear_time.size() == this->V){
    tmp_g.t_max = 0;
    tmp_g.appear_time.resize(tmp_g.V); 
    rep(i, tmp_g.V){  
      tmp_g.appear_time[i] = this->appear_time[old_index[i]];
      if(tmp_g.t_max < tmp_g.appear_time[i]){
        tmp_g.t_max = tmp_g.appear_time[i];
      }
    }
    tmp_g.t_max++;
  }

  // copy lambda_inf
  //tmp_g.lambda_inf = this->lambda_inf;

  // set number_community
  tmp_g.number_community = 1;

  return tmp_g;
}


/////////////////////////
long double
Graph::recalculate_tunedQ(int type,vector<long double> &s,long double &sum_skout,long double &sum_skin,int index){

  long double res = 0.0;

  //if type 0 or 1 we have to divide sum_sk-s[index]k[index] by E or 2E
  if(type == 0){

    sum_skin -= (k_in[index]*s[index]);

    //l'th row of A*s
    for(int j=row[index]; j<row[index+1]; j++){
      res += s[col[j]];
    }

    res -= (k_out[index]*sum_skin)/(2.0*original_E); 

  }else if(type == 1){

    sum_skout = (sum_skout - (k_out[index]*s[index]));
    sum_skin = (sum_skin - (k_in[index]*s[index]));

    for(int j=row[index]; j<row[index+1]; j++){
      res += s[col[j]];
    }

    res -= ((k_out[index]*sum_skin) + (k_in[index]*sum_skout))/(1.0*original_E);

  }else if(type == 2){

    rep(i,V){
      if(appear_time[i] >= appear_time[index]) sum_skin -= (k_in[i]*s[i]);
      if(appear_time[i] <= appear_time[index]) sum_skout -= (k_out[i]*s[i]);
    }

    for(int j=row[index]; j<row[index+1]; j++){
      res += s[col[j]];
    }
    res -= ((k_out[index]*sum_skin) + (k_in[index]*sum_skout));
  }
  //res /= 4.0*original_E;
  res *= 4.0*s[index];

  //if res is negative, we have an increase in Q if we change the community of i  
  return(res); 
}


/////////////////////////
void
Graph::fine_tuning(int type,vector<long double> &s){


  rep(i, s.size()){
    if(s[i] > 0){
      s[i] = 1.0;
    }
    else{
      s[i] = -1.0;
    }
  }

  long double sum_skout = 0, sum_skin = 0;

  rep(i,s.size()){
    sum_skout += k_out[i]*s[i];
    sum_skin += k_in[i]*s[i];
  }

  //long double tuned_Q = diff_modularity(type,s)/(4.0*original_E); //we want to do better than this
  long double tuned_Q = diff_modularity(type,s);
  long double res = 0.0;

  int changes;
  int it = 0;

  if(0){
    //in this tuning method, we are always changing the community of a node if it is beneficial
    //i.e. when res < 0  
    do{
      it++;
      changes = 0;
      rep(i,s.size()){
        //this operation changes sum_skout and sum_skin which is corrected below
        //recalculate_tunedQ only recalculates the term affected by s[i]
        //O(max degree) (i.e. usually O(log(n)))
        res = recalculate_tunedQ(type,s,sum_skout,sum_skin,i);
        if(res < 0){
          //if(it > 0) fprintf(stderr,"tuning/n");
          tuned_Q -= res;
          s[i] = -s[i];
          changes++;
        }
        if(type != 2){
          //sum_skout and sum_skin updated (because recalculate_tunedQ has changed the values above)
          sum_skout += k_out[i]*s[i];
          sum_skin += k_in[i]*s[i];
        }else{
          rep(j,V){
            if(appear_time[j] >= appear_time[i]) sum_skin += (k_in[j]*s[j]);
            if(appear_time[j] <= appear_time[i]) sum_skout += (k_out[j]*s[j]);
          }
        }
      }
      fprintf(stderr,"t_Q =%Lf, res=%Lf \n",tuned_Q,res);
      if( (V/2 < it) && (it < V/2+1) ) fprintf(stderr,"Already tuned %d times \n",it);
    }while(changes != 0 && it < 1000);
  }

  if(0==0){
    //in this tuning method, we are changing the community of the node, for which we have a maximal increase
    //in modularity. Every node is allowed to change the community change_bound times

    int change_bound = 1;

    vector<int> changed(V,0); //stores how often node i has changed communities
    long double res_max=0; //stores current maximum of res
    int change_index = V+1 ; //stores argmax of res_max
    do{
      it++; //to track the total number of changes

      //reinitalize
      change_index = -1.0;
      res_max = 0.0;

      rep(i,s.size()){
        //!!this operation changes sum_skout and sum_skin, which is therefore corrected below
        //the function recalculate_tunedQ only recalculates the term affected by s[i]
        //Complexity: O(max degree) (i.e. usually O(log(n)))
        long double foo = sum_skout, goo = sum_skin;
        res = recalculate_tunedQ(type,s,sum_skout,sum_skin,i);

        sum_skout = foo;
        sum_skin = goo;

        //if the change in modularity is larger than res_old, it is stored as the current maximum
        if((res < res_max) && (changed[i] < change_bound)){
          //fprintf(stderr,"%d, res = %Lf \n",it, res);
          res_max = res;
          change_index = i;
        } 
      }

      //if there is a change to make, we make the change here
      if((changed[change_index] < change_bound) && (change_index != -1)){   
        tuned_Q -= res_max;
        s[change_index] = -1.0 * s[change_index];
        sum_skout += 2.0 * k_out[change_index]*s[change_index];
        sum_skin += 2.0 * k_in[change_index]*s[change_index];
        changed[change_index]++;  
      }else change_index = V+1;    
      //fprintf(stderr,"t_Q =%Lf, change=%Lf \n",tuned_Q,change_index);

    }while(change_index != V+1);
  }


  if(0){
    //in this tuning method, we are changing the community of the node, for which we have a maximal increase
    //in modularity. Every node is allowed to change the community change_bound times

    int change_bound = 1;
    long double tuned_Q_tmp = tuned_Q;
    long double tuned_Q_max = tuned_Q;    
    do{
      vector<long double> s_tmp=s;
      tuned_Q = tuned_Q_max;
      tuned_Q_tmp = tuned_Q;
      vector<long double> best_partition=s;  //stores the changes made for the partition of tuned_Q_max
        
      long double sum_skout_best=sum_skout;
      long double sum_skin_best=sum_skin;
      vector<long double> changed(V,0); //stores how often node i has changed communities
      long double res_max; //stores current maximum of res
      int change_index; //stores argmax of res_old
      do{
        //reinitalize

        long double foo = sum_skout, goo = sum_skin;
        int tmp=0;
        while(changed[tmp] == change_bound){
          tmp++;
        } 
        assert(tmp < V);
        res_max = recalculate_tunedQ(type,s_tmp,sum_skout,sum_skin,tmp)+1.0; 
        sum_skout = foo;
        sum_skin = goo;
       
        rep(i,s_tmp.size()){
          //!!this operation changes sum_skout and sum_skin, which is therefore corrected below
          //the function recalculate_tunedQ only recalculates the term affected by s[i]
          //Complexity: O(max degree) (i.e. usually O(log(n)))
          foo = sum_skout, goo = sum_skin;
          res = recalculate_tunedQ(type,s_tmp,sum_skout,sum_skin,i);
          sum_skout = foo;
          sum_skin = goo;

          //if the change in modularity is larger than res_old, it is stored as the current maximum
          if((res < res_max) && (changed[i] < change_bound)){
           //fprintf(stderr,"%d, res = %Lf \n",it, res);
            res_max = res;
            change_index = i;
          } 
        }
        assert(changed[change_index] < change_bound);
        it++; //to track the total number of changes
        tuned_Q_tmp -= res_max;

        s_tmp[change_index] = - 1.0 * s_tmp[change_index];
        sum_skout += 2.0 * k_out[change_index]*s_tmp[change_index];
        sum_skin += 2.0 * k_in[change_index]*s_tmp[change_index];
        changed[change_index]++;  


        if(tuned_Q_tmp > tuned_Q_max){
          tuned_Q_max = tuned_Q_tmp;
          best_partition = s_tmp;
          sum_skout_best = sum_skout;
          sum_skin_best = sum_skin;  
        }

        //fprintf(stderr,"t_Q =%Lf, change=%Lf \n",tuned_Q,change_index);
      }while(it < V);

      if(tuned_Q_max - tuned_Q > 0) {
        s = best_partition;
        sum_skout = sum_skout_best;
        sum_skin = sum_skin_best;
      }
      //rep(i,s.size()){
      //  fprintf(stderr,"s[%d]=%Lf\n",i,best_partition[i]);
      //}
      //fprintf(stderr,"tuned_Q= %Lf \n",tuned_Q_max);
    }while((tuned_Q_max - tuned_Q)/original_E > 0);
  }






  //-disp
  //fprintf(stderr,"The output of the spectral method has been tuned %d times.\n",it);
  //fprintf(stderr,"tuned_Q= %Lf \n", tuned_Q);
}



////////////////////////////////////////////
void
Graph::divide_community(int comm, vector<long double> &s){

  // Community comm is divided into two communities:
  // comm and number_community+1
  // Then number_community is inclimented with unity.
  vector<int> comm_member_ids((int)s.size(),0);
  int count = 0;
  long double sum_s = 0.0;
  rep(i, V){
    if(community_id[i] == comm){
      comm_member_ids[count] = i;
      count++;
    }
  }
  assert(comm_member_ids.size() == s.size());

  rep(i,s.size()){
    assert(fabs(s[i])==1);
    sum_s += s[i];
  }

  // the larger part is always assigned with a smaller community id.
  int sign_s = (int)(sum_s/fabs(sum_s));
  //int sign_s = 1;


  rep(i, s.size()){
    if(sign_s*s[i] > 0.0){
      community_id[comm_member_ids[i]] = comm;
    }
    else{
      community_id[comm_member_ids[i]] = number_community;
    }
  }

  number_community++;
}



Graph
Graph::divide_components(int comm, int type, long double &Q, Graph tmp_g){

  //still todo: giant component

  //parent vector
  vector< pair<int,int> > parent;
  int parent_lcc;
  //get parent list
  bool is_connected = 0;
  is_connected = tmp_g.check_connectivity(comm,parent,parent_lcc);
  if(is_connected == 0) return(tmp_g);

  //relate indices of tmp_g to the indices of g in community comm
  vector<int> comm_member_ids;
  rep(i, V){
    if(community_id[i] == comm){
      comm_member_ids.push_back(i);
    }
  }

  int count_check=0;
  vector< long double > s(tmp_g.V,1.0);

  rep(i,tmp_g.V){
    int iter = parent[i].first;
    if((community_id[comm_member_ids[iter]] == comm) && (iter != parent_lcc)){
      count_check++;
      do{
        community_id[comm_member_ids[iter]] = number_community;
        iter = parent[iter].second;
        s[iter] = -1.0;
      }while(iter != parent[i].first);
      //number_community++;
    }  
  }
  number_community++;
  //fprintf(stderr,"count_check=%d\n",count_check);

  //assert(tmp_g.diff_modularity(type, s) > 0);  

  //Q += (tmp_g.diff_modularity(type, s) / (4.0*original_E));
  Q += tmp_g.diff_modularity(type, s);

  //extract tmp_g, which is now connected
  tmp_g = this->extract(comm);

  //fprintf(stderr,"number_community = %d\n",number_community);

  return(tmp_g);

}


//////////////////////////////////
long double
Graph::diff_modularity(int type, vector<long double> &s){

  //leo: this calculates the modularity
  //leo: s is the dominant eigenvector of B

  //leo: "approximate" the dominant eigenvector such that it only has entries in {-1,1}
  //debug
  //rep(i,V){
  //	fprintf(stderr,"s[i]=%Lf\n",s[i]);
  //

  int size=0;


  rep(i, V){
    //debug
    //if(-EPS < s[i] && s[i] < EPS) fprintf(stderr,"s[%d]=%Lf \n",i,s[i]);

    if(s[i] > 0){
      s[i] = 1.0;
      if(s[0]>0) size++;
    }
    else{
      s[i] = -1.0;
      if(s[0]<=0) size++;
    }
  }

  //fprintf(stderr,"size of community of node 0 = %d / %d\n",size,V);

  vector<long double> u((int)s.size(), 1.0);
  vector<long double> st((int)s.size()), ut((int)s.size());

  if(type == 0){
    // undirected
    //leo: st = B*s
    product_undirected(s, st);
    //leo: ut = B*u (=B*1)
    product_undirected(u, ut);
  }
  else if(type == 1){
    // directed
    product_directed(s, st);
    product_directed(u, ut);
  }
  else if(type == 2){
    // ordered
    product_ordered(s, st);
    product_ordered(u, ut);
  }

  //leo: res = s^T*B*s - s^2 * B * 1
  long double res = (innprod(s, st) - innprod(u, ut));
  //long double res = innprod(s, st);
  //fprintf(stderr, "innprod(s,st):%Lf innprod(u, ut):%Lf\n", innprod(s, st), innprod(u, ut));

  return res;
}


////////////////////////////////
long double
Graph::innprod(vector<long double> &x, vector<long double> &y){
  assert((int)x.size() == (int)y.size());
  long double res = 0.0;
  rep(i, (int)x.size()){
    res += x[i]*y[i];
  }
  return res;
}


//////////////////////////
void
Graph::community_assignment(vector<vector<int> > &community_membership,vector<int> &community_size, int &max_comm_size){

  max_comm_size = 0;

  rep(i,V){
    community_membership[community_id[i]].push_back(i);
    community_size[community_id[i]]++;
  }

  rep(i,community_size.size()){
    if(community_size[i] > max_comm_size){
      max_comm_size = community_size[i];
    }
  }

}


//////////////////////////
void
Graph::get_adjacent_communities(int k,vector<int> &community_member,vector<int> &community_size,vector<int> &adjacent_comms,int max_comm_index){

  int n_adjacent_comms = 0;

  vector<int> is_included(max_comm_index,1);

  rep(i,community_member.size()){  
    rep(l,adj[community_member[i]].size()){ 
      //if((is_included[community_id[adj[i][l]]] == 1) && (community_size[community_id[adj[i][l]]] != 0) && (community_id[adj[i][l]] != k))
      if((community_size[community_id[adj[community_member[i]][l]]] != 0) && (community_id[adj[community_member[i]][l]] != k) && (is_included[community_id[adj[community_member[i]][l]]] == 1)){    
        n_adjacent_comms++;
        adjacent_comms.resize(n_adjacent_comms);
        adjacent_comms[n_adjacent_comms-1] = community_id[adj[community_member[i]][l]];
        is_included[community_id[adj[community_member[i]][l]]] = 0;
      }
    }
  }

}


//////////////////////////
long double
Graph::isolate_node(int node_index,int type,int max_comm_index,int former_community,vector<int> &community_size,vector<int> &community_membership,int is_final){

  bool include; //if 0 then include if 1 then isolate

  if(former_community != max_comm_index){
    //we want to integrate node in some community
    community_id[node_index] = former_community;
    include = 0;
  }else{
    //we want to isolate node from former community
    former_community = community_id[node_index];
    include = 1;
  }
  // Node i in the original network has index j in the extracted network.  
  map<int, int> new_index;
  int count_new_index = 0;


  rep(i, this->V){
    if(community_id[i] == former_community){
      new_index[i] = count_new_index;
      count_new_index++;
    }
  }  
  //fprintf(stderr,"new index= %d \n", new_index[0]);
  //fprintf(stderr,"count new index = %d \n",count_new_index);

  Graph tmp_g;

  tmp_g = extract(former_community);

  assert(tmp_g.get_V() == count_new_index);

  vector<long double> s(tmp_g.get_V(),1.0);


  /*
  long double sum_skout = 0, sum_skin = 0;

  rep(i,s.size()){
    sum_skout += k_out[i]*s[i];
    sum_skin += k_in[i]*s[i];
  }
  long double foo = sum_skout, goo = sum_skin;
  long double res = tmp_g.recalculate_tunedQ(type,s,sum_skout,sum_skin,new_index[node_index]);
  sum_skout = foo;
  sum_skin = goo;
  */

  s[new_index[node_index]] = -1.0;
  //fprintf(stderr,"new_index = %d\n",new_index[indices_changed[i]]);
  long double Q_comp = tmp_g.diff_modularity(type,s);

  //fprintf(stderr,"Q_tmp = %Lf , res = %Lf \n",Q_comp,res);

  if(include == 1){
    community_size[former_community]--;
    community_id[node_index] = max_comm_index;
    community_membership[0] = node_index;
    assert((int)community_membership.size() == 1);
    //fprintf(stderr,"comm_id =%d \n",community_id[node_index]);
  }
  if((include == 0) && (is_final == 0)){
    community_id[node_index]=max_comm_index;
  }
  if((include == 0) && (is_final == 1)){
    community_size[former_community]++;
    //do we have to adjust community membership?
  }
  //Q_comp is negative for isolate and negative for integrate
  return(Q_comp);
}


//////////////////////////
long double
Graph::merge_communities(int comm_tobeincluded, int comm_toinclude, int type, vector<int> &community_size, vector<vector<int> > &community_membership,int is_final){

  //find better names for comm_tobeincluded, comm_toinclude
  vector<int> indices_changed;

  //assert(community_size[comm_tobeincluded] == (int)community_membership[comm_tobeincluded].size());
  //change community id of single to comm_include
  //fprintf(stderr,"comm size =%d \n",community_size[comm_toinclude]);
  //fprintf(stderr,"comm size =%d \n",(int)community_membership[comm_tobeincluded].size());
  assert(community_size[comm_tobeincluded] == (int)community_membership[comm_tobeincluded].size());
  rep(i,community_size[comm_tobeincluded]){
    community_id[community_membership[comm_tobeincluded][i]] = comm_toinclude;
    indices_changed.push_back(community_membership[comm_tobeincluded][i]);
  }

  // Node i in the original network has index j in the extracted network.  
  map<int, int> new_index;
  int count_new_index = 0;


  rep(i, this->V){
    if(community_id[i] == comm_toinclude){
      new_index[i] = count_new_index;
      count_new_index++;
    }
  }  
  //fprintf(stderr,"new index= %d \n", new_index[0]);
  //fprintf(stderr,"count new index = %d \n",count_new_index);

  Graph tmp_g;

  tmp_g = extract(comm_toinclude);

  assert(tmp_g.get_V() == count_new_index);

  vector<long double> s(tmp_g.get_V(),1.0);

  rep(i,indices_changed.size()){
    s[new_index[indices_changed[i]]] = -1.0;
    //fprintf(stderr,"new_index = %d\n",new_index[indices_changed[i]]);
  }

  long double Q_comp = tmp_g.diff_modularity(type,s);
  //fprintf(stderr,"Q_comp = %Lf \n", Q_comp);

  //if not final merge then undo change in community assignment
  //if final merge adjust community size
  if(is_final == 0){
    rep(i,indices_changed.size()){
      community_id[indices_changed[i]] = comm_tobeincluded;
    }
  }else{
    community_size[comm_toinclude] += community_size[comm_tobeincluded];
    community_size[comm_tobeincluded] = 0;
    rep(i,indices_changed.size()){
      community_membership[comm_toinclude].push_back(indices_changed[i]);
    }
  }

  //fprintf(stderr,"Q_comp= %Lf \n",Q_comp);

  return(Q_comp);

}


//////////////////////////
void
Graph::adjust_community_indices(int max_comm_index, vector<int> &community_size){

  //if community_size == 0 then decrease comm indices of all comms larger than this index

  rep(i,max_comm_index+1){
    int inverse_count = max_comm_index - i;
    if(community_size[inverse_count] == 0){
      rep(i,V){
        if(community_id[i] > inverse_count){
          community_id[i]--;
        }
      }
    }
  }

}



//////////////////////////
int
Graph::get_V(){
  return V;
}


//////////////////////////
int
Graph::get_original_E(){
  return original_E;
}


//////////////////////////
int
Graph::get_number_community(){
  return number_community;
}


///////////////////////////
void
Graph::show_community_id(){  
  assert(V ==(int)community_id.size());
  rep(i, V){
    printf("%d\n", community_id[i]);
    //fprintf(stderr,"%d\n", community_id[i]);
  }
}
