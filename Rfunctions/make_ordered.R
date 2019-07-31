
make_ordered <- function(el){
  
  el <- as.matrix(el) + 1
  
  library(igraph)
  
  g <- graph.edgelist(el)
  N <- length(V(g))
  
  if(is.dag(g)==F) return("not a DAG")
  
  #permutate node ids
  
  in_deg <- degree(g,mode="in")
  out_deg <- degree(g,mode="out")
  adj_list <- get.adjlist(g,mode="in")
  
  #print(head(in_deg))
  #print(head(out_deg))
  #return(adj_list[1:22])
  
  #take nodes with out_degree = 0 and decrease out_deg of nodes that are connected to these nodes
  
  count <- 1
  permutation <- rep(0,N)
  layer <- rep(0,N)
  current_layer <- 0
  #print(out_deg)
  
  while(max(out_deg) >= 0){
    out_deg_0 <- which(out_deg == 0)
    #print(out_deg_0)
    new_count <- count+length(out_deg_0) - 1
    permutation[out_deg_0] <- c(count:new_count)
    layer[count:new_count] <- current_layer
    count <- new_count + 1
    current_layer <- current_layer + 1
    out_deg[out_deg_0] <- out_deg[out_deg_0] - 1
    #decrease out_deg of nodes that are connected to these nodes
    for(i in out_deg_0){
      out_deg[adj_list[i][[1]]] <- out_deg[adj_list[i][[1]]] - 1
    }
  }
 
  #print(out_deg)
  if(any(out_deg == -1) == F) return("error")
  
  #print(permutation)
  #when finished do el[,1] <- permutation[el[,1]]
  
  el[,1] <- permutation[el[,1]]-1
  el[,2] <- permutation[el[,2]]-1
  
  write.table(layer,"layer.txt",row.names = F, col.names = F)
  
  return(el)
  
}