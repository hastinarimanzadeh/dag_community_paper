
pre_proc_graph <- function(file){

  source(here::here('Rfunctions/is_ordered.R')) 
  source(here::here('Rfunctions/make_ordered.R')) 
  
  C <- read.table(file)
  C <- as.matrix(C)[,1:2]
  if(min(C) == 1) C <- C - 1 
  print(is_ordered(C))
  el <- make_ordered(C)
  print(is_ordered(el))
  
  g <- graph.edgelist(el+1)
  plot(g)
  
  write.table(el,"edgelist.txt",col.names = F, row.names = F)
  
}
