plot_ordered_graph <- function(data,method,lyt){
  
  
  file_el <- paste(data,"lcc_edgelist.txt",sep="_")
  file_layer <- paste(data,"lcc_layer.txt",sep="_")
  file_comm <- paste(data,"community",method,sep="_")
  file_comm <- paste(file_comm,".txt",sep="")
  file_out <- paste("C6",method,sep="_")
  file_out <- paste(file_out,".pdf",sep="")

  #file_comm <- paste(data,"community_opt-dir.txt",sep="_")
  #file_out <- paste(data,"opt-dir.png",sep="_")
  
  library("igraph");
 
  
  community <- read.table(file_comm,skip=0,header=FALSE)
  community <- as.matrix(community)
  community <- community + 1
  n.comm <- max(community)
  cat("number of communities",n.comm,"\n")
  
  if(0){
  #to sort communities
  comm_size <- rep(0,n.comm)
  for(i in 1:n.comm){
    comm_size[i] <- length(which(community==i))
  }
  order_comm <- order(comm_size,decreasing=T)
  community <- order_comm[community]
  cat("community indices order by size\n")
  }
  
  #whole graph
  el <- read.table(file_el)
  el <- el + 1
  el <- as.matrix(el)
  #print(length(index_meta))
  g <- graph.edgelist(el,directed=T)
  
  #graph <- read.graph(str,format="edgelist",directed=TRUE)
  
  A <- get.adjacency(g,sparse=FALSE)
  N <- length(A[1,]) #number of nodes
  print(N)
  #read layer index
  layer <- as.matrix(read.table(file_layer,header=FALSE)) + 1
  
  #create a layout and color nodes according to community_meta
  layout <- matrix(0,N,2)
  if(lyt == 1){ 
    n.layer <- max(layer)
    current_layer <- 1
    max_length <- length(which(layer==1))
    
    while(current_layer <= n.layer){
  
      index <- which(layer==current_layer)  
      if(max_length < length(index)) max_length <- length(index)
      order_index <- order(community[index],decreasing=F)
      #x around center 0 layer 
      x.top <- length(index)/2 - 0.5
      layout[index[order_index],1] <- max_length/length(index)*(-x.top:x.top)
      #layout[index[order_index],1] <- community[index[order_index]]
     
      #for(k in 1:n.comm){
      #  index_comm <- which(community[index[order_index]]==k)
      #  x.top <- length(index_comm)/2 - 0.5
      #  layout[index[order_index][index_comm],1] <- 20*k+ 1.3* max_length/length(index)*(-x.top:x.top)
      #}
      
      layout[index[order_index],2] <- 0.8*rep(current_layer,length(index))
      current_layer <- current_layer + 1
      
    }
    #layout <- layout.fruchterman.reingold(g)
    write.table(layout,"layout.txt",col.names=F,row.names=F)
  }else{
    layout <- as.matrix(read.table("layout.txt"))  
  }

  col <- rainbow(n.comm)
  #col <- c("black")
  
  for(k in 1:n.comm){
    comm_index <- which(community == k)
    V(g)[comm_index]$color <- col[k]
  }
  
  #E(g)$color <- "gray14"
  E(g)$color <- "black"
  el <- get.edgelist(g)
  for(k in 1:length(el[,1])){
    if(community[el[k,1]] != community[el[k,2]]){
      #E(g)[k]$color <- "grey"
      E(g)[k]$color <- "black"
    }
  }
  
  pdf(file_out,
      width     = 12,
      height    = 6.5,
      #units     = "in",
      #res       = 700,
      pointsize = 5
  )
  
  par(mfrow=c(1,1)#,las=1,yaxs="i",xaxs="i",cex=2.5
      #mar=c(4,4.5,4.5,1),
      #mgp=c(2.5,1,0)
  );
  
  #C6
  #layout[52,1] <- layout[52,1]+5
  #layout[54,1] <- layout[54,1]-5
  #layout[55,1] <- layout[55,1]+10
  #layout[56,1] <- layout[56,1]-10
  
  #C1
  #layout[19,1] <- layout[19,1] - 2.3
  layout[19,1] <- layout[19,1] - 3
  #layout[18,1] <- layout[4,1] - 0.3
  layout[18,1] <- layout[4,1] - 0.48
  #layout[17,1] <- layout[17,1] + 3.5
  layout[17,1] <- layout[17,1] + 3.3
  #layout[16,1] <- layout[16,1] - 0.25
  layout[16,1] <- layout[16,1] + 0.75
  
  #C2
  #layout[32,1] <- layout[32,1] - 2.5
  
  #layout[32,1] <- layout[20,1]
  #layout[31,1] <- layout[31,1] + 6.5
  #layout[26,1] <- layout[26,1] + 4
  #layout[29,1] <- layout[26,1] 
  #layout[28,1] <- layout[14,1] 
  #layout[25,1] <- layout[18,1] 
  #layout[27,1] <- layout[25,1]
  
  #layout[31,1] <- layout[31,1] + 2
  #layout[29,1] <- layout[29,1] + 2
  #layout[26,1] <- layout[26,1] + 2
  #layout[30,1] <- layout[30,1] - 2
  #layout[28,1] <- layout[28,1] - 2
  #layout[27,1] <- layout[27,1] - 2
  
  #layout <- lyt
  
  plot(g,
       vertex.size=13,
       vertex.frame.color=NA,
       edge.width=2.5,
       vertex.label="",
       #vertex.label=V(g)$index,
       #vertex.label=community,
       vertex.label.cex = 1,
       vertex.label.color=NA,
       layout = layout,
       #layout = layout.fruchterman.reingold,
       edge.arrow.size=3,
       #vertex$color <- V(g)$color
       asp=0
  )
  
  dev.off()
  
}