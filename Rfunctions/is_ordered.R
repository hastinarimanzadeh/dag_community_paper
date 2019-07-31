
is_ordered <- function(el){
  print(head(el))
  
  layer <- max(el[,1],el[,2]):min(el[,1],el[,2]);
  
  #check wether each entry in el satisfies layer[el[,1]]>layer[el[,2]]
  
  comp <- c(layer[el[,1]+1] > layer[el[,2]+1])
  
  checkF <- any(comp == F)
  checkT <- any(comp == T)
  
  foo <- which(comp == F)
  goo <- which(comp == T)
  if(length(foo) < length(goo)) print(foo)
  if(length(foo) >= length(goo)) print(goo)

  if((checkF == T) && (checkT == T)) return(F);
  return(T)
  
}