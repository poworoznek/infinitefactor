# mean of a p times k matrix with samples stored as a p times k times m array

amean = function(ar){
  return(apply(ar, c(1, 2), mean))
}
