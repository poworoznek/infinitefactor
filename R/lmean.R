# mean of a list of elements for which the + operator
# and scalar division is well defined

lmean = function(list){
  l = length(list)
  o = list[[1]]/l
  for(i in 2:length(list)){
    o = o + list[[i]]/l
  }
  return(o)
}




