# plot a matrix with a nice color scheme

plotmat = function(mat, color = "green", title = NULL){
  mat = apply(mat, 2, rev)
  longmat = melt(mat)
  
  p = ggplot(longmat, aes(x = Var2, y = Var1)) + 
    geom_tile(aes(fill=value), colour="grey20") 
  
  if(color == "green") p = p + scale_fill_gradient2(low = "#3d52bf", high = "#33961b", mid = "white")
  if(color == "red") p = p + scale_fill_gradient2(low = "#191970", high = "#800000", mid = "white")
  
  p = p + theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.title = element_text(),
          plot.title = element_text(hjust = 0.5)) + 
    labs(fill = " ")
    if(! is.null(title)) p = p + ggtitle(title)
  p
}