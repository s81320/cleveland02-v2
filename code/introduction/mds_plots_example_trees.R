# file code / introduction / mds_plots_example_trees.R
# plotting example trees in context of forest
# figure fig:intro:mds:example
# 10.5.2022

load('data/nursery02/nursery02_01.rda')
DM <- doc[[1]]$DM
par(mar=c(2,2,1,1)+0.5)
#par(mar=rep(1.1,4))
for(metric in c('d0','d1','d2','sb')){
  dm <-  DM[[metric]]
  attr(dm,'metric') <- metric

  mp <-  cmdscale(d=dm, k=2) # for plotting we don't need eigenvalues, just the points
  central_tree <- which.min(lapply(1:nrow(dm), function(i) sum(dm[i,-i])))
  plot(mp, type='n' , axes=F)
  if(max(mp[,1]<1)){
    axis(1, at=c(0,round(range(mp[,1]),1)), cex.axis=1.5)
    axis(2, at =c(0,round(range(mp[,2]),1)) , cex.axis=1.5)
  }else{
    axis(1, at=c(0,round(range(mp[,1]))), cex.axis=1.5)
    axis(2, at =c(0,round(range(mp[,2]))) , cex.axis=1.5)
  }
  text(mp[1:3,] , col=1 , labels = 1:3 , cex=1.5)
  points(x=mp[central_tree,1] , y=mp[central_tree,2] , pch=0 , col='red', cex=1.5)
  legend('topright', legend=metric, cex=1.5)
}
