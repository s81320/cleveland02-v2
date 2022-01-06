library(dplyr)
library(ranger)
library(igraph)



plotTree <- function(rg,tri,...){
  #' plot a single ranger tree, the split nodes only
  #' 
  #' throws away the terminal nodes and info on them
  #' additional arguments can be passed to the plot. like a main (title)

  ti <- treeInfo(rg,tri)
  
  # remove terminal nodes
  ti <- ti[!ti$terminal,]

  names(ti) <-  NULL

  el <- rbind(as.matrix(ti[,c(1,2)]) , as.matrix(ti[,c(1,3)]))
  el <-  el+1 # originally nodeID in ranger starts at 0 but we need positive numbers for the nodes

  graph_from_edgelist(el) -> g1

  lo <-  layout_as_tree(g1,root=1)
  
  g1 %>% plot(layout=lo
              , vertex.label.cex=0.7
            , vertex.label=ti[,5] # ti$splitvarName
            , vertex.color=NA
            , vertex.frame.color=NA
            #, vertex.size=50
            #, vertex.shape='none'
            #, vertex.label.dist=0.5 # dist from the center of the vertex
            , edge.color ='red'
          , rescale = T
          # , ylim=0.5*range(lo[,2])
          # ,xlim=0.5*range(lo[,1])
          , asp = 0
          , ...
          )
}



rg <- ranger(Species~. , num.trees=3 , data=iris)
plotTree(rg, 1, main='tree 1')
plotTree(rg, 2, main='tree 2')
plotTree(rg, 3, main='tree 3')
createDMd0(rg$forest)

# d1, d2 kann man nicht visuell zu überprüfen
df <- iris[sample(1:nrow(iris),5,F),]
createDMd1(rg$forest, iris)
createDMd1(rg$forest, df)

createDMsb(rg$forest)

load('data/data_SupMat4.rda')

names(Cleve)
Cleve$CAD_fac <-  NULL
rg <-  ranger(CAD~. , num.trees=3,  data=Cleve)            

plotTree(rg,1, main='blubs')
