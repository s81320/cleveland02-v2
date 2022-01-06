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
  ti
  ti$splitvarName %>% unique -> d0.info
  
  sn <- unique(ti[,c('nodeID','splitvarID','splitvarName')]) # the split nodes by ID
  sn
  el <- matrix(NA,2*nrow(ti),ncol=2)
  ct <-  1
  for(i in 1:nrow(ti)){
    if(ti[i,'leftChild'] %in% sn$nodeID ){
      #print(paste('add row' , i, ct))
      el[ct,] <- ti[i,c('nodeID','leftChild')] %>% unlist # add edge
      ct <- ct+1
    }
    if(ti[i,'rightChild'] %in% sn$nodeID ){
      el[ct,] <- ti[i,c('nodeID','rightChild')] %>% unlist # add edge
      ct <- ct+1
    }
  }
  
  el <- el[!is.na(el[,1]),]
  #el %>% print
  
  el <- el+1 # originally nodeID in ranger starts at 0 but we need positive numbers for the nodes
  sn$nodeID <-  sn$nodeID+1
  el # nodes need to be numbered 1,2,3.. no gaps allowed
  # splitvarName has to be tied to the numbering of nodes
  
  ist <- unique(as.vector(el))
  ist <-  ist[order(ist)]
  soll <- 1:length(ist)
  
  for(i in soll){
    el[el==ist[i]] <- soll[i]
  }
  
  sn$nodeID <- soll
  
  graph_from_edgelist(el) -> g1

  lo <-  layout_as_tree(g1,root=1)
  
  # g1 %>% plot(layout=lo, main=tri, mode='out') # simpler plot
  g1 %>% plot(layout=lo
            , vertex.label.cex=0.7
            , vertex.label=sn$splitvarName # needs to be in sync with edge list el
            #, vertex.color=NA
            #, vertex.frame.color=NA
            #, vertex.size=50
            , vertex.shape='none'
            #, vertex.label.dist=0.5 # dist from the center of the vertex
            , edge.color ='red'
            , edge.arrow.size=.4
            , rescale = T
          # , ylim=0.5*range(lo[,2])
          # ,xlim=0.5*range(lo[,1])
          , asp = 0
          , ...
          )
  
  return(d0.info)
}

rg <- ranger(Species~. , num.trees=3 , data=iris , max.depth = 4)

treeInfo(rg,2)

plotTree(rg, 1, main='tree 1')
plotTree(rg, 2, main='tree 2')
plotTree(rg, 3, main='tree 3')
createDMd0(rg$forest)

# d1, d2 kann man nicht visuell zu überprüfen
df <- iris[sample(1:nrow(iris),5,F),]
createDMd1(rg$forest, iris)
createDMd1(rg$forest, df)

createDMsb(rg$forest)

rt2fbt(rg$forest,1)

# full binary trees built from ranger trees
lapply(1:rg$forest$num.trees, function(i) rt2fbt(rg$forest,i) ) ->
  A
A

treeInfo(rg,1)

load('data/data_SupMat4.rda')

names(Cleve)
Cleve$CAD_fac <-  NULL
rg <-  ranger(CAD~. , num.trees=3,  data=Cleve)            

plotTree(rg,1, main='blubs')
