# plotting trees with igraph
# plotTree1 : creating edge list from ranger tree , plot with layout_as_tree, plot L or R on edges
# plotTreeFb : plot a full binary tree and make the nodes missing in tree vanish. Typical layout, no need to indicate L or R.

library(dplyr)
library(ranger)
library(igraph)

plotTree1 <- function(rg,tri,...){
  #' plot a single ranger tree, the split nodes only
  #' 
  #' throws away the terminal nodes and info on them
  #' additional arguments can be passed to the plot. like a main (title)

  ti <- treeInfo(rg,tri)
  
  # remove terminal nodes
  ti <- ti[!ti$terminal,]
  
  ti$splitvarName %>% unique -> d0.info
  
  # the split nodes by ID : s_plit n_odes
  sn <- ti[,c('nodeID','splitvarID','splitvarName')]
  # sn %>% print
  # el : e_dge l_ist
  # final number of rows (=edges) unknown , maximal number of edges is number of nodes -1
  # columns : nodeID of split node , nodeID of child node , info if it is a left or right child
  el <- matrix(NA,2*nrow(ti)-1,ncol=3) 
  
  ct <-  1
  for(i in 1:nrow(ti)){
    if(ti[i,'leftChild'] %in% sn$nodeID ){
      #print(paste('add row' , i, ct))
      el[ct,1:2] <- ti[i,c('nodeID','leftChild')] %>% unlist # add edge
      ct <- ct+1
    }
    if(ti[i,'rightChild'] %in% sn$nodeID ){
      el[ct,1:2] <- ti[i,c('nodeID','rightChild')] %>% unlist # add edge
      ct <- ct+1
    }
  }
  
  el <- el[!is.na(el[,1]),] # remove rows with NA nodeID -> rows with nodeID NA whenever binary tree is not full binary tree
  #el %>% print
  
  # el%>% print # nodes need to be numbered 1,2,3.. no 0 and no gaps allowed
  # splitvarName has to be tied to the numbering of nodes
  
  ## current node numbering : in sn$nodeID
  #for(i in 1:length(sn$nodeID)){
  #  el[el==sn[i,'nodeID']] <- i
  #  sn[i,'newNodeID'] <- i
  #}
  
  done <-  matrix(0,nrow=nrow(el),ncol=2) # same layout as edgelist el
  ct <- 1 # new IDs will be sequential, starting at 1
  for(i in sn$nodeID){
    act <-  (!done & el[,1:2]==i) # act on positions in el where nothing has been done before and entries = i
    el[act] <- ct # change the nodeID from i to ct
    done[act] <-  1 # document the positions where the nodeID has been changed
    ct <- ct+1
    }
  el
  el[,3] <- sn[-1,'nodeID'] %%2 # remove 1st entry , it is for the root , other entries are for outgoing edges
  
  #sn$nodeID <- 1:length(nrow(sn))
  
  graph_from_edgelist(el[,1:2]) -> g1

  lo <-  layout_as_tree(g1,root=1)
  
  #print(sn)
  #print(el) # el has a row less than sn (nodes vs edges)
  # left or right child is coded in sn : nodeID even = left child , nodeID odd = rightChild
  
  # g1 %>% plot(layout=lo, main=tri, mode='out') # simpler plot
  g1 %>% plot(layout=lo
            #, vertex.label.cex=0.7
            , vertex.label=sn$splitvarName # needs to be in sync with edge list el
            , edge.label= ifelse(el[,3] == 0,'R','L') 
          #  , edge.label.dist = 5 # not working , would like to put the label above the edge
           # , edge.label.cex=0.7
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
  
  return('d0.info'=d0.info)
}

plotTreeFb <-  function(rg,tri,...){
  #' plot a tree on the layout of a full binary tree
  
  # turn into a full binary tree
  fbt <- rt2fbt(rg$forest,tri)
  #fbt
  N <- length(fbt) # nodes in the full binary tree
  # E <-  N-1 # number of edges in full binary tree
  D <- log(N+1 , base=2)-1 # D for depth of tree , only root D equals 0
  #el <- matrix(0,nrow=N-1,ncol=2)
  el <- matrix(0,nrow=N-1,ncol=3)
  #el
  ct <- 1
  for(k in 1:(2^D-1)){
    el[ct,] <- c(k,2*k,fbt[2*k]>0) # draw edge only if the child node is a split node (value !=0)
    el[ct+1,] <- c(k,2*k+1,fbt[2*k+1]>0) # same
    ct <- ct+2
  }
  el
  graph_from_edgelist(el[,1:2]) -> g1
  lo <-  layout_as_tree(g1,root=1)
  plot(g1
       , layout=lo
       , vertex.label=fbt # needs to be in sync with edge list el
       , vertex.color='white'
       , vertex.frame.color= ifelse(fbt==0,'white','black')
       , vertex.label.color= ifelse(fbt==0,'white','black')
       , edge.color= ifelse(el[,3]==1,'black','white')
       , edge.arrow.size=.4
       , rescale=T
       , ...)
  d0.info <- unique(fbt)
  d0.info <- d0.info[order(d0.info)][-1]
  return(list('d0.info'=d0.info,'sb.info'=fbt))
}

