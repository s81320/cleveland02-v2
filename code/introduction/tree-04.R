# introduction , methods , trees and random forests

# what is a single tree in a random forest
# based on first 3 trees in first simulation

rm(list=ls())

file_name_script <-  'code/introduction/tree-04.R'
list.files(path = 'code/introduction' , recursive = T)
assertthat::assert_that(file_name_script %in% list.files(recursive = T))

library(ranger)
library(xtable)
library(caret)
library(dplyr)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest , tNodes
source('code/source/plotTree.R') # plotTree1 , plotTreeFb
source('code/source/sim-prep-v2.R') # loads functions create_prob , new_bootstrap


{ # block copied from code/nursery/01.R
  # may not be changed! Creates exactly the same (first 3) trees as in first simulation
  # but with added keep.inbag=T
  Cleve_enriched <- create_prob(Cleve)
  nBs <- 50
  
  seed <- 13
  set.seed(seed)
  cr<-createResample(Cleve_enriched$CAD
                   , times = nBs
                   , list = T)

  ct <- 1
  i <- 1
  data.train <- new_bootstrap(Cleve_enriched , cr[[i]])[,-12] # no probabilities
  
  # grow small forest (size nT)
  rg <- ranger(CAD~.
               , data = data.train 
               , num.trees = 500
               , replace = F 
               , mtry= 3 
               , importance = 'impurity'
               , probability = T 
               , min.node.size = 13 
               , keep.inbag = T # NEW
    )
}

tri <-  1
ti <- treeInfo(rg,tri) ; ti
plotTree1(rg,tri)

# too large : data for all 500 trees , we only need 1 at a time
# rg$inbag.counts %>% lapply(function(x) which(x==1)) -> inbag.obs # for each tree the inbag obs of the simulated Cleveland training data

rg$inbag.counts[[tri]] %>% (function(x) which(x==1)) -> inbag.obs.tri # for tree tri in forest rg$forest
# rownumbers of data.train , the training data of the ranger object rg
rg$call$data
get(x = rg$call$data %>% as.character)

ID2pos <-  function(x) x+1 # IDs start at 0, positions at 1

f6 <-  function(data, nodeID){
  #' based in the ranger package
  #' given a split node by nodeID and some data
  #' (needs d.split in parent environment, builds d.split)
  #' 
  #' return a list with the row numbers of observations in each node.
  #' tested by comparing the counts to predictions of type='terminalNodes'
  
  # old , relying on treeInfo
  # svN <- ti[ti$nodeID==nodeID, 'splitvarName'] # svN : splitvarName
  nodeID %>% 
    ID2pos %>% 
    rg$forest$split.varIDs[[tri]][.] %>% # returns a nodeID
    ID2pos %>% 
    rg$forest$independent.variable.names[.] -> svN
 
  nodeID %>%
    ID2pos %>%
    rg$forest$split.values[[tri]][.] -> sVal
  # old , based on treeInfo
  #sVal <- ti[ti$nodeID==nodeID, 'splitval'] # sVal : splitValue
  
  if(nodeID==0){
    # idcs <- 1:nrow(data) # old , wrong when working with individual tree's training data , might be OK on a full data set
    idcs <- inbag.obs.tri # refers to rows of rg object training data
  }else{
    nodeID %>%
      ID2pos %>%
      d.split[[.]] -> idcs
  }
  
  #print(idcs)
  
  idcs %>% 
    (function(x)
      # as.numeric is needed only for factors - but does not hurt when applying to already numeric features
      as.numeric(data[x,svN]) <= sVal
    ) %>% # mask for Cleveland rows of input (not original Cleveland rows)
    
    which %>% 
    idcs[.] %>% # transform back to original Cleveland rows
    return
}

inbag.obs.tri # indices relate to data.train

f6(data=data.train ,nodeID = 0) # works only for nodeID=0 , else we need d.split which is not yet created

ti <-  treeInfo(rg, tri)

d.split <- list()
d.split[[1]] <- inbag.obs.tri
all.nodes <- unique(rg$forest$child.nodeIDs[[tri]] %>% unlist)
#all.nodes <- base::setdiff(all.nodes,0)
all.nodes <- all.nodes[order(all.nodes)]
for(nodeID in all.nodes){
  #nodeID %>% paste('nodeID') %>% print
  #ti[ti$nodeID == nodeID , 'terminal'] %>% paste('terminal') %>% print
  
  nodeID %>% 
    ID2pos %>% 
    rg$forest$child.nodeIDs[[tri]][[1]][.] %>% 
    (function(x) x==0) -> terminal

  if(!terminal){
    
    nodeID %>%
      ID2pos %>%
      rg$forest$child.nodeIDs[[tri]][[1]][.] %>%
      ID2pos -> pos1
    
    #pos1 %>% paste('pos for left child, pos1')  %>% print
    
    nodeID %>%
      ID2pos %>%
      rg$forest$child.nodeIDs[[tri]][[2]][.] %>%
      ID2pos -> pos2
    
    #pos2 %>% paste('pos for right child, pos2')  %>% print
    
    d.split[[pos1]] <- f6(data=data.train, nodeID=nodeID)
    #d.split[[pos1]] %>% length %>% print
    d.split[[pos2]] <- base::setdiff(d.split[[ID2pos(nodeID)]], d.split[[pos1]])
    #d.split[[pos2]] %>% length %>% print
  }
}

# nodeIDs appear multiple times, in each node along its way to the terminal node
d.split %>% lengths
d.split %>% unlist %>% length


# compare terminal class counts of ranger (in terminal nodes) 
# with empirical distributions of observations in d.split
# ranger : 
rg$forest$terminal.class.counts[[tri]][1+ti[ti$terminal, 'nodeID']] %>% # 1+... nodeID makes it correct : ID2pos!
  unlist %>% 
  matrix(ncol=2, byrow=T) %>%
  .[,rg$forest$class.values[2] ] -> a1 ; a1

# silke :
f7 <- function(tri,nodeID){
  #' for a node in a tree get the positive class (emp.) probability 
  #' at the node (split or terminal)
  #' needs d.split and data.train in parent environment

  data.train[d.split[[ID2pos(nodeID)]],'CAD'] %>% table %>% prop.table %>% .['Yes']
}

nodeID <- ti[ti$terminal,'nodeID']
df1 <- cbind(nodeID 
      , size=nodeID %>% ID2pos %>% d.split[.] %>% lengths
      , p=lapply(nodeID , FUN=f7 , tri=1) %>% unlist %>% unname
) %>% data.frame

df1$Gini <-  with(df1,2*p*(1-p))

df1$p 

# check
assertthat::assert_that(df1$size %>% sum==191)

df1[,c('nodeID', 'size', 'p')] %>% t %>% xtable -> xt
digits(xt) <-  c(0,0,0,2)
xt
################################################################################
#### TEST #### should be TRUE #### TEST should be TRUE #### TEST should be TRUE 
################################################################################
############## a1 is the probability for the positive class in the terminal nodes
############## calculated by ranger terminal.class.count #######################
############## df1$p is my calculation #########################################
################################################################################
a1 == df1$p
all(a1 == df1$p)
################################################################################

# voting of forest in probability estimation

pred <- predict(rg, data = data.train[1,] , predict.all=T )
#pred$predictions[1,,1:3] 
pred$predictions[1,'Yes',1:3]
pred$predictions[1,'Yes',1:3] %>% mean
#predict(rg$forest, data = data.train[1,])$predictions
predict(subforest(rg$forest,1:3), data = data.train[1,])$predictions[1,'Yes']

# this 'proves' that a tree's predictions are made from the empirical distribution at the terminal node
px <- function(obs1, trindcs){
  #' calculate predictions in ranger based on empirical distributions at terminal nodes
  #'
  #' uses the ranger object rg in parent environment
  
  pred.tn <- predict(rg, data = obs1 , type = 'terminalNodes')
  tn1 <- pred.tn$predictions[1,trindcs] # nodeIDs - no probabilities!
  # we need the empirical distributions of training data CAD at the terminal nodes

  print('observation falls to these terminal nodes (one nodeID per tree):')
  print(tn1)
  
  # rg$forest$terminal.class.counts[[1]][tn1[1]+1]

  print('function returns the empirical probabilities at these terminal nodes (one emp. dist. per tree):')
  ((function(x) x %>% 
    tn1[.] %>%
    ID2pos %>%
    rg$forest$terminal.class.counts[[x]][.]
  ) %>%
      Vectorize)(trindcs)
}

px(data.train[1,],1:3)
px(Cleve[3,],1:3)

#### careful : do not do this (type='terminalNodes' is overruled by predict.all, but the predictions change)
predict(rg, data = data.train[1,] , type = 'terminalNodes' , predict.all=T)$predictions[1,,1:3]
# because you do not know which probability is for Yes / the positive class or for No 
################################################################################

# prune trees to height 2


f8 <- function(tri,nodeID){
  #' for a node in a tree get the positive class (emp.) probability 
  #' at the node (split or terminal)
  #' needs d.split and data.train in parent environment
  
  print(nodeID)
  
  # checks if nodeID is terminal in forest sf - not in ranger object rg
  nodeID %>% 
    ID2pos %>% 
    sf$child.nodeIDs[[tri]][[1]][.] %>% 
    (function(x) x==0) -> terminal
  
  print(terminal)
  print(d.split[[ID2pos(nodeID)]])
  
  if(terminal){
    data.train[d.split[[ID2pos(nodeID)]],'CAD'] %>% table %>% prop.table %>% as.vector
  }else{
    return( vector("numeric", 0))
  }
}

sf <- subforest(rg$forest,1:3)
tri <- 2
terminal.nodes <- tNodes(forest = rg$forest , tri=tri)

ti <- treeInfo(rg,tri)
ti
print('Select a split node:')
print(base::setdiff(1:(terminal.nodes%>% (function(x) x[length(x)])) , terminal.nodes)) # the largest nodeID will be fo a terminal node
# terminal nodes are already in ascending order
last.nodeID.as.split <- 6 # not generally tested. # for trees 1,2,3 put 6,5,6 as last.nodeID.split

if(last.nodeID.as.split %in% terminal.nodes){
  simpleError(paste('this node is not a split node, it is a terminal node in the full tree: ', last.nodeID.as.split))
}else{
  last.nodeID.as.split %>%
    ID2pos %>%
    sf$child.nodeIDs[[tri]][[2]][.] -> last.nodeID
}

last.node.pos <-  ID2pos(last.nodeID)
# cutoff at last tree 
sf$child.nodeIDs[[tri]][[1]] <- sf$child.nodeIDs[[tri]][[1]][1:last.node.pos]
sf$child.nodeIDs[[tri]][[2]] <- sf$child.nodeIDs[[tri]][[2]][1:last.node.pos]
# new terminal nodes need to have children as 0
# terminal node = has no child
(sf$child.nodeIDs[[tri]][[1]] > last.nodeID ) -> mask
sf$child.nodeIDs[[tri]][[1]][mask] <- 0
(sf$child.nodeIDs[[tri]][[2]] > last.nodeID ) -> mask
sf$child.nodeIDs[[tri]][[2]][mask] <- 0


sf$split.varIDs[[tri]] <- sf$split.varIDs[[tri]][1:last.node.pos]
sf$split.values[[tri]] <- sf$split.values[[tri]][1:last.node.pos]

# replace class.counts
sf$terminal.class.counts[[1]]

#lapply(0:last.nodeID , FUN=f8 , tri=tri)
sf$terminal.class.counts[[tri]] <- lapply(0:last.nodeID , FUN=f8 , tri=tri)

forestHull(sf) -> sfh

sfh$num.trees
sfh$treetype <- rg$treetype # 'Probability estimation'

treeInfo(object = forestHull(sf), tree = tri)

plotTree1(forestHull(sf),tri)
plotTree1(forestHull(sf),3)

sfh$importance.mode <- rg$importance.mode
predict(sfh,data = data.train[inbag.obs.tri,], type='terminalNodes')$prediction
