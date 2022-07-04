# introduction , methods , trees and random forests

# what is a single tree in a random forest
# based on first 3 trees in first simulation

rm(list=ls())

file_name_script <-  'code/introduction/tree-03.R'
list.files(path = 'code/introduction' , recursive = T)
assertthat::assert_that(file_name_script %in% list.files(recursive = T))

library(ranger)
library(xtable)
# library(caret)
library(dplyr)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/plotTree.R') # plotTree1 , plotTreeFb

load('data/nursery02/nursery02_01.rda')

rg <- doc[[1]]$rg
data.train <- doc[[1]]$`bootstrapped training data`
rownames(data.train)

f5 <-  function(nodeID){
  #' based in the ranger package
  #' given a split node by nodeID
  #' given treeInfo and d.split in parent environment
  #' 
  #' return a list with the row numbers of observations in each node.
  #' tested by comparing the counts to predictions of type='terminalNodes'

  svN <- ti[ti$nodeID==nodeID, 'splitvarName'] # svN : splitvarName
  #print(svN)
  #print(d.split[[nodeID]])
  sVal <- ti[ti$nodeID==nodeID, 'splitval'] # sVal : splitValue
  #print(sVal)
  
  if(nodeID==0){
    idcs <- 1:nrow(data.train)
  }else{
    idcs <- d.split[[nodeID]]
  }
  
  #print(idcs)
  
  idcs %>% 
    (function(x)
      # as.numeric is needed only for factors - but does not hurt when applying to already numeric featurs
      as.numeric(data.train[x,svN]) <= sVal
    ) %>% # mask for Cleveland rows of input (not original Cleveland rows)
    
    which %>% 
    idcs[.] %>% # transform back to original Cleveland rows
    return
}

pred.tn <- predict(rg, type='terminalNodes' , data=data.train)

tri <-  1
ti <- treeInfo(rg,tri) ; ti
plotTree1(rg,tri)
#plotTreeFb(rg,tri)

# which observations are in individual trees training data? This will affect their empirical distributions at terminal nodes!
# we need keep.inbag==TRUE and we do not have that in forests from Cleveland simulations in nursery or nursery02
rg$inbag.counts %>% lapply(function(x) which(x==1)) -> inbag.obs


# terminal nodes
ti[ti$terminal,'nodeID']

d.split <- NULL
for(nodeID in ti[!ti$terminal, 'nodeID']){ # going through the split nodes / non-terminal nodes

  # initial split , root split
  if(nodeID==0){
    d.split[[1]] <- f5(0)
    d.split[[2]] <- base::setdiff(1:nrow(data.train),d.split[[1]])
  }else{
    lc.nodeID <- ti[ti$nodeID==nodeID,'leftChild']
    d.split[[lc.nodeID]] <- f5(nodeID)
    rc.nodeID <- ti[ti$nodeID==nodeID,'rightChild']
    d.split[[rc.nodeID]] <- base::setdiff(d.split[[nodeID]],d.split[[lc.nodeID]])
  }
}

# terminal nodes 
ti[ti$terminal,'nodeID']
# observations (of forest's training data, initial training data) in each terminal node
d.split[ti[ti$terminal,'nodeID']] %>% lengths
# check against ranger predictions for terminal nodes
pred.tn$predictions[,tri] %>% table 

# terminal node
tN <- ti[ti$terminal,'nodeID'][1]
# forest's training data observations in terminal node
which(pred.tn$predictions[,tri]==tN)
inbag.obs[[tri]] ; inbag.obs[[tri]] %>% length
# inbag obs in terminal node
inbag.obs.tN <- base::intersect(which(pred.tn$predictions[,tri]==tN) , inbag.obs[[tri]])

# inbag obs in terminal node (individual tree's training data in terminal)
d.split[[tN]] 
d.split.inbag.obs <- base::intersect(d.split[[tN]] , inbag.obs[[tri]])

# remove: includes observations not in trees training data
# empirical distribution in 1st terminal node
# data.train[d.split[[tN]],'CAD'] %>% table # abs. numbers
# data.train[d.split[[tN]],'CAD'] %>% table %>% prop.table # proportions

# empirical distribution in terminal node tN
data.train[d.split.inbag.obs,'CAD'] %>% table # abs. numbers
data.train[d.split.inbag.obs,'CAD'] %>% table %>% prop.table -> emp.dist.inbag ; emp.dist.inbag # proportions

# check : empirical distribution of training data generates predictions
pred.c <- predict(rg, predict.all = T ,type = 'response', data=data.train) # 3 dim : obs x (N/Y) x tree
# predictions by fixed tree tri for observations falling to terminal node tN
pred.c$predictions[d.split.inbag.obs,'Yes',tri] # all observations equal, makes sense!

# TEST #
# test how predictions for the individual tree are generated
assertthat::assert_that(all(pred.c$predictions[d.split.inbag.obs,'Yes',tri]==emp.dist.inbag['Yes'])) # should be TRUE

# test majority vote in forest #
# move to forests training data , move away from trees training data

# fix three observations : 1:3
# prediction mean over all trees
predict(rg, type='response', data=data.train[1:3,])$predictions[,'Yes'] 
# single tree predictions , probability for positive class 'Yes'
pred.c$predictions[1:3,'Yes',]
pred.c$predictions[1:3,'Yes',] %>% apply(1,mean)

assertthat::assert_that(all(predict(rg, type='response', data=data.train)$predictions[,'Yes']-(pred.c$predictions[,'Yes',] %>% apply(1,mean))<0.5/nrow(data.train)), msg='too bad!')
#####
# Check Gini index on the terminal nodes and forest's training data
doc1 <- matrix(0, nrow=sum(ti$terminal), ncol=4)
ct <- 1
for(i in ti[ti$terminal, 'nodeID']){
  print(i)
  if(i>0){
  (d.split[[i]] %>% data.train[.,'CAD'] %>% as.numeric %>% mean -1) -> p 
    #d.split[[i]] %>% print
    data.train[d.split[[i]],'CAD'] %>% as.numeric %>% print
    2*p*(1-p) %>% print
    doc1[ct,] <-  c(i, length(d.split[[i]]), 1-p, 2*p*(1-p) )
    ct <- ct+1
  }
}
doc1 <- as.data.frame(doc1)
colnames(doc1) <-  c('nodeID','obs.in.node','p','gini')
View(doc1)
sum(doc1$obs.in.node) # check : should be the number of observations in training
assertthat::assert_that(nrow(data.train)==sum(doc1$obs.in.node))

hist(doc1$gini, breaks=10, main='a fixed tree\'s Gini index \nfor forests training data at terminal nodes')
mean(doc1$gini)

# pe trees document the empirical distributions in their terminal nodes in
# rg $forest$terminal.class.counts

# mean of Gini index over the terminal nodes
lapply(rg$forest$terminal.class.counts[[1]] # list elements empty or an empirical probability (over CAD in terminal node)
       , function(x) {if(length(x)>0) 2*x[2]*(1-x[2])} # if list element is a probability , calculate the Gini Index
       ) %>% unlist %>% mean

# check for ensemble magic
if('probability' %in% names(rg$call))
  if(rg$call$probability=='T'){
    print(c(
      calcLogloss(forest=subforest(rg$forest,1), df = data.train)
    , calcLogloss(forest=subforest(rg$forest,2), df = data.train)
    , calcLogloss(forest=subforest(rg$forest,3), df = data.train)
    , calcLogloss(forest=subforest(rg$forest,1:3), df = data.train)
    , calcLogloss(forest=rg$forest, df = data.train)
    ))
  }

# evolution of Gini index
# needs keep.inbag=T during the forest's creation
# initial Gini index
rg$inbag.counts %>% lapply(function(x) which(x==1)) -> inbag.obs # rg$numbtrees lists , 1 list per tree 

# inbag obs in terminal node (individual tree's training data in terminal)
d.split[[tN]] 
d.split.inbag.obs <- base::intersect(d.split[[tN]] , inbag.obs[[tri]])

ti[ti$terminal, 'nodeID'] # terminal nodes

data.train[inbag.obs[[tri]],'CAD'] %>% table
data.train[inbag.obs[[tri]],'CAD'] %>% table %>% prop.table -> emp.dist ; emp.dist

# inspect single node (nodeID) of a tree (tri, variable set above)
(function(nodeID) { # 
  obs <- base::intersect(d.split[[nodeID]], inbag.obs[[tri]])
  data.train[obs,'CAD'] %>% table %>% prop.table -> emp.dist
  emp.dist %>% (function(x) unname(2*x[1]*x[2])) -> gi # Gini index
  terminal <- ti[ti$nodeID==nodeID, 'terminal']
  return.split <- if(terminal){NULL}else{ti[ti$nodeID == nodeID, c('splitvarName','splitval')]}
  return(list('nodeID'=nodeID
              , 'tree'=tri
              , 'terminal'= terminal
              , 'split' = return.split
              , 'emp.dist'=emp.dist
              , 'Gini'=gi 
              , 'obs'=obs 
              , 'num.obs'=length(obs)))
  })(3)
