# introduction , methods , trees and random forests

# what is a single tree in a random forest

rm(list=ls())

file_name_script <-  'code/introduction/tree.R'
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

data.train <- Cleve[,1:11] 

set.seed(1)
rg <- ranger(CAD~.
             , data = data.train
             , num.trees = 500
             , replace = F # neu nach Absprache mit AZ
             , mtry= 3 # default : 3
             , importance = 'impurity'
             # the following (in this script) needs classical classification , so the next line should be commented out
             , probability = T # this makes it a random forest of type 'Probability Estimation'
             , min.node.size = 13 # optimized in A.Z. paper: BiomJ 56, 2014, 
            # , max.depth=3
             #, respect.unordered.factors = 'partition' # this option turns splitvals into characters # factor selection 1,3 cannot be done with < >= but needs to split the string 1,3 into a lost of integers and compare with %in%
             , keep.inbag = TRUE
)



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
    idcs <-  1:nrow(data.train)
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

pred.tn <- predict(rg, predict.all = T ,type = 'terminalNodes', data=data.train)

tri <-  3
ti <- treeInfo(rg,tri) ; ti
plotTree1(rg,tri)
#plotTreeFb(rg,tri)
pred.tn$predictions[,tri] %>% table # observations per terminal node

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

pred.c <- predict(rg, predict.all = T ,type = 'response', data=data.train) # 2 dim : obs x tree

doc1 <-  matrix(0, nrow=sum(ti$terminal), ncol=4)
ct <- 1
for(i in ti[ti$terminal, 'nodeID']){
  print(i)
  if(i>0){
  (d.split[[i]] %>% data.train[.,'CAD'] %>% as.numeric %>% mean -1) -> p 
    #d.split[[i]] %>% print
    data.train[d.split[[i]],'CAD'] %>% as.numeric %>% print
    2*p*(1-p) %>% print
    doc1[ct,] <-  c(i, length(d.split[[i]]), p, 2*p*(1-p) )
    ct <- ct+1
  }
}
doc1 <- as.data.frame(doc1)
colnames(doc1) <-  c('nodeID','obs.in.node','p','gini')
View(doc1)
sum(doc1$obs.in.node) # check : should be the number of observations in training
assertthat::assert_that(nrow(data.train)==sum(doc1$obs.in.node))

hist(doc1$gini, breaks=10)
mean(doc1$gini)

