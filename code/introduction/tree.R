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

data.train <- Cleve[1:50,1:11] 

set.seed(1)
rg <- ranger(CAD~.
             , data = data.train
             , num.trees = 3
             , replace = F # neu nach Absprache mit AZ
             , mtry= 3 # default : 3
             , importance = 'impurity'
             #, probability = T # this makes it a random forest of type 'Probability Estimation'
            # , min.node.size = 13 # optimized in A.Z. paper: BiomJ 56, 2014, 
             , max.depth=3
             #, respect.unordered.factors = 'partition' # this option turns splitvals into characters # factor selection 1,3 cannot be done with < >= but needs to split the string 1,3 into a lost of integers and compare with %in%
             , keep.inbag = TRUE
)

pred <- predict(rg, predict.all = T ,type = 'terminalNodes', data=data.train)

tri <-  2
ti <- treeInfo(rg,tri) ; ti
plotTree1(rg,tri)
pred$predictions[,tri] %>% table # observations per terminal node

# total data
data.train %>% nrow
data.train$CAD %>% table
data.train$CAD %>% table %>% prop.table

# purity , Gini index
# ModelMetrics::gini(data.train[d.split[[1]],'CAD'], as.factor(rep('Yes', length(d.split[[1]]))))

# tree data structure 
#tri <- NULL
#tri.level <- NULL
d.split <-  NULL
# root split : rownumbers of observations in each group after first split 

# as.numeric is required when we have respect.unordered.factors=T
# d.split[[1]] <- which(Cleve[,ti$splitvarName[1]] < as.numeric(ti$splitval[1]))

# when handling factor variables for respect.unordered.factors=T
# splitval[4] could be a character, to be interpreted as a list of values like "2,3".
# needs to be transformed to "2,3" -> c(2,3) with strsplit(',') %>% unlist
# Cleve[,ti$splitvarName[4]] %>% as.integer %in% (ti$splitval[4] %>% strsplit(',') %>% unlist %>% as.integer)

d.split[[1]] <- which(data.train[,ti$splitvarName[1]] <= ti$splitval[1])
f3(1:nrow(data.train),1)
d.split[[2]] <- base::setdiff(1:nrow(data.train), d.split[[1]])

# test that conditions are met in each node
data.train[d.split[[1]], ti$splitvarName[1]] %>% max
ti$splitval[[1]]
data.train[d.split[[2]], ti$splitvarName[1]] %>% min

# observations at each node
data.train[d.split[[1]],]  %>% nrow
data.train[d.split[[2]],]  %>% nrow
# same as
# Cleve[Cleve[,ti$splitvarID[1]+1] >= ti$splitval[1],] %>% nrow

# distribution of classes (No , Yes) in each node
data.train[d.split[[1]],'CAD'] %>% table %>% prop.table
data.train[d.split[[2]],'CAD'] %>% table %>% prop.table
# same as
# Cleve[Cleve[,ti$splitvarID[1]+1] >= ti$splitval[1],'CAD'] %>% table %>% prop.table

# first level split
# original , chaining conditions
# which(Cleve[,ti$splitvarName[1]] < ti$splitval[1] 
#      & Cleve[,ti$splitvarName[2]] < ti$splitval[2]) 

# recursive
d.split[[3]] <- which(data.train[d.split[[1]],ti$splitvarName[2]] <= ti$splitval[2]) %>%d.split[[1]][.]
d.split[[4]] <- base::setdiff(d.split[[1]], d.split[[3]])

data.train[d.split[[3]],] %>% nrow
data.train[d.split[[4]],] %>% nrow
# same as 
# Cleve[Cleve[,ti$splitvarID[1]+1] < ti$splitval[1] 
#      & Cleve[,ti$splitvarID[2]+1] >= ti$splitval[2],] %>% nrow

data.train[d.split[[3]],'CAD'] %>% table %>% prop.table
data.train[d.split[[4]],'CAD'] %>% table %>% prop.table

data.train[d.split[[3]], ti$splitvarName[2]] %>% max
ti$splitval[[2]]
data.train[d.split[[4]], ti$splitvarName[2]] %>% min

# first level , 2nd node on this level
d.split[[5]] <-  which(data.train[d.split[[2]],ti$splitvarName[3]] <= ti$splitval[3]) %>%d.split[[2]][.]
d.split[[6]] <- base::setdiff(d.split[[2]], d.split[[5]])

data.train[d.split[[5]],] %>% nrow
data.train[d.split[[6]],] %>% nrow

data.train[d.split[[5]],'CAD'] %>% table %>% prop.table
data.train[d.split[[6]],'CAD'] %>% table %>% prop.table

# test that conditions are met in each node
data.train[d.split[[5]], ti$splitvarName[3]] %>% max
ti$splitval[[3]]
data.train[d.split[[6]], ti$splitvarName[3]] %>% min


f3 <- function(row.idx , nodeID){
  #' return what is dropped to the left child
  #' needs the row numbers of observations in the splitting node and the nodeID of the splitting node
    
      row.idx %>% 
    
        (function(x)
          # as.numeric is needed only for factors - but does not hurt when applying to already numeric featurs
         as.numeric(data.train[x,ti$splitvarName[nodeID]]) <= ti$splitval[nodeID]
          ) %>% # mask for Cleveland rows of input (not original Cleveland rows)
  
        which %>% 
        row.idx[.] %>% # transform back to original Cleveland rows
        return
}


f4 <-  function(nodeID){
  svN <- ti[ti$nodeID==nodeID, 'splitvarName'] # svN : splitvarName
  #print(svN)
  #print(d.split[[nodeID]])
  sVal <- ti[ti$nodeID==nodeID, 'splitval']
  # print(sVal)
  
  d.split[[nodeID]] %>% 
    
    (function(x)
      # as.numeric is needed only for factors - but does not hurt when applying to already numeric featurs
      as.numeric(data.train[x,svN]) <= sVal
    ) %>% # mask for Cleveland rows of input (not original Cleveland rows)
    
    which %>% 
    d.split[[nodeID]][.] %>% # transform back to original Cleveland rows
    return
}


f5 <-  function(nodeID){
  #' based in the ranger package
  #' given a split node by nodeID
  #' given treeInfo and d.split in parent environment
  #' 
  #' return 

  svN <- ti[ti$nodeID==nodeID, 'splitvarName'] # svN : splitvarName
  #print(svN)
  #print(d.split[[nodeID]])
  sVal <- ti[ti$nodeID==nodeID, 'splitval']
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

pred <- predict(rg, predict.all = T ,type = 'terminalNodes', data=data.train)

tri <-  3
ti <- treeInfo(rg,tri) ; ti
plotTree1(rg,tri)
pred$predictions[,tri] %>% table # observations per terminal node

 d.split <- NULL
for(nodeID in ti[!ti$terminal, 'nodeID']){ # going through the split nodes / non-terminal nodes

  # initial split , root split
  if(nodeID==0){
    d.split[[1]] <- f5(0)
    d.split[[2]] <-  base::setdiff(1:nrow(data.train),d.split[[1]])
  }else{
  lc.nodeID <- ti[ti$nodeID==nodeID,'leftChild']
  d.split[[lc.nodeID]] <-  f5(nodeID)
  rc.nodeID <- ti[ti$nodeID==nodeID,'rightChild']
  d.split[[rc.nodeID]] <-  base::setdiff(d.split[[nodeID]],d.split[[lc.nodeID]])
  }
}

rg$inbag.counts %>% bind_rows()
rg$inbag.counts %>% unlist %>% matrix(nrow=500, byrow=F) %>% apply(1,sum) %>% summary

