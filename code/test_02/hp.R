# code/test_02/hp.R
# 31.5.2022

# script to test models with repetitions
# test like a data scientist : mix data and split it (homogenization),
# split it many different ways and take the mean over those splits: 
# Monte Carlo cross validation

rm(list=ls())

file_name_script <-  'code/test_02/hp.R'
list.files(path = 'code/test_02' , recursive = T)
assertthat::assert_that(file_name_script %in% list.files(recursive = T))


library(ranger)
library(xtable)
library(caret)
library(dplyr)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection
source('code/source/distance-matrices.R')
source('code/source/prep.R')

# create the data differently: mix Cleve, Swiss, Hung
# and partition it into test, dev, train

# find differences in how to gain the devset: 
# data.dev 'OOB' results in worse results than the default forest
# data.dev 'split' is and high performers outperform the default forest.
# -> OOB tests the splitting into test and train, if these sets are close, a tree build on train performs good on train and will perform good on test (as test is similar to train)
# if train and test are different, the the tree built on train performs good on train but not on test. 
# selecting trees that were grown on training sets with a similar test set is no quality measure for the tree.

# load data from beginning of script (Cleve, Swiss, Hung)

data.all <-  rbind(Cleve[,1:11], Hung[,1:11])
dim(data.all)

test_hp_forest <- function(data.all, data.val, sz=50 , nLoops=100, returnLL = F){
  
  print(data.val)
  if(!(data.val %in% c('OOB','split'))){ print('data.val should be OOB or split.') ; return() }
  
  set.seed(1)
  dp <- createDataPartition(data.all$CAD ,times = nLoops , p= 0.5) # dp : data partition
  doc <- data.frame(matrix(NA,nrow=nLoops, ncol=4))
  for(i in 1:nLoops){
    data.train <- data.all[dp[[i]],]
    data.test <- data.all[-dp[[i]],] # actual testing will be on data.test.set which is set depending on data.val
    
    rg <- ranger(CAD~.
                 , data = data.train 
                 , num.trees = 500
                 , replace = T # replace = T is good for doing predictions later on? would not have guessed from the name!
                 , mtry= 3 # default : 3
                 , importance = 'impurity'
                 , probability = T 
                 , min.node.size = 13 
                 , keep.inbag = T # need inbag data only for data.val =='OOB'
    )
    forest.full <- rg$forest
    
    if(data.val=='OOB'){
      data.test.set <- data.test
      
      # on individual OOB observations , (may be) different for each tree
      # prediction probabilities : pp
      pp <- predict(forest.full 
                    , data=data.train 
                    , predict.all = T)$predictions[,2,] # dim 303 x 500
      
      lapply(1:forest.full$num.trees
             , function(t){ 
               OOB <- which(rg$inbag.counts[[t]]==0)
               pp[OOB,t] %>% 
                 calcLogloss2( df=data.train[OOB,] ) %>% 
                 unlist}
      ) %>% 
        unlist -> LL # LL created , no data.val.set needed
    }
    
    if(data.val=='split'){
      fold1 <- createDataPartition(data.test$CAD, 1, 0.5) %>% unlist # split the ORIGINAL data.test which is a df
      data.val.set <-  data.test[fold1,]
      data.test.set <- data.test[-fold1,] # overwriting setting data.test.set at beginning of function
      LL <- calc_LL(forest.full, data.val.set) # LL for all trees in forest , cf Chipman.R
    }
    
    hp <- order(LL)[1:sz]
    
    doc[i,] <- c( calcLogloss(forest.full, data.test.set) # full forest
                  , sz # size of (selected sub-) forest
                  , calcLogloss(subforest(forest.full, hp), data.test.set) # logloss for (selecte sub-) forest
                  , calcLogloss(subforest(forest.full, 1:sz), data.test.set)) # logloss for regular small forest of same size as (selected sub-) forest
  }
  
  names(doc) <- c('default','size','high performers','regular small')
  return(doc)
}

####################
#### first test ####
####################

test_hp_forest(data.all , 'split' , nLoops = 10)

# give options on models to test directly
doc.hp <- test_hp_forest(data.all, data.val='split' ,sz=5 , nLoops=1000 )
doc.hp %>% apply(2,function(x) c(mean(x),sd(x)))
# doc.hp 

oversh <- doc.hp[,4]-doc.hp[,3]
par(mar=c(4,4,2,1)+0.1)
hist(oversh, breaks=100) # looks multimodal?
hist(oversh, breaks=10) # looks unimodal, bit skewed?

# loop over all / some models to test
# use expand grids to generate all combinations of parameters to be tested
parameter <- expand.grid(c('split','OOB') 
                         , c(5,50)
                         , stringsAsFactors = F 
                         , KEEP.OUT.ATTRS = F)
colnames(parameter) <- c('data.val','sz') 
parameter[1,] # is now a test case

results <-  list()
for(i in 1:nrow(parameter)){
  print(Sys.time())
  print(parameter[i,])
  results[[i]] <- test_hp_forest(data.all=data.all 
                                 , data.val=parameter[i,1] 
                                 , sz =parameter[i,2] 
                                 , nLoops=10000)
  results[[i]] %>% apply(2,function(x) c(mean(x),sd(x))) %>% t %>% xtable -> xt
  digits(xt) <-  4
  print(xt)
}

res.doc <- list(parameter=parameter
                , results=results
                , info=paste('training data is 50% of mixed Cleve and Hung. Code in'
                             , file_name_script)
)
# save(res.doc , file='data/test_02/results-hp-50-50split-of-Cleve-and-Hung-10000loops.rda')


################################################################################


# t.test

t.test(results[[1]]$`high performers` - results[[2]]$`high performers`, alternative = 'l')
(results[[1]]$`high performers` - results[[2]]$`high performers`) %>% hist(breaks=50)

t.test(results[[1]]$`high performers` - results[[1]]$`regular small`, alternative='l')

t.test(results[[2]]$`high performers` - results[[2]]$`regular small`, alternative='l')
