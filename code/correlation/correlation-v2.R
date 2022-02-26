######################################################################
#### simplification : calc pred metric, logloss and mae ##############
#### get correlation (and document it) ###############################
######################################################################

#### relevance : #####################################################
#### 1) d2 carries over between the training and validation / test set 
#### between cleveland and hungary data
#### 2) d2 and logloss are indeed connected. measuring d2 on training or validation data helps
#### helps with logloss on test data (our true aim)
#######################################################################

#### method : looking at correlation mattix only. Simple and straight.
#######################################################################

rm(list=ls())

library(caret)
library(ranger)
library(dplyr)

source('code/source/subforest.R')
source('code/source/prep.R')

load('data/data_SupMat4.rda')

Cleve$CAD_fac <- NULL
data.train <- Cleve

calcPValue <- function(a) pf(a[[1]],a[[2]],a[[3]],lower.tail=FALSE)

# check it is working
lm(logloss~mDiss, data=et.s) %>% summary -> s1
s1
calcPValue(s1$fstatistic)
# end of check

nLoops <- 10

#doc <- vector('list', nLoops)
doc <- array(0,dim=c(nLoops,5,5))
set.seed(1999)

createDataPartition(Hung$CAD, times=nLoops, p=0.5) -> val.set.partition

# random indices: 5 indices , 100 times
idcsList <- lapply(seq(1,500,by=5), function(i) i:(i+4))
#idcsList <- lapply(rep(5,100), function(i) sample(1:rg$num.trees,i,F)) 

for(i in 1:nLoops){
  
  rg <- ranger(CAD~.
               , data = data.train
               , num.trees = 500
               , replace = F # neu nach Absprache mit AZ
               , mtry= 3 # default : 3
               , importance = 'impurity'
               , probability = T # this makes it a random forest of type 'Probability Estimation'
               , min.node.size = 13 # optimized in AZ paper: BiomJ 56, 2014
  )
  
  print(paste(i , Sys.time()))
  
  val <-  val.set.partition[[i]]
  
  #data.val <- Hung[val,]
  #data.test <- Hung[-val,]
  
  data.val <- VA
  data.test <- Hung
  
  # 100 sub-forests of 5 trees
  forestList <- lapply(idcsList, function(vec) subforest(rg$forest, vec))
  
  f1 <-  function(vec, data){
    (vec - predForest(rg$forest,data)) %>% abs %>% mean
  } # prediction dissimilarity based on mean absolute error of predictions
  
  # predictions of the sub-forests on training data
  ppList <-  lapply(forestList, predForest, data=data.train) # passed additional argument to predForest
 
  # comparing predictions with the prediction of the full forest : that's the prediction metric
  lapply(ppList, FUN=function(vec) f1(vec,data.train)) %>%
    unlist -> d2.train
  
  # predictions of the sub-forests on training data
  ppList <-  lapply(forestList, predForest, data=data.val) # passed additional argument to predForest
  lapply(ppList, FUN=function(vec) f1(vec,data.val)) %>%
    unlist -> d2.val
  
  # predictions of the sub-forests on training data
  ppList <-  lapply(forestList, predForest, data=data.test) # passed additional argument to predForest
  lapply(ppList, FUN=function(vec) f1(vec,data.test)) %>%
    unlist -> d2.test

  lapply(forestList, function(for1) calcLogloss(for1,df=data.val)) %>% 
    unlist -> LL.val # logloss on validation set
  
  lapply(forestList, function(for1) calcLogloss(for1,df=data.test)) %>% 
    unlist -> LL.test # logloss on test set
  
  df1 <- data.frame(cbind(d2.train, d2.val, d2.test, LL.val,LL.test)) # d2 dissimilarities on 3 different data sets
  cm <- cor(df1)
  
  loglosses <- data.frame(cbind(LL.val,LL.test))
  print(cor(loglosses))
  
  doc[i,,] <- cm
  
  #idcsList <-  NULL
  forestList <- NULL
  ppList <-  NULL
  d2.train <- NULL
  d2.val <-  NULL
  d2.test <-  NULL
  LL.val <-  NULL
  LL.test <-  NULL
  df1 <-  NULL
  cm <-  NULL
  
}


View(doc)

doc[1,,]
apply(doc,c(2,3), mean)
apply(doc,c(2,3), sd)

df1 <- doc[1,,]
colnames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')
rownames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')

par(mar=c(5,5,5,2)+0.1)
corrplot::corrplot(df1
                   , method='pie'
                   , type='lower'
                   , title='correlation over a list of forests'
                   , sub='a single correlation matrix'
                   , diag=TRUE
                   , outline='white'
                   , order='original')

df1 <- apply(doc,c(2,3), mean)
colnames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')
rownames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')
corrplot::corrplot( df1
                   , method='pie'
                   , type='lower'
                   , title='correlation over a list of forests\nmean over partitions of validation-test split'
                   , diag=TRUE
                   , outline='white'
                   , order='original')


df1 <- apply(doc,c(2,3), sd)
colnames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')
rownames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')
corrplot::corrplot( df1
                    , method='pie'
                    , type='lower'
                    , title='sd of correlation over a list of forests\nsd over partitions of validation-test split'
                    , diag=TRUE
                    , outline='white'
                    , order='original')
