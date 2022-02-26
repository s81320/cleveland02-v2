######################################################################
#### simplification : calc pred metric, logloss and mae ##############
#### get correlation (and document it) ###############################
######################################################################

#### relevance : #####################################################
#### 1) d2 only slightly carries over between to the validation data 
#### i.e. from Cleveland to VA data (correlation of about 25%)
#### 2) the bond between the splits of the validation data is stronger (about 70%)
#### which tells us there is a considerable statistic difference in the training and validation data.
#### 3) d2 and logloss are quite independent of each other
#### measuring d2 on training or validation data does not help
#### with logloss on test data (our true aim)
#######################################################################
#### actually we can stop here. 
#### No, but we have to hope that clustering does the magic all by itself
#######################################################################

#### method : looking at correlation matrix only. Simple and straight.
#######################################################################

#### there is probably an error hidden (15.2.2022)
#######################################################################

#### d2 is the distance to the default forest. 
#### So it belongs to cleveland04 not cleveland02 project (15.2.2022)
#######################################################################

rm(list=ls())

library(caret)
library(ranger)
library(dplyr)

source('code/source/subforest.R')
source('code/source/prep.R')

load('data/data_SupMat4.rda') # loads Hung, VA , Swiss , and Cleve (we do not work with Cleve directly!)
data.val <- VA

calcPValue <- function(a) pf(a[[1]],a[[2]],a[[3]],lower.tail=FALSE)

folder <- 'data/nursery'
files <- list.files(folder)
load(paste(folder,files[[1]], sep='/')) # loads doc
nLoops <- length(doc)
#nLoops <- 10

#doc <- vector('list', nLoops)
doc2 <- array(0,dim=c(nLoops,5,5))

sizeSF <- 50
# indices: 5 indices , 100 times , not random
idcsList <- lapply(seq(1,500,by=sizeSF), function(i) i:(i+4))

# createDataPartition(Hung$CAD, times=nLoops, p=0.5) -> val.set.partition
createDataPartition(data.val$CAD, times=nLoops, p=0.5) -> val.set.partition

for(i in 1:nLoops){
  #print(paste(i/nBs , Sys.time()))
  if(i%%10 ==0) print(paste(i/nLoops , Sys.time()))
  data.train <- doc[[i]]$`bootstapped training data`

  val <-  val.set.partition[[i]]
  data.v <- data.val[val,]
  data.t <- data.val[-val,] # t for test , we pretend half of the validation data was our test data
  
  # grow small forest (size nT)
  rg <- doc[[i]]$rg
  forest<-rg$forest
  
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
  ppList <-  lapply(forestList, predForest, data=data.v) # passed additional argument to predForest
  lapply(ppList, FUN=function(vec) f1(vec,data.v)) %>%
    unlist -> d2.val
  
  # predictions of the sub-forests on training data
  ppList <-  lapply(forestList, predForest, data=data.t) # passed additional argument to predForest
  lapply(ppList, FUN=function(vec) f1(vec,data.t)) %>%
    unlist -> d2.test

  lapply(forestList, function(for1) calcLogloss(for1,df=data.v)) %>% 
    unlist -> LL.val # logloss on validation set
  
  lapply(forestList, function(for1) calcLogloss(for1,df=data.t)) %>% 
    unlist -> LL.test # logloss on test set
  
  df1 <- data.frame(cbind(d2.train, d2.val, d2.test, LL.val,LL.test)) # d2 dissimilarities on 3 different data sets
  cm <- cor(df1)
  
  doc2[i,,] <- cm
}

View(doc2)

doc2[10,,]
apply(doc2,c(2,3), mean)
apply(doc2,c(2,3), sd)

df1 <- doc2[1,,]
colnames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')
rownames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')

corrplot::corrplot(df1
                   , method='pie'
                   , type='lower'
                   , title='correlation over a list of forests'
                   , diag=TRUE
                   , outline='white'
                   , order='original')

df1 <- apply(doc2,c(2,3), mean)
colnames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')
rownames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')
corrplot::corrplot( df1
                   , method='pie'
                   , type='lower'
                   , title='correlation over a list of forests\nmean over partitioned validation set'
                   , diag=TRUE
                   , outline='white'
                   , order='original'
                   , tl.srt = 45
                  # , sub='blubs' # very low below the plot, needs the margin (x,x,5,x)
                   , mar=c(1,1,2,1))

df1 <- apply(doc2,c(2,3), sd)
colnames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')
rownames(df1) <- c('d2.train', 'd2.val', 'd2.test', 'LL.val','LL.test')
corrplot::corrplot( df1
                    , method='pie'
                    , type='lower'
                    , title='sd of correlation over a list of forests\nsd over partitions of validation-test split'
                    , diag=TRUE
                    , outline='white'
                    , order='original')

