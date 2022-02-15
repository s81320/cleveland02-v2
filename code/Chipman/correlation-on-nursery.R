#################################################
#### simplification : calc logloss ##############
#### get correlation (and document it) ##########
#################################################

#### relevance : #####################################################
#### Chipman works with the best models first then adds models in 
#### decreasing performance if they give new insights (or some benefits)
#### OK for descriptive statistics
#### to use it for predictive statistics, we would like to see correlation 
#### of the performance of trees/forests on different data sets
#### i.e. : a good tree / good forest is good on all / most data sets
#### a low performing tree / forest does so on all / most data sets
############################################################################

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
doc2 <- array(0,dim=c(nLoops,3,3))

sizeSF <- 50
# indices: 5 indices , 100 times , not random
idcsList <- lapply(seq(1,500,by=sizeSF), function(i) i:(i+(sizeSF-1)))

# createDataPartition(Hung$CAD, times=nLoops, p=0.5) -> val.set.partition
createDataPartition(data.val$CAD, times=nLoops, p=0.5) -> val.set.partition

for(i in 1:nLoops){
  #print(paste(i/nBs , Sys.time()))
  if(i%%10 ==0) print(paste(i/nLoops , Sys.time()))
  
  data.train <- doc[[i]]$`bootstapped training data`
  
  val <-  val.set.partition[[i]]
  data.v <- data.val[val,]
  data.t <- data.val[-val,] # t for test , we pretend half of the validation data was our test data
  
  # out of bag observations from sampling
  # observations not in the training set data.train the current forest was built on
  OOB <- doc[[i]]$resample %>% unique %>% setdiff(1:nrow(Cleve),.)
  
  # use forest from nursery
  rg <- doc[[i]]$rg
  forest<-rg$forest
  
  pp.OOB <- predict(forest , data=Cleve[OOB,] , predict.all = T)$predictions 
  pp.all.v <- predict(forest , data=data.v , predict.all = T)$predictions 
  pp.all.t <- predict(forest , data=data.t , predict.all = T)$predictions 
  
  if(sizeSF >1){
  lapply(idcsList , function(i) calcLogloss2(pp=apply(pp.OOB[,2,i],1,mean) , df=Cleve[OOB,])) %>% 
    unlist -> LL.OOB # logloss on OOB of simulation
  
  lapply(idcsList , function(i) calcLogloss2(pp=apply(pp.all.v[,2,i],1,mean) , df=data.v)) %>% 
    unlist -> LL.val # logloss on validation set
  
  lapply(idcsList , function(i) calcLogloss2(pp=apply(pp.all.t[,2,i],1,mean) , df=data.t)) %>% 
    unlist -> LL.test # logloss on test set
  }else{
    lapply(idcsList , function(i) calcLogloss2(pp=pp.OOB[,2,i] , df=Cleve[OOB,])) %>% 
      unlist -> LL.OOB # logloss on OOB of simulation
    
    lapply(idcsList , function(i) calcLogloss2(pp=pp.all.v[,2,i] , df=data.v)) %>% 
      unlist -> LL.val # logloss on validation set
    
    lapply(idcsList , function(i) calcLogloss2(pp=pp.all.t[,2,i] , df=data.t)) %>% 
      unlist -> LL.test # logloss on test set
  }
  df1 <- data.frame(cbind( LL.OOB, LL.val , LL.test)) # d2 dissimilarities on 3 different data sets
  par(mar=c(5,4,4,4)+0.1)
  if(i%%10 ==0) boxplot(df1
          , main='Logloss on different sets'
          , sub=paste(sizeSF , 'tree(s) per forest')
          , ylab='logloss')
  
  cm <- cor(df1)
  
  doc2[i,,] <- cm
}


#View(doc2)

#doc2[10,,]
apply(doc2,c(2,3), mean)
apply(doc2,c(2,3), sd)

df1 <- apply(doc2,c(2,3), mean)
colnames(df1) <- c('LL.OOB','LL.val','LL.test')
rownames(df1) <- c('LL.OOB','LL.val','LL.test')
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
colnames(df1) <- c('LL.val','LL.test')
rownames(df1) <- c('LL.val','LL.test')
corrplot::corrplot( df1
                    , method='pie'
                    , type='lower'
                    , title='sd of correlation over a list of forests\nsd over partitions of validation-test split'
                    , diag=TRUE
                    , outline='white'
                    , order='original')
