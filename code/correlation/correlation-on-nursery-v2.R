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

#### ERROR : sometimes I get good values, at least 0.5 , sometimes small ones (max 0.2)

#### ERROR (?): I get higher correlations (more than 0.5) for forests with 1 tree
#### than for forests with 5 trees (0.05 to 0.15)
#### but it changes continuously as the size moves 1,2,3,4,5
#### is there an explanation : more stability for fewer trees? Sounds nonsense
#### covariances drop as number of trees grows, and variances drop too.
#### since correlation drops, too (as number of trees increases), the covariances drop faster than the variances (we use stdevs in the formula)?

rm(list=ls())

library(caret)
library(ranger)
library(dplyr)

source('code/source/subforest.R')
source('code/source/prep.R')

load('data/data_SupMat4.rda') # loads Hung, VA , Swiss , and Cleve (we do not work with Cleve directly!)
data.val <- VA

folder <- 'data/nursery'
files <- list.files(folder)
load(paste(folder,files[[1]], sep='/')) # loads doc
nLoops <- length(doc)
#nLoops <- 3

#doc <- vector('list', nLoops)
doc2 <- array(0,dim=c(nLoops,3,3))

sizeSF <- 1
# indices: 5 indices , 100 times , not random
idcsList <- lapply(seq(1,500,by=sizeSF), function(i) i:(i+(sizeSF-1)))
if(max(unlist(idcsList))>500){
  len <-  length(idcsList)
  idcsList <-  idcsList[1:(len-1)]
}

# createDataPartition(Hung$CAD, times=nLoops, p=0.5) -> val.set.partition
createDataPartition(data.val$CAD, times=nLoops, p=0.5) -> val.set.partition

for(i in 1:nLoops){
  #print(paste('i' , i))
  #print(paste(i/nBs , Sys.time()))
  if(i%%10 ==0) print(paste(i/nLoops , Sys.time()))
  
  # data.train <- doc[[i]]$`bootstapped training data`
  
  val <-  val.set.partition[[i]]
  # out of bag observations from sampling
  # observations not in the training set data.train the current forest was built on
  OOB <- doc[[i]]$resample %>% unique %>% setdiff(1:nrow(Cleve),.)
  data.sets <- list(Cleve[OOB,], data.val[val,], data.val[-val,])
  
  # use forest from nursery
  rg <- doc[[i]]$rg
  forest<-rg$forest
  
  #print('read forest from nursery')
  
  # predictions for positive class on data.sets
  pp1 <- lapply(1:3, 
              function(k) predict(forest , data=data.sets[[k]] , predict.all = T)$predictions[,2,])
  pp1 <- simplify2array(pp1, higher=FALSE)
  
  #print('predictions done')
  #print(length(pp1))
  #print(dim(pp1[[1]]))
  
  # use predictions (mean prediction for forests (sizeSF>1)) to calculate logloss
  if(sizeSF>1){
    # take the mean over the predictions of trees in forest
    LL <- lapply(1:3,
               function(j)
                 lapply(idcsList 
                        , function(k){
                          #print(paste('j=',j))
                          #print(paste('k=',k))
                          unlist(calcLogloss2(pp=apply(pp1[[j]][,k],1,mean) 
                                                        , df=data.sets[[j]] ) )
                          })
               )}else{
      # no need for a mean with only one prediction by one tree (sizeSF==1 : one tree in forest)
      LL <- lapply(1:3,
                   function(j)
                     lapply(idcsList 
                            , function(k) unlist(calcLogloss2(pp=pp1[[j]][,k] 
                                                              , df=data.sets[[j]] ) ) )
      )}
  
  LL <- cbind(unlist(LL[[1]]),unlist(LL[[2]]),unlist(LL[[3]]))
  
  #if(i%%10 ==0){
  #  boxplot(LL, main='Logloss on different sets'
  #               , sub=paste(sizeSF , 'tree(s) per forest')
  #                , ylab='logloss')
  #}
  
  doc2[i,,] <- cor(LL) 
}


#View(doc2)

doc2[10,,]
apply(doc2,c(2,3), mean)
apply(doc2,c(2,3), sd)

df1 <- apply(doc2,c(2,3), mean)
colnames(df1) <- c('LL.OOB','LL.val.1','LL.val.2')
rownames(df1) <- c('LL.OOB','LL.val.1','LL.val.2')
corrplot::corrplot( df1
                   , method='pie'
                   , type='lower'
                  , title=ifelse(sizeSF==1, 
                                 paste('correlation of logloss over a list of trees\nmean over out of bag observations from simulations\nand partitioned validation set\n',length(idcsList) ,' trees per simulation, ' , nLoops , 'simulations'),
                                 paste('correlation of logloss over a list of forests\nmean over out of bag observations from simulations\nand partitioned validation set\n',sizeSF,'tree(s) per forest, ' ,length(idcsList) ,'forests per simulation, ' , nLoops , 'simulations'))
                   , diag=TRUE
                   , outline='white'
                   , order='original'
                   , tl.srt = 45
                  # , sub='blubs' # very low below the plot, needs the margin (5,x,x,x)
                   , mar=c(0,0,4,0)+0.5)

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
