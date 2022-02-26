# 17.2.2022
#### correlation #####################################################
######################################################################
#### trees / forests performing well on one set ######################
#### will they also perform well on another (similar) set?  ##########
######################################################################

#### relevance : #####################################################
#### Chipman works with the best models first then adds models in 
#### decreasing performance if they give new insights (or some benefits)
#### OK for descriptive statistics
#### to use it for predictive statistics, we would like to see correlation 
#### of the performance of trees/forests on different data sets
#### i.e. : a good tree / good forest is good on all / most data sets
#### a low performing tree / forest does so on all / most data sets
######################################################################

#### not sure how this works as we look for correlation of small forests
#### does it indicate an error if correlation for forests is smaller than for single trees?
#### I'd expect higher correlation (-> stability!) from forests rather than from trees ??

#### HOW CAN I TEST THIS CODE ?? ####
#####################################

rm(list=ls())

library(caret)
library(ranger)
library(dplyr)

source('code/source/subforest.R')
source('code/source/prep.R')

load('data/data_SupMat4.rda') # loads Hung, VA , Swiss , and Cleve 
data.train <- Cleve
data.train$CAD_fac <-  NULL
data.val.1 <- VA
data.val.2 <- Swiss

sizeSF <- 1
num.forests <-  500
rg <- ranger(CAD~.
             , data = data.train
             , num.trees = sizeSF*num.forests
             , replace = F # neu nach Absprache mit AZ
             , mtry= 3 # default : 3
             , importance = 'impurity'
             , probability = T # this makes it a random forest of type 'Probability Estimation'
             , min.node.size = 13 # optimized in A.Z. paper: BiomJ 56, 2014, 
)
forest<-rg$forest

#doc <- vector('list', nLoops)
doc2 <- array(0,dim=c(3,3))


# indices: 5 indices , 100 times , not random
idcsList <- lapply(seq(1,rg$num.trees,by=sizeSF), function(i) i:(i+(sizeSF-1)))
if(max(unlist(idcsList))>rg$num.trees){
  len <-  length(idcsList)
  idcsList <-  idcsList[1:(len-1)]
}

createDataPartition(data.val.1$CAD, times=1, p=0.5) -> val.1.partition
createDataPartition(data.val.2$CAD, times=1, p=0.5) -> val.2.partition

val1 <-  val.1.partition[[1]]
val2 <-  val.2.partition[[1]]

data.sets <- list(data.val.1[val1,], data.val.1[-val1,], data.val.2[val2,],data.val.2[-val2,])
  
# predictions for positive class on data.sets
pp1 <- lapply(1:4, 
              function(k) predict(forest , data=data.sets[[k]] , predict.all = T)$predictions[,2,])
pp1 <- simplify2array(pp1, higher=FALSE)
  
# use predictions (mean prediction for forests (sizeSF>1)) to calculate logloss
if(sizeSF>1){
  # take the mean over the predictions of trees in forest
  LL <- lapply(1:4,
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
               LL <- lapply(1:4,
                   function(j)
                     lapply(idcsList 
                            , function(k) unlist(calcLogloss2(pp=pp1[[j]][,k] 
                                                              , df=data.sets[[j]] ) ) )
                   )}
LL <- cbind(unlist(LL[[1]]),unlist(LL[[2]]),unlist(LL[[3]]),unlist(LL[[4]]))
colnames(LL) <- c('LL.VA.1','LL.VA.2','LL.Swiss.1','LL.Swiss.2')

#boxplot(LL, main='Logloss on different sets'
#                 , sub=paste(sizeSF , 'tree(s) per forest')
#                  , ylab='logloss')
  
  
#View(doc2)

corrplot::corrplot( cor(LL)
                   , method='pie'
                   , type='lower'
                  , title=ifelse(sizeSF==1, 
                                 paste('correlation of logloss over a list of trees\nover two partitioned validation sets (VA, Swiss)\n', rg$num.trees ,' trees'),
                                 paste('correlation of logloss over a list of forests\nover two partitioned validation sets (VA, Swiss)\n',sizeSF, 'tree(s) per forest, ' ,length(idcsList) ,'forests' ))
                   , diag=TRUE
                   , outline='white'
                   , order='original'
                   , tl.srt = 45
                  # , sub='blubs' # very low below the plot, needs the margin (5,x,x,x)
                   , mar=c(0,0,3,0)+0.5)

