# 12.5.2022
#### correlation #####################################################
######################################################################
#### trees / forests performing well on one set ######################
#### will they also perform well on another (similar) set?  ##########
######################################################################

rm(list=ls())

library(caret)
library(ranger)
library(dplyr)

source('code/source/subforest.R')
source('code/source/prep.R')

load('data/data_SupMat4.rda') # loads Hung, VA , Swiss , and Cleve 
# use Cleve as training data in ranger call , use 10 predictors, one target, columns 1:11

opt1 <- FALSE 

if(opt1){
  # option 1
  data.sets <- list('Cleve'=Cleve, 'VA'=VA, 'Swiss'=Swiss, 'Hung'=Hung)
}else{
  # option 2
  createDataPartition(VA$CAD, times=1, p=0.5)[[1]]  -> VA.partition
  createDataPartition(Swiss$CAD, times=1, p=0.5)[[1]] -> Swiss.partition
  data.sets <- list('Cleve'=Cleve
                  , 'VA.1'=VA[VA.partition,]
                  , 'VA.2'=VA[-VA.partition,] 
                  , 'Swiss.1'=Swiss[Swiss.partition,]
                  , 'Swiss.2'=Swiss[-Swiss.partition,])
}


sizeSF <- 1 # sizeSF <- 1 looks at the individual trees and the correlation of their performance. We get a correlation around 0.3, 0.5
# sizeSF <- 5 would look at small forests of 5 trees and their correlated performance. Strangely here correlation is basically 0
num.forests <-  500 # if sizeSF==1 this is the num.trees in the ranger call
rg <- ranger(CAD~.
             , data = Cleve[,1:11]
             , num.trees = sizeSF*num.forests
             , replace = F # neu nach Absprache mit AZ
             , mtry= 3 # default : 3
             , importance = 'impurity'
             , probability = T # this makes it a random forest of type 'Probability Estimation'
             , min.node.size = 13 # optimized in A.Z. paper: BiomJ 56, 2014, 
)
forest<-rg$forest

# indices: sizeSF==5 : chopping the large forest 1,2,...500 into many small forests of size 5 , sequential (1,..,5), (6,.., 10) nothing random
idcsList <- lapply(seq(1,rg$num.trees,by=sizeSF), function(i) i:(i+(sizeSF-1)))
if(max(unlist(idcsList))>rg$num.trees){
  len <-  length(idcsList)
  idcsList <-  idcsList[1:(len-1)]
}

# predictions for positive class on data.sets
pp1 <- lapply(1:length(data.sets), 
              function(k) predict(forest , data=data.sets[[k]] , predict.all = T)$predictions[,2,])
pp1 <- simplify2array(pp1, higher=FALSE)
  
# use predictions (mean prediction for forests (sizeSF>1)) to calculate logloss
if(sizeSF>1){
  # take the mean over the predictions of trees in forest
  LL <- lapply(1:length(data.sets),
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
               LL <- lapply(1:length(data.sets),
                   function(j)
                     lapply(idcsList 
                            , function(k) unlist(calcLogloss2(pp=pp1[[j]][,k] 
                                                              , df=data.sets[[j]] ) ) )
                   )}

LL <- matrix(unlist(LL), nrow=length(LL), byrow=TRUE) %>% t
colnames(LL) <- names(data.sets)

par(mar=c(0,0,0,0)+0.1)
corrplot::corrplot( cor(LL)
                   , method='pie'
                   , type='lower'
                 # , title=ifelse(sizeSF==1, 
                #                 paste('correlation of logloss over a list of trees\nover two partitioned validation sets (VA, Swiss)\n', rg$num.trees ,' trees'),
                #                 paste('correlation of logloss over a list of forests\nover two partitioned validation sets (VA, Swiss)\n',sizeSF, 'tree(s) per forest, ' ,length(idcsList) ,'forests' ))
                   , diag=TRUE
                   , outline='white'
                   , order='original'
                   , tl.srt = 45
                  # , sub='blubs' # very low below the plot, needs the margin (5,x,x,x)
                   # , mar=c(0,0,0,0)+0.5
                )

