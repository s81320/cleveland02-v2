# ranger_num_trees_v3.R
# look at different sizes in a ranger random forest 
# and compare their performance in logloss , especially in mean and sd (median and IQR)

# produces as boxplot and table(s)

# can be used with simulation and without
# and for different sets as training and test sets

# 3.2.2022

# Use Cleve as training data, Swiss as validation data, Hung as test

rm(list=ls())

library(ranger)
library(dplyr)
library(xtable)

source('code/source/prep.R') # calcLogloss

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

#data.test <-  VA
data.test <- Hung

# nt <- seq(5, 5000, 5)
nt <- c(5,10,50,500)
length(nt)

load('data/nursery/nursery01_02_2x500.rda') # loads doc and info
info

nLoops <-  length(doc)
LL <- matrix(0,nLoops, length(nt))

# build forest up , starting with nt[1] trees , ending with max(nt) trees

for(j in 1:nLoops){
  set.seed(2*j)
  if(j%%10 ==0) print(paste(j/nLoops , Sys.time()))
  
  data.train <- doc[[j]]$`bootstapped training data`
  
  rg <- doc[[j]]$rg
  
  # predictions for all observations and all trees
  preds <- predict(rg, data.test, predict.all = TRUE)$predictions[,2,]
  dim(preds) # obs , trees
  
  # build forest up , starting with nt[1] trees , ending with max(nt) trees
  for(i in 1:length(nt)){
  
    pp <- preds[,1:nt[i]] %>% rowMeans
    correctedpp <- ifelse(data.test$CAD=='Yes',pp,1-pp) # problematic when this returns 0
    correctedpp <- winsorize_probs(correctedpp) 
    
    LL[j,i] <- -mean(log(correctedpp))
  }
}

LL <-  as.data.frame(LL)
names(LL) <- nt

boxplot(LL
        , main=paste('logloss for forests of different sizes\n(trained on simulated Cleveland, tested on Hung)') 
        , sub=paste('N=',nrow(LL),'(nr of observations per size)')
        , xlab= 'forest size , number of trees in forest'
        , ylab='logloss')

xt <- apply(LL,2,function(x)c(mean(x),sd(x)))
rownames(xt) <- c('mean','sd')

t(xt) %>% xtable -> xt
digits(xt) <- 4
xt

xt <- apply(LL,2,function(x)c(median(x),IQR(x)))
rownames(xt) <- c('median','IQR')
t(xt) %>% xtable -> xt
digits(xt) <- 4
xt
