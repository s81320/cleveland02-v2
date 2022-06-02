# ranger_num_trees_v3-nursery.R
# look at different sizes in a ranger random forest 
# and compare their performance in logloss , especially in mean and sd (median and IQR)

# produces as boxplot and table(s)

# can be used with simulation and without
# and for different sets as training and test sets

# 11.3.2022

# Use Cleve as training data, Swiss as validation data, Hung as test

rm(list=ls())

library(ranger)
library(dplyr)
library(xtable)

source('code/source/prep.R') # calcLogloss

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

calc_LL_matrix <- function(doc , nt, data.test){
  #' uses random forests from the nursery to calculate loglosses for forests of different sizes (as subforests of the pre-grown forest)
  #' forests are nested, smaller forests are always subforests of any larger forest.
  #' 
  #' @param doc should contain a list with each element containing a ranger random forest
  #' @param nt for the number of trees in each smaller forest
  #' @param data.test test data to calculate logloss for each forest on
  
  nLoops <-  length(doc)
  LL <- matrix(0,nLoops, length(nt))
  
  # build forest up , starting with nt[1] trees , ending with max(nt) trees
  
  for(j in 1:nLoops){
    set.seed(2*j)
    if(j%%10 ==0) print(paste(j/nLoops , Sys.time()))
    
    # data.train <- doc[[j]]$`bootstapped training data`# typo in nursery
    data.train <- doc[[j]]$`bootstrapped training data` # typo corrected in nursery02
    
    rg <- doc[[j]]$rg
    
    # predictions for all observations and all trees
    # for test data
    preds <- predict(rg, data.test, predict.all = TRUE)$predictions[,2,]
    # dim(preds) # obs , trees
    
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
  
  return(LL)
}

# to base the result on more bootstraps
folder <- 'data/nursery02'
files <- list.files(folder)[1:2] ; files

collector <-  list()
ct <-  1 # counter for the above collector

data.test.name <-  'Swiss'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name

nt <- c(5,10,50,500)
#nt <- c(2,3,4,5,6,7,8,9, seq(10,500,10))

for(file in files){
  # run loops over doc loaded from file
  load(paste(folder,file, sep='/'))
  LL <- calc_LL_matrix(doc, nt=nt, data.test=data.test)
  # keep result from loop and merge them all together
  collector[[ct]] <- LL
  ct <-  ct+1
}

LL <- bind_rows(collector)
attr(LL,'data.test.name') <- attr(data.test,'data.test.name') 

#save(LL,file='data/LL_nested_regular_forests_size_2_to_500.rda')

boxplot(LL
        , main=paste('logloss for forests of different sizes\n(trained on simulated Cleveland, tested on', attr(data.test,'data.test.name'),')') 
        , sub=paste('N=',nrow(LL),'(nr of observations per size)')
        , xlab= 'forest size , number of trees in forest'
        , ylab='logloss')

xt <- apply(LL,2,function(x)c(mean(x),sd(x)))
rownames(xt) <- c('mean','sd')
t(xt) %>% xtable -> xxt
digits(xxt) <- 7
xxt

if(ncol(xt>9)){
  plot(colnames(xt)[1:9]
     , xt['mean',1:9]
     , type='b'
     , main='sucess over size, regular forests'
     , xlab='size: number of trees in forest (2,..,10)'
     , ylab='sucess: mean logloss')

  plot(colnames(xt)[9:18]
     , xt['mean',9:18]
     , type='b'
     , main='sucess over size, regular forests'
     , xlab='size: number of trees in forest (seq(10,100,10))'
     , ylab='sucess: mean logloss')

  plot(colnames(xt)[18:ncol(xt)]
     , xt['mean',18:ncol(xt)]
     , type='l'
     , main='sucess over size, regular forests'
     , xlab='size: number of trees in forest (seq(100,500,10))'
     , ylab='sucess: mean logloss')
}

xt <- apply(LL,2,function(x)c(median(x),IQR(x)))
rownames(xt) <- c('median','IQR')
t(xt) %>% xtable -> xxt
digits(xxt) <- 4
xxt

# tradeoff
(LL[,'50'] %>% mean)/((LL[,'500'] %>% mean))
(LL[,'5'] %>% mean)/(LL[,'500'] %>% mean)
(LL[,'5'] %>% mean-LL[,'500'] %>% mean)/((LL[,'500'] %>% mean))

(LL[,'50'] %>% sd)/((LL[,'500'] %>% sd))
(LL[,'5'] %>% sd)/((LL[,'500'] %>% sd))

