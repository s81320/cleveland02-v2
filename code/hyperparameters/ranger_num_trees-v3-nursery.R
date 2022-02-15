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
  
    data.train <- doc[[j]]$`bootstapped training data`
  
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
folder <- 'data/nursery'
files <- list.files(folder)
# dir(folder)
collector <-  list()
ct <-  1 # counter for the above collector

#data.test <- Hung
#attr(data.test,'data.test.name') <- 'Hung'

data.test.name <-  'VA'
# data.test <-  VA
data.test <-  get(data.test.nayme)
attr(data.test,'data.test.name') <- data.test.name

nt <- c(5,10,50,500)

for(file in files){
  # run loops over doc loaded from file
  load(paste(folder,file, sep='/'))
  LL <- calc_LL_matrix(doc, nt=nt, data.test=data.test)
  # keep result from loop and merge them all together
  collector[[ct]] <- LL
  ct <-  LL_ct+1
}

LL <- bind_rows(collector)
attr(LL,'data.test.name') <- attr(data.test,'data.test.name') 

boxplot(LL[,c(3,4)]
        , main=paste('logloss for forests of different sizes\n(trained on simulated Cleveland, tested on', attr(data.test,'set.name'),')') 
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

