
rm(list=ls())

library(dplyr)
library(caret)
library(ranger)
library(effects)
library(cluster)

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

# set test data by name : VA, Swiss or Hung
data.test.name <-  'Swiss'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name

sizeSF <- 5

calc_LL_for_selection <- function(doc, sizeSF){
 
  nT <- doc[[1]]$rg$num.trees
  sizeDefault <- nT
  
  nBs <- length(doc)
  #nBs <- 5
  
  evalT <- matrix(NA,nrow=nBs,ncol=4) # table of evaluations
  
  ct <- 1
  for(i in 1:nBs){
    #print(paste(i/nBs , Sys.time()))
    if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
    
    data.train <- doc[[i]]$`bootstapped training data`

    OOB <-  base::setdiff(1:nrow(Cleve), unique(doc[[i]]$resample))
    data.set.val <- Cleve[OOB,] # goes back to original Cleveland data

    # grow small forest (size nT)
    rg <- doc[[i]]$rg
    forest<-rg$forest
    
    pp <- predict(forest 
              , data=data.set.val 
              , predict.all = T)$predictions[,2,]
    #pp <- simplify2array(pp, higher=FALSE)

    # is this the same as Vectorize ??!!
    lapply(1:rg$num.trees
       , function(k){ 
         pp[,k] %>% 
           calcLogloss2( df=data.set.val ) %>% 
           unlist
         # unlist(calcLogloss2(pp=pp[,k] , df=data.set.val ) ) 
       }) %>% 
    unlist -> LL
    
    hpst <- order(LL)[1:sizeSF] # high performing single trees for second basecase
    
  # I could do this as above , pp , then calcLogloss2 (maybe with Vectorize) for the full forest
  # and then calculate logloss with subsets of pp
    
  evalT[i,] <- c(i
                , calcLogloss(subforest(forest,sample(1:rg$num.trees, sizeSF)), data.test)
                , calcLogloss( subforest(forest, hpst ), data.test)
                , calcLogloss(forest, data.test)
                )
  }
  evalT <-  data.frame(evalT)
  names(evalT) <- c('sim','LL.test.random','LL.test.hp', 'LL.test.default')
  return(evalT)
}

# to base the result on more bootstraps
folder <- 'data/nursery'
files <- list.files(folder)[1:5]
# dir(folder)
collector <-  list()
ct <-  1 # counter for the above collector

for(file in files){
  # run loops over doc loaded from file
  load(paste(folder,file, sep='/'))
  collector[[ct]] <- calc_LL_for_selection(doc , sizeSF)
  ct <-  ct+1
}
et <- bind_rows(collector)

names(et)[2:4] <- c('random','unimod hp','default')
boxplot(et[,2:4]
        , ylab='logloss'
        #, main=paste(attr(et,'metric'), 'metric ,', num.clusters , 'clusters ,', sizeSF, 'trees in forest\nfor non-default')
        )

#names(et)
et %>% 
  select(random, 'unimod hp', default) %>%
  summarize_all(function(x) c(mean(x), sd(x))) %>%
  t %>% 
  xtable -> et.xt

et.xt
digits(et.xt) <- 4
et.xt

et %>% 
  select(random, 'unimod hp', default) %>%
  summarize_all(function(x) c(median(x), IQR(x)))
