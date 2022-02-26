
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
metric <- 'd0'
cutoff.parameter <-  0.7

calc_LL_for_selection <- function(doc, sizeSF){
 
  nT <- doc[[1]]$rg$num.trees
  sizeDefault <- nT
  
  nBs <- length(doc)
  #nBs <- 5
  
  evalT <- matrix(NA,nrow=nBs,ncol=7) # table of evaluations
  
  ct <- 1
  for(i in 1:nBs){
    #print(paste(i/nBs , Sys.time()))
    if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
    
    #DM <-  doc[[i]]$DM
    dm <-  doc[[i]]$DM[[metric]]
    
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
       }) %>% 
    unlist -> LL
    
    # LL can be calculated with vectorize :
    #((function(k){ 
    #  pp[,k] %>% 
    #    calcLogloss2( df=data.set.val ) %>% 
    #    unlist
    #} )%>%
    #Vectorize(.))(1:rg$num.trees) -> LL
    
    oLL <- order(LL) # tree indices , ordered by logloss on OOB
    
    chipForest <- oLL[1] # start forest with the best tree
    
    cutoff <- quantile(dm,cutoff.parameter) # rough starting point
    #print(cutoff)

    for(j in 2:rg$num.trees){
      if(length(chipForest)==sizeSF){
        break
        }else{
        trindx <- oLL[j]
        #print(paste(j,trindx))
        #print(paste(length(chipForest)))
        #print(paste('distances to current candidate:',dm[chipForest,trindx]))
        if(min(dm[chipForest,trindx])>cutoff){
          #print(paste('logloss single tree now added: ' , LL[trindx]))
          chipForest <- c(chipForest, trindx)
        }
        }
    } # exits with j the first index after all trees have been collected
    
  evalT[i,] <- c(i
                 , length(chipForest)
                 , j
                , calcLogloss(subforest(forest,sample(1:rg$num.trees, sizeSF)), data.test)
                , calcLogloss( subforest(forest, oLL[1:sizeSF] ), data.test)
                , calcLogloss( subforest(forest, chipForest ), data.test)
                , calcLogloss(forest, data.test)
                )
  }
  evalT <-  data.frame(evalT)
  names(evalT) <- c('sim','sizeChipForest','j','LL.test.random','LL.test.hp', 'LL.test.chip', 'LL.test.default')
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
et01 <- bind_rows(collector)

mean(et01$j)
mean(et01$LL.test.chip)

names(et01)[4:7] <- c('random','unimod hp','chip','default')
boxplot(et01[,4:7]
        , ylab='logloss'
        , main=paste(sizeSF,'trees in small forest\n cutoff.parameter: ', cutoff.parameter, ' metric:' , metric)
        #, main=paste(attr(et,'metric'), 'metric ,', num.clusters , 'clusters ,', sizeSF, 'trees in forest\nfor non-default')
        )

#names(et)
et01 %>% 
  select(random, 'unimod hp', chip,  default) %>%
  summarize_all(function(x) c(mean(x), sd(x))) %>%
  t #%>% 
  xtable -> et.xt
  
et.xt
digits(et.xt) <- 4
et.xt

et01 %>% 
  select(random, 'unimod hp', chip, default) %>%
  summarize_all(function(x) c(median(x), IQR(x)))
