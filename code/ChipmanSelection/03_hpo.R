# 23.2.2022
# code/chipmanSelection/03_hpo.R
# parameter optimization for each metric

# validation data Swiss -> check if still valid !

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

calc_chipForest <- function(dm , forest, oLL, parameter){
  #' calculates the Chipman forest
  #' 
  #' @param parameter is a list with items cutoff and sizeSF 
  #'
  #' returns its logloss on data.test (from higher environment)
  #' returns the rank of the last inspected tree (or the index +1 , like 6 when collecting the first 5 trees)
  #' returns the number of trees in the Chipman forest (enforced to be parameters$sizeSF or smaller)

  chipForest <- oLL[1] # start chipForest with the best tree
  cutoff <- quantile(dm,parameter$cutoff) 
  #print(cutoff)
  
  # go through trees (ordered by performance) until sizeSF trees are collected / selected
  for(j in 2:length(oLL)){
    if(length(chipForest)==parameter$sizeSF){
      break
    }else{
      trindx <- oLL[j]
      # if new tree far from (current) chipForest , add it to chipForest
      # if cutoff =0 then each tree is added with certainty -> same as unimodal!
      if(min(dm[chipForest,trindx]) >= parameter$cutoff){ 
        chipForest <- c(chipForest, trindx)
      }
    }
  } # for exits with j the first index after all trees have been collected
  return(list(calcLogloss( subforest(forest, chipForest ), data.test)
              ,j
              , length(chipForest)))
}

calc_meinForest <- function(dm , forest, oLL, parameter){
  #' calculates the Meiner forest
  #' 
  #' @param parameter is a list with items cutoff and sizeSF 
  #'
  #' returns its logloss on data.test (from higher environment)
  #' returns the rank of the last inspected tree (or the index +1 , like 6 when collecting the first 5 trees)
  #' returns the number of trees in the Meiner forest (enforced to be parameters$sizeSF or smaller)
  
  meinForest <- oLL[1] # start chipForest with the best tree
  cutoff <- quantile(dm,parameter$cutoff) 
  #print(cutoff)
  
  # go through trees (ordered by performance) until sizeSF trees are collected / selected
  for(j in 2:length(oLL)){
    if(length(meinForest)==parameter$sizeSF){
      break
    }else{
      trindx <- oLL[j]
      # if new tree far from (current) meinForest , add the best tree (lowest logloss) to the meinForest that represents trindx.
      # if cutoff =0 then each tree is added with certainty -> same as unimodal!
      if(min(dm[meinForest,trindx]) >= parameter$cutoff){ 
        base::setdiff(oLL[1:j], meinForest) %>%
          dm[.,trindx] %>%
          (function(x) x<parameter$cutoff) %>%
          which.min %>%
          base::setdiff(oLL[1:j], meinForest)[.] %>%
          c(meinForest) -> meinForest
          # alternative to calc_chipForest :
         # chipForest <- c(chipForest, trindx)
      }
    }
  } # for exits with j the first index after all trees have been collected
  return(list(calcLogloss( subforest(forest, meinForest ), data.test)
              , j
              , length(meinForest)))
}

calc_LL_for_selection <- function(doc, parameter){
  
  nT <- doc[[1]]$rg$num.trees
  sizeDefault <- nT
  
  nBs <- length(doc)
  #nBs <- 5
  
  evalT <- matrix(NA,nrow=nBs,ncol=12) # table of evaluations
  
  ct <- 1
  for(i in 1:nBs){
    #print(paste(i/nBs , Sys.time()))
    #if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
    
    DM <-  doc[[i]]$DM
    
    data.train <- doc[[i]]$`bootstapped training data`
    OOB <-  base::setdiff(1:nrow(Cleve), unique(doc[[i]]$resample))
    data.set.val <- Cleve[OOB,] # goes back to original Cleveland data
    
    forest<-doc[[i]]$rg$forest
    
    pp <- predict(forest 
                  , data=data.set.val 
                  , predict.all = T)$predictions[,2,]
    
    lapply(1:forest$num.trees
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
    
    metrices <- names(DM)
    res <- list()
    for(metric in metrices){
      DM[[metric]] %>%
        #calc_chipForest(forest, oLL, parameter) -> res[[metric]]
        calc_meinForest(forest, oLL, parameter) -> res[[metric]]
    }
    evalT[i,] <- unlist(res)
  }
  
  evalT <-  data.frame(evalT)
  names(evalT) <- c(paste(rep(c('LL.test.chip.','j.','size.'),4),rep(c('d0','d1','d2','sb'),each=3),sep=''))
  return(evalT)
}

#paste(rep(c('LL.test.chip.','j.'),4),rep(c('d0','d1','d2','sb'),each=2),sep='')

# to base the result on more bootstraps
folder <- 'data/nursery'
files <- list.files(folder)[1:2]
# dir(folder)

cutoff <- seq(0, 0.5,by=0.05)
# setting sizeSF to the number of trees (500)  generates a Chipman forest 
# that only stops when all trees are represented (at the specified level)
# cutoff 0 means only trees with dissim 0 represent each other. Happens only for d0.
sizeSF <- rep(5,length(cutoff))

parameter <-  data.frame(cbind(cutoff, sizeSF))
# parameter$cutoff
# parameter <- list('sizeSF' = 5
#            , 'cutoff'=0.5) # old
collector.p <-  list()
ct.p <-  1

for(p in 1:nrow(parameter)){
  print(paste('parameter ', p ,'at', Sys.time()))

  collector.f <-  list()
  ct.f <-  1 # counter for the above collector 
  for(file in files){
    # run loops over doc loaded from file
    load(paste(folder,file, sep='/')) # loads doc
    collector.f[[ct.f]] <- calc_LL_for_selection(doc , parameter[p,])
    ct.f <-  ct.f+1
  }
  collector.p[[ct.p]] <- c(parameter[p,], apply(bind_rows(collector.f),2,mean)) %>% unlist
  ct.p <- ct.p+1
}

et03 <- data.frame(bind_rows(collector.p))

#save(et03 , file='data/chipman/hyperparameter_cutoff_5trees*.rda')
#load('data/chipman/hyperparameter_cutoff_50trees.rda')

# makes no sense. et03 is already aggregated, averaged over...
#et03 %>% apply(2,function(x) c(mean(x),sd(x))) %>% t

et03.s <-  et03 %>% select(tidyr::starts_with('LL.'))
ylim <- range(et03.s)
# run only with loaded data that has not been created in this session
# cutoff <- unique(et03$cutoff)
# sizeSF <-  unique(et03$sizeSF)
plot(cutoff
     , et03.s[,1]
     , type='b'
     , main=paste('loglosses for cutoff parameters, ', 50*length(files) ,' simulations\n(filled dots for minimum values)')
     , ylim=ylim
     , xlab='cutoff parameter (quantile of dissimilarities)'
     , ylab='mean logloss')
wm <- which.min(et03.s[,1])
points(cutoff[wm], et03.s[wm,1],col=1, pch=19 )
for(k in 2:4){
  points(cutoff
         , et03.s[,k]
         , type='b'
         , col=k)
  wm <- which.min(et03.s[,k])
  points(cutoff[wm], et03.s[wm,k],col=k, pch=19 )
}
legend('topleft', legend=c('d0','d1','d2','sb') , pch=1, col=1:4 , cex=0.8)

et03.s <- et03 %>% select(tidyr::starts_with('j'))
ylim <- range(et03.s)
plot(cutoff
     , et03.s[,1]
     , type='b'
     , main=paste('mean number of trees in selection process\n(',50*length(files),' simulations, ', sizeSF[[1]], ' trees)', sep='')
     , ylim=ylim
     , xlab='cutoff parameter (quantile of dissimilarities)'
     , ylab='trees oredered by OOB performance')
#points(cutoff[wm], et03.s[wm,1],col=1, pch=19 )
for(k in 2:4){
  points(cutoff
         , et03.s[,k]
         , type='b'
         , col=k)
 # points(cutoff[wm], et03.s[wm,k],col=k, pch=19 )
}
legend( #'topleft'
  'bottomright'
       , legend=c('d0','d1','d2','sb') 
       , pch=1
       , col=1:4
       , cex=0.8)

et03.s <-  et03 %>% select(tidyr::starts_with('size.'))
ylim <- range(et03.s)
# run only with loaded data that has not been created in this session
# cutoff <- unique(et03$cutoff)
# sizeSF <-  unique(et03$sizeSF)
plot(cutoff
     , et03.s[,1]
     , type='b'
     , main=paste('mean sizes of Chipman forests,', 50*length(files) , ' simulations')
     , ylim=ylim
     , xlab='cutoff parameter (quantile of dissimilarities)'
     , ylab='size')
wm <- which.min(et03.s[,1])
points(cutoff[wm], et03.s[wm,1],col=1, pch=19 )
for(k in 2:4){
  points(cutoff
         , et03.s[,k]
         , type='b'
         , col=k)
  wm <- which.min(et03.s[,k])
  points(cutoff[wm], et03.s[wm,k],col=k, pch=19 )
}
legend('topright', legend=c('d0','d1','d2','sb') , pch=1, col=1:4 , cex=0.8)

