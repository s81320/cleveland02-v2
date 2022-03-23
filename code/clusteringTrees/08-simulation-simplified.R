
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
data.test.name <-  'Hung'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name

num.clusters <-  20
sizeSF <- 20
treesSelected <- vector('list',num.clusters)

calc_LL_for_selection <- function(doc, sizeSF, num.clusters){
  #metrices<-c('d0','d1','d2','sb')
  #metrices <- 'd0'
  
  nT <- doc[[1]]$rg$num.trees
  sizeDefault <- nT
  
  nBs <- length(doc)
  #nBs <- 5
  
  evalT <- matrix(NA,nrow=nBs,ncol=6) # table of evaluations
  
  ct <- 1
  for(i in 1:nBs){
    #print(paste(i/nBs , Sys.time()))
    if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
    
    DM <-  doc[[i]]$DM
    metric <- 'd2'
    dm <-  DM[[metric]]

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

    pam.obj <- cluster::pam(dm
                        , k= num.clusters
                        , diss=TRUE
                        , medoids='random'
                        , nstart = 5)

    # NAMES is important here!!
    # for pam.obj$silinfo$width
    # rows are ordered in descending sil width - not by indices of trees in the dissim matrix or forest

    (pam.obj$silinfo$width[,'sil_width']>0) %>% 
    which %>% # which returns row numbers
    names %>%  # names gets the original indices (as strings)
    as.integer -> goodTrees # tree indices from forest

    #goodTrees[order(goodTrees)]

    which(pam.obj$clusinfo[,'size'] >1) -> clustersSelected 

    for(nClus in clustersSelected){
      A1 <- which(pam.obj$clustering==nClus) # tree indices in dm , built on forest
      S1 <- base::intersect(A1,goodTrees) 
      # S1 : trees with pos sil width belonging to cluster nClus
      #length(S1) %>% print 
      
      if(sizeSF==num.clusters){
        howMany <- 1
      }else{
        howMany <- round(sizeSF*pam.obj$clusinfo[nClus,1] / sum(pam.obj$clusinfo[,1]),0)
      }
      
      LL.cluster <-  LL[S1] # logloss on validation set (OOB of simulation)
      
      # select trees with highest rank(LL.cluster)
      # print(paste('select ', howMany,' from cluster ', nClus))
      # which(rank(LL.cluster)<=howMany) %>% S1[.] %>% print

      treesSelected[[nClus]] <- which(rank(LL.cluster)<=howMany) %>% S1[.]
      }

  #print('trees selected')
  #treesSelected %>% unlist %>% length %>% print

    hpst <- which(rank(LL)<=sizeSF) # high performing single trees for second basecase
   
  # I could do this as above , pp , then calcLogloss2 (maybe with Vectorize) for the full forest
  # and then calculate logloss with subsets of pp
    
  evalT[i,] <- c(i
                , length(treesSelected)
                , calcLogloss(subforest(forest,unlist(treesSelected)), data.test)
                , calcLogloss(subforest(forest,1:sizeSF), data.test)
                , calcLogloss( subforest(forest, hpst ), data.test)
                , calcLogloss(forest, data.test)
                )
  }
  evalT <-  data.frame(evalT)
  names(evalT) <- c('sim','sizeSF','LL.test.select','LL.test.random','LL.test.hp', 'LL.test.default')
  attr(evalT,'metric') <- metric
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
  collector[[ct]] <- calc_LL_for_selection(doc , sizeSF , num.clusters)
  ct <-  ct+1
}
boxplot(collector[[1]][,3:6])

boxplot(collector[[3]][,4]- collector[[3]][,3]
        , main='logloss overshoot for random vs selected subforest')
abline(0,0, col='grey')

boxplot(collector[[3]][,5]- collector[[3]][,3]
        , main='logloss overshoot for default vs selected subforest')
abline(0,0, col='grey')

et <- bind_rows(collector)
names(et)[3:6] <- c('select','random','unimod hp','default')
boxplot(et[,3:6]
        , ylab='logloss'
        , main=paste(attr(et,'metric'), 'metric ,', num.clusters , 'clusters ,', sizeSF, 'trees in forest\nfor non-default')
        )

#names(et)
et %>% 
  select(select, random, 'unimod hp', default) %>%
  summarize_all(function(x) c(mean(x), sd(x)))
