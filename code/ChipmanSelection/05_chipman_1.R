# 23.2.2022
# code/chipmanSelection/02.R
# trying to do all metrices in one loop
# cutoff parameter is set individually for each metric . cutoff parameter should have been selected with 03_hyperparamter_selection.R

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

calc_chipForest_1 <- function(dm, sizeSF, oLL, forest, alpha, metric){
  #print(pco)
  cutoff <- quantile(dm,alpha)
  #print(cutoff)
  
  # dm should be ordered
  dm2 <-  dm[oLL,oLL]
  #oLL[1]
  #LL[oLL]
  # N <- forest$num.trees 
  N <-  nrow(dm2)
  
  # everything is represented at a cutoff larger than (or equal to the) maximum dissimilarity
  represented <-  function(I){
    #' for I >1 , 
    #' representation by one tree alone it only at level of largest dissimilarity
    if(I %in% c(1,N-1)){
      dm2[1:I,(I+1):N] %>% max %>%
        (function(x){x<=cutoff})
    }else{
    dm2[1:I,(I+1):N] %>%
      apply(2,min) %>% # not working for I = 1 or I+1=N
      max %>%
      (function(x){x<=cutoff})
    }
  }
  I <- 1
  
  if(cutoff<max(dm2)){
    I <- 2
    while(!represented(I) && (I<N)){
      I <- I+1
      }
  }else{I <- 1}
  # exits with smallest I with represented(I) TRUE
  represented(I)
  print(paste('representing with 1..', I))
  
  #m <- c( "average", "single", "complete", "ward")
  #names(m) <- c( "average", "single", "complete", "ward")
  # function to compute coefficient
  #ac <- function(x) {
  #  agnes(dm2[1:I,1:I], method = x, diss=T)$ac
  #}
  #hpo <-  lapply(m, ac)
  
  # make this a basecase: Large Forest, smaller than default. Better?
  #calcLogloss( subforest(forest, 1:I ), data.test) %>% print
  
  hc <- cluster::agnes(dm2[1:I,1:I], 
                       method="ward", # always optimal - whenever I did hpo
                       #method = m[which.max(hpo)] , 
                       diss=T)
  #pltree(hc, cex = 0.6, hang = -1, main = paste(m[which.max(hpo)], ',', metric))
  
  hc.clus <- cutree(hc, k = sizeSF)
  hc.clus %>% table
  # cluster at a specified height , but we do not need that
  # ch6 <- cutree(as.hclust(hc), h = 0.3)
  # ch6 %>% table
  # optimal number of clusters?
  
  # from each cluster select the tree with the smallest logloss (on OOB data)
  lapply(1:5,function(i) which(hc.clus==i)[1]) %>% # we take the first, because indices are ordered, first is smallest
    unlist -> R # representing subset
  # if indices were not ordered, we'd do 
  #lapply(1:5,function(i) which(hc.clus==i) %>% min ) %>% unlist
  
  # that is the smallest index under the oLL ordering (as used for dm2)
  # then go back to the original tree indices!
  #oLL # original indices
  #oLL[R]
  #LL[oLL[R]] # check that the loglosses are small
  
  return(list(calcLogloss( subforest(forest, R ), data.test),I))
}

calc_LL_for_selection <- function(doc, sizeSF){
  
  nT <- doc[[1]]$rg$num.trees
  sizeDefault <- nT
  
  nBs <- length(doc)
  nBs <- 5
  
  evalT <- matrix(NA,nrow=nBs,ncol=12) # table of evaluations
  
  ct <- 1
  for(i in 1:nBs){
    #print(paste(i/nBs , Sys.time()))
    if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
    
    DM <-  doc[[i]]$DM
    
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
    
    metrices <- names(DM)
    res <- list()
    
    sizeSF <- 5
    
    if(sizeSF==5){
    # for 5 trees
    parameter.cutoff <-  list('d0'=0.2, 'd1'=0.2 , 'd2'=0.2, 'sb'=0.2) 
    }else{
      if(sizeSF==50){
        # for 50 trees
        parameter.cutoff <-  list('d0'=0.02, 'd1'=0.04 , 'd2'=0.02 , 'sb'=0.02)
      }else{
        print('no known cutoff parameter for given size of Chipman foorest')
      }
    }
    
    for(metric in metrices){
      #print(parameter.cutoff[metric])
      DM[[metric]] %>%
        calc_chipForest_1(sizeSF, oLL, forest, unlist(parameter.cutoff[metric]), metric) -> res[[metric]]
    }
    
    evalT[i,] <- c(i
                   , calcLogloss( forest, data.test)
                   , calcLogloss( subforest(forest, 1:sizeSF), data.test)
                   , calcLogloss( subforest(forest, oLL[1:sizeSF] ), data.test)
                   , unlist(res)
    )
  }
  
  evalT <-  data.frame(evalT)
  names(evalT) <- c('sim','LL.test.default', 'LL.test.random','LL.test.hp', 
                    paste(rep(c('LL.test.chip.','I.'),4),rep(c('d0','d1','d2','sb'),each=2),sep=''))
  return(evalT)
}

#paste(rep(c('LL.test.chip.','j.'),4),rep(c('d0','d1','d2','sb'),each=2),sep='')

# to base the result on more bootstraps
folder <- 'data/nursery'
files <- list.files(folder)[1:2]
# dir(folder)
collector <-  list()
ct <-  1 # counter for the above collector

for(file in files){
  # run loops over doc loaded from file
  load(paste(folder,file, sep='/'))
  collector[[ct]] <- calc_LL_for_selection(doc , sizeSF)
  ct <-  ct+1
}
et02 <- bind_rows(collector)

names(et02)[c(2:5,7,9,11)] <- c('default','random','not diverse','d0','d1','d2','sb')

#[4:7] <- c('random','unimod hp','chip','default')
boxplot(et02[,c(2:5,7,9,11)]
        , ylab='logloss'
        , main=paste('basecases and Chipman forests wrt dissimilarities\n(',50*length(files),' simulations)', sep='')
       # , main=paste(attr(et,'metric'), 'metric ,', num.clusters , 'clusters ,', sizeSF, 'trees in forest\nfor non-default')
)
abline(median(et02$random),0,col='red')
abline(median(et02$default),0, col='green')

#names(et)
et02[,c(2:5,7,9,11)] %>%
  summarize_all(function(x) c(mean(x), sd(x))) %>% 
  t %>% xtable -> et.xt

#et.xt
digits(et.xt) <- 4
et.xt

et02 %>% 
  select(default, random, 'not diverse', metrices) %>%
  summarize_all(function(x) c(median(x), IQR(x)))

