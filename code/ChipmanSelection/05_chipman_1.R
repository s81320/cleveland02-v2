# 23.2.2022
# code/chipmanSelection/05_chipman_1.R
# representation parameter alpha

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

sizeSF <-  5
# parameter.alpha <-  list('d0'=0.3, 'd1'=0.6 , 'd2'=0.45 , 'sb'=0.4) # I around 150 , 75% above 100
parameter.alpha <-  list('d0'=0.4, 'd1'=0.8 , 'd2'=0.6 , 'sb'=0.55) # I around 50 , 75% below 100

#calc_chipForest_1 <- function(dm, oLL, forest, parameter.alpha, metric){
calc_chipForest_1 <- function(dm, oLL, forest, pa){ # pa : parameter alpha
  #print(dim(dm))
  #print(pa)
  
  alpha <- quantile(dm,pa)
  #print(alpha)
  
  represented <-  function(I){
    #' for I >1 , 
    #' representation by one tree alone it only at level of largest dissimilarity
    #' no representation by almost all trees (I==N-1) only if last tree, 
    #' index N, is far away from all others, further than cutoff
    #assertthat::assert_that(class(I)=='numeric', 'argument I has to be numeric') %>% print
    #assertthat::assert_that( (I >0) , 'I must be at least 1 and at most one less than trees in forest') %>% print
    
    if(I==1){ # max
      dm2[1,2:N] %>% max %>%
        (function(x){x<=alpha}) # cutoff.a  used to be cutoff , will become alpha
    }else{
      if(I==N-1){ # min
        dm2[1:(N-1),N] %>% min %>%
          (function(x){x<=alpha})
      }else{ # min and max
        dm2[1:I,(I+1):N] %>%
          apply(2,min) %>% # not working for I = 1 or I=N-1
          max %>%
          (function(x){x<=alpha})
      }
    }
  }
  
  # dm should be ordered
  dm2 <-  dm[oLL,oLL]
  #oLL[1]
  #LL[oLL]
  # N <- forest$num.trees 
  N <-  nrow(dm2)

  I <- 1

  if(alpha<max(dm2)){
    I <- 2
    while((I<N) && !represented(I)){
      I <- I+1
      }
  }else{I <- 1}
  # exits with smallest I with represented(I) TRUE
  #represented(I)
  #print(paste('representing with 1..', I))
  
  #m <- c( "average", "single", "complete", "ward")
  #names(m) <- c( "average", "single", "complete", "ward")
  # function to compute coefficient
  #ac <- function(x) {
  #  agnes(dm2[1:I,1:I], method = x, diss=T)$ac
  #}
  #hpo <-  lapply(m, ac)
  
  # make this a basecase: Large Forest, smaller than default. Better?
  #calcLogloss( subforest(forest, 1:I ), data.test) %>% print
  
  if(I>sizeSF){
    hc <- cluster::agnes(dm2[1:I,1:I], 
                         method="ward", # always optimal - whenever I did hpo
                         #method = m[which.max(hpo)] , 
                         diss=T)
    #pltree(hc, cex = 0.6, hang = -1, main = paste(m[which.max(hpo)], ',', metric))
    
    hc.clus <- cutree(hc, k = sizeSF)
    #hc.clus %>% table
    # cluster at a specified height , but we do not need that
    # ch6 <- cutree(as.hclust(hc), h = 0.3)
    # ch6 %>% table
    # optimal number of clusters?
  
    # from each cluster select the tree with the smallest logloss (on OOB data)
    lapply(1:sizeSF,function(i) which(hc.clus==i)[1]) %>% # we take the first, because indices are ordered, first is smallest
      unlist -> R # representing subset
    # if indices were not ordered, we'd do 
    #lapply(1:sizeSF,function(i) which(hc.clus==i) %>% min ) %>% unlist
    }else{
      R=1:I
      }
  
    # that is the smallest index under the oLL ordering (as used for dm2)
    # then go back to the original tree indices!
    #oLL # original indices
    #oLL[R]
    #LL[oLL[R]] # check that the loglosses are small
  
    return(list(calcLogloss( subforest(forest, R ), data.test)
                , calcLogloss( subforest(forest, 1:I ), data.test)
                ,I)
           )
}

calc_LL_for_selection <- function(doc){
  
  nBs <- length(doc)
  #nBs <- 5
  
  evalT <- matrix(NA,nrow=nBs,ncol=16) # table of evaluations
  
  ct <- 1
  for(i in 1:nBs){
    #print(paste(i/nBs , Sys.time()))
    if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
    
    DM <-  doc[[i]]$DM
    
    data.train <- doc[[i]]$`bootstapped training data`
    
    OOB <-  base::setdiff(1:nrow(Cleve), unique(doc[[i]]$resample))
    data.set.val <- Cleve[OOB,] # goes back to original Cleveland data
    
    forest<- doc[[i]]$rg$forest
    
    pp <- predict(forest 
                  , data=data.set.val 
                  , predict.all = T)$predictions[,2,]
    #pp <- simplify2array(pp, higher=FALSE)
    
    lapply(1:forest$num.trees
           , function(k){ 
             pp[,k] %>% 
               calcLogloss2( df=data.set.val ) %>% 
               unlist
           }) %>% 
      unlist -> LL
    
    oLL <- order(LL) # tree indices , ordered by logloss on OOB
    
    metrices <- names(DM)
    res <- list()
    
    for(metric in metrices){
      DM[[metric]] %>%
        calc_chipForest_1(oLL=oLL
                          , forest=forest
                          , pa=parameter.alpha[metric] %>% as.numeric) -> res[[metric]]
        #calc_chipForest_1(oLL, forest, unlist(parameter.alpha[metric]), metric) -> res[[metric]]
    }
    
    evalT[i,] <- c(i
                   , calcLogloss( forest, data.test)
                   , calcLogloss( subforest(forest, 1:sizeSF), data.test)
                   , calcLogloss( subforest(forest, oLL[1:sizeSF] ), data.test)
                   , unlist(res)
    )
  }
  
  evalT <-  data.frame(evalT)
  names(evalT) <- c('sim','LL.default', 'LL.random','LL.hp',
                    paste(rep(c('LL.chip.','LL.chip.I.','I.'),4),rep(c('d0','d1','d2','sb'),each=3),sep=''))
  return(evalT)
}

# to base the result on more bootstraps
folder <- 'data/nursery'
files <- list.files(folder)#[1:2]
# dir(folder)
collector <-  list()
ct <-  1 # counter for the above collector

for(file in files){
  # run loops over doc loaded from file
  load(paste(folder,file, sep='/'))
  collector[[ct]] <- calc_LL_for_selection(doc)
  ct <-  ct+1
}
et02 <- bind_rows(collector)

et02 %>% 
  select(starts_with('LL.')) %>%
  boxplot( ylab='logloss'
        , main=paste('basecases and Chipman forests wrt dissimilarities\n(',50*length(files),' simulations)', sep='')
       # , main=paste(attr(et,'metric'), 'metric ,', num.clusters , 'clusters ,', sizeSF, 'trees in forest\nfor non-default')
)
abline(median(et02$LL.random),0,col='red')
abline(median(et02$LL.default),0, col='green')

#names(et)
et02 %>% 
  select(starts_with('LL.')) %>%
  summarize_all(function(x) c(mean(x), sd(x))) %>% 
  t %>% xtable -> et.xt

#et.xt
digits(et.xt) <- 4
et.xt

et02 %>% 
  select(starts_with('I.')) %>%
  summarize_all(function(x) c(mean(x), sd(x))) %>% 
  t %>% xtable -> et.xt

#et.xt
digits(et.xt) <- 4
et.xt

et02 %>% 
  select(starts_with('I.')) %>%
  boxplot(main=paste('sizes of representing subforests, Chipman selection before clustering\n(sizeSF', sizeSF,', representation parameter', paste(unlist(parameter.alpha),collapse=',') ,')'))

# we need to know:
# do the small representing Chipman forests perform badly , and the large ones well?
plot(et02$LL.chip.I.d0~et02$I.d0)
lm(et02$LL.chip.I.d0~et02$I.d0) %>% summary
plot(et02$LL.chip.I.d1~et02$I.d1)
lm(et02$LL.chip.I.d1~et02$I.d1) %>% summary
plot(et02$LL.chip.I.d2~et02$I.d2)
lm(et02$LL.chip.I.d2~et02$I.d2) %>% summary
plot(et02$LL.chip.I.sb~et02$I.sb)
lm(et02$LL.chip.I.sb~et02$I.sb) %>% summary

