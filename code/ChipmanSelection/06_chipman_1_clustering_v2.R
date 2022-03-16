# 23.2.2022
# code/chipmanSelection/05_chipman_1.R
# representation parameter alpha

rm(list=ls())

library(dplyr)
library(ranger)
library(cluster)
library(xtable)

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # calc_oLL

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

# set test data by name : VA, Swiss or Hung
data.test.name <-  'Swiss'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name

#parameter.alpha <-  list('d0'=0.111, 'd1'=0.465 , 'd2'=0.4125 , 'sb'=0.36) # I around 250
parameter.alpha <-  list('d0'=0.5, 'd1'=0.5 , 'd2'=0.5 , 'sb'=0.5)

calc_chipForest_1 <- function(dm, oLL, forest, pa){ # pa : parameter alpha
  #print(dim(dm))
  #print(pa)
  
  selectCentralTree <-  !TRUE
  selectBestTree <-  TRUE
  assertthat::assert_that(!selectCentralTree==selectBestTree, msg='Either select best or central tree per cluster, not both, not either.')
  
  # level of representation
  alpha <- quantile(dm,pa)
  #print(alpha)
  
  represented <-  function(I){
    #' testing for representation of best I trees (R= dense representing forest = {1,..,I}) at level alpha
    #' 
    #' uses dm2 and alpha from parent environment. i.e. function calc_chipForest_1
    #' for I >1 , 
    #' representation by one tree alone can happen only at level of largest dissimilarity
    #' no representation by almost all trees (I==N-1) only if last tree, 
    #' index N, is far away from all others, further than cutoff
    #assertthat::assert_that(class(I)=='numeric', 'argument I has to be numeric') %>% print
    #assertthat::assert_that( (I >0) , 'I must be at least 1 and at most one less than trees in forest') %>% print
    
    if(I==1){ # max
      dm2[1,2:N] %>% max %>%
        (function(x){x<=alpha}) 
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
  N <-  nrow(dm2)
  
  if(alpha<max(dm2)){
    I <- 2
    while((I<N) && !represented(I)){
      I <- I+1
    }
  }else{I <- 1}
  # exits with smallest I with represented(I) TRUE
  
  # dissimilarity matrix needed only for R
  dm2 <-  dm[oLL[1:I],oLL[1:I]] # careful :  dm[oLL[1:I],oLL[1:I]] != dm2[1:I,1:I]
  # which 2 trees are meant in dm2[1,2] ?? the best and the second best , indexed oLL[1], oLL[2]
  
  # dip test for all trees in R (dense representing forest)
  lapply(1:I , function(i) dip.test(dm2[i,-i])$p.value) %>%  # changed to 1:I , used to be 1:nrow(dm) , 8.3.2022
    unlist %>% 
    min %>% 
    (function(x) x<0.05) -> multimod # logical , TRUE if pvalue of at least one dip test less than 0.05
  
  # cluster only if unimodality is rejected
  if(multimod){ # if hypothesis of unimodality is rejected : cluster
  #print('multimod')
  hc <- cluster::agnes(dm2, 
                         method="ward", # always optimal - whenever I did hpo
                         diss=T)

  col1 <- list()
  # alle möglichen Anzahlen von clustern , nicht alle machen Sinn
  # lapply(min(I-1,10):min(I-1,100), # nur zwischen 10 und 100 cluster bzw representative Bäume
  lapply(2:(I-1), 
           function(k){
             cutree(hc, k = k) %>% 
               silhouette(dmatrix=dm2)%>%
               .[,'sil_width'] %>%
               mean
           } ) -> col1
  unlist(col1) %>% which.max +  1 -> kOpt # which liefert eine Zahl in 1 bis I-2, das muss umgerechnet werden auf die ursprünglichen range , hier 2:(I-1)
  #unlist(col1) %>% which.max + min(I-1,10) - 1 -> kOpt # umrechnen von 1,...min(I-1,100) - min(I-1,10) +1 auf min(I-1,10) ... min(I-1,100)
    
  hc.clus <- cutree(hc, k = kOpt)
  #print(paste(I, kOpt))
    #hc.clus %>% table
    # cluster at a specified height , but we do not need that
    # ch6 <- cutree(as.hclust(hc), h = 0.3)
    # ch6 %>% table
    # optimal number of clusters?
    
    # from each cluster select the tree with the smallest logloss (on OOB data)
  if(selectBestTree){
  lapply(1:kOpt,function(i) which(hc.clus==i)[1] %>% oLL[.]) %>% # we take the first, because indices are ordered, first is smallest
      # oLL to go back to original index
      unlist -> R # representing subset
    # if indices were not ordered, we'd do 
    #lapply(1:sizeSF,function(i) which(hc.clus==i) %>% min ) %>% unlist
  }
 
  
  if(selectCentralTree){
  lapply(1:kOpt,
         function(i){
           x <-  which(hc.clus==i) # trees of cluster i
           lapply(1:length(x), 
                  function(j){sum(dm[x[j],x])} %>% # sum of dissimilarities to all other trees in the same cluster (cluster i)
                    unlist ) %>% 
             which.min %>% # smallest sum of dissimilaritites
             x[.] %>% # zugehöriges Element in x
             oLL[.] # original index in default forest
         }
  ) %>% unlist -> R
  }
  
  # Alternative: cluster into fewer clusters and draw multiple trees from each cluser. Always the best (lowest logloss on OOB) possible
  # cluster into 5 clusters and draw 10 from each cluster
  # may retund NA if some clusters do not have 10 trees
  #if(I>sizeSF){
  #  hc <- cluster::agnes(dm2[1:I,1:I], 
  #                       method="ward", # always optimal - whenever I did hpo
  #                       #method = m[which.max(hpo)] , 
  #                       diss=T)
  #  hc.clus <- cutree(hc, k = 5)
  #  # from each cluster select the tree with the smallest logloss (on OOB data)
  #  lapply(1:sizeSF,function(i) which(hc.clus==i)[1:10]) %>% # we take the first, because indices are ordered, first is smallest
  #    unlist -> R # representing subset
  #  # print(R) # why is R of length 500 ??
  #  R <- R[!is.na(R)] # drop NA , happens when some clusters are too small
  #}else{
  #  R=1:I
  #}
  }else{ # if hypothesis of unimodality is not rejected : return unclustered representing subforest
    R <- 1:I
    kOpt <-  NA
    #print(paste('unclustered with I=',I))
  }
  
  return(list(calcLogloss( subforest(forest, R ), data.test)
              , I
              , kOpt
              )
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
    
    oLL <- calc_oLL(forest, data.set.val)
    
    metrices <- names(DM)
    res <- list()
    
    for(metric in metrices){
      #print(metric)
      DM[[metric]] %>%
        calc_chipForest_1(oLL=oLL
                          , forest=forest
                          , pa=parameter.alpha[metric] %>% as.numeric) -> res[[metric]]
      #calc_chipForest_1(oLL, forest, unlist(parameter.alpha[metric]), metric) -> res[[metric]]
    }
    
    evalT[i,] <- c(i
                   , calcLogloss( forest, data.test)
                   , calcLogloss( subforest(forest, 1:50), data.test)
                   , calcLogloss( subforest(forest, oLL[1:50] ), data.test)
                   , unlist(res)
    )
  }
  
  evalT <-  data.frame(evalT)
  names(evalT) <- c('sim','LL.default', 'LL.random','LL.hp',
                    paste(rep(c('LL.chip.','I.','sizeR.'),4),rep(c('d0','d1','d2','sb'),each=3),sep=''))
  return(evalT)
}

# to base the result on more bootstraps
folder <- 'data/nursery'
files <- list.files(folder)[1:2]
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

# too much information
#et02 %>% 
#  select(starts_with('LL.')) %>%
#  boxplot( ylab='logloss'
#        , main=paste('basecases and Chipman forests wrt dissimilarities\n(',50*length(files),' simulations)', sep='')
#       # , main=paste(attr(et,'metric'), 'metric ,', num.clusters , 'clusters ,', sizeSF, 'trees in forest\nfor non-default')
#)
#abline(median(et02$LL.random),0,col='red')
#abline(median(et02$LL.default),0, col='green')

# reduce information
par(mar=c(4,4,3,2)+0.1)
et02 %>% 
  select(starts_with('LL.chip.')) %>%
  select(!ends_with(paste('.I.',c('d0','d1','d2','sb'),sep=''))) %>%
  rename(c('d0'='LL.chip.d0','d1'='LL.chip.d1','d2'='LL.chip.d2','sb'='LL.chip.sb')) %>%
  boxplot( ylab='logloss'
           , xlab='dissimilarity for representation in Chipman I forest'
           , main=paste('Chipman I forests wrt dissimilarities\n(',50*length(files),' simulations)', sep='')
           # , main=paste(attr(et,'metric'), 'metric ,', num.clusters , 'clusters ,', sizeSF, 'trees in forest\nfor non-default')
  )
abline(median(et02$LL.random),0,col='red')
abline(median(et02$LL.default),0, col='green')
legend('topleft' 
       , legend=c('mean logloss default forest, size 500', 'mean logloss small forest, size 5')
       , cex =0.8
       , col=c('green','red')
       , pch = '-')

# count NA : unimodal I has sizeR = NA
# calculate the size of the final Chipman forest
# if multimodal it is the size of R , if unimodal it is I
metrices <- c('d0','d1','d2','sb')
for(metric in metrices){
  et02 %>% 
    select(ends_with(metric)) %>%
    select(starts_with(c('I.','size'))) %>%
    apply(1,min, na.rm=T) -> et02[,paste('sizeCh.',metric,sep='')] 
}


et02 %>% 
  summarize_all(function(x) c(mean(x, na.rm=T), sd(x, na.rm=T), sum(1*is.na(x)))) %>% 
  t %>% xtable -> et.xt
#et.xt
digits(et.xt) <- 4
et.xt

et02 %>% 
  summarize_all(function(x) c(mean(x, na.rm=T), sd(x, na.rm=T))) %>% 
  t %>% xtable -> et.xt
#et.xt
digits(et.xt) <- 4
et.xt

#names(et)
et02 %>% 
  select(starts_with('LL.')) %>%
  summarize_all(function(x) c(mean(x), sd(x), median(x), IQR(x))) %>% 
  t %>% xtable -> et.xt

#et.xt
digits(et.xt) <- 4
et.xt

et02 %>% 
  select(starts_with('sizeR.')) %>%
  summarize_all(function(x) c(mean(x), sd(x))) %>% 
  t %>% xtable -> et.xt

#et.xt
digits(et.xt) <- 4
et.xt

#sizes of representing subforests before clustering
et02 %>% 
  select(starts_with('I.')) %>%
  rename(c('d0'='I.d0','d1'='I.d1','d2'='I.d2','sb'='I.sb')) %>%
  boxplot(main=paste('Chipman I: sizes of representing dense subforests before clustering\n(representation parameter', paste(unlist(parameter.alpha),collapse=',') ,')')
          , xlab='dissimilarity for representation in Chipman I forest'
          , ylab=' number of trees')

#sizes after clustering
et02 %>% 
  select(starts_with('sizeR.')) %>%
  rename(c('d0'='sizeR.d0','d1'='sizeR.d1','d2'='sizeR.d2','sb'='sizeR.sb')) %>%
  boxplot(main=paste('Chipman I: sizes of final forest\n(representation parameter', paste(unlist(parameter.alpha),collapse=',') ,')')
          , xlab='dissimilarity for representation in Chipman I forest'
          , ylab=' number of trees')


metric <- 'd0'
plot(et02[,paste('LL.chip.', metric,sep='')]~ et02[,paste('sizeR.', metric , sep='')])
abline(lm(et02$LL.chip.d2~ et02$sizeR.d2)$coeff)
lm(et02$LL.chip.d2~ et02$sizeR.d2) %>% summary
 