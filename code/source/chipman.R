# test the function calc_chipForest 
# seems to be not working on the sb metric, 03.hpo does not give results for sb


calc_chipForest_1 <- function(dm, forest , oLL, pa){ # pa : parameter alpha
  #print(dim(dm))
  #print(pa)
  
  selectCentralTree <-  !TRUE
  selectBestTree <-  TRUE
  assertthat::assert_that(!selectCentralTree==selectBestTree, msg='Either select best or central tree per cluster, not both, not either.')
  
  # level of representation , remove the diagonal, it has only 0 values
  alpha <- quantile(dm[upper.tri(dm)],pa$cutoff) # diagonal with all values 0 not included
  
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
  N <-  min(nrow(dm2), sizeSF)
  
  if(alpha<max(dm2)){
    I <- 2
    while((I<N) && !represented(I)){
      I <- I+1
    }
  }else{I <- 1}
  # exits with smallest I with represented(I) TRUE
  
  R <-  1:I
  kOpt <-  NA
  
  # dissimilarity matrix needed only for R
  dm2 <-  dm2[1:I,1:I]
  # which 2 trees are meant in dm2[1,2] ?? the best and the second best , indexed oLL[1], oLL[2]
  
  # dip test for all trees in R (dense representing forest)
  lapply(1:I , function(i) dip.test(dm2[i,-i], simulate=T)$p.value) %>% 
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
           } ) %>% unlist -> col1
    col1 %>% which.max +  1 -> kOpt # which liefert eine Zahl in 1 bis I-2, das muss umgerechnet werden auf die ursprünglichen range , hier 2:(I-1)
    #col1 %>% which.max + min(I-1,10) - 1 -> kOpt # umrechnen von 1,...min(I-1,100) - min(I-1,10) +1 auf min(I-1,10) ... min(I-1,100)
    
    hc.clus <- cutree(hc, k = kOpt)
   
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
    
  }
  
  return(list(calcLogloss( subforest(forest, R ), data.test)
              , I
              , ifelse(multimod,kOpt,I) # when multimodal, (opt.) number of clusters is the size
              #, kOpt # produces NA for unimodal R , allows to count unimodal vs multimodal cases
  )
  )
}

calc_chipForest_2 <- function(dm , forest, oLL, parameter){
  #' calculates the Chipman forest by adding diverse trees , Chipman 2
  #' 
  #' @param parameter is a list with items cutoff and sizeSF 
  #'
  #' returns its logloss on data.test (from higher environment)
  #' returns the rank of the last inspected tree (or the index +1 , like 6 when collecting the first 5 trees)
  #' returns the number of trees in the Chipman forest (enforced to be parameters$sizeSF or smaller)
  
  chipForest <- oLL[1] # start chipForest with the best tree
  # loDiv <- quantile(dm,parameter$cutoff)  # old
  loDiv<- quantile(dm[upper.tri(dm)],parameter$cutoff)

  #print(paste('level of diversity=',loDiv))
  
  # go through trees (ordered by performance) until sizeSF trees are collected / selected
  for(I in 2:length(oLL)){
    if(length(chipForest)==parameter$sizeSF){
      I<-I-1
      break
    }else{
      trindx <- oLL[I]
      # if new tree far from (current) chipForest , add it to chipForest
      # if cutoff =0 then each tree is added with certainty -> same as unimodal!
      if(min(dm[chipForest,trindx]) > loDiv){ # if adding trindex to R keeps it diverse: then add it
        chipForest <- c(chipForest, trindx)
      }
    }
  } # for exits with j the first index after all trees have been collected
  return(list(calcLogloss( subforest(forest, chipForest ), data.test)
              , I
              , length(chipForest)))
}

calc_meinForest <- function(dm , forest, LL, parameter){
  #' calculates the Meiner forest
  #' 
  #' @param parameter is a list with items cutoff and sizeSF 
  #'
  #' returns its logloss on data.test (from higher environment)
  #' returns the rank of the last inspected tree (or the index +1 , like 6 when collecting the first 5 trees)
  #' returns the number of trees in the Meiner forest (enforced to be parameters$sizeSF or smaller)
  
  oLL <- order(LL)
  #print(parameter)
  meinForest <- oLL[1] # start chipForest with the best tree
  
  alpha <- quantile(dm[upper.tri(dm)], parameter$cutoff) # level of representation
  #print(cutoff)
  
  # go through trees (ordered by performance) until sizeSF trees are collected / selected
  for(I in 2:length(oLL)){
    if(length(meinForest)==parameter$sizeSF){
      I <- I-1
      break
    }else{
      trindx <- oLL[I]
      # if new tree far from (current) meinForest , add the best tree (lowest logloss) to the meinForest that represents trindx.
      # if cutoff =0 then each tree is added with certainty -> same as unimodal!
      
      if(min(dm[meinForest,trindx]) > alpha){ # if meinForest does not represent trindx
        allCandTrees <- base::setdiff(oLL[1:I], meinForest) # all candidate trees
        allCandTrees %>% # low logloss, trees not yet selected into meinForest
          dm[.,trindx] %>%
          (function(x) which(x<=alpha)) %>% # a mask , for representation
          allCandTrees[.] -> closeCandTrees # candidates close to trindex
        # if this is empty , add trindex!
        if(length(closeCandTrees)==0){
          meinForest <- c(meinForest,trindx) # if nothing better can be found : add trindx
          }else{
            closeCandTrees %>%
            LL[.] %>%
            which.min %>%
            closeCandTrees[.] %>%
            c(meinForest) -> meinForest
          }
        #print(I)
        #print(round(oLL[meinForest],4)[1:10])
        #print(round(oLL[order(oLL[-meinForest])],4)[1:10])
        #write.csv(c(I,round(oLL[meinForest],4)), file="./log.csv")
        #write.csv(c(I,round(oLL[meinForest],4)), file="./log.csv")
        # alternative to calc_chipForest :
        # chipForest <- c(chipForest, trindx)
      }
    }
  } # for exits with j the first index after all trees have been collected
  return(list(calcLogloss( subforest(forest, meinForest ), data.test)
              , I
              , length(meinForest)))
}


calc_LL <- function(forest , data){
  #' calculate ordered logloss for trees in forest,
  #' logloss evaluated on predictions on on data
  
  pp <- predict(forest 
                , data=data 
                , predict.all = T)$predictions[,2,]
  
  lapply(1:forest$num.trees
         , function(k){ 
           pp[,k] %>% 
             calcLogloss2( df=data ) %>% 
             unlist
         }) %>% 
    unlist -> LL
  
  return(LL) # logloss for tree indices order as in forest
}

calc_oLL <- function(forest , data){
  #' calculate ordered logloss for trees in forest,
  #' logloss evaluated on predictions on on data
  
  pp <- predict(forest 
                , data=data 
                , predict.all = T)$predictions[,2,]
  
  lapply(1:forest$num.trees
         , function(k){ 
           pp[,k] %>% 
             calcLogloss2( df=data ) %>% 
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
  
  return(order(LL)) # tree indices , ordered by logloss on OOB
}


calc_LL_for_selection <- function(doc, parameter){
  #' gets arguments and calculates arguments for calc_chipForest_2 from doc
  #' then runs loops over calc_chipForest_2 for all metrices and all simulations in doc
  #' 
  #' parameter is (only and directly) passed to calc_chipForest_2 
  
  nBs <- length(doc)
  #nBs <- 5
  
  evalT <- matrix(NA,nrow=nBs,ncol=12) # table of evaluations
  
  ct <- 1
  for(i in 1:nBs){
    #print(paste(i/nBs , Sys.time()))
    #if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
    
    DM <-  doc[[i]]$DM
    
    forest<-doc[[i]]$rg$forest
    
    data.train <- doc[[i]]$`bootstapped training data`
    OOB <-  base::setdiff(1:nrow(Cleve), unique(doc[[i]]$resample))
    data.set.val <- Cleve[OOB,] # goes back to original Cleveland data
    
    LL <- calc_LL(forest, data.set.val)
    oLL <- order(LL)
    
    metrices <- names(DM)
    res <- list() # result from following loop
    
    # parameter will be used as a named list
    # with a sub-list for each metric , 
    # like this : parameter <- list('d0'=list('cutoff'=0.4, 'sizeSF'=500), 'd1'= list('cutoff'=0.3, 'sizeSF'=500) ...))
    
    # if parameter is passed as only a single list
    # parameter = list('cutoff'=0.4, 'sizeSF'=500)
    # then this list is repeated for each metric to become of the required format
    
    if(!is.null(parameter)){ # no parameter if we use Chipman 1 enforce
    if(all(c('cutoff','sizeSF') %in% names(parameter))){
      P1 <-  list()
      for(metric in metrices){
        P1[[metric]] <- parameter
      }
      parameter <- P1
    }
    }
    
    # check for required format of named list with same names as the dissimilarity matrices
    assertthat::assert_that(all(names(parameter) == names(DM))
                            , msg = 'names for parameter and names for dissimilarity matrices do not match.')
    
    for(metric in metrices){
      #print(parameter[[metric]])
      DM[[metric]] %>%
       calc_chipForest_2(forest, oLL, parameter[[metric]]) -> res[[metric]]
       #calc_meinForest(forest, LL, parameter[[metric]]) -> res[[metric]]
        #calc_chipForest_1(forest, oLL, parameter[[metric]]) -> res[[metric]]
      #calc_chipForest_1_enforce(oLL=oLL, forest=forest) -> res[[metric]] # needs no parameter
    }
    evalT[i,] <- unlist(res)
  }
  
  evalT <-  data.frame(evalT)
  names(evalT) <- c(paste(rep(c('LL.test.chip.','I.','size.'),4),rep(c('d0','d1','d2','sb'),each=3),sep=''))
  return(evalT)
}
# returns 12 columns

"
# load arguments (dm , forest, oLL, parameter) to feed the function
{load('data/nursery/nursery01_10_50x500.rda')
metric <- 'd1'
DM <- doc[[1]]$DM
dm <-  doc[[1]]$DM[[metric]]

data.train <- doc[[1]]$`bootstapped training data`
OOB <-  base::setdiff(1:nrow(Cleve), unique(doc[[1]]$resample))
data.set.val <- Cleve[OOB,] # goes back to original Cleveland data

forest<-doc[[1]]$rg$forest

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

oLL <- order(LL) # tree indices , ordered by logloss on OOB

parameter <-  list('cutoff'=0.5 , 'sizeSF'=500)
pa <- parameter
}

## call function
calc_chipForest_1(dm , forest, oLL, pa)

mean(dm[1,])
summary(dm[1,])
quantile(dm[1,],0.5)
quantile(dm,0.5)

# 1)
dm.test <- matrix(0.4+rnorm(25,0,0.05),5,5)
# 2)
dm.test <- matrix(rnorm(25,0,1),5,5)
dm.test <- 3+dm.test
all(dm.test>0) # check
#
diag(dm.test) <- 0
dm.test[1,2] <- 0.1
dm.test[2,1] <- 0.1
dm.test[3,4] <- 0.1
dm.test[4,3] <- 0.1
dm2 <-  dm.test
LL.test <- c(0.5,0.7,0.45,0.8,0.6)
LL <-  LL.test

dip.test(c(-0.1,-0.2,0,2,1.5,1.6,1.7)) # very multimodal
dip.test(c(-4,-4.1,-4.21,-4.11,-0.1,-0.2,0,2,1.5,1.6,1.7)) # very multimodal
# construct as multimodal
data <- c(rnorm(5,-4,0.5) , rnorm(7,0,0.5) , rnorm(9,4,0.5))
plot(density(data))
plot(data)
dip.test(data) # recognise as multimodal, small p.value
#### overlapping , unimodlity rejected
data <- c(rnorm(50,-2,0.5) , rnorm(70,0,0.5) , rnorm(90,2,0.5))
plot(density(data))
plot(data)
dip.test(data)
#### more overlapping , unimodality not rejected
data <- c(rnorm(50,-1.7,0.5) , rnorm(70,0,0.5) , rnorm(90,1.7,0.5))
plot(density(data))
plot(data)
dip.test(data)

#### more overlapping and more (!) data (!), unimodality rejected
data <- c(rnorm(150,-1.7,0.5) , rnorm(150,0,0.5) , rnorm(150,1.7,0.5))
plot(density(data))
plot(data)
dip.test(data)
"