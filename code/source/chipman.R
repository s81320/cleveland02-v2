# test the function calc_chipForest 
# seems to be not working on the sb metric, 03.hpo does not give results for sb


calc_chipForest_2 <- function(dm , forest, oLL, parameter){
  #' calculates the Chipman forest by adding diverse trees , Chipman 2
  #' 
  #' @param parameter is a list with items cutoff and sizeSF 
  #'
  #' returns its logloss on data.test (from higher environment)
  #' returns the rank of the last inspected tree (or the index +1 , like 6 when collecting the first 5 trees)
  #' returns the number of trees in the Chipman forest (enforced to be parameters$sizeSF or smaller)
  
  chipForest <- oLL[1] # start chipForest with the best tree
  loDiv <- quantile(dm,parameter$cutoff) 
  #print(paste('level of diversity=',loDiv))
  
  # go through trees (ordered by performance) until sizeSF trees are collected / selected
  for(I in 2:length(oLL)){
    if(length(chipForest)==parameter$sizeSF){
      break
    }else{
      trindx <- oLL[I]
      # if new tree far from (current) chipForest , add it to chipForest
      # if cutoff =0 then each tree is added with certainty -> same as unimodal!
      if(min(dm[chipForest,trindx]) >= loDiv){ 
        chipForest <- c(chipForest, trindx)
      }
    }
  } # for exits with j the first index after all trees have been collected
  return(list(calcLogloss( subforest(forest, chipForest ), data.test)
              , I
              , length(chipForest)))
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
    
    oLL <-  calc_oLL(forest, data.set.val)
    
    metrices <- names(DM)
    res <- list()
    
    # parameter will be used as a named list
    # with a sub-list for each metric , 
    # like this : parameter <- list('d0'=list('cutoff'=0.4, 'sizeSF'=500), 'd1'= list('cutoff'=0.3, 'sizeSF'=500) ...))
    
    # if parameter is passed as only a single list
    # parameter = list('cutoff'=0.4, 'sizeSF'=500)
    # then this list is repeated for each metric to become of the required format
    
    if(all(c('cutoff','sizeSF') %in% names(parameter))){
      P1 <-  list()
      for(metric in metrices){
        P1[[metric]] <- parameter
      }
      parameter <- P1
    }
    
    # check for required format of named list with same names as the dissimilarity matrices
    assertthat::assert_that(all(names(parameter) == names(DM))
                            , msg = 'names for parameter and names for dissimilarity matrices do not match.')
    
    for(metric in metrices){
      DM[[metric]] %>%
        calc_chipForest_2(forest, oLL, parameter[[metric]]) -> res[[metric]]
      #calc_meinForest(forest, oLL, parameter) -> res[[metric]]
    }
    evalT[i,] <- unlist(res)
  }
  
  evalT <-  data.frame(evalT)
  names(evalT) <- c(paste(rep(c('LL.test.chip.','I.','size.'),4),rep(c('d0','d1','d2','sb'),each=3),sep=''))
  return(evalT)
}


# load arguments (dm , forest, oLL, parameter) to feed the function
"{load('data/nursery/nursery01_10_50x500.rda')
metric <- 'sb'
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

parameter <-  list('cutoff'=0.9 , 'sizeSF'=40)
}

## call function
calc_chipForest_2(dm , forest, oLL, parameter)

mean(dm[1,])
summary(dm[1,])
quantile(dm[1,],0.5)
quantile(dm,0.5)
"
