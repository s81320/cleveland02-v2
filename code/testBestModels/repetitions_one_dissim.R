rm(list=ls())

library(ranger)
library(xtable)
library(caret)
library(dplyr)
library(diptest)
library(cluster)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection
source('code/source/distance-matrices.R')
source('code/source/prep.R')

#### training set
{data.set.name <- 'Cleve'
  data.train <- get(data.set.name)[,1:11]
  attr(data.train,'data.set.name') <- data.set.name
  data.set.name <- NULL}
# attribute data.set.name will be used in documentation (end of script)


#### test set
{data.set.name <- 'Hung'
  data.test <- get(data.set.name)[,1:11]
  attr(data.test,'data.set.name') <- data.set.name
  data.set.name <- NULL}
# attribute data.set.name will be used in documentation 
# end of script, returned list of function f1

################################################################################
#### testing Chipman and Meiner forests with parameters from hpo 
################################################################################


f3_one_dissimilarity <- function(method, metric, parameter, nLoops=100){
  #' build forest for given method and evaluate its performance
  #' similar to f1 but:
  #' 1 metric (=> one dissimilarity matrix) and several methods and parameters
  #' 
  #' @param method : list of methods chip1, chip2, meiner, chip1_simplified
  #' @param metric : one(!) of d0, d1, d2, sb
  #' @param parameter list of lists with cutoff, sizeSF (500 indicating original Chipman methods, smaller: stopped algorithms to return forests of (at most the) indicated size). chip1 and chip1_simplified needs selection method best or central. Same length as method
  #' @param nLoops number of repetitions / loops to run. default is 100.
  #' uses data.train from parent environment to built the random forest
  #' uses attribute data.set.name of data sets data.train and data.test from parent environment to document results
  
  # parameter may be a list of parameters, 1 for each method. Parameter and method then have the same length
  # parameter may be a list of lists, with the top-layer list of same length as method
  # parameter and method should have the same length

  assertthat::assert_that(length(method)==length(parameter),
                          msg='parameter and method have to have the same length')

  doc <- matrix(rep(NA,4*nLoops*length(method)), ncol=5)
  
  for(i in 1:nLoops){
    if(i %% 10 == 0){ print(paste(i/nLoops, Sys.time())) }
    rg <- ranger(CAD~.
                 , data = data.train # from parent environment
                 , num.trees = 500
                 , replace = F # T required for predictions ? Of course we want predictions!
                 , mtry= 3 # default : 3
                 #, importance = 'impurity'
                 , probability = T 
                 , min.node.size = 13 
                 , keep.inbag = T # needed to generate OOB obs as validation data
    )
    
    forest <- rg$forest
    
    dm <- createDM(forest=forest , type=metric , dft=data.train)
    
    # use OOB observations for each tree to calculate the tree's logloss
    # prediction probabilities : pp
    predict(forest
            , data=data.train
            , predict.all = T)$predictions[,2,] -> pp # dim 303 x 500
    
    # List of Loglosses : LL
    lapply(1:forest$num.trees
           , function(t){ 
             OOB <- which(rg$inbag.counts[[t]]==0)
             pp[OOB,t] %>% 
               calcLogloss2( df=data.train[OOB,] ) %>% 
               unlist
           }
    ) %>% 
      unlist -> LL
    
    
    mthd.ct <- 1
    for(mthd in method){
      
      if(mthd=='chip1'){
      grow_chipForest_1(dm = dm
                        , oLL= order(LL)
                        #, oLL = calc_oLL(forest=forest, data=data.val)
                        , parameter = parameter[[mthd.ct]]) -> cf1 
      #print(cf1$forest)
      
      sz <- length(cf1$forest)
      
      c(calcLogloss(forest, data.test)
        , calcLogloss(subforest(forest, 1:sz), data.test)
        , mthd
        , sz
        , calcLogloss(subforest(forest,cf1$forest), data.test)
        ) -> doc[i,] 
      }
    
      if(method=='chip1_simplified'){
        grow_chipForest_1_simplified(dm = dm
                                   , oLL= order(LL)
                                   #, oLL = calc_oLL(forest=forest, data=data.val)
                                   , parameter = parameter[[mthd.ct]]) -> cf1 
      #print(cf1$forest)
      
      sz <- length(cf1$forest)
      
      c(calcLogloss(forest, data.test)
        , calcLogloss(subforest(forest, 1:sz), data.test)
        , mthd
        , sz
        , calcLogloss(subforest(forest,cf1$forest), data.test)
        ) -> doc[i,] 
    }
    
      if(method =='meiner'){
        mf <- grow_meinForest(dm=dm
                            , LL=LL
                            , parameter=parameter[[mthd.ct]]) # meiner forest
        
        sz <- ifelse(parameter[[mthd.ct]]$sizeSF==500,length(mf$forest),parameter[[mthd.ct]]$sizeSF) # size
        # why so complicated ?? 7.7.2022
        
        doc[i,] <- c( calcLogloss(forest, data.test) # full forest
                      , calcLogloss(subforest(forest,sample(1:forest$num.trees, sz)), data.test) # make the regular small forest random
                      # , calcLogloss(subforest(forest, 1:sz), data.test) # logloss for regular small forest of same size as (selected sub-) forest
                      , mthd
                      , sz # size of (selected sub-) forest
                      , calcLogloss(subforest(forest, mf$forest), data.test) # logloss for (selecte sub-) forest
                      )
    }
    
      if(method=='chip2'){
      cf2 <- grow_chipForest_2(dm=dm
                               , oLL= order(LL)
                               , parameter=parameter[[mthd.ct]]) # chipman 2 forest
      
      sz <- ifelse(parameter[[mthd.ct]]$sizeSF==500,length(cf2$forest),parameter[[mthd.ct]]$sizeSF) # size
      
      doc[i,] <- c( calcLogloss(forest, data.test) # full forest
                    # , calcLogloss(subforest(forest, 1:sz), data.test) # logloss for regular small forest of same size as (selected sub-) forest
                    , calcLogloss(subforest(forest,sample(1:forest$num.trees, sz)), data.test) # make the regular small forest random
                    , mthd
                    , sz # size of (selected sub-) forest
                    , calcLogloss(subforest(forest, cf2$forest), data.test) # logloss for (selecte sub-) forest
                    )
      }
      
      mthd.ct <- mthd.ct+1
    } # end of for all methods
  }
  
  doc <- data.frame(doc)
  names(doc) <- c('logloss.full.forest','logloss.base.case', method,'size',paste('logloss',method,sep='.'))
 
  return(list('res'=doc, 'call'=list(method=method
                                     , metric=metric
                                     , parameter=parameter
                                     , data.train=attr(data.train,'data.set.name')
                                     , data.test=attr(data.test,'data.set.name'))))
}

# calling f1 on the models in test.models 
# keeping results in list results
results <-  list()
# i <- 3
for(i in 1:length(test.models)){
  print(i)
  print(Sys.time())
  results[[i]] <- f1(method=test.models[[i]]$method
                     , metric=test.models[[i]]$metric 
                     , parameter=test.models[[i]]$parameter
                     , nLoops=100)
}