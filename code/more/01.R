# build Meiner forest for multiple parameters
# This forms a forest list, compare to default forest,
# select Meiner forest closest to default forest.

library(ranger)
library(dplyr)
library(caret)

calc_d0 <- function(fL , defRF , pow=1){
  #' calculates d0 dissimilarity (based on variable importance) for a list of forests
  #' comparing to a given ranger random forest
  
  #' @param fL forest list , a list of forests
  #' @param defF a ranger random forest 
  
  irg <- importance(defRF) # importance vector of the default forest
  
  # calculate weights from the variable importance vi
  wgt <- (irg^pow)/sum(irg^pow)
  
  # calculate variable importance for all forests in the forest list
  sapply(fL, importance) %>% t  -> vImp
  
  # calculate difference in variable importance for each variable
  scImp <- lapply(1:length(irg),function(i){abs(vImp[,i]-irg[i])*wgt[i]})
  scImp <- simplify2array(scImp)
  
  scImp %>% apply(1,sum) %>% return
}

calc_d2 <- function(fL, tD , defRF){
  # predictions for the subforests
  ppList <-  lapply(fL, predForest, data=tD)
  
  # mean distance (mae) to predictions by full forest
  f1 <- function(vec,myData){ (vec - predForest(rg$forest,myData) )  %>% abs %>% mean}
  
  lapply(ppList,f1,myData=data.train) %>% unlist %>% return
}

# d0 : 0, 0.1, 0.2
# d1 : 0.1 , 0.2, 0.3, 0.4, 0.5
# d2 : 0.2, 0.3, 0.4
# sb : 0.1 , 0.2, 0.3

# maybe add parameters in between these parameter values to increase the number of forests to select from

# create the forest list, then pass it to code in cleveland04

# parameters per metric
# named list , names must be names of metrices
L1 <- list( 'd0'=c(0, 0.1, 0.2)
            , 'd1'=c(0.1 , 0.2, 0.3, 0.4, 0.5)
            , 'd2'=c(0.2, 0.3, 0.4)
            , 'sb'=c(0.1 , 0.2, 0.3))
# , 'd0'=c(0, 0.05,0.1, 0.15, 0.2)
# , 'd1'= c(0.1 , 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 , 0.5)

folder <- 'data/nursery02'
files <- list.files(folder)[1:2]
# for(file in files){
  file <- files[[1]]
  # load doc
  load(paste(folder,file,sep='/')) # loads doc , info
  
  # nBs <- info$nBs
  nBs <-  length(doc)
  
  #for(i in 1:nBs){
  i <- 1
    DM <- doc[[i]]$DM
    rg <- doc[[i]]$rg
    forest <-  rg$forest
    fL <-  list() # forestList
    
    LL <-  calc_LL(forest, Swiss)
    ct.a <- 1
    for(metric in names(L1)){
      for(j in 1:length(L1[[metric]])){
        res1 <-  grow_meinForest(dm=DM[[metric]]
                                , LL=LL
                                , parameter=list('cutoff'=L1[[metric]][j], 'sizeSF'=500))
        res1$parameter$'metric' <-  metric
        fL[[ct.a]] <- res1
        ct.a <- ct.a+1
        }
      }
    calc_d2(fL=fL , tD = Swiss, defRF = doc[[i]]$rg) %>% # calc_d2 also needs data to calculate each forest's predictions and compare them 
      which.min -> bestForest.idx
    # fL[[bestForest.idx]] -> bestForest
    fL.best[ct.b, ] <- c(logloss(fL[bestForest.idx]) # same as LL[bestForest.idx] ??
                           , metric
                           , L1[[metric]][bestForest.idx])
 # } # for i
#} # for file

pf1 <- function(model, data){
      predict(object=model , data=data)$predictions[,2]
    }

my_metric <-  function(actual, predicted){
  
  print(actual[1:10])
  print(predicted[1:10])
  corrected_predicted <- ifelse(actual=='Yes',predicted,1-predicted) # problematic when this returns 0
  
  correctedpp <- winsorize_probs(corrected_predicted) 
  
  return( -mean(log(corrected_predicted)))
}

data <- Swiss
#data$CAD <- as.integer(data$CAD)-1
library(flashlight)
# create flashlight object
fl1 <- flashlight(model=forest 
           , data=Swiss[,1:11] 
           , y='CAD' 
           , label='ols'
           , predict_function = pf1
           , metric = list('my_metric'=metric1)
           )
light_importance(fl1,type='permutation')
    # pass a prediction function
    # calc variable importance

my_metric(Swiss$CAD, pf1(forest,Swiss))

metric1 <- function (actual, predicted, w = NULL, ...) 
{
  1
}

metric1(as.numeric(Swiss$CAD)-1 , pf1(forest,Swiss))
