#function winsorize copied from Calibration_Sup4_V01.R

print('sourcing p r e p .R : functions winsorize_probs , predForest(forest , data) , calclosloss(forest, data) , calclogloss2(predictedProbabilities , data) , f2(forest, data , vec)')

winsorize_probs <- function(x) {
  eps <- .Machine$double.eps
  pmax(pmin(x, 1 - eps), eps)
}

predForest <-  function(forest1, data) {
  #' get the probabilities predicted for data by forest
  #' 
  #' @param forest a ranger pe forest
  #' @param data to be predicted on
  
  if('treetype' %in% names(forest1)){
    if(forest1$treetype=='Probability estimation'){
      #print('doing great')
    }else{
      #print('Problem alert')
    }
  }else{
    #print('no treetype for forest1')
  }
  
  predict(forest1 
          , data=data 
  )$predictions[,2]
}

calcLogloss <-  function(forest, df){
  
  if('CAD' %in% names(df)){
    # we're good
  }else{
    # we're not good
  }
  
  pp <- predForest(forest,df)
  
  correctedpp <- ifelse(df$CAD=='Yes',pp,1-pp) # problematic when this returns 0
  
  correctedpp <- winsorize_probs(correctedpp) 
  
  return( -mean(log(correctedpp)))
}

calcLogloss2 <-  function(pp , df){
  #' calc logloss when you already have the predicted probabilities
  #' no need to calculate them again
  
  if('CAD' %in% names(df)){
    # we're good
  }else{
    # we're not good
  }
  
  # pp should be predictions on the df data
  if(nrow(df)==length(pp)){
    # we're fine
  }else{
    # error
  }
  
  correctedpp <- ifelse(df$CAD=='Yes',pp,1-pp) # problematic when this returns 0
  
  correctedpp <- winsorize_probs(correctedpp) 
  
  return( -mean(log(correctedpp)))
}

f2 <- function(forest2 , data, vec){
  #' mean absolute error for predictions by forest2 on data with vec the true values
  #' or: the mean difference of vec and the predictions by forest2 on data
  (predForest(forest1=forest2,data=data) - vec) %>% 
    abs %>% 
    mean
}
