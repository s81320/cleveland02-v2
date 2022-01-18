# optimal number of trees in ranger random forest by OOB error 

rm(list=ls())

library(ranger)
library(dplyr)

source('code/source/prep.R') # calcLogloss
load('data/data_SupMat4.rda')

Cleve$CAD_fac <- NULL

plot1 <- function(restrict=1:rg$num.trees , ...){
  plot(x = nt[restrict]
       , y = oob_error[restrict]
       , col = "red"
       #, type = "b"
       , main=paste('Cleveland data : OOB error over size of forest\n', rg$treetype)
       , ylab=ifelse(rg$treetype=='Classification','classification error (1-accuracy)','logloss')
       , xlab='number of trees in forest'
       , ...
  )
}

###########################################################
#### Cleveland - probability estimation ###########################
###########################################################

nt <- seq(5, 1000, 5)
length(nt)
oob_error <- vector("numeric", length(nt))

#i <- 20
for(i in 1:length(nt)){
  rg <- ranger(CAD ~ .
               , data=Cleve
               , num.trees = nt[i] # 6 trees for i=2
               , write.forest = T
               , replace=F
               , probability=T
               , keep.inbag=T
               )
  
  # Identify OOB cases 
  oob_idx <- ifelse(simplify2array(rg$inbag.counts) == 0, TRUE, NA) # matrix , mask
  # some observations may be OOB to no tree, so we cannot get an OOB prediction for it
  sometimes_oob_idx <- which(rowMeans(oob_idx, na.rm =TRUE)==1) # rowMean oob_idx can only be 1 , is NA if there is no singel TRUE entry in the row
  c1 <- Cleve[sometimes_oob_idx,]
  # Average over trees for just OOB cases
  preds <- predict(rg
                   ,Cleve[sometimes_oob_idx,]
                   , predict.all = TRUE)$predictions[,2,] # probabilities for CAD=='Yes'
  pp <- rowMeans(oob_idx[sometimes_oob_idx,] *preds, na.rm = TRUE) 
  # pp contains NaN if some observations are not oob to any tree
  correctedpp <- ifelse(c1$CAD=='Yes',pp,1-pp) # problematic when this returns 0
  correctedpp <- winsorize_probs(correctedpp) 

  oob_error[i] <- -mean(log(correctedpp))
}

plot1(sub='new forest for each size', type='l')
plot1(1:20, sub='new forest for each size', type='b')
plot1(20:200, sub='new forest for each size', type='l')
#plot1(200:length(nt), sub='new forest for each size')

##########################################################
#### Cleveland - classification ###########################
##########################################################


oob_error <- vector("numeric", length(nt))

for(i in 1:length(nt)){
  rg <- ranger(CAD ~ .
               , data=Cleve[,-12]
               , num.trees = nt[i]
               , write.forest = TRUE
               , replace=FALSE
               , probability=FALSE
               , keep.inbag=TRUE)
  # Identify OOB cases 
  oob_idx <- ifelse(simplify2array(rg$inbag.counts) == 0, TRUE, NA) # matrix , mask
  # some observations may be OOB to no tree, so we cannot get an OOB prediction for it
  sometimes_oob_idx <- which(rowMeans(oob_idx, na.rm =TRUE)==1) # rowMean oob_idx can only be 1 , is NA if there is no singel TRUE entry in the row
  c1 <- Cleve[sometimes_oob_idx,]
  # Average over trees for just OOB cases
  preds <- predict(rg, c1, predict.all = TRUE)$predictions # values 1 or 2
  vote <- rowMeans(oob_idx[sometimes_oob_idx,] *preds, na.rm = TRUE) 
  vote <- ifelse(vote>1.5,'No','Yes')
  vote <- factor(vote, levels=c('No','Yes'))
  oob_error[i] <- (vote==c1$CAD) %>% which %>% length/length(sometimes_oob_idx)
}

plot1(sub='new forest for each size', type='l')
plot1(1:20, sub='new forest for each size', type='b')
plot1(20:200, sub='new forest for each size', type='l')
#plot1(200:length(nt), sub='new forest for each size')



#### build only one forest of maximal size #####################################
#### classification ############################################################

# one large fixed forest
rg <- ranger(CAD ~ .
             , data=Cleve
             , num.trees = max(nt)
             , write.forest = TRUE
             , replace=FALSE
             , probability=FALSE
             , keep.inbag=TRUE)

oob_error <- vector("numeric", length(nt))

# Identify OOB cases for all trees
oob_idx <- ifelse(simplify2array(rg$inbag.counts) == 0, TRUE, NA) # matrix , mask
# indexing : observation , tree

# predictions for all observations and all trees
preds <- predict(rg, Cleve, predict.all = TRUE)$predictions
dim(preds) # obs , trees

# build forest up , starting with nt[1] trees , ending with max(nt) trees
for(i in 1:length(nt)){
  
  # oob_idx indexing : obs , trees
  # reduce to trees of smaller forest of size nt[i]
  oob_idx_1 <- oob_idx[,1:nt[i]]
  
  # some observations may be OOB to no tree, so we cannot get an OOB prediction for it
  # is calculated inside the loop as with growing number of trees more and more observations are oob for some tree
  sometimes_oob_idx <- which(rowMeans(oob_idx_1, na.rm =TRUE)==1) # rowMean oob_idx can only be 1 , is NA if there is no singel TRUE entry in the row
  #c1 <- Cleve[sometimes_oob_idx,] # data remove observations that are never oob
  #c1 %>% dim
  # OOB error based on how many observations? should increase to max number
  # nrow(c1) %>% print
  
  # Average over trees for just OOB cases
  # calculate again as c1 may have changed (happens only early on for small forests ...)
  
  # maybe ok, but in general error was too low
 # values 1 or 2
  preds1 <- preds[sometimes_oob_idx,1:nt[i]] # reduce to current set of trees , current forest size
  preds1 %>% dim
  
  dim(oob_idx[sometimes_oob_idx,1:nt[i]])
  vote <- rowMeans(oob_idx[sometimes_oob_idx,1:nt[i]]*preds1, na.rm = TRUE)
  vote <- ifelse(vote>1.5,'No','Yes')
  vote <- factor(vote, levels=c('No','Yes'))
  oob_error[i] <- (vote==Cleve[sometimes_oob_idx,]$CAD) %>% which %>% length/length(sometimes_oob_idx)
}

plot1(sub='one forest sequentially increasing', type='l')
plot1(1:20, sub='one forest sequentially increasing', type='b')
plot1(20:200, sub='one forest sequentially increasing', type='l')
# plot1(200:length(nt), sub='one forest sequentially increasing')

#### build only one forest of maximal size #####################################
#### probability estimation ####################################################

# one large fixed forest
rg <- ranger(CAD ~ .
             , data=Cleve
             , num.trees = max(nt)
             , write.forest = TRUE
             , replace=FALSE
             , probability=TRUE
             , keep.inbag=TRUE)

oob_error <- vector("numeric", length(nt))

# Identify OOB cases for all trees
oob_idx <- ifelse(simplify2array(rg$inbag.counts) == 0, TRUE, NA) # matrix , mask
# indexing : observation , tree

# predictions for all observations and all trees
preds <- predict(rg, Cleve, predict.all = TRUE)$predictions[,2,]
dim(preds) # obs , trees

# build forest up , starting with nt[1] trees , ending with max(nt) trees
for(i in 1:length(nt)){
  
  # oob_idx indexing : obs , trees
  # reduce to trees of smaller forest of size nt[i]
  oob_idx_1 <- oob_idx[,1:nt[i]]
  
  # some observations may be OOB to no tree, so we cannot get an OOB prediction for it
  # is calculated inside the loop as with growing number of trees more and more observations are oob for some tree
  sometimes_oob_idx <- which(rowMeans(oob_idx_1, na.rm =TRUE)==1) # rowMean oob_idx can only be 1 , is NA if there is no singel TRUE entry in the row

  preds1 <- preds[sometimes_oob_idx,1:nt[i]] # reduce to current set of trees , current forest size
  preds1 %>% dim # indexing obs, prob for each class , tree
  pp <- rowMeans(oob_idx[sometimes_oob_idx,1:nt[i]] *preds1, na.rm = TRUE) 
  # pp contains NaN if some observations are not oob to any tree
  correctedpp <- ifelse(c1$CAD=='Yes',pp,1-pp) # problematic when this returns 0
  correctedpp <- winsorize_probs(correctedpp) 
  
  oob_error[i] <- -mean(log(correctedpp))
}

plot1(sub='one forest sequentially increasing', type='l')
plot1(1:20, sub='one forest sequentially increasing', type='b')
plot1(20:200, sub='one forest sequentially increasing', type='l')
#plot1(200:length(nt), sub='one forest sequentially increasing')

