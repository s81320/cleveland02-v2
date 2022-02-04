# optimal number of trees in ranger random forest by OOB error 

rm(list=ls())

library(ranger)
library(dplyr)

source('code/source/prep.R') # calcLogloss
load('data/data_SupMat4.rda')

Cleve$CAD_fac <- NULL

# nt <- seq(5, 5000, 5)
nt <- seq(50, 5000, 100)
length(nt)
oob_error <- matrix(0,5, length(nt))

#### build only one forest of maximal size #####################################
#### classification ############################################################

oob_error <- matrix(0,5, length(nt))

for(j in 1:5){

# one large fixed forest
rg <- ranger(CAD ~ .
             , data=Cleve
             , num.trees = max(nt)
             , write.forest = TRUE
             , replace=FALSE
             , probability=FALSE
             , keep.inbag=TRUE)

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

  preds1 <- preds[sometimes_oob_idx,1:nt[i]] # reduce to current set of trees , current forest size
  preds1 %>% dim
  
  dim(oob_idx[sometimes_oob_idx,1:nt[i]])
  vote <- rowMeans(oob_idx[sometimes_oob_idx,1:nt[i]]*preds1, na.rm = TRUE)
  vote <- ifelse(vote>1.5,'No','Yes')
  vote <- factor(vote, levels=c('No','Yes'))
  
  (vote==Cleve[sometimes_oob_idx,]$CAD) %>% 
    which %>% 
    length/length(sometimes_oob_idx) -> oob_error[j,i]
}
}


#### build only one forest of maximal size #####################################
#### probability estimation ####################################################

oob_error <- matrix(0,20, length(nt))

for(j in 1:20){
  set.seed(2*j)
rg <- ranger(CAD ~ .
             , data=Cleve
             , num.trees = max(nt)
             , write.forest = TRUE
             , replace=FALSE
             , probability=TRUE
             , keep.inbag=TRUE)

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
  correctedpp <- ifelse(Cleve[sometimes_oob_idx,]$CAD=='Yes',pp,1-pp) # problematic when this returns 0
  correctedpp <- winsorize_probs(correctedpp) 
  
  oob_error[j,i] <- -mean(log(correctedpp))
}
}

# tails , plots
par(mar=c(5,4,1,1) + 0.1)
for(r in list(1:5,10:30,50:length(nt))){
restrict <- r
matplot(nt[restrict]
        ,t(oob_error[,restrict])
        , type='l'
        , lty=1
        #, main=paste('Cleveland data : OOB error over size of forest\n', rg$treetype)
        , sub=paste('5 forests, each sequentially increasing, starting with ',min(nt[restrict]),' trees', sep='')
        , ylab=ifelse(rg$treetype=='Classification','classification error (1-accuracy)','logloss')
        , xlab='number of trees in forest')
legend('topright',legend=1:5, pch='-', col=1:5, cex=0.8)
}

par(mar=c(5,4,2,2) + 0.1)


dim(oob_error)
boxplot(oob_error)
# table
