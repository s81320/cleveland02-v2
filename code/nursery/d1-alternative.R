# new d1 formula

rm(list=ls())

source('code/source/distance-matrices.R')

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

library(ranger)
library(dplyr)
library(xtable)

set.seed(1)
rg <- ranger(CAD~.
                , data = Cleve[,1:11] 
                , num.trees = 5
                , replace = F 
                , mtry= 3 
                , importance = 'impurity'
                , probability = T 
                , min.node.size = 13 
)


treeInfo(rg,1)
plotTree1(rg,1)

nObs <- 10
pred.tn <- predict(rg, data=Cleve[1:nObs,], type='terminalNodes')
pred.tn$predictions

# pretty much the definition verbatim
d1.alt.pred <- function(pred1, pred2){
  nObs <-  length(pred1)
  assertthat::assert_that(length(pred1)==length(pred2), msg='No can do. Need predictions of same length.')
  
  # masks / indicator matrix for observations in the same node for predictions 1
  I1 <- outer(1:nObs , 1:nObs , (function(i,j) pred1[i]==pred1[j])) 
  I1 <- I1[upper.tri(I1)]
  # ... and predictions 2
  I2 <- outer(1:nObs , 1:nObs , (function(i,j) pred2[i]==pred2[j]))
  I2 <- I2[upper.tri(I2)]
  
  # do predictions agree in putting observations in the same terminal node or not?
  # which pairs are treated the differently in predictions 1 and 2
  (I1!=I2) %>% factor(levels=c(TRUE,FALSE)) %>% table %>% prop.table %>% .['TRUE']
}

d1.alt.tri <- function(forest,data,t1,t2){
  pred.tn <- predict(subforest(forest,c(t1,t2))
                     , data=data
                     , type='terminalNodes')$predictions
  d1.alt.pred(pred.tn[,1],pred.tn[,2])
}

d1.alt.tri(forest=rg$forest, data=Cleve[1:nObs,], 1, 2)
d1.alt.tri(rg$forest, Cleve[1:nObs,], 1, 3)
d1.alt.tri(rg$forest, Cleve[1:nObs,], 1, 4)
d1.alt.tri(rg$forest, Cleve[1:nObs,], 1, 5)

# not working ...
# outer(1:rg$num.trees,1:rg$num.trees, (function(i,j) d1.alt.tri(forest=rg$forest, data=Cleve[1:10,], i, j)))

createDMd1(forest=rg$forest, dft=Cleve[1:nObs,])

################################################################################
#### ALTERNATIVE 

# the definition, quite verbatim , working on the predictions (type terminalNodes)
d1.def.1 <- function(pred1, pred2){
  nObs <-  length(pred1)
  assertthat::assert_that(length(pred1)==length(pred2), msg='No can do. Need predictions of same length.')
  
  # masks / indicator matrix for observations in the same node for predictions 1
  I1 <- outer(1:nObs , 1:nObs , (function(i,j) pred1[i]==pred1[j])) 
  I1 <- I1[upper.tri(I1)]
  # ... and predictions 2
  I2 <- outer(1:nObs , 1:nObs , (function(i,j) pred2[i]==pred2[j]))
  I2 <- I2[upper.tri(I2)]
  
  # do predictions agree in putting observations in the same terminal node or not?
  # which pairs are treated the differently in predictions 1 and 2
  (I1!=I2) %>% factor(levels=c(TRUE,FALSE)) %>% table %>% prop.table %>% .['TRUE']
}

#f2 <- function(i,j) d1.def.1(preds[,i],preds[,j])%>% unname
#f2(1,2)

d1.def <-  function(rg , data){
  nObs <-  nrow(data)
  preds <- predict(rg, data, type='terminalNodes')$predictions
  outer(1:rg$num.trees 
        , 1:rg$num.trees 
        , Vectorize(function(i,j) d1.def.1(preds[,i],preds[,j])))
}
d1.def(rg, Cleve[1:10,])

################################################################################

# pt counts the pairs give predictions of terminal nodes
pt <- function(pred.tn){
  # predict the terminal nodes the obs in df fall to
  table(pred.tn) %>% unlist %>% unname %>% (function(x) x*(x-1)/2) %>% sum
}

pt(pred.tn$predictions[,2])

pt.idx <- function(pred.tn){
  nObs <-  length(pred.tn)
  I1 <- outer(1:nObs , 1:nObs , (function(i,j) pred.tn[i]==pred.tn[j])) 
  r1 <- which(I1==TRUE, arr.ind = TRUE)
  sameTn.idx.1 <- r1[(r1[,'row'] < r1[,'col']),]
  return(sameTn.idx.1)
}

pt.idx(pred.tn$predictions[,1])

f.alt <- function(pred1, pred2){
  nObs <-  length(pred1)
  assertthat::assert_that(length(pred1)==length(pred2), msg='No can do. Need predictions of same length.')
  I1 <- outer(1:nObs , 1:nObs , (function(i,j) pred1[i]==pred1[j])) 
  r1 <- which(I1==TRUE, arr.ind = TRUE)
  sameTn.idx.1 <- r1[(r1[,'row'] < r1[,'col']),]
  I2 <- outer(1:nObs , 1:nObs , (function(i,j) pred2[i]==pred2[j]))
  I2[sameTn.idx.1] %>% sum
}

f.alt(pred.tn$predictions[,2],pred.tn$predictions[,1])

# check this :
i <- 4
f.alt(pred.tn$predictions[,i],pred.tn$predictions[,i])
pt(pred.tn$predictions[,i])
# always equal


d1.alt <-  function(ptnp, i,j){
  nObs <-  length(ptnp[,i])
  (pt(ptnp[,i]) 
   + pt(ptnp[,j]) 
   - 2*f.alt(ptnp[,i],ptnp[,j]))*2/(nObs*(nObs-1))
  }
  
# compare to results of createDMd1

d1.alt(1,2)
d1.alt(1,3)
d1.alt(1,4)
d1.alt(1,5)



################################################################################

nObs <- 100
rg$call$num.trees

bench.2 <- function(rg , data){
  nObs <- nrow(data)
  ptnp <- predict(rg, data=data, type='terminalNodes')$predictions # P_red(... type = T_erminal N_odes)$ P_redictions
  A <- matrix(0,nrow=rg$num.trees, ncol=rg$num.trees)
  for(i in 2:rg$num.trees)
    for(j in 1:(i-1)){
      A[i,j] <- d1.alt(ptnp,i,j) 
      A[j,i] <- d1.alt(ptnp,i,j) # make it symmetric
    }
  return(A)
  }

bench.3 <- function(rg , data){
  nObs <- nrow(data)
  ptnp <- predict(rg, data=data, type='terminalNodes')$predictions # P_red(... type = T_erminal N_odes)$ P_redictions
  # outer(1:rg$num.trees , 1:rg$num.trees , function(i,j) d1.alt(ptnp,i,j))
  outer(1:rg$num.trees , 1:rg$num.trees , Vectorize(function(i,j) d1.alt(ptnp,i,j)))
}



bench.2(rg, Cleve[1:nObs,])
bench.3(rg, Cleve[1:nObs,])


# compare runtimes

library(microbenchmark)

d1.def(rg, Cleve[1:nObs,])
bench.2(rg, Cleve[1:nObs,])
createDMd1(rg$forest, Cleve[1:nObs,])

microbenchmark(createDMd1(rg$forest, Cleve[1:nObs,]) 
              , d1.def(rg, data=Cleve[1:nObs,]) 
               , bench.2(rg, data=Cleve[1:nObs,])
              # , bench.3(rg, data=Cleve[1:nObs,])
) -> mb.tab

mb.tab
autoplot(mb.tab)
mb.tab %>% summary %>% xtable



