# 30.3.2022
# something for testing

# Testing should be on the full Cleveland data, not again with simulations.

# we have to work with the full Cleveland data,
# train forests on them, calculate dissimilarity matrices
# set parameters
# select the forests by the new method(s)

# test on Hungary

results <-  list()
#rm(list=ls())

library(ranger)
library(xtable)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection
source('code/source/distance-matrices.R')
source('code/source/prep.R')

data.train <-  Cleve[,1:11]
data.val <- Swiss[,1:11]
data.test <- Hung[,1:11]

DM <- list()
# base cases

nLoops <- 100
doc.bc <- matrix(rep(NA,10*nLoops), ncol=10)
#fL <- list()
for(i in 1:nLoops){
  if(i %% 10 == 0){ print(paste(i/nLoops, Sys.time())) }
  rg <- ranger(CAD~.
               , data = data.train
               , num.trees = 500
               , replace = F 
               , mtry= 3 # default : 3
               , importance = 'impurity'
               , probability = T 
               , min.node.size = 13 
  )
  
  forest <-  rg$forest
  
  szs <-  c(100, 150, 200) # szs : sIzEs
  doc.bc[i,] <- c( Vectorize(function(n) calcLogloss(subforest(forest, 1:n), data.test))(szs)
                   , calcLogloss(forest, data.test)
  )
  
}
doc.bc <- data.frame(doc.bc)
names(doc.bc) <-  c(szs , 500)

rbind(
  apply(doc.bc,2,function(x) c(mean(x),sd(x)))
  , apply(doc.bc,2 ,sd) %>% (function(x) x/x[10])
)  %>% xtable -> xdoc
digits(xdoc) <- 4
xdoc

plot(sqrt(as.numeric(names(doc.bc)))[-10] , apply(doc.bc,2 ,sd) %>% 
       unlist %>%
       .[-10])


# build model on full Cleveland data (no bootstrap or simulation)

# stopped Meiner
metric <- 'sb'
sizeSF <-  5
cutoff <-  0.1
parameter <- list('cutoff'=cutoff, sizeSF=sizeSF)

# unstopped Meiner
metric <- 'd1'
sizeSF <-  500
cutoff <-  0.25
parameter <- list('cutoff'=cutoff, sizeSF=sizeSF)

# simplified Chipman 1
metric <- 'd1'
parameter <- list('sizeSF'=50, 'selection'='central')

nLoops <- 100
doc <- matrix(rep(NA,3*nLoops), ncol=3)
#fL <- list()
for(i in 1:nLoops){
  if(i %% 10 == 0){ print(paste(i/nLoops, Sys.time())) }
  rg <- ranger(CAD~.
               , data = data.train
               , num.trees = 500
               , replace = F 
               , mtry= 3 # default : 3
               , importance = 'impurity'
               , probability = T 
               , min.node.size = 13 
  )
  
  forest <-  rg$forest
  
  dm <- createDM(forest=forest , type=metric , dft=data.train)
  "
  mf <- grow_meinForest(dm=dm 
  , LL=calc_LL(forest, data.val) 
  , parameter=parameter)
  
  if(parameter$sizeSF==500){
    # unstopped
    doc[i,] <- c(calcLogloss(forest, data.test)
                 , length(mf$forest)
                 , calcLogloss(subforest(forest, mf$forest), data.test))
  }else{
    # stopped forest
    doc[i,] <- c( calcLogloss(forest, data.test)
                ,calcLogloss(subforest(forest, mf$forest), data.test)
                , calcLogloss(subforest(forest, 1:sizeSF), data.test))
  }
"
  # uses Swiss for validation -> not a good idea ??!!
  calc_chipForest_1_enforce(forest=forest 
                            , dm = dm
                            , oLL = calc_oLL(forest=forest, data=data.val)
                            , parameter = parameter) %>% 
    unlist -> doc[i,]
  
}
doc <- data.frame(doc)

# stopped Meiner forest, size deterministic
names(doc) <- c('default','meiner','regular small')

# unstopped Meiner forest , size is part of the outcome
names(doc) <- c('default','size','meiner')

# simplified Chipman 1
doc <- doc[,c(1,3)]
names(doc) <- c('logloss','size')

# results <-  list()
doc <-  data.frame(doc)
names(doc) <- c(5,40,50,60,70,80,90,100,500) %>% as.character
results[['regular small']] <-  list(doc=doc.bc)
results[['stopped Meiner 2']] <-  list(doc=doc, 'type'='stopped Meiner' , metric=metric ,  parameter = parameter)
results[['simplified Chipman 1 d']] <-  list(doc=doc, 'type'='simplified Chipman 1' , metric=metric ,  parameter = parameter)

apply(doc,2,mean)
apply(doc,2,function(x) c(mean(x),sd(x)))
(doc[,3] - doc[,2]) %>% 
  hist(main='overshoot for regular small forest vs meiner forest'
       , breaks=nLoops/2)
# test model on Hungary

#save(results, file='data/results-01*.rda')

