# 31.3.2022
# how to test a model without simulation, with repetition

# We need a validation set to rank individual trees.
# when doing simulations, this was the set of OOB observations of the Cleveland data,
# observations not used to grow the default forest.

# here we try different validation sets: 
# Swiss, corrected Swiss, Cleve , VA, and a combination of corrected Swiss and VA

# when there is no validation set (is.null(data.val)) then the OOB observations per tree are used.

library(ranger)
library(xtable)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection
source('code/source/distance-matrices.R')
source('code/source/prep.R')

data.train <-  Cleve[,1:11]

#data.val <- rbind(Swiss[,1:11], VA[,1:11])
#data.val <- VA[,1:11]

data.val <-  Swiss[,1:11]
#range(Swiss$STDepression)

# this correction does not seem correct any more, as of April 6th 2022
#corSwiss <- Swiss[,1:11]
#corSwiss$STDepression <-  corSwiss$STDepression - min(corSwiss$STDepression)
#range(corSwiss$STDepression)
#data.val <-  corSwiss
#data.val <- NULL

#data.val <- Hung[,1:11]
#data.val <-  rbind(VA[,1:11], Cleve[,1:11], Swiss[,1:11])

data.test <- Hung[,1:11]

nLoops <- 1000
sizeSF <-  50 # 125 was optimal in simulations , on (corrected??) Swiss ?
doc <- matrix(rep(NA,3*nLoops), ncol=3)
#fL <- list()
for(i in 1:nLoops){
  #for(i in 1:1){
  if(i %% 10 == 0){ print(paste(i/nLoops, Sys.time())) }
  rg <- ranger(CAD~.
               , data = data.train
               , num.trees = 500
               , replace = F 
               , mtry= 3 # default : 3
               #, importance = 'impurity'
               , probability = T 
               , min.node.size = 13 
               , keep.inbag=T
  )
  
  forest <-  rg$forest
  
  # create validation data from OOB Cleveland
  pp <- predict(forest 
                , data=data.train 
                , predict.all = T)$predictions[,2,] # dim 303 x 500
  
  if(is.null(data.val)){
  # OOB <- which(rg$inbag.counts[[t]]==0)
  lapply(1:forest$num.trees
         , function(t){ 
           OOB <- which(rg$inbag.counts[[t]]==0)
           pp[OOB,t] %>% 
             calcLogloss2( df=data.train[OOB,] ) %>% 
             unlist
         }) %>% 
    unlist -> LL
  }else{
    calc_LL(forest=forest, data=data.val) -> LL
  }

  LL %>% order -> oLL 
  
  doc[i,] <- c( calcLogloss(forest, data.test)
                ,calcLogloss(subforest(forest, oLL[1:sizeSF]), data.test)
                , calcLogloss(subforest(forest, 1:sizeSF), data.test))
}

r1 <- apply(doc,2,function(x) c(mean(x),sd(x))) %>% t
r1 <-  cbind(c(500,sizeSF,sizeSF),r1)
colnames(r1) <-  c('size','mean','sd')
rownames(r1) <- c('default','high performers','regular small')
r1 %>% xtable -> xt
digits(xt) <-  4
xt

(doc[,2]-doc[,1]) %>% 
  hist(breaks=20
       , main='histogram, trees ranked by OOB-Cleveland'
       , xlab='logloss overshoot for high performers vs default')

doc.LL <- list()
doc.LL[['Cleve']] <- calc_LL(forest=forest, data= Cleve[,1:11])
doc.LL[['Swiss']] <- calc_LL(forest=forest, data=Swiss[,1:11])
doc.LL[['VA']] <- calc_LL(forest=forest, data= VA[,1:11])
doc.LL[['Hung']] <- calc_LL(forest=forest, data= Hung[,1:11])
doc.LL[['correctedSwiss']] <- calc_LL(forest=forest, data=corSwiss)
doc.LL[['correctedSwiss_VA']] <- calc_LL(forest=forest, data= rbind(VA[,1:11], corSwiss))

m3 <- cbind(doc.LL[['Cleve']] , doc.LL[['Swiss']] 
            , doc.LL[['VA']] , doc.LL[['Hung']] 
           # , doc.LL[['correctedSwiss']] 
          #  , doc.LL[['correctedSwiss_VA']]
          )
colnames(m3) <- c('Cleve','Swiss','VA' , 'Hung'
                 # ,'correctedSwiss','correctedSwiss_VA'
                  )
corrplot::corrplot(cor(m3), method='pie')


