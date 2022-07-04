# introduction , methods , trees and random forests

# logloss performance for the single tree and 
# logloss for 3 individual trees vs the forest of these 3 trees 

# based on first 3 trees in first simulation

rm(list=ls())

file_name_script <-  'code/introduction/ensemble-magic.R'
list.files(path = 'code/introduction' , recursive = T)
assertthat::assert_that(file_name_script %in% list.files(recursive = T))

library(ranger)
library(xtable)
# library(caret)
library(dplyr)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/plotTree.R') # plotTree1 , plotTreeFb

{ # block copied from code/nursery/01.R
  # may not be changed! Creates exactly the same (first 3) trees as in first simulation
  # but with added keep.inbag=T
  Cleve_enriched <- create_prob(Cleve)
  nBs <- 50
  
  seed <- 13
  set.seed(seed)
  cr<-createResample(Cleve_enriched$CAD
                     , times = nBs
                     , list = T)
  
  ct <- 1
  i <- 1
  data.train <- new_bootstrap(Cleve_enriched , cr[[i]])[,-12] # no probabilities
  
  # grow small forest (size nT)
  rg <- ranger(CAD~.
               , data = data.train 
               , num.trees = 500
               , replace = F 
               , mtry= 3 
               , importance = 'impurity'
               , probability = T 
               , min.node.size = 13 
               , keep.inbag = T # NEW
  )
}

df <- data.train
rbind(
  Inbag = c(calcLogloss(forest=subforest(rg$forest,1), data.train[rg$inbag.counts[[1]],]),
            calcLogloss(forest=subforest(rg$forest,2), data.train[rg$inbag.counts[[2]],]),
            calcLogloss(forest=subforest(rg$forest,3), data.train[rg$inbag.counts[[3]],]),
            NA , 
            NA
  ),
training=c(
calcLogloss(forest=subforest(rg$forest,1), df),
calcLogloss(forest=subforest(rg$forest,2), df),
calcLogloss(forest=subforest(rg$forest,3), df),
calcLogloss(forest=subforest(rg$forest,1:3), df),
calcLogloss(forest=rg$forest, df)
) ,

unseen=c(
calcLogloss(forest=subforest(rg$forest,1), df = Swiss),
calcLogloss(forest=subforest(rg$forest,2), df = Swiss),
calcLogloss(forest=subforest(rg$forest,3), df = Swiss),
calcLogloss(forest=subforest(rg$forest,1:3), df = Swiss),
calcLogloss(forest=rg$forest, df = Swiss)
)
) %>% xtable 

# compare to best tree , CART

#install.packages('rpart')
library(rpart)
library(rpart.plot)

tri <- 1
df <- data.train[which(rg$inbag.counts[[tri]]==1) ,]
#df <- data.train
#tree <- rpart::rpart(CAD ~ ., data = df, control = rpart::rpart.control(cp = 0.0001))
tree <- rpart::rpart(CAD ~ ., data = df)
rpart.plot::rpart.plot(tree)


# Step3: Prune the tree using the best cp.
bestcp <- tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"]
tree.pruned <- rpart::prune(tree, cp = bestcp)

# predictions do not help when they are for a classification tree
# cannot calculate logloss for classification
# predict(tree.pruned, Swiss, type="class")

# we can compare trees directly: See how the trees in the random forest are not optimal, 
# not on training data for the forest (different data!) and not on training data for the individual tree 

rpart.plot::rpart.plot(tree.pruned)

