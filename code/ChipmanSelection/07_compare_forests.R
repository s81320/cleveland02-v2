# code / ChipmanSelection / 07 compare forests.R
# 31.5.2022

# grow different forests on a simulation set and compare them
# it is also kind of a test for how the selection of trees work

# it is troubling that the Chipman 2 and Meiner forest for the d0 dissimilarity 
# at representation 0 give exactly the same mean logloss

rm(list=ls())

file_name_script <-  'code/**.R'
list.files(path = 'code/' , recursive = T)
# assertthat::assert_that(file_name_script %in% list.files(recursive = T))


library(ranger)
library(xtable)
library(caret)
library(dplyr)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection , calc_LL (for all trees)
source('code/source/distance-matrices.R')
source('code/source/prep.R')

load('data/nursery02/nursery02_01.rda')

i.sim <-  1
metric <-  'd0'
dm <- doc[[i.sim]]$DM[[metric]]
forest <- doc[[i.sim]]$rg$forest
LL <- calc_LL(forest=forest, data=doc[[i.sim]]$`bootstrapped training data`)

mf.2 <- grow_meinForest(dm = dm , LL = LL , parameter = list(cutoff=0.2, sizeSF=500) )
cf1.best <- grow_chipForest_1(dm=dm[1:200,1:200], forest=subforest(forest,1:200) ,oLL=order(LL[1:100])
                         , parameter=list(cutoff=0.2, selection='best') 
                         , output=T)
cf1.central <- grow_chipForest_1(dm=dm, forest=forest ,oLL=order(LL) 
                              ,parameter=list(cutoff=0.2, selection='central')
                              , output=T)
# parameter selection = central should not require oLL at all , but should be invariant under permutations of oLL
# of course - there is variation in the clustering that might depend on the ordering of elements ...

cf2 <- grow_chipForest_2(dm=dm, forest=forest, oLL = order(LL)
                         , parameter = list(cutoff=0.2, sizeSF=500))

cf1.best$forest
cf1.central$forest
cf2$forest
mf$forest
mf.2$forest

setequal(cf1.best$forest, cf1.central$forest)
setequal(cf1.best$forest, cf2$forest)
setequal(cf1.best$forest, mf$forest)
setequal(cf2$forest, mf$forest)

intersect(cf1.best$forest, mf$forest)
intersect(cf2$forest, mf$forest)

calcLogloss(forest, Swiss)
calcLogloss(subforest(forest,cf1.best$forest), Swiss)
calcLogloss(subforest(forest,cf1.central$forest), Swiss)
calcLogloss(subforest(forest,cf2$forest), Swiss)
calcLogloss(subforest(forest,mf$forest), Swiss)

data.test <-  Swiss
eval_meinForest(dm, forest, LL, list(cutoff=0.2, sizeSF=500))
mf$progress




