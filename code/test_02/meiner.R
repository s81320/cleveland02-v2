# code/test_02/meiner.R
# 31.5.2022

# script to test models with repetitions
# test like a data scientist : mix data and split it (homogenization),
# split it many different ways and take the mean over those splits: 
# Monte Carlo cross validation

rm(list=ls())

file_name_script <-  'code/test_02/meiner.R'
list.files(path = 'code/test_02' , recursive = T)
assertthat::assert_that(file_name_script %in% list.files(recursive = T))


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

# create the data differently: mix Cleve and Hung , leave ot Swiss 
# and partition it into test, dev, train
# partition as 50 / 25 / 25 so the sizes of training sets are similar to simulation and repetitions done before.

# find differences in how to gain the devset: 
# data.dev 'OOB' results in worse results than the default forest
# data.dev 'split' is and high performers outperform the default forest.
# -> OOB tests the splitting into test and train, if these sets are close, a tree build on train performs good on train and will perform good on test (as test is similar to train)
# if train and test are different, the the tree built on train performs good on train but not on test. 
# selecting trees that were grown on training sets with a similar test set is no quality measure for the tree.

# load data from beginning of script (Cleve, Swiss, Hung)

# data.all <-  rbind(Cleve[,1:11], Swiss[,1:11], Hung[,1:11])
data.all <-  rbind(Cleve[,1:11], Hung[,1:11])
dim(data.all)

# find optimal representation parameter using a validation set
test_meiner_forest <- function(data, metric, nLoops=5){
  #' Generate a 50/25/25 split, grow a full forest on the training data.
  #' Generate a Meiner forest (tree ranking based on OOB training data)
  #' for each representation parameter 0,0.1,...,0.9 
  #' and pick the Meiner forest that performs best (by logloss) on the validation set.
  #' Test the Meiner forest with the best parameter on the test set (logloss)
  #' 
  
  dp <- createDataPartition(data$CAD ,times = nLoops , p= 0.5) # dp : data partition
  doc <- data.frame(matrix(NA,nrow=nLoops, ncol=6))
  for(i in 1:nLoops){
    if(i%%100==0){
      print(paste(i/nLoops , Sys.time()))
    }
    
    data.train <- data[dp[[i]],]
    data.test <- data[-dp[[i]],] # will soon be overwritten
    
    fold1 <- createDataPartition(data.test$CAD, 1, 0.5) %>% unlist
    data.val <-  data.test[fold1,]
    data.test <- data.test[-fold1,]
    
    rg <- ranger(CAD~.
                 , data = data.train 
                 , num.trees = 500
                 , mtry= 3 # default : 3
                 #, importance = 'impurity' # we do not need variable importance
                 , probability = T 
                 , min.node.size = 13 
                 , keep.inbag = T # need inbag data for ranking trees
    )
    forest.full <- rg$forest
    
    #dm <- createDMd2(forest.full, dft = data.val)
    #dm <- createDMd0(forest.full)
    dm <-  createDM(forest=forest.full, type=metric, dft = data.train)
    # summary(as.dist(dm)) %>% print

    # ranking trees 
    # on individual OOB observations , (may be) different for each tree
    # prediction probabilities : pp
    pp <- predict(forest.full 
                    , data=data.train 
                    , predict.all = T)$predictions[,2,] # indexing : obs x trees
      
    lapply(1:forest.full$num.trees
             , function(t){ 
               OOB <- which(rg$inbag.counts[[t]]==0)
               pp[OOB,t] %>% 
                 calcLogloss2( df=data.train[OOB,] ) %>% 
                 unlist}
      ) %>% 
        unlist -> LL # LL created , no validation data needed
    
    # validation
    # find optimal parameter
    
    # criterion for optimality: smallest overall logloss?? 
    doc.forest <-  list()
    ct <- 1
    doc.val <-  matrix(0,nrow=10, ncol=3) # 
    for(r in seq(0,0.9,0.1)){
      mf <- grow_meinForest(dm, LL, list(cutoff=r , sizeSF=500), output = F)
      #View(mf$progress)
      doc.forest[[ct]] <- mf$forest
      doc.val[ct,] <- c(r , length(mf$forest), calcLogloss(forest = subforest(forest.full, mf$forest), data.val))
      ct <- ct+1
    }
    
    ct.opt <-  which.min(doc.val[,3] %>% unlist)
    "print(ct.opt)
    print(paste('optimal for' , doc.val[ct.opt,3] ))
    print(paste('at size' , doc.val[ct.opt,2]))
    print(doc.val)
    print('***************')"
    
    forest.selected <- doc.forest[[ct.opt]]
    "print('selected forest')
    print(forest.selected)
    print('***************')"
    
    sz <-  length(forest.selected)
    doc[i,] <- c( calcLogloss(forest.full, data.test) # full forest
                  , doc.val[ct.opt,1] # optimal parameter for level of representation
                  , doc.val[ct.opt,3] # logloss at optimal parameter
                  , sz # size of (selected sub-) forest
                  , calcLogloss(subforest(forest.full, forest.selected), data.test) # logloss for (selecte sub-) forest
                  , calcLogloss(subforest(forest.full, 1:sz), data.test)) # logloss for regular small forest of same size as (selected sub-) forest
  }
  
  names(doc) <- c('logloss.default','rep.parameter','validation logloss','size','logloss.meiner','logloss.regular.small')
  return(doc)
}

# fix a representation parameter r
test_meiner_forest_2 <- function(data, metric, r , nLoops=5){
  
  dp <- createDataPartition(data$CAD ,times = nLoops , p= 0.5) # dp : data partition
  doc <- data.frame(matrix(NA,nrow=nLoops, ncol=4))
  for(i in 1:nLoops){
    if(i%%10==0){
      print(paste(i/nLoops , Sys.time()))
    }
    
    data.train <- data[dp[[i]],]
    data.test <- data[-dp[[i]],] 

    rg <- ranger(CAD~.
                 , data = data.train 
                 , num.trees = 500
                 , mtry= 3 # default : 3
                 #, importance = 'impurity' # we do not need variable importance
                 , probability = T 
                 , min.node.size = 13 
                 , keep.inbag = T # need inbag data for ranking trees
    )
    forest.full <- rg$forest
    
    #dm <- createDMd2(forest.full, dft = data.val)
    #dm <- createDMd0(forest.full)
    dm <-  createDM(forest=forest.full, type=metric, dft = data.train)
    # summary(as.dist(dm)) %>% print
    
    # ranking trees 
    # on individual OOB observations , (may be) different for each tree
    # prediction probabilities : pp
    pp <- predict(forest.full 
                  , data=data.train 
                  , predict.all = T)$predictions[,2,] # indexing : obs x trees
    
    lapply(1:forest.full$num.trees
           , function(t){ 
             OOB <- which(rg$inbag.counts[[t]]==0)
             pp[OOB,t] %>% 
               calcLogloss2( df=data.train[OOB,] ) %>% 
               unlist}
    ) %>% 
      unlist -> LL # LL created , no validation data needed
    
    mf <-  grow_meinForest(dm, LL, parameter=list(cutoff=r, sizeSF=500), output = F)
    forest.selected <- mf$forest

    sz <-  length(forest.selected)
    doc[i,] <- c( calcLogloss(forest.full, data.test) # full forest
                  , sz # size of (selected sub-) forest
                  , calcLogloss(subforest(forest.full, forest.selected), data.test) # logloss for (selecte sub-) forest
                  , calcLogloss(subforest(forest.full, 1:sz), data.test)) # logloss for regular small forest of same size as (selected sub-) forest
  }
  
  names(doc) <- c('logloss.default','size','logloss.meiner','logloss.regular.small')
  return(doc)
}

####################
#### first test ####
####################

doc.t <- test_meiner_forest(data=data.all , metric='d1', nLoops = 100)
res.doc <- list(parameter=list(metric='d1'
                               , type='Meiner'
                               , criterion='lowest logloss on validation set')
                , results=doc.t
                , info=paste('training data is 50% of mixed Cleve, Hung. Code in'
                             , file_name_script)
)
#save(res.doc , file='data/test_02/results-meiner-d1-10**.rda') # encode the metric!

doc.t <-  res.doc$results[[1]]

# give options on models to test directly
doc.t %>% apply(2,function(x) c(mean(x),sd(x)))
# doc.hp 
doc.t[doc.t$rep.parameter>0,] %>% apply(2,function(x) c(mean(x),sd(x)))

doc.t$rep.parameter %>% table

doc.t %>% group_by(rep.parameter) %>% summarise(m=mean(size))


oversh <- doc.t[,6]-doc.t[,5] # test with looking for optimal parameter
oversh <- doc.t[,4]-doc.t[,3] # result of test_2 with fixed parameter
par(mar=c(4,4,2,1)+0.1)
hist(oversh, breaks=100) # looks multimodal?
# hist(oversh, breaks=10) # looks unimodal, bit skewed?

# same as next line of code
# with(doc.t[doc.t$size>40,] , doc.t$logloss.regular.small - doc.t$logloss.meiner)  %>% hist

doc.t[doc.t$rep.parameter>0,] %>%
  with(hist(logloss.regular.small - logloss.meiner, breaks=20))

# loop over all / some models to test
# use expand grids to generate all combinations of parameters to be tested
parameter <- as.data.frame(c('d0','d2','sb'))
colnames(parameter) <- c('metric') 
parameter[1,] # is now a test case

results <-  list()
for(i in 1:nrow(parameter)){
  print(Sys.time())
  print(parameter[i,])
  results[[i]] <- test_meiner_forest(data=data.all 
                                 , metric =parameter[i,1] 
                                 , nLoops=1000)
  results[[i]] %>% apply(2,function(x) c(mean(x),sd(x))) %>% t %>% xtable -> xt
  digits(xt) <-  4
  print(xt)
}

res.doc <- list(parameter=parameter
                , results=results
                , info=paste('training data is 50% of mixed Cleve and Hung. Code in'
                             , file_name_script)
)
#save(res.doc , file='data/test_02/results-meiner-d0-d2-sb-50-25-25split-1000loops**.rda') # encode the metric!

# strangely high correlation for d2 metric and meiner and regular small forest. How can that be??
corrplot::corrplot(cor(results[[1]][c(1,3,5,6)]))
corrplot::corrplot(cor(results[[2]][c(1,3,5,6)]))
corrplot::corrplot(cor(results[[3]][c(1,3,5,6)]))

# how representation parameter explains size . Very reasonable!
j <- 3
plot(results[[j]]$size~results[[j]]$rep.parameter, xlab='representation parameter', ylab='size')
legend('topright', legend=c('d0','d2','sb')[j])

# which representation parameters gave good results on the validation set?
# when parameter > 0?!
results[[j]]$rep.parameter[results[[j]]$rep.parameter>0]
results[[2]]$rep.parameter %>% table


#################

test_and_document <-  function(metric , rep.parameter, nLoops , filename='***.rda'){
  
  doc.t <- test_meiner_forest_2(data.all, metric , rep.parameter , nLoops)
  res.doc <- list(parameter=list(metric=metric
                               , type='Meiner'
                               , rep.parameter=rep.parameter)
                , results=doc.t
                , info=paste('training data is 50% of mixed Cleveland and Hungary. Code in'
                             , file_name_script)
  )
  save(res.doc , file=paste('data/test_02/results-meiner',metric,rep.parameter,nLoops,'loops*.rda', sep='-')) # encode the metric!
  
  doc.t %>% apply(2,function(x) c(mean(x),sd(x)))
  
  return(res.doc)
}

test_and_document(metric='d2'
                  , rep.parameter=0.2
                  , nLoops=10000
                  , filename=paste('data/test_02/results-meiner',metric,rep.parameter,nLoops,'loops*x*.rda', sep='-'))
