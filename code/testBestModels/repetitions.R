# script to test models with parameters found in hpo
# 18.5.2022

rm(list=ls())

library(ranger)
library(xtable)
library(caret)
library(dplyr)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection
source('code/source/distance-matrices.R')
source('code/source/prep.R')

#### training set
{data.set.name <- 'Cleve'
data.train <- get(data.set.name)[,1:11]
attr(data.train,'data.set.name') <- data.set.name
data.set.name <- NULL}
# attribute data.set.name will be used in documentation (end of script)


#### test set
{data.set.name <- 'Hung'
data.test <- get(data.set.name)[,1:11]
attr(data.test,'data.set.name') <- data.set.name
data.set.name <- NULL}
# attribute data.set.name will be used in documentation 
# end of script, returned list of function f1

################################################################################

# base cases : bc
nLoops <- 100
doc.bc <- matrix(rep(NA,10*nLoops), ncol=10) # bc for base case
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
  
  # documenting the base case : logloss for regular small forests and the full forest
  szs <-  c(5,20,40,50,60,70,80,90,100) # szs : sIzEs
  doc.bc[i,] <- c( Vectorize(function(n) calcLogloss(subforest(forest, 1:n), data.test))(szs)
                  , calcLogloss(forest, data.test)
                  )
  
}
doc.bc <- data.frame(doc.bc)
names(doc.bc) <-  c(szs , 500)

rbind(
apply(doc.bc,2,function(x) c(mean(x),sd(x)))
, apply(doc.bc,2 ,mean) %>% (function(x) x/x[10]) 
, apply(doc.bc,2 ,sd) %>% (function(x) x/x[10])
)  %>% xtable -> xdoc
digits(xdoc) <- 4
rownames(xdoc) <- c('mean logloss','sd logloss','mean ratio', 'sd ratio')
xdoc

# plot(sqrt(as.numeric(names(doc.bc)))[-10] , apply(doc.bc,2 ,sd) %>% unlist %>% .[-10])
par(mar=c(4,4,2,1)+0.1)
plot(names(doc.bc)[-10] 
     , type='l'
     , apply(doc.bc,2 ,sd) %>% unlist %>% .[-10]
     , main='sd of regular (small) forests'
     , ylab='sd logloss'
     , xlab='size regular small forest')

################################################################################
#### testing Chipman and Meiner forests with parameters from hpo 
################################################################################

test.models <-  list()
test.models[[1]] <- list(method='meiner'
                         , metric='sb'
                         , parameter=list('cutoff'=0.1, 'sizeSF'=500))
test.models[[2]] <- list(method='meiner'
                         , metric='d0'
                         , parameter=list('cutoff'=0, 'sizeSF'=500))
test.models[[1]] <- list(method='meiner'
                         , metric='d1' # takes 10 hours to run!
                         , parameter=list('cutoff'=0.1, 'sizeSF'=500))
                         
test.models[[1]] <- list(method='meiner'
                         , metric='d0' 
                         , parameter=list('cutoff'=0.5, 'sizeSF'=500))
test.models[[2]] <- list(method='meiner'
                         , metric='d2' 
                         , parameter=list('cutoff'=0.7, 'sizeSF'=500))
test.models[[1]] <- list(method='meiner'
                         , metric='sb' 
                         , parameter=list('cutoff'=0.6, 'sizeSF'=500))
test.models[[2]] <- list(method='meiner'
                         , metric='d1' # takes 10 hours to run!
                         , parameter=list('cutoff'=0.8, 'sizeSF'=500))

f1 <- function(method, metric, parameter, nLoops=100){
  #' build forest for given method and evaluate its performance
  #' 
  #' @param method : chip1, chip2, meiner
  #' @param metric : d0, d1, d2, sb
  #' @param parameter list with cutoff, sizeSF (500 indicating original Chipman methods, smaller: stopped algorithms to return forests of (at most the) indicated size). chip1 needs selection method best or central
  #' @param nLoops number of repetitions / loops to run. default is 100.
  #' uses data.train from parent environment to built the random forest
  #' uses attribute data.set.name of data sets data.train and data.test from parent environment to document results
  
  if(method=='meiner'){
    doc <- matrix(rep(NA,4*nLoops), ncol=4)
  }else{
    doc <- matrix(rep(NA,3*nLoops), ncol=3)
  }
  
  for(i in 1:nLoops){
    if(i %% 10 == 0){ print(paste(i/nLoops, Sys.time())) }
    rg <- ranger(CAD~.
               , data = data.train # from parent environment
               , num.trees = 500
               , replace = F # T required for predictions ? Of cource we want predictions!
               , mtry= 3 # default : 3
               , importance = 'impurity'
               , probability = T 
               , min.node.size = 13 
               , keep.inbag = T # needed to generate OOB obs as validation data
    )
  
    forest <- rg$forest
  
    dm <- createDM(forest=forest , type=metric , dft=data.train)

    # use OOB observations for each tree to calculate the tree's logloss
    # prediction probabilities : pp
    pp <- predict(forest 
    , data=data.train 
    , predict.all = T)$predictions[,2,] # dim 303 x 500
    
    # List of Loglosses : LL
    lapply(1:forest$num.trees
           , function(t){ 
           OOB <- which(rg$inbag.counts[[t]]==0)
           pp[OOB,t] %>% 
             calcLogloss2( df=data.train[OOB,] ) %>% 
             unlist
           }
         ) %>% 
      unlist -> LL
    
    if(method=='chip1'){
    calc_chipForest_1_enforce(forest=forest 
                              , dm = dm
                              , oLL= order(LL)
                                #, oLL = calc_oLL(forest=forest, data=data.val)
                              , parameter = parameter) %>% 
      unlist -> doc[i,]
    }
    
    if(method=='meiner'){
      mf <- grow_meinForest(dm=dm
                            , LL=LL
                            , parameter=parameter) # meiner forest
      
      sz <- ifelse(parameter$sizeSF==500,length(mf$forest),parameter$sizeSF) # size
      
      doc[i,] <- c( calcLogloss(forest, data.test) # full forest
                    , sz # size of (selected sub-) forest
                    , calcLogloss(subforest(forest, mf$forest), data.test) # logloss for (selecte sub-) forest
                    # , calcLogloss(subforest(forest, 1:sz), data.test) # logloss for regular small forest of same size as (selected sub-) forest
                    , calcLogloss(subforest(forest,sample(1:forest$num.trees, sz)), data.test) # make the regular small forest random
      )
    }
    
    if(method=='chip2'){
      print('not implemented')
      return()
    }
  }
 
  doc <- data.frame(doc)
  if(method=='meiner'){
    names(doc) <- c('logloss.full.forest','size','logloss.meiner','logloss.base.case')
  }else{
    doc <- doc[,c(1,3)]
    names(doc) <- c('logloss.chip','size')
  }
  
  return(list('res'=doc, 'call'=list(method=method
                                     , metric=metric
                                     , parameter=parameter
                                     , data.train=attr(data.train,'data.set.name')
                                     , data.test=attr(data.test,'data.set.name'))))
    }

# calling f1 on the models in test.models 
# keeping results in list results
results <-  list()
for(i in 1:length(test.models)){
  results[[i]] <- f1(method=test.models[[i]]$method
     , metric=test.models[[i]]$metric 
     , parameter=test.models[[i]]$parameter
     , nLoops=1000)
}

test.case <- 1
res.cal <- results[[test.case]]$call ; res.cal
res.res <- results[[test.case]]$res

apply(res.res,2,mean)
apply(res.res,2,function(x) c(mean(x),sd(x)))

{with.main <-  F
if(with.main){
  par(mar=c(4,4,2,1)+0.1)
}else{
  par(mar=c(4,4,0,1)+0.1)
}
(res.res[,'logloss.base.case'] - res.res[,'logloss.meiner']) %>% 
  hist(xlab='logloss overshoot of base case'
       , main=ifelse(with.main,'overshoot for regular small forest vs meiner forest','')
       , breaks=100)
legend('topright',legend=c(res.cal$method, res.cal$metric , res.cal$parameter$cutoff))
box()
}

dip.test(res.res[,'logloss.base.case'] - res.res[,'logloss.meiner'])
t.test(res.res[,'logloss.base.case'] - res.res[,'logloss.meiner'], alternative='greater')

doc.tests <-  matrix(NA, nrow=length(results), ncol=4)
for(r in 1:length(results)){
  oversh <- results[[r]]$res[,'logloss.base.case'] - results[[r]]$res[,'logloss.meiner']
  tt <- t.test(oversh, alternative='greater')
  doc.tests[r,] <-  c(dip.test(oversh)$p.value 
                      , tt$p.value
                      , tt$conf.int[[1]]
                      , attr(tt$conf.int,'conf.level'))
}
doc.tests <-  as.data.frame(doc.tests)
colnames(doc.tests) <-  c('dip test p-value','t-test p-value','lower bound of mean overshoot','t-test confidence level')
doc.tests %>% xtable -> dt.xt
digits(dt.xt) <- c(0,4,4,4,2)
dt.xt
# I need a normal distribution (of the difference) to apply the t test

# I do NOT need a normal distribution for the guys I do the difference with, 
# still, it might be interesting:
hist(res.res[,'logloss.base.case'], breaks=10 , main='looks exponential')
hist( res.res[,'logloss.meiner'], breaks=10, main='looks exponential')

hist(res.res[,'logloss.base.case'], breaks=40, main='looks strange and multimodal')
hist( res.res[,'logloss.meiner'], breaks=40,main='looks strange and multimodal')

#save(results, file='data/testNewMethods/results-meiner-*.rda')
#save(test.models , results , file='data/testNewMethods/***.rda ')

#results[[length(results)+1]] <-  res1

################################################################################
# subforest of high performers


"{data.set.name <- 'Hung'
data.set <- get(data.set.name)[,1:11]
fold1 <- createDataPartition(data.set$CAD, 1, 0.5) %>% unlist
data.val <-  data.set[fold1,]
attr(data.val,'data.set.name') <- data.set.name
data.test <- data.set[-fold1,]
attr(data.test,'data.set.name') <- data.set.name
data.set.name <- NULL
data.set <-  NULL
}"

# attribute data.set.name will be used in documentation (end of script)

f2 <- function(data.train, data.val , data.test, sz=50 , nLoops=100, returnLL = F){
  #' grow a forest on data.train, select sz trees that perform best on data.val into a subforest and test this subforest on data.test
  #' 
  #' @param data.val : strings 'OOB' or 'split' indicating how to obtain validation data from training or test data or a data frame of validation data, needed to rank trees (lowest logloss on validation data is highest rank)
  #' @param data.test : a data frame for testing
  #' @param sz : size of the sub-forest to build
  #' @param nLoops number of repetitions / loops to run. default is 100.
  #' uses attribute data.set.name of data sets data.train and data.test from parent environment to document results
  
  LL_all_trees <- function(forest,data){  
    pp <- predict(forest 
                  , data=data
                  , predict.all = T)$predictions[,2,]
    lapply(1:forest$num.trees
           , function(t){ 
             pp[,t] %>% 
               calcLogloss2( df=data ) %>% 
               unlist}
    ) %>% 
      unlist 
  }
  
  if(is.character(data.val)){
    data.val.set <- NA
    }else{
    data.val.set <- data.val # we might test for data.val to be a data frame or a matrix?!
    keep.inbag <- FALSE
    } 
  
  data.test.set <- data.test 
  
  if(returnLL){
    # 3 columns, for correlations of logloss on training, validation, and test data
    doc.LL <- matrix(0,nrow=nLoops,3) #%>% data.frame
  }
  
  LL <-  NULL
  
  {print('data sets')
  print(paste('data.train', data.train %>% dim %>% as.vector))
  print(paste('data.val.set', data.val.set %>% dim))
  print(paste('data.test.set', data.test.set %>% dim))
  }
  
  doc <- matrix(rep(NA,4*nLoops), ncol=4)
  
  for(i in 1:nLoops){
    if(i %% 10 == 0){ print(paste(i/nLoops, Sys.time())) }
    rg <- ranger(CAD~.
                 , data = data.train 
                 , num.trees = 500
                 , replace = T # replace = T is good for doing predictions later on? would not have guessed from the name!
                 , mtry= 3 # default : 3
                 , importance = 'impurity'
                 , probability = T 
                 , min.node.size = 13 
                 , keep.inbag = T # need inbag data only for data.val =='OOB'
    )
    
    forest.full <- rg$forest
     
    # rank trees : get List of Loglosses : LL
  
    if(is.character(data.val)){
        if(data.val=='OOB'){
        # on individual OOB observations , (may be) different for each tree
        # prediction probabilities : pp
        pp <- predict(forest.full 
                      , data=data.train 
                      , predict.all = T)$predictions[,2,] # dim 303 x 500
        
        lapply(1:forest.full$num.trees
               , function(t){ 
                 OOB <- which(rg$inbag.counts[[t]]==0)
                 pp[OOB,t] %>% 
                   calcLogloss2( df=data.train[OOB,] ) %>% 
                   unlist}
        ) %>% 
          unlist -> LL # LL created , no data.val.set needed
      }
        if(data.val=='split'){
          fold1 <- createDataPartition(data.test$CAD, 1, 0.5) %>% unlist # split the ORIGINAL data.test which is a df
          data.val.set <-  data.test[fold1,]
          data.test.set <- data.test[-fold1,] # overwriting setting data.test.set at beginning of function
        }
      }
    
    if(is.null(LL)){ # LL has already been created for data.val 'OOB' , the others have created data.val.set
      # prediction probabilities : pp
      
      LL <- LL_all_trees(forest.full, data.val.set)
    }
    
    hp <- order(LL)[1:sz]
    
    if(returnLL){
      cor(cbind(LL_all_trees(forest.full, data.train)
                , LL 
                , LL_all_trees(forest.full, data.test.set))) %>%
        as.dist %>%
        as.vector -> doc.LL[i,]
    }
    doc[i,] <- c( calcLogloss(forest.full, data.test.set) # full forest
                    , sz # size of (selected sub-) forest
                    , calcLogloss(subforest(forest.full, hp), data.test.set) # logloss for (selecte sub-) forest
                    , calcLogloss(subforest(forest.full, 1:sz), data.test.set)) # logloss for regular small forest of same size as (selected sub-) forest
    LL <-  NULL # needed for if(is.null(LL)){ ...} above
    } # end of for loop
  
  doc <- data.frame(doc)
  names(doc) <- c('logloss.full.forest','size','logloss.high.performer','logloss.base.case')
  
  returnList <- list('res'=doc, 'call'=list(method='high performers'
                                            , data.train=attr(data.train,'data.set.name')
                                            , data.val=ifelse(is.character(data.val), data.val, attr(data.val,'data.set.name'))
                                            , data.test=attr(data.test,'data.set.name')
                                            , size=sz
                                            , repetitions=nLoops
  ))
  
  if(returnLL){
    doc.LL <- data.frame(doc.LL)
    names(doc.LL) <-  c('cor.train.val','cor.train.test','cor.val.test')
    returnList$corLL <- doc.LL
  }
  
  return(returnList)

}

set.seed(1)
res1 <- f2(data.train, data.val='OOB' , data.test, sz=50 , nLoops=1000, returnLL=T) # no validation data, just a strategy for creating the validation data from training or test data
#res1 <- f2(data.train=data.train , data.val=data.train , data.test=data.test, sz=50 , nLoops=100, returnLL=T) 
res1$call
res1$res %>% apply(2,function(x) c(mean(x),sd(x)))
res1$corLL %>% apply(2,function(x) c(mean(x),sd(x)))

