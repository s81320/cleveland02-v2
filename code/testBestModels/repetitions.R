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

# data.val <-  Swiss[,1:11]
# data.val <- NULL
# data.test <- Hung[,1:11]
data.test <- Swiss[,1:11]

# fold1 <- createDataPartition(Hung$CAD, 1, 0.5) %>% unlist
# data.val <-  Hung[fold1,]
# data.test <- Hung[-fold1,]

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
  
  szs <-  c(5,20,40,50,60,70,80,90,100) # szs : sIzEs
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

plot(sqrt(as.numeric(names(doc.bc)))[-10] , apply(doc.bc,2 ,sd) %>% unlist %>% .[-10])

# build model on Cleveland data
method <- 'meiner'

# set other parameters required to test a model
if(method=='meiner'){
  metric <- 'sb'
  parameter <- list('cutoff'=0.25, 'sizeSF'=500)
}else{
  if(method=='chip1'){
    metric <- 'd1'
    parameter <- list('sizeSF'=50, 'selection'='central')
  }else{
      if(method=='chip2'){
        #pass
      }
    }
}

f1 <- function(method, metric, parameter, nLoops=100){

  if(method=='meiner'){
    doc <- matrix(rep(NA,4*nLoops), ncol=4)
  }else{
    doc <- matrix(rep(NA,3*nLoops), ncol=3)
  }
  
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
               , keep.inbag = T # needed to generate OOB obs as validation data
    )
  
    forest <-  rg$forest
  
    dm <- createDM(forest=forest , type=metric , dft=data.train)

    # use OOB observations for each tree to calculate the tree's logloss
    pp <- predict(forest 
    , data=data.train 
    , predict.all = T)$predictions[,2,] # dim 303 x 500
    
    lapply(1:forest$num.trees
           , function(t){ 
           OOB <- which(rg$inbag.counts[[t]]==0)
           pp[OOB,t] %>% 
             calcLogloss2( df=data.train[OOB,] ) %>% 
             unlist
           }
         ) %>% 
      unlist -> LL
  
  " (function(tri){
      (rg$inbag.counts[[tri]]==0) %>%
        which %>% 
        data.train[.,] %>%
        calcLogloss(subforest(forest,tri),.)
    }) %>% 
      Vectorize %>%
      lapply(1:rg$num.trees,.) %>%
      unlist -> LL"
    
    if(method=='meiner'){
      mf <- grow_meinForest(dm=dm
                            , LL=LL
                            , parameter=parameter)
    
      sz <- ifelse(parameter$sizeSF==500,length(mf$forest),parameter$sizeSF) # size
      
      doc[i,] <- c( calcLogloss(forest, data.test)
                    , sz
                    , calcLogloss(subforest(forest, mf$forest), data.test)
                    , calcLogloss(subforest(forest, 1:sz), data.test))
    }
      
    if(method=='chip1'){
    calc_chipForest_1_enforce(forest=forest 
                              , dm = dm
                              , oLL= order(LL)
                                #, oLL = calc_oLL(forest=forest, data=data.val)
                              , parameter = parameter) %>% 
      unlist -> doc[i,]
    }
  }
 
  doc <- data.frame(doc)
  if(method=='meiner'){
    names(doc) <- c('default','size','logloss.meiner','logloss.bc')
  }else{
    doc <- doc[,c(1,3)]
    names(doc) <- c('logloss.chip','size')
  }
  
  return(list('res'=doc, 'call'=list(method=method
                                     ,metric=metric
                                     , parameter=parameter
                                     , data.train='Cleve'
                                     , data.test='Swiss')))
    }

res1 <- f1(method='meiner',metric=metric, parameter=parameter, nLoops=100)

res <- res1$res

apply(res,2,mean)
apply(res,2,function(x) c(mean(x),sd(x)))
(res[,'logloss.bc'] - res[,'logloss.meiner']) %>% 
  hist(main='overshoot for regular small forest vs meiner forest'
       , breaks=nLoops/2)

#save(results, file='data/results-03*.rda')
results[[length(results)+1]] <-  res1
