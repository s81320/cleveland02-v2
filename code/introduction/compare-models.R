# show differences in Chipman 1, Chipman 2 and Meiner forests
# introduction
# 10.7.2022

# should also compare and visualize d0 for best vs central

rm(list=ls())

library(ranger)
library(xtable)
library(caret)
library(dplyr)
library(diptest)
library(cluster)

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


rg <- ranger(CAD~.
             , data = data.train
             , num.trees = 500
             , replace = F # T required for predictions ? Of course we want predictions!
             , mtry= 3 # default : 3
             , importance = 'impurity'
             , probability = T 
             , min.node.size = 13 
             , keep.inbag = T # needed to generate OOB obs as validation data
)

forest <- rg$forest

metric <- 'd0'
dm <- createDM(forest=forest , type=metric , dft=data.train)

# to load from
{#load('data/compare-models.rda') # loads doc1
  k <- 5
  metric <- doc1[[k]]$metric
  dm <- doc1[[k]]$dm
  parameter <- doc1[[k]]$mf$parameter
  cf1 <-  doc1[[k]]$cf1
  cf2 <-  doc1[[k]]$cf2
  mf <-  doc1[[k]]$mf
  LL <-  doc1[[k]]$LL
}

# use OOB observations for each tree to calculate the tree's logloss
# prediction probabilities : pp
predict(forest
        , data=data.train
        , predict.all = T)$predictions[,2,] -> pp # dim 303 x 500

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

parameter <- list('cutoff'=0.2 , 'sizeSF'=500)

{ parameter$selection <- 'best'
  cf1.best <- grow_chipForest_1(dm = dm
                                , oLL= order(LL)
                                #, oLL = calc_oLL(forest=forest, data=data.val)
                                , parameter = parameter
                                , output=T)
}

{ parameter$selection <- 'central'
  cf1.central <- grow_chipForest_1(dm = dm
                    , oLL= order(LL)
                    #, oLL = calc_oLL(forest=forest, data=data.val)
                    , parameter = parameter
                    , output=T)
}

cf1$forest <-  cf1.best


cf2 <- grow_chipForest_2(dm=dm
                         , oLL= order(LL)
                         , parameter=parameter) # chipman 2 forest
cf2$forest

mf <- grow_meinForest(dm=dm
                      , LL=LL
                      , parameter=parameter
                      , output=F) # meiner forest
mf$forest

# compare Chipman 1 and Chipman 2
# d0 at parameter 0 : forests are equal
length(cf1$forest) ; length(cf2$forest)
base::intersect(cf1$forest, cf2$forest) %>% length


# compare Chipman 2 and Meiner
# d0 at parameter 0 : forests are equal
length(cf2$forest) ; length(mf$forest)
base::intersect(cf2$forest, mf$forest) %>% length

LL[cf1$forest] %>% plot(type='l', main='Chipman 1', xlim=c(0,500) , ylim=range(LL))
LL[cf1.central$forest] %>% plot(type='l', main='Chipman 1', xlim=c(0,500) , ylim=range(LL))
LL[cf2$forest] %>% plot(type='l', main='Chipman 2' , xlim=c(0,500), ylim=range(LL))
LL[mf$forest] %>% plot(type='l', main='Meiner' , xlim=c(0,500), ylim=range(LL))

oLLmf <- order(LL[mf$forest])
LL[mf$forest[oLLmf]] %>% plot(type='l', main='Meiner, ordered' , xlim=c(0,500), ylim=range(LL))

# plot line for all trees in forest (500)
# plot dots for trees in representing forests 
{par(mar=c(4,4,0,1)+0.2)
  # par(mar=c(4,4,2,1)+0.2) # when plotting a title
plot(y=LL[order(LL)]
     , x=1:forest$num.trees
     , type='l'
     #, main='Logloss for all trees in forest, ordered'
     , xlim=c(0,500)
     , ylim=range(LL)
     , xlab='trees ordered by performance'
     , ylab='logloss (of individual / single tree)')

oLL <-  order(LL)
lapply(cf1.central$forest, function(i) which(i==oLL)) %>% unlist -> xXcf1
points(x=xXcf1 , y=LL[oLL[xXcf1]], col='yellow')

lapply(cf2$forest, function(i) which(i==oLL)) %>% unlist -> xXcf2
points(x=xXcf2 , y=LL[oLL[xXcf2]]+0.2, col='orange')

lapply(mf$forest, function(i) which(i==oLL)) %>% unlist -> xXmf
points(x=xXmf , y=LL[oLL[xXmf]]-0.2, col='red')

legend('topleft', legend = c('Chipman 1','Chipman 2', 'Meiner') , col=c('yellow','orange','red'), pch=1)
legend('bottomright' , legend=paste(metric, parameter$cutoff,sep=', '))
}

length(cf1$forest) ; calcLogloss(forest=subforest(rg$forest, cf1$forest), df=data.test)
length(cf2$forest) ; calcLogloss(forest=subforest(rg$forest, cf2$forest), df=data.test)
length(mf$forest) ; calcLogloss(forest=subforest(rg$forest, mf$forest), df=data.test)

#doc1 <- list()
#doc1[[5]] <- list(rg=rg , metric=metric , dm=dm , LL=LL , cf1=cf1 , cf2=cf2 , mf=mf)

#################################################################
#### compare Chipman 1 forest with best and central tree per cluster

{par(mar=c(4,4,0,1)+0.2)
  # par(mar=c(4,4,2,1)+0.2) # when plotting a title
  plot(y=LL[order(LL)]
       , x=1:forest$num.trees
       , type='l'
       #, main='Logloss for all trees in forest, ordered'
       , xlim=c(0,500)
       , ylim=range(LL)
       , xlab='trees ordered by performance'
       , ylab='logloss (of individual / single tree)')
  
  oLL <-  order(LL)
  lapply(cf1.best$forest, function(i) which(i==oLL)) %>% unlist -> xXcf1.b
  points(x=xXcf1.b , y=LL[oLL[xXcf1.b]]+0.4, col='yellow', pch='x' )
  
  lapply(cf1.central$forest, function(i) which(i==oLL)) %>% unlist -> xXcf1.c
  points(x=xXcf1.c , y=LL[oLL[xXcf1.c]]-0.4, col='yellow')
  
#  lapply(cf1.central$denseR, function(i) which(i==oLL)) %>% unlist -> xXcf1.d
#  points(x=xXcf1.d , y=LL[oLL[xXcf1.d]], pch=0)
  
  legend('topleft', legend = c('Chipman 1, best tree per cluster','Chipman 1, central tree per cluster') 
         , col=c('yellow','yellow'), pch=c('x','o'))
  legend('bottomright' , legend=paste(metric, parameter$cutoff,sep=', '))
}

LL[cf1.best$forest] %>% mean ; calcLogloss(forest=subforest(rg$forest, cf1.best$forest), df=data.test)
LL[cf1.central$forest] %>% mean ; calcLogloss(forest=subforest(rg$forest, cf1.central$forest), df=data.test)

length(cf1.best$forest) ; length(cf1.best$denseR)
cf1.best$forest %in% cf1.central$forest %>% table
cf1.central$forest %in% cf1.best$forest %>% table
base::intersect(cf1.best$forest,cf1.central$forest) # which cluster do they belong to??

# clusters with only 1 tree => best is always central
cf1.best$hc.clus %>% table %>% (function(x) which(x==1)) %>% length
# clusters with 2 trees, for (around) half of them central should be best
cf1.best$hc.clus %>% table %>% (function(x) which(x==2)) %>% length

##########################################################################
#### compare models by plotting the subforests (mds)
#### probably not insightful as mds in k=2 dimensions has too much stress
#### and cannot represent original distances / dissimilarities

mp <- cmdscale(d=dm, k=2) # for plotting we don't need eigenvalues, just the points
nT <- nrow(dm)
par(mar=c(1,1,1,1))
#col <- c('yellow','orange','red')
ct <- 1

plot(mp , col='black' , pch=0 , axes=F , xlab='' , ylab='')
box(which = "plot", lty = "solid")
legend('bottomleft', legend=paste(metric, parameter$cutoff, sep=', '))

for(f in list(cf1$forest, cf2$forest, mf$forest)){
#col <-  rep(0,nT)
#col[f] <- ct
#pch <- rep(0, nT)
plot(mp[f,] , col=ct , pch=0 , axes=F , xlab='' , ylab='')
box(which = "plot", lty = "solid")
ct <- ct+1
}

