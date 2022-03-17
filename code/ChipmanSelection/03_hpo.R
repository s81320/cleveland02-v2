# 23.2.2022
# code/chipmanSelection/03_hpo.R
# parameter optimization for each metric for the cipman 2 algorithm, adding diverse trees

# loops calc_LL_for_selection over parameter as a matrix

# validation data Swiss -> check if still valid !

rm(list=ls())

library(dplyr)
library(ranger)
library(cluster)
library(diptest)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection

# set test data by name : VA, Swiss or Hung
data.test.name <-  'Swiss'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name

# to base the result on more bootstraps
folder <- 'data/nursery'
files <- list.files(folder)#[1:2]
# dir(folder)

cutoff <- seq(0.0, 0.9, by=0.1)
# setting sizeSF to the number of trees (500)  generates a Chipman forest 
# that only stops when all trees are represented (at the specified level)
# cutoff 0 means only trees with dissim 0 represent each other. Happens only for d0.
sizeSF <- rep(500,length(cutoff))

parameter <-  data.frame(cbind(cutoff, sizeSF))
# parameter$cutoff
# parameter <- list('sizeSF' = 5
#            , 'cutoff'=0.5) # old
collector.p <-  list()
ct.p <-  1

for(p in 1:nrow(parameter)){
  print(paste('parameter ', p ,'at', Sys.time()))

  collector.f <-  list()
  ct.f <-  1 # counter for the above collector 
  for(file in files){
    # run loops over doc loaded from file
    load(paste(folder,file, sep='/')) # loads doc
    #print( parameter[p,])
    collector.f[[ct.f]] <- calc_LL_for_selection(doc , parameter[p,])
    ct.f <-  ct.f+1
  }
  
  apply(bind_rows(collector.f),2,mean) %>% # remove NA before taking the mean??
    c(parameter[p,]) %>% 
    unlist -> 
    collector.p[[ct.p]]
  
  ct.p <- ct.p+1
}

et03 <- data.frame(bind_rows(collector.p))

et03 %>% t %>% xtable -> xtb.ch2
digits(xtb.ch2) <-  4
xtb.ch2

# save(et03 , file='data/chipman/hyperparameter_cutoff_5trees*.rda')
# save(hpo_mei_stopped_5 , file='data/hpo_mei_stopped_5trees*.rda')
# load('data/chipman/hyperparameter_cutoff_50trees.rda')

# makes no sense. et03 is already aggregated, averaged over...
#et03 %>% apply(2,function(x) c(mean(x),sd(x))) %>% t

{
  et03.s <-  et03 %>% select(tidyr::starts_with('LL.'))
  ylim <- range(et03.s)
  # run only with loaded data that has not been created in this session
  # cutoff <- unique(et03$cutoff)
  # sizeSF <-  unique(et03$sizeSF)
  par(mar=c(4,4,3,1)+0.2)
  plot(cutoff
     , et03.s[,1]
     , type='l'
     , main=paste('logloss for cutoff parameters\n(', 50*length(files) ,' simulations)' , sep='')
     , ylim=ylim
     , xlab='cutoff parameter (quantile of dissimilarities)'
     , ylab='mean logloss')
  #wm <- which.min(et03.s[,1])
  #points(cutoff[wm], et03.s[wm,1],col=1, pch=19 )
  for(k in 2:4){
    points(cutoff
         , et03.s[,k]
         , type='l'
         , col=k)
    #wm <- which.min(et03.s[,k])
    #points(cutoff[wm], et03.s[wm,k],col=k, pch=19 )
  }
  legend('topleft', legend=c('d0','d1','d2','sb') , pch=1, col=1:4 , cex=0.8)
}

et03.s <- et03 %>% select(tidyr::starts_with('I'))
ylim <- range(et03.s)
plot(cutoff
     , et03.s[,1]
     , type='b'
     , main=paste('mean number of trees in selection process\n(',50*length(files),' simulations, ', sizeSF[[1]], ' trees)', sep='')
     , ylim=ylim
     , xlab='cutoff parameter (quantile of dissimilarities)'
     , ylab='trees oredered by OOB performance')
#points(cutoff[wm], et03.s[wm,1],col=1, pch=19 )
for(k in 2:4){
  points(cutoff
         , et03.s[,k]
         , type='b'
         , col=k)
 # points(cutoff[wm], et03.s[wm,k],col=k, pch=19 )
}
legend( 'topleft'
  #'bottomright'
       , legend=c('d0','d1','d2','sb') 
       , pch=1
       , col=1:4
       , cex=0.8)

# plot for sizes over parameter
{
  par(mar=c(4,4,3,1)+0.2)
  et03.s <-  et03 %>% select(tidyr::starts_with('size.'))
  ylim <- range(et03.s)
  # run only with loaded data that has not been created in this session
  # cutoff <- unique(et03$cutoff)
  # sizeSF <-  unique(et03$sizeSF)
  plot(cutoff
     , et03.s[,1]
     , type='l'
     , main=paste('size of Chipman forests\n(', 50*length(files) , ' simulations)', sep='')
     , ylim=ylim
     , xlab='cutoff parameter (quantile of dissimilarities)'
     , ylab='mean size')
  #wm <- which.min(et03.s[,1])
  #points(cutoff[wm], et03.s[wm,1],col=1, pch=19 )
  for(k in 2:4){
    points(cutoff
         , et03.s[,k]
         , type='l'
         , col=k)
  #wm <- which.min(et03.s[,k])
  #points(cutoff[wm], et03.s[wm,k],col=k, pch=19 )
  }
  abline(50,0, col='grey')
  abline(5,0, col='grey')
  legend('topright', legend=c('d0','d1','d2','sb') , pch=1, col=1:4 , cex=0.8)
}


# for selected metric : size and logloss in one plot
{
  par(mar=c(4,4,3,4)+0.2)
  metric <- 'd1'
  et03 %>% select(cutoff , ends_with(metric)) -> etA
  x <-  etA$cutoff
  y <-  etA %>% select(starts_with('size')) %>% unlist
  z <-  etA %>% select(starts_with('LL')) %>% unlist
  plot(y~x 
     , main=paste('tradeoff for size and success in Chipman 2 forests\n(dissim',metric, ', mean values for ', 50*length(files) , ' simulations)')
     , type='l' 
     , ylab='size: number of trees in Chipman forest' 
     , xlab='parameter')
  par(new = TRUE)
  plot(z~x , axes = FALSE, bty = "n", xlab = "", ylab = "", type='l' , col='blue')
  axis(side=4, at = pretty(range(z)), col='blue')
  mtext("success: logloss", side=4, line=3 , col='blue')
}

# plot logloss over size
# must be some kind of smoothing of the many pairs of (size, logloss) we generate over the 100 or 1000 simulations

par(mar=c(4,4,2,1)+0.1)
ylim <-  range(et03 %>% select(starts_with('LL')))
plot(et03$cutoff
     , et03$LL.test.chip.d0
     , type='l'
     , ylim=ylim 
     , xlab='cutoff parameter' 
     , ylab='logloss'
     , main='logloss for Meiner forest stopped at 5 trees')
points(et03$cutoff, et03$LL.test.chip.d1, type='l', col=2)
points(et03$cutoff, et03$LL.test.chip.d2, type='l', col=3)
points(et03$cutoff, et03$LL.test.chip.sb, type='l', col=4)
legend('topleft', legend=c('d0','d1','d2','sb') , pch=1, col=1:4 , cex=0.8)
