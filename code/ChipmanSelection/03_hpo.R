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
library(xtable)

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA
Swiss <- Swiss[,1:11]

# corrected Swiss
#corSwiss <- Swiss[,1:11]
#corSwiss$STDepression <-  corSwiss$STDepression - min(corSwiss$STDepression)
#range(corSwiss$STDepression)

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection


hpo <- function(method, parameter, data){
  
  data.test <<- get(data) # data should be a variable name , character string 
  # get data from parent environment
  # is not passed to calc_LL_for_selection
  # is accessed from within the function (no good programming style!!)
  
  # to base the result on more bootstraps
  folder <- 'data/nursery02'
  files <- list.files(folder)[1:2]
  # dir(folder)
  
  # for the parameter : collector and counter
  collector.p <-  list()
  ct.p <-  1
  
  for(p in 1:nrow(parameter)){
    #print(paste('parameter ', p ,'at', Sys.time()))
    
    collector.f <-  list()
    ct.f <-  1 # counter for the above collector 
    for(file in files){
      # run loops over doc loaded from file
      load(paste(folder,file, sep='/')) # loads doc
      #print( parameter[p,])
      collector.f[[ct.f]] <- calc_LL_for_selection(method=method
                                                   , doc=doc 
                                                   , parameter=parameter[p,])
      ct.f <-  ct.f+1
    }
    
    apply(bind_rows(collector.f),2,mean) %>% # remove NA before taking the mean??
      c(parameter[p,]) %>% 
      unlist -> 
      collector.p[[ct.p]]
    
    ct.p <- ct.p+1
  }
  return(list('et'=data.frame(bind_rows(collector.p))
              , 'call'= list('method'=method, 'parameter'=parameter, 'data'=data)
              , 'simulations'=files)
         )
}

method='chip1'

#### build parameter
####################
# chipman 1 simplified does not need a parameter. 
# Multiple cutoff parameters just repeat the same code and give the same results.
cutoff <- seq(0, 0.9, by=0.1)
# setting sizeSF to the number of trees (500)  generates a Chipman forest 
# that only stops when all trees are represented (at the specified level)
# cutoff 0 means only trees with minimal dissim represent each other. 
# Happens for some trees under d0, happens for exactly 2 trees (0r 3,4?? few!) under d1,d2,sb
sizeSF <- rep(500,length(cutoff))
parameter <-  data.frame(cbind(cutoff, sizeSF))

if(method %in% c('chip1','chip1_simplified')){
  parameter$selection <- rep('central',length(cutoff))
}

#### call function hpo
######################
res1 <- hpo(method=method
            , parameter=parameter
            , data='Swiss'
            ) # result

et03 <- res1$et
files <- res1$simulations

# next run:
#save(res1 , file='data/chipman/hpo_method*.rda')
# save(hpo_mei_stopped_5 , file='data/hpo_mei_stopped_5trees*.rda')
# load('data/chipman/hyperparameter_cutoff_50trees.rda')

et03 %>% t %>% xtable -> xtb.ch2
digits(xtb.ch2) <-  4
xtb.ch2

# if neccessary , for chipman 1 with selection as character -> all columns are character
# exclude selection column in names(et03)
for(nam in names(et03)){
  et03[,nam] <- as.numeric(et03[,nam])
}

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
     , main=paste('Chipman forests \n(method', res1$call$method , ',', 50*length(files) , ' simulations)')
     , ylim=ylim
     , xlab='representation parameter (quantile of dissimilarity)' 
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
     , type='l'
     , main=paste('Chipman forests \n(method', res1$call$method , ' , ' , 50*length(files), 'simulations)')
     , ylim=ylim
     , xlab='representation parameter (quantile of dissimilarity)' 
     , ylab='representing trees')
#points(cutoff[wm], et03.s[wm,1],col=1, pch=19 )
for(k in 2:4){
  points(cutoff
         , et03.s[,k]
         , type='l'
         , col=k)
 # points(cutoff[wm], et03.s[wm,k],col=k, pch=19 )
}
legend( 'topright'
  #'bottomright'
       , legend=c('d0','d1','d2','sb') 
       , pch=1
       , col=1:4
       , cex=0.8)

# plot for sizes over parameter
metrices <- c('d0','d1','d2','sb')
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
     , main=paste('Chipman forests \n(method', res1$call$method , ', ', 50*length(files) , ' simulations)')
     , ylim=ylim
     , xlab='representation parameter (quantile of dissimilarity)' 
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

# free lunch for d0
# for selected metric : size and logloss in one plot
{
  par(mar=c(4,4,1,4)+0.2)
  #par(mar=c(4,4,3,4)+0.2)
  metric <- 'd0'
  et03 %>% select(cutoff , ends_with(metric)) -> etA
  x <- etA$cutoff
  y <- etA %>% select(starts_with('size')) %>% unlist
  z <- etA %>% select(starts_with('LL')) %>% unlist
  plot(y~x 
     #, main=paste('Chipman forest \n(method', res1$call$method , ' , dissimilarity' , metric , ' , ' , 50*length(files) , ' simulations)')
     , type='l' 
     , ylab='size: number of trees in Chipman forest' 
     , xlab='representation parameter (quantile of dissimilarity)' 
     , cex.lab=1.2)
  par(new = TRUE)
  plot(z~x , axes = FALSE, bty = "n", xlab = "", ylab = "", type='l' , col='blue')
  axis(side=4, at = pretty(range(z)), col='blue' , col.axis='blue',cex=1.2)
  mtext("success: logloss", side=4, line=3 , col='blue', cex=1.2)
}

# plot logloss over size
# must be some kind of smoothing of the many pairs of (size, logloss) we generate over the 100 or 1000 simulations

par(mar=c(4,4,3,1)+0.1)
ylim <-  range(et03 %>% select(starts_with('LL')))
plot(et03$cutoff
     , et03$LL.test.chip.d0
     , type='l'
     , ylim=ylim 
     , xlab='representation parameter (quantile of dissimilarity)' 
     , ylab='logloss'
     , main=paste('Chipman forest\n(method', res1$call$method , ')', sep=' ')
)
points(et03$cutoff, et03$LL.test.chip.d1, type='l', col=2)
points(et03$cutoff, et03$LL.test.chip.d2, type='l', col=3)
points(et03$cutoff, et03$LL.test.chip.sb, type='l', col=4)
legend('topleft', legend=c('d0','d1','d2','sb') , pch=1, col=1:4 , cex=0.8)

# zoom in
par(mar=c(4,4,1,1)+0.1)
# par(mar=c(4,4,3,1)+0.1)
N <- 6 # number of parameters in plot
ylim <-  range(et03[1:N,] %>% select(starts_with('LL')))
x <- et03$cutoff[1:N]
y <- et03$LL.test.chip.d0[1:N]
plot(x
     , y
     , type='l'
     , ylim=ylim 
     # , main=paste('Chipman forest\n(method', res1$call$method , ')', sep=' ')
     , xlab='representation parameter (quantile of dissimilarity)'
     , ylab='logloss'
     , axes =F )
axis(1)
axis(2 , at=seq(round(min(ylim),2),round(max(ylim),2),by=0.01))
box()
points(x, et03$LL.test.chip.d1[1:N], type='l', col=2)
points(x, et03$LL.test.chip.d2[1:N], type='l', col=3)
points(x, et03$LL.test.chip.sb[1:N], type='l', col=4)
legend('topleft', legend=c('d0','d1','d2','sb') , pch=1, col=1:4 , cex=0.8)

et03$LL.test.chip.d0 %>% hist(breaks=50)
