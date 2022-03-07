# code / ChipmanSelection / 04_vis_single_simulations.R

# loads a single file from the nursery : 50 simulation only
# looks at each of the single simulation 
# and how sizes and logloss changes with the cutoff paramter (diversity)

# for moderate parameters 0.1 .. 0.7 logloss is stable and decline in size is convex

rm(list=ls())

library(dplyr)
library(ranger)
library(cluster)

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # load calc_LL_for_selection

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

# set test data by name : VA, Swiss or Hung
data.test.name <-  'Swiss'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name

folder <- 'data/nursery'
files <- list.files(folder)
file <- files[1]
load(paste(folder, file,sep='/'))

cutoff <- seq(0, 1,by=0.01)
# for Chipman to stop early , set sizeSF
# if sizeSF=500 the algorithm finishes with an Chipman forest of unknown size
sizeSF <- rep(500,length(cutoff)) 

parameter <-  data.frame(cbind(cutoff, sizeSF))

# where to load colc_LL_for_selection from ??
# it should return a / many columns named size.d0, size.d1, ..
# load from 03_hpo.R
X1 <-  calc_LL_for_selection(doc , list('cutoff'=0.5, 'sizeSF'=500))

# no linear connection for size and logloss , picked one parameter
names(X1)
apply(X1,2,function(x) c(mean(x),sd(x))) %>% t

metrices <-  c('d0','d1','d2','sb')
for(metric in metrices){
  x <- X1[,paste('size.',metric,sep='')]
  y <- X1[,paste('LL.test.chip.',metric,sep='')]
  plot(y~x 
       , main=paste('Chipman forest wrt ', metric ,', parameter 0.8', sep='')
       , xlab='size of Chipman forest'
       , ylab='logloss')
  lm(y~x) %>% summary %>% print
}

collector.p <-  list()
ct.p <-  1

for(p in 1:nrow(parameter)){
  res <-  calc_LL_for_selection(doc , parameter[p,])[,c(4,6)] #c(4,6) selects LL and size for d1 metric
  collector.p[[ct.p]] <- cbind(res , 'sim'=1:50, 'cutoff'=parameter[p,'cutoff'])
  ct.p <- ct.p+1
}

Y1 <- data.frame(bind_rows(collector.p))

Y1 # rows 1:50 have smallest cutoff 


#### plots, natural order of simulations
dim(Y1)
Y1 %>% 
  ggplot(aes(x=cutoff,y=size.d1,colour=sim)) +
  geom_line(aes(group=sim))

Y1 %>% 
  ggplot(aes(x=cutoff,y=LL.test.chip.d1,colour=sim)) +
  geom_line(aes(group=sim))

######  

nParam <- nrow(parameter)

# logloss over cutoff
# order by initial logloss / logloss at smallest cutoff
Y2 <-  Y1
oLL <- order(Y1[1:50,'LL.test.chip.d1'])
Y1[oLL,'LL.test.chip.d1']
#order by value at other cutoff point
oLL <- order(Y1[Y1$cutoff==0.1,'LL.test.chip.d1'])
Y2  %>% 
  filter(cutoff==0.1) %>%
  select(starts_with('LL')) %>%
  unlist %>%
  (function(x) x[oLL]) %>%
  plot # to check the ascending order induced by oLL

for(i in 0:49){ 
#  i <- 2
  Y2[i*nParam+(1:nParam), ]  <- Y1[Y1[,'sim']==oLL[i+1],]
  Y2[i*nParam+(1:nParam),'sim'] <- i+1
}
Y2 %>% 
  #filter(cutoff >0) %>% 
  filter(cutoff<0.8) %>%
  ggplot(aes(x=cutoff,y=LL.test.chip.d1,colour=sim)) +
  geom_line(aes(group=sim)) +
  labs(title="fluctuations of logloss over cutoff param. for Chipman 2 forest\n(dissim d1, 50 simulations)",
       x= "cutoff parameter",
       y= 'logloss',
       color='simulation\n(reordered)')

# size over cutoff
#####
# order by initial size / size at smallest cutoff
Y2 <-  Y1
# simulations ordered by initial size
oS <- order(Y1[1:50,'size.d1'])
Y1[oS,'size.d1']
#order by value at other cutoff point
oS <- order(Y1[Y1$cutoff==0.1,'size.d1']) # works for some 0.1 , 0.5 but not for others 0.3 ?

for(i in 0:49){ 
  #  i <- 2
  Y2[i*nParam+(1:nParam), ]  <- Y1[Y1[,'sim']==oS[i+1],]
  Y2[i*nParam+(1:nParam),'sim'] <- i+1
}
Y2 %>%
  filter(cutoff >0.1) %>% 
  filter(cutoff<0.8) %>%
  ggplot(aes(x=cutoff,y=size.d1,colour=sim)) +
  geom_line(aes(group=sim)) +
  labs(title=paste('size mostly decreasing in cutoff parameter for Chipman 2 forest\n(dissim d1 , 50 simulations)'),
       x="cutoff parameter",
       y= 'size of Chipman forest',
       color='simulation\n(reordered)')

(function(co) Y1 %>% 
    filter(cutoff==co) %>% 
    select(starts_with('LL')) %>% 
    unlist %>% 
    (function(x){c(mean(x),sd(x), min(x), max(x))})) %>% 
  lapply(X = seq(0,1,0.1)) %>%
  simplify2array -> XA
XA

par(mar=c(4,4,3,1)+0.1)
x <- unique(Y1$cutoff)
plot(x,XA[1,], ylim=c(0,max(XA[4,])), ylab='logloss', xlab='parameter', type='l')
points(x,XA[1,]+ XA[2,], type='l', col='yellow')
points(x,XA[1,]- XA[2,], type='l', col='yellow')  
points(x,XA[3,], type='l' , col='orange')
points(x,XA[4,], type='l' , col='orange')
legend('topleft', legend=c('mean','mean+/-sd','min, max')
       , col=c('black','yellow','orange')
       , pch='--'
       , cex =0.8)

