# code / ChipmanSelection / 04_vis_single_simulations.R

rm(list=ls())

library(dplyr)
library(caret)
library(ranger)
library(effects)
library(cluster)

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

# set test data by name : VA, Swiss or Hung
data.test.name <-  'Swiss'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name

folder <- 'data/nursery'
files <- list.files(folder)
file <- files[2]
load(paste(folder, file,sep='/'))

cutoff <- seq(0.1, 0.7,by=0.05)
sizeSF <- rep(500,length(cutoff))

parameter <-  data.frame(cbind(cutoff, sizeSF))

X1 <-  calc_LL_for_selection(doc , parameter[3,])

names(X1)
apply(X1,2,function(x) c(mean(x),sd(x))) %>% t

plot(X1$LL.test.chip.d0~X1$size.d0, main=parameter[3,])


collector.p <-  list()
ct.p <-  1

for(p in 1:nrow(parameter)){
  res <-  calc_LL_for_selection(doc , parameter[p,])[,c(4,6)]
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
# order by initial logloss / logloss at smallest cutoff
Y2 <-  Y1
oLL <- order(Y1[1:50,'LL.test.chip.d1'])
Y1[oLL,'LL.test.chip.d1']
for(i in 0:49){ 
#  i <- 2
  Y2[i*13+(1:13), ]  <- Y1[Y1[,'sim']==oLL[i+1],]
  Y2[i*13+(1:13),'sim'] <- i+1
}
Y2 %>% 
  ggplot(aes(x=cutoff,y=LL.test.chip.d1,colour=sim)) +
  geom_line(aes(group=sim)) +
  labs(title="fluctuations of logloss over cutoff param. for Chipman forest (dissim d1)",
       x="cutoff parameter",
       y= 'logloss',
       color='simulation\n(reordered)')

#####
# order by initial size / size at smallest cutoff
Y2 <-  Y1
oS <- order(Y1[1:50,'size.d1'])
Y1[oS,'size.d1']
for(i in 0:49){ 
  #  i <- 2
  Y2[i*13+(1:13), ]  <- Y1[Y1[,'sim']==oS[i+1],]
  Y2[i*13+(1:13),'sim'] <- i+1
}
Y2 %>% 
  ggplot(aes(x=cutoff,y=size.d1,colour=sim)) +
  geom_line(aes(group=sim)) +
  labs(title="size mostly decreasing in cutoff parameter for Chipman forest (dissim d1)",
       x="cutoff parameter",
       y= 'size of Chipman forest',
       color='simulation\n(reordered)')

