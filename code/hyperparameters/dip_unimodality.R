# calculate the dip statistic

rm(list=ls())

#install.packages('diptest')
library(diptest)
library(xtable)

# option A
#load('data/dms_for_a_forest_500_trees.rda') # loads DM

# option B
load('data/nursery02/nursery02_01.rda') # loads DM with 4 dissimilarity matrices, one for each dissimilarity

names(DM)
metric <- names(DM)[1] ; metric
dm <-  DM[[metric]]

calcDip <- function(i) dip(dm[i,-i]) # uses dm from parent environment
Vectorize(calcDip)(1:100) %>% hist(xlim=c(0,1))

# names(DM)
if(!exists('metrices')){
  if(is.null(names(DM))){
    metrices <-  c('d0','d1','d2','sb')
  }else{
    metrices <- names(DM)
  }
}

par(mar=c(4,4,3,2)+0.1)

for(i in 1:4){
  dm <- DM[[i]]
  Vectorize(calcDip)(1:500) -> vdip 
  vdip %>% 
    hist(xlim=c(0,1)
         , main=paste('histogram of dip statistic for',metrices[[i]])
    )
  c(min(vdip), mean(vdip),max(vdip)) %>% print
  ecdf(vdip)%>%plot(main=paste('empirical distribution function of dip statistic\non all trees of forest with distances measured in',metrices[[i]]))
  density(vdip)%>%plot(main=paste('density of dip statistic on all trees of forest\ndissimilarity',metrices[[i]]))
}

# fair enough, but we also need the p-values of the dip test

# the dip test is a test for multimodality
# the hypothesis is unimodality, if it is rejected (on whatever level) we claim we found multimodality
# a small value of the statistic implies the regular behaviour, strengthens the hypothesis
# a larger value for the statistic is a vote against the hypothesis.
# the p-value tells us, when the vote is strong enough to reject the hypothesis (p-value relative to level)

metric <- metrices[2]
dm <- DM[[metric]]
doc.dip <- data.frame(matrix(0,nrow(dm),2))

for(i in 1:nrow(dm)){
  # for(i in 1:5){
  #print(i)
  dt <-   dip.test(dm[i,-i])
  doc.dip[i,] <-  c(dt$statistic, dt$p.value)
}
doc.dip %>% apply(2,function(x) c(min(x), mean(x), max(x))) %>% t

table(doc.dip[,2]<0.05)
min(doc.dip[,2])
# option A results
####### reject unimodality:
# d0 dissim : p-value 0 -> hypothesis can be rejected on any level. Trees fall into natural clusters
# sb : p-value 0 -> hypothesis can be rejected on any level. Trees fall into natural clusters
####### no way to be unimodal
# d1 : p-values larger than 0.34 (min p-value) always accept the Hypothesis : really unimodal!
# d2 : p-values larger than 0.17 (min p-value) always accept the Hypothesis : really unimodal!

# It is very clear: either reject on any level, or quite large p-values, 0.17 , 0.34. No fuzzy inbetween at 6%, 7%, or 4.5%

#### visualize ####

mds<-function(D, xylim=FALSE ,main=NULL , col=NULL, pch=NULL){
  m <- cmdscale(D, eig = TRUE, k = 2)
  x <- m$points[, 1]
  y <- m$points[, 2]
  
  if(xylim){
    xlim<- c(-1,1)
    ylim<- c(-1,1)
  }else{
    xlim<-range(x)
    ylim<-range(y)
  }
  
  if(is.null(main)) main<-'mds'
  if(is.null(col)) col <- 1
  if(is.null(pch)) pch <- 1
  
  plot(x, y,  xlim = xlim, ylim=ylim , main=main, col=col, pch=pch)
  #text(x, y, pos = 4, labels =1:nrow(D))
  idx.min <- which.min(apply(DM[[nm]],2,sum))
  plot(x=x[c(1:3,idx.min)], 
       y=y[c(1:3,idx.min)],  
       xlim = xlim, 
       ylim=ylim , 
       main=main, 
       col=c(1,1,1,2), 
       pch=c(1,1,1,2))
}

for(m in 1:4){
  dm <- DM[[m]]
  
  doc <- data.frame(matrix(0,nrow(dm),2))
  
  for(i in 1:nrow(dm)){
    dt <-   dip.test(dm[i,-i])
    doc[i,] <-  c(dt$statistic, dt$p.value)
  }
  names(doc) <- c('dip','pval')
  
  # falling curve, not linear
  # plot(doc$dip , doc$p)
  
  plevel <-  0.01
  
  mds(dm 
    , main=paste('forest: mds for',metrices[[m]],'\nred : reject hypothesis of unimodality, p-level',plevel)
    , col=1+1*(doc$pval<plevel))
}

#### create a LaTex table ####

# create doc in long format to work well with dplyr from the tidyverse

doc <- data.frame(matrix(0,length(DM)*nrow(DM[[1]]),3))
ct <- 1 
for(m in 1:4){
  dm <- DM[[m]]
  
  for(i in 1:nrow(dm)){
    dt <-   dip.test(dm[i,-i])
    doc[(ct-1)*nrow(dm)+i,] <-  c(m,dt$statistic, dt$p.value)
  }
  ct <- ct+1
}
names(doc) <- c('metric','dip','pval')

doc$metric <- metrices[doc$metric]
doc %>% 
  group_by(metric) %>% 
  #filter(pval<0.05) %>%
  #filter(dip>0.1) %>%
  #summarise(mean(dip),mean(pval), length(which(pval<0.05)) , min(pval)) %>%
  summarise(length(which(pval<0.05)) , length(which(pval<0.01))) %>% 
  xtable -> xt

#digits(xt) <- 3
xt

dm <- DM[[3]]
dt <- dip(dm[i,-i], full.result=TRUE)
dt$mj %>% hist
dt$dip 
dip.test(dm)
plot(dt, main='')

