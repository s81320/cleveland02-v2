# calculate the dip statistic

rm(list=ls())

load('data/dms_for_a_forest.rda')

#install.packages('diptest')
library(diptest)
library(xtable)

calcDip <- function(i) dip(dm[i,-i]) # uses dm from parent environment

dm <-  DM[[1]]
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
  min(vdip) %>% print
  ecdf(vdip)%>%plot(main=paste('empirical distribution function of dip statistic\non all trees of forest with distances measured in',metrices[[i]]))
}

# fair enough, but we also need the p-values of the dip test

# the dip test is a test for multimodality
# the hypothesis is unimodality, if it is rejected (on whatever level) we claim we found multimodality
# a small value of the statistic implies the regular behaviour, strengthens the hypothesis
# a larger value for the statistic is a vote against the hypothesis.
# the p-value tells us, when the vote is strong enough to reject the hypothesis (p-value relative to level)

dm <- DM[[2]]
doc <- data.frame(matrix(0,nrow(dm),2))

for(i in 1:nrow(dm)){
  # for(i in 1:5){
  #print(i)
  dt <-   dip.test(dm[i,-i])
  doc[i,] <-  c(dt$statistic, dt$p.value)
}


# d0 dissim : p-value 0 -> hypothesis can be rejected on any level. Trees fall into natural clusters
####### no way to be unimodal
# d1 : always accept the Hypothesis : really unimodal!
# d2 : reject the H for 58 , accept the H for 442
####### it is mulimodal , there are 58 trees that see the mulimodality
# sb : always accept the H , this forest is really unimodal, cannot be divided

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

dm <- DM[[4]]
dt <- dip(dm[i,-i], full.result=TRUE)
dt$mj %>% hist
dt$dip 
dip.test(dm)
plot(dt)

