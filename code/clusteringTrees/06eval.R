# mail Ulrike 6.11.2021
# Vielleicht ist es sinnvoller, getrennte Modelle für die verschiedenen metrics 
# aufzustellen ? Dann könntest Du auch die Dissimilarities unskaliert belassen 
# und erst einmal mit einem einfachen Streudiagramm logloss gegen dissimilarity 
# beginnen. 
# Wenn da nichts zu sehen ist, ist eigentlich schon fertig, 
# sonst könntest Du versuchen, das Gesehene zu modellieren, 
# z.B. mittels metric-spezifischen Quartil- oder Quintilgruppen oder polynomisch 
# (z.B. mit der poly()-Funktion). 
# Das müsste doch dem logloss wurscht sein, wie Du letztlich zu den 
# dissimilarities gekommen bist, d.h. die clustering Methode könnte man einfach 
# außen vor lassen (hast Du ja auch im linearen Modell gemacht, mit dieser 
# Begründung finde ich das auch plausibel :-))?
  

rm(list=ls())

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
library(dplyr)
library(ggplot2)

load('data/cluster/06_49.rda')
#load('03_02.rda')
info
moreInfo %>% names
names(et)
dim(et)

# no infinite logloss
(et$logloss==Inf) %>% table

# some subforests are not of the desired size (50) at all
# should we remove them? will give unequal length when doing paired t-test
# should we use a smaller number of clusters, esp. for metric d2?
# for metric d2 the best number of clusters was 2 (and maybe no clustering would be best -> unimodal?!)
et[,c('size','metric')] %>% table %>% addmargins

# number of bootstrap samples
moreInfo$nBs

# number of (sub-) forest per type - metric- combination
et[,c('type','metric')] %>% table


et %>% filter(type %in% c('random','default')) -> et.dr # default and random 
## including the random forests and their mDiss in all 4 dissimilarity metrices, i.e. each randomly selected subforest is measured in each of the 4 dissimilarities
dim(et.dr)
et.dr %>% group_by(type) %>% summarize(mLL=mean(logloss)%>% round(3),sdLL=sd(logloss)%>% round(3)) -> LL
LL

et %>% filter(type %in% c('random',paste('clustering',1:4,sep=''))) -> et.s # selected : subforests built on clusters and randomly selected subforests
dim(et.s)

metrices <-  moreInfo$metrices # moreInfo is in loaded data
calcPValue <- function(a) pf(a[[1]],a[[2]],a[[3]],lower.tail=FALSE)

# check it is working, calcPValue should give the same p-value as the summary
lm(logloss~mDiss, data=et.s) %>% summary -> s1
s1

calcPValue(s1$fstatistic)
#### end of check
f1 <-  function(a,b){
  data <- data.frame(cbind('a'=a,'b'=b))
  lm(a~b , data=data)%>%
    summary %>%
    .$fstatistic %>%
    calcPValue()
}  

f2 <-  function(a,b){
  data <- data.frame(cbind('a'=a,'b'=b))
  lm1 <-  lm(a~b , data=data)
  lm1 %>%
    summary %>%
    .$fstatistic %>%
    calcPValue() -> pv
  return(list('pv'=pv , 'coeff'=lm1$coefficients))
}

# check
x <- et.s[(et.s$metric=='d0')&(et.s$type=='clustering1'),'mDiss']
y <- et.s[(et.s$metric=='d0')&(et.s$type=='clustering1'),'logloss']

lm(y~x , data=data.frame(cbind(x,y))) %>% summary -> s1
# is the same as
# lm(logloss~mDiss, data=et.s[(et.s$metric=='d0')&(et.s$type=='clustering1'),]) %>% summary -> s1
s1

f2(y,x)
#### end of check

####################################
#### this is the relevant table ####
####################################
et.s %>% 
  group_by(type,metric) %>% 
  summarize('lm p.value'=f2(logloss,mDiss)$pv
            , 'lm intercept'=f2(logloss,mDiss)$coeff[[1]]
            , 'lm slope'=f2(logloss,mDiss)$coeff[[2]] # cannot use a vector valued function in summarise yet
            , 'cor'=cor(logloss,mDiss))

# not working
#et.s %>% 
#  group_by(type,metric) %>% 
#  ggplot()+
#  geom_point(aes(x=mDiss, y=logloss))


######################################
#### these are the relevant plots ####
######################################

par(mar=c(4,4,2,1)+0.1)
for(m in metrices){
  for(t in unique(et.s$type)){
    x <- et.s[(et.s$type==t)&(et.s$metric==m), 'mDiss']
    y <- et.s[(et.s$type==t)&(et.s$metric==m),'logloss']
    plot(x=x
         , y=y
         , main=paste(m,t,', model p-value ', round(f1(x,y),3))
         , xlab=paste('mean dissimilarity in',m)
         , ylab='logloss'
        # , ylim=c(0.35,0.55)
        # , col=round(et.s$mSplits,0)
         )
    abline(lm(y~x)$coefficients)
    }
}

