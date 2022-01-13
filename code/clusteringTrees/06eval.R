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

load('data/cluster/06_75.rda')
# load('data/cluster/06_87.rda')

# combine several et from different simulations (each around 50 loops, which is a small number...)

x <- list()
folder <- 'data/cluster'
lf <- list.files(folder)
for(file in lf){
  load(paste(folder,file,sep='/'))
  print(paste('loaded et with dim(et)=',dim(et)))
  print(paste('names(et):',names(et)))
  if(length(x)==0){
    x <- et
  }else{
      x <- rbind(x,et)
      }
  }
et <- x

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

et %>% 
  filter(type %in% c('random',paste('clustering',1:4,sep=''))) -> et.s # selected : subforests built on clusters and randomly selected subforests
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
et %>% 
  group_by(metric, type) %>%
  summarize(mean(logloss), sd(logloss)) %>% 
  data.frame() -> et.xt
dim(et.xt)
names(et.xt) <- c('metric','clustering strategy','mean logloss' , 'sd logloss')

et.xt %>%
  xtable -> et.xt
digits(et.xt) <- 3
et.xt
# remove multiple rows of randomly selected , without clustering
# start with default forest and randomly selected forest


#### problem : evrything is too close. 
#### t-tests: can't even say that default is better than random?
et %>% 
  filter(metric=='none' & type=='default') %>%
  select(logloss) ->Ld
et %>% 
  filter(metric=='d0' & type=='random') %>%
  select(logloss) ->Lr

(Lr-Ld) %>% unlist %>% hist(breaks=20
                            , xlab='logloss overshoot'
                            , main='histogram for logloss overshoot \nfor randomly sampled subforest (compare to default forest)')
# Hypothesis: randomly selected forest has larger logloss than default
t.test(Lr-Ld, alternative='g') # reject on a 5% level / p value 0.025
# this is what we wanted , what we expected

et %>% 
  filter(metric=='d1' & type=='clustering3') %>%
  select(logloss) ->Ls
# we would like to see our selection as better than random
# Hypo: randomly selected forest has larger logloss than carefully selected forest
t.test(Lr-Ls, alternative='g')


doc2 <- data.frame(matrix(NA,16,4))
ct <- 1
for(m in metrices){
  for(t in unique(et$type)){
    print(paste(m,t))
    et %>% 
      filter(metric==m & type==t) %>%
      select(logloss) -> Ls
    if(nrow(Ls)>0){
      print(paste(m,t))
      doc2[ct,] <- c(m
                     , t
                     , t.test(Lr-Ls, alternative='g')$p.value # would like to reject
                     , t.test(Ls-Ld, alternative='g')$p.value # expect to reject
                     )
      ct <- ct+1
    }
  }
}

names(doc2) <- c('metric','type','p.value.r','p.value.d','p.value.d + 0.01')
doc2[,3] %>% as.numeric %>% round(4) -> doc2[,3]
doc2[,4] %>% as.numeric %>% round(4) -> doc2[,4]
#doc2[,5] %>% as.numeric %>% round(4) -> doc2[,5]
#doc2[,6] %>% as.numeric %>% round(4) -> doc2[,6]
doc2 %>% xtable -> doc2.xt
digits(doc2.xt) <- 5
doc2.xt

#### this is to check if in the selected forests there is a connection of mean dissimilarity and logloss
#### but we should check it on regular randomly selected forests, not on the results of clustering

et.s %>% 
  group_by(type,metric) %>% 
  summarize('lm p.value'=f2(logloss,mDiss)$pv
            , 'lm intercept'=f2(logloss,mDiss)$coeff[[1]]
            , 'lm slope'=f2(logloss,mDiss)$coeff[[2]] # cannot use a vector valued function in summarise yet
            , 'cor'=cor(logloss,mDiss)) -> et.s.t # the table for the evaluation table of selected
et.s.t %>% xtable -> et.s.xt

align(et.s.xt) <- xalign(et.s.xt)
digits(et.s.xt) <- 3
display(et.s.xt) <- xdisplay(et.s.xt)
et.s.xt

# looking into 16 p-values, probably some are small , some are large ...
et.s.t$`lm p.value` %>% hist(breaks=10, main='p values for linear models for\ndifferent metrices and cluster strategies')
runif(16,0,1) %>% hist(breaks=10)
rnorm(16) %>% 
  (function(x) (x-min(x))/(max(x)-min(x))) %>% 
  hist(breaks=10, main='scaled normal random numbers')
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

