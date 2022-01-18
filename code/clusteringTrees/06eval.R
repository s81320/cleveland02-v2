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

#load('data/cluster/06_75.rda')
# load('data/cluster/06_87.rda')

# combine several et from different simulations (each around 50 loops, which is a small number...)

x <- list()
folder <- 'data/cluster'
lf <- list.files(folder)
for(file in lf){
  load(paste(folder,file,sep='/'))
  print(paste('loaded et with dim(et)=',dim(et)))
  #print(paste('names(et):',names(et)))
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
  return(list('pv'=pv 
              , 'coeff'=lm1$coefficients 
              , 'rse'=lm1$residuals %>% sd %>% round(4)))
}

# check
x <- et.s[(et.s$metric=='d0')&(et.s$type=='clustering1'),'mDiss']
y <- et.s[(et.s$metric=='d0')&(et.s$type=='clustering1'),'logloss']

print('combare coefficients in summary and f2 , should be the same')
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


#### problem : everything is too close. 
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

names(doc2) <- c('metric','type','p.value.r','p.value.d')
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
  summarize('lm p.value'=f2(logloss,mDiss)$pv ,
             'lm intercept'=f2(logloss,mDiss)$coeff[[1]] , 
             'lm slope'=f2(logloss,mDiss)$coeff[[2]] ,  # cannot use a vector valued function in summarise yet
             'rse'=f2(logloss,mDiss)$rse) -> et.s.t # the table for the evaluation table of selected
et.s.t %>% xtable -> xt

align(xt) <- xalign(xt)
digits(xt) <- 3
display(xt) <- xdisplay(xt)
xt

et.s.t %>% 
  select(metric, type, 'lm p.value') %>%
  tidyr::pivot_wider(names_from=c(type) , values_from = 'lm p.value') %>%
  xtable -> xt

align(xt) <- xalign(xt)
digits(xt) <- 3
display(xt) <- xdisplay(xt)
xt


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
    
    lm1 <- lm(y ~ x)
    plot(x[order(x)],lm1$residuals[order(x)] , main='residuals in simple linear model')
    plot(lm1$fitted.values,lm1$residuals , main='residuals over fitted values simple linear model')
    }
}

######################################################
#### randomly selected forests under all 4 dissim ####
######################################################

# not working
et %>%
  filter(type=='random') %>%
  select(c(bootstrap,metric,mDiss,logloss)) %>%
  tidyr::pivot_wider(names_from=metric, values_from=mDiss) -> x

et %>% 
  filter(metric=='d0' & type=='random') %>%
  select(logloss, mDiss) -> x1
names(x1) <- c('logloss','d0')

for(m in c('d1','d2','sb')){
  et %>% 
    filter(metric==m & type=='random') %>%
    select(logloss, mDiss) -> x2
  names(x2) <- c('logloss',m)
  x1 <- cbind(x1,x2)
}

# check that all loglosses are the same in each row
x1[,c(1,3,5,7)] %>% apply(1,function(y){all(y==y[1])}) %>% table
# table should have all entries True

# the relevant data we will work with
x <-  x1[,c(2,4,6,8)] %>% 
  scale %>% 
  cbind(x1$logloss) %>%
  as.data.frame
names(x)[5] <- 'logloss'

myLM(formula=logloss~ d0*d1*d2*sb , x) # all interactions
myLM(formula=logloss~ d0+sb+d0:d2+d0:d1:sb , x) # interactions that were significant in first model
myLM(formula=logloss~ d0+sb+d0:d1:sb , x) # interactions that were significant in previous model
myLM(formula=logloss~ d0:d1:sb + sb , x)
myLM(formula=logloss~  d0*sb , x)
myLM(formula=logloss~ d0:sb +sb , x)

myLM(logloss~ d0:d1 , data=x)

myLM(logloss~ d2*d0 , data=x) # ok model p-value but no significant coefficients

myLM(logloss~ d0 , data=x)
myLM(logloss~ sb , data=x)
myLM(logloss~ d0+sb ,x)
myLM(logloss~ d2:d0, x)


myLM <- function(formula , data){
  par(mar=c(5,4,3,1)+0.1)
  lm1 <- lm(formula , data) 
  lm1 %>% 
    summary %>% 
#    xtable %>% 
    print
  plot(lm1$residuals ~ lm1$fitted.values 
       , main=paste('residuals in linear model, p-value: '
                    , lm1 %>% summary %>% .$fstatistic %>% calcPValue %>% round(4)
                    , '\n', as.character(formula)[[2]] , as.character(formula)[[1]], as.character(formula)[[3]] , sep='')
       , sub=paste('residuals absolute mean: ', lm1$residuals %>% abs %>% mean %>% round(4)
                     , ' and sd: ', lm1$residuals %>% sd %>% round(4))
         )
  lm1$residuals %>% hist(main='residuals')
  
}

#############################################################
#### mean dissimilarities under clustering ##################
#############################################################

et %>%
  group_by(metric, type) %>%
  summarise(md=mean(mDiss)) %>% 
  filter(type!='default') %>%
  tidyr::pivot_wider(names_from=type, values_from = md) %>%
  xtable



  
