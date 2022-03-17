# 23.2.2022
# code/chipmanSelection/06_chipman_1 clustering enforce size.R
# no representation
# no parameter needed to set alpha

# take better half of trees, cluster into 5 or 50 clusters, take best tree per cluster

rm(list=ls())

library(dplyr)
library(ranger)
library(cluster)
library(xtable)

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # calc_oLL , calc_chipForest_1_enforce , calc_LL_for_selection

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

# set test data by name : VA, Swiss or Hung
data.test.name <-  'Swiss'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name

sizeSF <- 5

# to base the result on more bootstraps
folder <- 'data/nursery'
files <- list.files(folder)[1:2]
# dir(folder)
collector <-  list()
ct <-  1 # counter for the above collector

for(file in files){
  # run loops over doc loaded from file
  load(paste(folder,file, sep='/'))
  collector[[ct]] <- calc_LL_for_selection(doc, parameter=NULL)
  ct <-  ct+1
}
et02 <- bind_rows(collector)

# too much information
#et02 %>% 
#  select(starts_with('LL.')) %>%
#  boxplot( ylab='logloss'
#        , main=paste('basecases and Chipman forests wrt dissimilarities\n(',50*length(files),' simulations)', sep='')
#       # , main=paste(attr(et,'metric'), 'metric ,', num.clusters , 'clusters ,', sizeSF, 'trees in forest\nfor non-default')
#)
#abline(median(et02$LL.random),0,col='red')
#abline(median(et02$LL.default),0, col='green')

# reduce information
par(mar=c(4,4,3,2)+0.1)
et02 %>% 
  select(starts_with('LL.test.chip.')) %>%
  select(!ends_with(paste('.I.',c('d0','d1','d2','sb'),sep=''))) %>%
  rename(c('d0'='LL.test.chip.d0','d1'='LL.test.chip.d1','d2'='LL.test.chip.d2','sb'='LL.test.chip.sb')) %>%
  boxplot( ylab='logloss'
           , xlab='dissimilarity for representation in Chipman I forest'
           , main=paste('modified Chipman I forests wrt dissimilarities\n(',50*length(files),' simulations)', sep='')
           # , main=paste(attr(et,'metric'), 'metric ,', num.clusters , 'clusters ,', sizeSF, 'trees in forest\nfor non-default')
  )
abline(median(et02$LL.random),0,col='red')
abline(median(et02$LL.default),0, col='green')
legend('topleft' 
       , legend=c('mean logloss default forest, size 500', paste('mean logloss small forest, size' , sizeSF))
       , cex =0.8
       , col=c('green','red')
       , pch = '-')

# count NA : unimodal I has sizeR = NA
# calculate the size of the final Chipman forest
# if multimodal it is the size of R , if unimodal it is I
metrices <- c('d0','d1','d2','sb')
for(metric in metrices){
  et02 %>% 
    select(ends_with(metric)) %>%
    select(starts_with(c('I.','size'))) %>%
    apply(1,min, na.rm=T) -> et02[,paste('sizeCh.',metric,sep='')] 
}


et02 %>% 
  summarize_all(function(x) c(mean(x, na.rm=T), sd(x, na.rm=T))) %>% 
  t %>% xtable -> et.xt
#et.xt
digits(et.xt) <- 4
et.xt

et02 %>% 
  summarize_all(function(x) c(mean(x, na.rm=T), sd(x, na.rm=T))) %>% 
  t %>% xtable -> et.xt
#et.xt
digits(et.xt) <- 4
et.xt

#names(et)
et02 %>% 
  select(starts_with('LL.')) %>%
  summarize_all(function(x) c(mean(x), sd(x), median(x), IQR(x))) %>% 
  t %>% xtable -> et.xt

#et.xt
digits(et.xt) <- 4
et.xt

et02 %>% 
  select(starts_with('sizeR.')) %>%
  summarize_all(function(x) c(mean(x), sd(x))) %>% 
  t %>% xtable -> et.xt

#et.xt
digits(et.xt) <- 4
et.xt

#sizes of representing subforests before clustering
et02 %>% 
  select(starts_with('I.')) %>%
  rename(c('d0'='I.d0','d1'='I.d1','d2'='I.d2','sb'='I.sb')) %>%
  boxplot(main=paste('Chipman I: sizes of representing dense subforests before clustering\n(representation parameter', paste(unlist(parameter.alpha),collapse=',') ,')')
          , xlab='dissimilarity for representation in Chipman I forest'
          , ylab=' number of trees')

#sizes after clustering
et02 %>% 
  select(starts_with('sizeR.')) %>%
  rename(c('d0'='sizeR.d0','d1'='sizeR.d1','d2'='sizeR.d2','sb'='sizeR.sb')) %>%
  boxplot(main=paste('Chipman I: sizes of final forest\n(representation parameter', paste(unlist(parameter.alpha),collapse=',') ,')')
          , xlab='dissimilarity for representation in Chipman I forest'
          , ylab=' number of trees')


metric <- 'd0'
plot(et02[,paste('LL.chip.', metric,sep='')]~ et02[,paste('sizeR.', metric , sep='')])
abline(lm(et02$LL.chip.d2~ et02$sizeR.d2)$coeff)
lm(et02$LL.chip.d2~ et02$sizeR.d2) %>% summary
 
