# 23.2.2022
# code/chipmanSelection/02_single_method.R
# parameter optimization for each metric for the cipman 2 algorithm, adding diverse trees

# loops calc_LL_for_selection over parameter as a matrix

# validation data Swiss -> check if still valid !

rm(list=ls())

library(dplyr)
library(ranger)
library(cluster)
library(xtable)

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

# set test data by name : VA, Swiss or Hung
data.test.name <-  'Hung'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name

#paste(rep(c('LL.test.chip.','I.'),4),rep(c('d0','d1','d2','sb'),each=2),sep='')

# to base the result on more bootstraps
folder <- 'data/nursery'
files <- list.files(folder)#[1:2]
# dir(folder)

# give parameter for each dissimilarity
parameter <- list('d0'=list('cutoff'=0.1, 'sizeSF'=5)
, 'd1'=list('cutoff'=0.5, 'sizeSF'=5)
,'d2'=list('cutoff'=0.5, 'sizeSF'=1)
,'sb'=list('cutoff'=0.3, 'sizeSF'=5)
)

# names(parameter)

collector.f <-  list()
ct.f <-  1 # counter for the above collector 
for(file in files){
  # run loops over doc loaded from file
  load(paste(folder,file, sep='/')) # loads doc
  collector.f[[ct.f]] <- calc_LL_for_selection(doc , parameter)
  ct.f <-  ct.f+1
}
et02 <- data.frame(bind_rows(collector.f))

# save(et02 , file='data/chipman/chipman_2_5trees*.rda')
# load('data/chipman/chipman_2_50trees*.rda')

#et02 <-  et02 %>% select(!ends_with('d2'))

#et02 %>% filter(starts_with('size'))
names(et02)
par(mar=c(2,4,2,1)+0.2)
et02 %>% 
  select(starts_with('LL')) %>%
  rename('d0'=ends_with('d0') , 'd1'=ends_with('d1') , 'd2'=ends_with('d2') , 'sb'=ends_with('sb')) %>%
  boxplot(ylab='logloss'
        , main=paste('Chipman 2 forests wrt dissimilarities\n(',50*length(files),' simulations)', sep='')
        # , main=paste(attr(et,'metric'), 'metric ,', num.clusters , 'clusters ,', sizeSF, 'trees in forest\nfor non-default')
)

et02 %>%
  select(starts_with('size.')) %>%
  rename('d0'=ends_with('d0') , 'd1'=ends_with('d1') , 'd2'=ends_with('d2') , 'sb'=ends_with('sb')) %>%
  boxplot( ylab='number of trees'
        , main=paste('Sizes of Chipman 2 forests wrt dissimilarities\n(',50*length(files),' simulations)', sep='')
        # , main=paste(attr(et,'metric'), 'metric ,', num.clusters , 'clusters ,', sizeSF, 'trees in forest\nfor non-default')
)


#names(et)
et02 %>%
  #select(starts_with('LL')) %>%
  summarize_all(function(x) c(mean(x), sd(x))) %>% 
  t %>% xtable -> et.xt

#et.xt
digits(et.xt) <- 4
et.xt

