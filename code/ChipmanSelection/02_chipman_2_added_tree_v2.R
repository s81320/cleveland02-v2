# 23.2.2022
# code/chipmanSelection/03_hpo.R
# parameter optimization for each metric for the cipman 2 algorithm, adding diverse trees

# loops calc_LL_for_selection over parameter as a matrix

# validation data Swiss -> check if still valid !

rm(list=ls())

library(dplyr)
library(ranger)
library(cluster)

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)
source('code/source/subforest.R') # subforest
source('code/source/chipman.R') # for calc_chipForest_2 , calc_LL_for_selection

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

# set test data by name : VA, Swiss or Hung
data.test.name <-  'Swiss'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name

calc_meinForest <- function(dm , forest, oLL, parameter){
  #' calculates the Meiner forest
  #' 
  #' @param parameter is a list with items cutoff and sizeSF 
  #'
  #' returns its logloss on data.test (from higher environment)
  #' returns the rank of the last inspected tree (or the index +1 , like 6 when collecting the first 5 trees)
  #' returns the number of trees in the Meiner forest (enforced to be parameters$sizeSF or smaller)
  
  meinForest <- oLL[1] # start chipForest with the best tree
  loDiv <- quantile(dm,parameter$cutoff) # level of diversity
  #print(cutoff)
  
  # go through trees (ordered by performance) until sizeSF trees are collected / selected
  for(I in 2:length(oLL)){
    if(length(meinForest)==parameter$sizeSF){
      break
    }else{
      trindx <- oLL[I]
      # if new tree far from (current) meinForest , add the best tree (lowest logloss) to the meinForest that represents trindx.
      # if cutoff =0 then each tree is added with certainty -> same as unimodal!
      if(min(dm[meinForest,trindx]) >= loDiv){ 
        base::setdiff(oLL[1:I], meinForest) %>%
          dm[.,trindx] %>%
          (function(x) x<parameter$cutoff) %>%
          which.min %>%
          base::setdiff(oLL[1:I], meinForest)[.] %>%
          c(meinForest) -> meinForest
          # alternative to calc_chipForest :
         # chipForest <- c(chipForest, trindx)
      }
    }
  } # for exits with j the first index after all trees have been collected
  return(list(calcLogloss( subforest(forest, meinForest ), data.test)
              , I
              , length(meinForest)))
}


#paste(rep(c('LL.test.chip.','I.'),4),rep(c('d0','d1','d2','sb'),each=2),sep='')

# to base the result on more bootstraps
folder <- 'data/nursery'
files <- list.files(folder)#[1:5]
# dir(folder)

# give parameter for each dissimilarity
parameter <- list('d0'=list('cutoff'=0.2, 'sizeSF'=500)
, 'd1'=list('cutoff'=0.3, 'sizeSF'=500)
,'d2'=list('cutoff'=0.4, 'sizeSF'=500)
,'sb'=list('cutoff'=0.5, 'sizeSF'=500)
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

names(et02)

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

