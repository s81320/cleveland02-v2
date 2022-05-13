# chipman.R test

library(dplyr)
library(cluster)

load('data/data_SupMat4.rda')
data.test <- Swiss 

source('code/source/chipman.R')
source('code/source/prep.R') # needed in chipman.R ... maybe I should load it there?
source('code/source/subforest.R')

# load arguments (dm , forest, oLL, parameter) to feed the function
{
  #load('data/nursery/nursery01_10_50x500.rda')
  load('data/nursery02/nursery02_01.rda')
  
{
  metric <- 'd1'
  DM <- doc[[1]]$DM
  dm <-  doc[[1]]$DM[[metric]]
  attr(dm,'metric') <- metric
}

  #data.train <- doc[[1]]$`bootstapped training data` # typo, missing r in bootst apped
  data.train <- doc[[1]]$`bootstrapped training data` # typo corrected in nursery02
  OOB <-  base::setdiff(1:nrow(Cleve), unique(doc[[1]]$resample))
  data.set.val <- Cleve[OOB,] # goes back to original Cleveland data

  forest<-doc[[1]]$rg$forest

  pp <- predict(forest 
              , data=data.set.val 
              , predict.all = T)$predictions[,2,]

  lapply(1:forest$num.trees
       , function(k){ 
         pp[,k] %>% 
           calcLogloss2( df=data.set.val ) %>% 
           unlist
       }) %>% 
  unlist -> LL

  oLL <- order(LL) # tree indices , ordered by logloss on OOB

  parameter <-  list('cutoff'=0.2 , selection='best')
  pa <- parameter
}

# grow Chipman forest in original code in chipman.R till you get R
# stop before clustering!
R %>% length
# move to multidimensional scaling for the whole forest
m <-  cmdscale(dm,k=50, eig=T) # scaling for all trees
# dip test for representing subforest R
mp <-  m$points[R,]

# function my_diptest is in diptest_mds.R
my_diptest(mp)
my_diptest(m$points)
metric

# within the grow_chipForest_1 function we are 
# NOT APPLYING THE DIP TEST CORRECTLY
# the diptest on d0 and sb dissim will always find multimodality (for larger sizes)
# as d0 and d_sb are discreet
# should be first move to R^k and then test for multimodality? That would be similar to the Chipman paper

gch1 <- grow_chipForest_1(dm=dm , forest= rg$forest, oLL=oLL, parameter=pa, output=T)
gch1$metric <-  metric

gch1$sil %>% plot

calc_chipForest_1(dm=dm , forest= forest, oLL=oLL, parameter=pa)

lapply(1:length(R) , function(i) dip.test(dm2[i,-i], simulate=T)$p.value) %>% 
  unlist%>% 
  min %>% 
  (function(x) x<0.05) -> multimod # logical , TRUE if p.value of at least one dip test less than 0.05

