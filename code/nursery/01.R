# nursery : baumschule : build a collection of forests
# 4.2.2022

# 1) simulation : bootstrapped datasets
# 2) grow a ranger forest on it
# 3) calculate all dissimilarity matrices
# 4) save resample, bootstrap training data, forest, dissimilarity matrices

rm(list=ls())

library(caret)
library(ranger)
library(dplyr)

metrices <- c('d0','d1','d2','sb')
nT <- 500 # to set num.trees=nT in ranger

load('data/data_SupMat4.rda')
# original Cleve is needed in createDMd1

source('code/source/sim-prep-v2.R') # loads functions create_prob , new_bootstrap
Cleve_enriched <- create_prob(Cleve)

source('code/source/distance-matrices.R')

# should be 1000
nBs <- 50

doc <-  list()

seed <- 20
set.seed(seed)
cr<-createResample(Cleve_enriched$CAD
                   , times = nBs
                   , list = T)

ct <- 1
for(i in 1:nBs){
  #print(paste(i/nBs , Sys.time()))
  if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
  # data frame, bootstrapped
  data.train <- new_bootstrap(Cleve_enriched , cr[[i]])[,-12] # no probabilities
  
  # grow small forest (size nT)
  rg <- ranger(CAD~.
               , data = data.train 
               , num.trees = nT
               , replace = F # neu nach Absprache mit AZ
               , mtry= 3 # default : 3
               , importance = 'impurity'
               , probability = T # this makes it a random forest of type 'Probability Estimation'
               , min.node.size = 13 # optimized in A.Z. paper: BiomJ 56, 2014, 
  )
  
  DM <- list()

  DM$d0 <- createDMd0(forest=forest)
  DM$d1 <- createDMd1(forest=forest, dft=Cleve[cr[[i]],]) # working with original Cleveland data , no prob, no CAD needed
  DM$d2 <- createDMd2(forest=forest, dft=data.train)
  DM$sb <- createDMsb(forest=forest)
  
  # old code, changed 22.03.2022
  #for(metric in metrices){
  #    # print(paste('    ' , metric , Sys.time()))
  #  DM[[metric]] <- createDM(forest=rg$forest , type=metric, dft=data.train)
  #}
  
  doc[[i]] <-  list('resample'=cr[[i]] 
                    , 'bootstrapped training data'=data.train 
                    , 'rg'=rg 
                    , 'DM'=DM)
}
View(doc)    

info <-  list('version'='2, training forests on bootstraps of Cleveland data set.'
                  , 'seed'=seed
                  , 'nBs'=nBs 
                  , 'metrices'=metrices 
                  , 'date'=Sys.time() 
                  , 'created with script file'='code/nursery/01.R')
save(info, doc, file='data/nursery/nursery01_16_50x500*.rda')

# how to load data

contentRData <- function(file) {
  #' Function for getting the names of objects in an .rda file created by R's save() command
  #' allows to check if loading will overwrite objects in the global environment
  #' Inputs: RData file 
  E <- new.env()
  load(file=file, envir=E)
  return(ls(envir=E))
}

extractorRData <- function(file, object) {
  #' Function for extracting an object from a .RData file created by R's save() command
  #' Inputs: RData file, object name
  E <- new.env()
  load(file=file, envir=E)
  return(get(object, envir=E, inherits=F))
}

filename <- 'data/nursery/nursery01_01.rda'
contentRData(filename)
contentRData(filename) %in% ls()
extractorRData(filename, 'info')    

  

