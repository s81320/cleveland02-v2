rm(list=ls())

library(dplyr)
library(caret)
library(ranger)
library(effects)
library(cluster)

source('code/source/sim-prep.R')

Cleve <-  Cleve[,1:11]
Hung <-  Hung[,1:11]
VA <-  VA[,1:11]
Swiss <- Swiss[,1:11]

rg <- ranger(CAD~.
             , data = Cleve # exclude the probabilities
             , num.trees = 500
             , replace = F # neu nach Absprache mit AZ
             , mtry= 3 # default : 3
             , importance = 'impurity'
            # , probability = T # this makes it a random forest of type 'Probability Estimation'
             , min.node.size = 13 # optimized in A.Z. paper: BiomJ 56, 2014, 
)

# accuracy for classification 
# needs probability=F , or remove probability=T
# forest trained on Cleve, tested on all available sets
for(nameDataSet in c('Cleve','VA','Swiss','Hung')){
  teDa <-  get(nameDataSet)
  print(paste('accuracy:', nameDataSet))
  (predict(rg, teDa)$predictions == teDa$CAD) %>%
    which %>% 
    length %>%
    (function(x){x/nrow(teDa)}) %>% 
    print
}





