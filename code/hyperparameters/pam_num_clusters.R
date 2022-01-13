# Best number of clusters for CLeve and each metric

rm(list=ls())

library(dplyr)
library(caret)
library(ranger)
library(effects)

load('data/data_SupMat4.rda')
df <- Cleve
df$CAD_fac <-  NULL

rg <-  ranger(CAD~.
              , data=df
              , num.trees = 2000
              , min.node.size = 13
              , replace = F
              , probability = T
              )

source('code/source/distance-matrices.R')

#### clustering a forest is very different for the different dissimilarity measures
#### the different dissimilarity measures create / see very different structures

#df2 <- df[sample(1:nrow(df),10,F),] # for debugging, a smaller set , will run faster
DM <- list()

metrices <-  c('d0','d1','d2','sb')
for(m in metrices){
  print(paste(m , Sys.time()))
  #dm <- createDM(rg$forest, type=m, dft=df2) # faster ...
  dm <- createDM(rg$forest, type=m, dft=df)
  doc <- rep(NA,19)
  for(k in 2:20){
    #print(k)
    pam.obj <-  cluster::pam(x=dm
                             , k=k
                             , diss=T
                             , nstart = 5
                             , stand = T
                             , keep.diss = F
                             , keep.data = F)
    doc[k-1] <- pam.obj$silinfo$avg.width
  }
  plot(x=1+(1:length(doc))
       , y=doc
       , type='b'
       , main=paste('clustering based on metric' , m)
       , xlab='number of clusters'
       , ylab='average silhouette width'
       )
  (which.max(doc) +1 ) %>% print
  DM[[metric]] <- dm
}

# d0 : the more clusters the better, big jump / gain at 4 clusters
# d1 : 2,3 or 4 (argmax) clusters give a high (pos) average silhouette width
# d2 : always negative  avg. silhouette, maybe no cluster is best?
# sb : 2 clusters is best (pos avg. silhouette) , maybe no clusters is even better?

# save the generated dissimilarity matrices
#info <- paste('ranger forest with num.trees=',rg$num.trees,'2000, dissimilarity matrices for d0,d1,d,sb ')
#save(info, rg, DM,file=paste('data/dms_for_a_forest_',rg$num.trees,'_trees*.rda', sep=''))



