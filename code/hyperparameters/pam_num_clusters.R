# cleveland02-v2/code/hyperparameters/pam_num_clusters.R
# Best number of clusters for dissimilarity matrices (based on each metric)
# needs to load dissimilarity matrices (does not create / calculate them)
# can load dissimilarity matrices based on Cleveland data or simulations thereof

rm(list=ls())

library(dplyr)
library(caret)
library(ranger)
library(effects)

withSim <- T
if(withSim){
  #load('data/nursery/nursery01_01_50x500.rda') # loads doc[[x]]$DM for x = 1:50
  load('data/nursery02/nursery02_01.rda')
  DM <-  doc[[1]]$DM
  }else{
    load('data/dms_for_a_forest_500_trees.rda') # loads DM
  }

metrices <-  c('d0','d1','d2','sb')

maxNumClus <-  50

for(m in 1:length(metrices)){
  dm <- DM[[m]]
  doc <- rep(NA,maxNumClus-1)
  for(k in 2:maxNumClus){
    #print(k)
    pam.obj <-  cluster::pam(x=dm
                             , k=k
                             , diss=T
                             , nstart = 5
                             , keep.diss = F # faster
                             , keep.data = F #faster
                             )
    doc[k-1] <- pam.obj$silinfo$avg.width
  }
  par(mar=c(2,2,1,1)+0.2)
  plot(x=1+(1:length(doc))
       , y=doc
       , type='b'
       #, main=paste('clustering based on metric' , metrices[[m]])
       #, xlab='number of clusters'
       , xlab=''
       #, ylab='average silhouette width'
       , ylab=''
       ,cex.axis=1.2
       )
  legend('topright'
         , legend=metrices[[m]]
         , cex=1.2
         #, bty = 'n'
         )
  (which.max(doc) +1 ) %>% print
}

#### Hartigans' dip test on dissimilarity matrices

for(m in 1:length(metrices)){
  dm <- DM[[m]]
  doc.dip <-  rep(NA,nrow(dm))
  for(i in 1:nrow(dm)){
    doc.dip[i] <- dip.test(dm[i,-i])$p.value
  }
  print(paste('p values for trees under', metrices[m], 'dissimilarity, directly on dissimilarity matrix'))
  c( min(doc.dip), quantile(doc.dip, 0.01) ) %>% print
}

#### Hartigans' dip test on spacial data after mds

k <- 220 # dimensions in mds

for(m in 1:length(metrices)){
  dm <- DM[[m]]
  mp <- cmdscale(dm,k)
  doc.dip <-  rep(NA,nrow(mp))
  dm <- as.matrix(dist(mp))
  for(i in 1:nrow(dm)){
    doc.dip[i] <- dip.test(dm[i,-i])$p.value
  }
  print(paste('p values for trees under', metrices[m], 'dissimilarity, mds with',k,'dimensions'))
  c( min(doc.dip), quantile(doc.dip, 0.01) ) %>% print
}

# results depend on k in mds
