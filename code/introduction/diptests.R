# visualise trees of a forest with a mds (multidimensional scaling) or similar (t-SNE)

library(MASS)
library(ranger)
library(dplyr)
#install.packages('tsne')
library(tsne)

# load manually : mds from hyperparameters/dip_unimodality.R
# load manually : calcLogloss2, winsorize_probs from code/source/prep.R
source('code/source/prep.R')
# load manually : calc_LL from code/source/chipman.R
source('code/source/chipman.R')

# open first forest
load('data/nursery02/nursery02_01.rda')
rg <- doc[[1]]$rg
DM <- doc[[1]]$DM

LL <-  calc_LL(rg$forest, doc[[1]]$`bootstrapped training data`)

for(nm in names(DM)){
  idx.min <- which.min(apply(DM[[nm]],2,sum))
  col <- rep(1,500)
  col[idx.min] <- 2
  pch <- rep(1,500)
  pch[idx.min] <-  2
  mds(DM[[nm]], col=col, pch=pch , main=nm)
  col <- rep(0,500)
  print(paste(nm, idx.min))
  #plot(DM[[nm]][idx.min,], LL) # doesn't look good, not generating insights
  }

m <- cmdscale(DM$d0, eig = TRUE, k = 2)
print(names(m))
m$points[,1:2]

# sammon # did that before, right?!
metric <- 'd1'
sammon(DM[[metric]],k=2) # not working for d0 with dissimilarities 0
# stress ist too high for the plot to show representation and / or diversity
sammon(DM[[metric]],k=2)$stress

# goodness of fit
# distance to tree 1
metric
m <- cmdscale(DM[[metric]], eig = TRUE, k = 2)
new_dist <-  (m$points - m$points[1,])^2 %>% apply(1,sum) %>% sqrt
assertthat::assert_that(length(new_dist)==500)
plot(DM[[metric]][1,], new_dist)

abline(0,1) # is slope 1 the ideal?
abline(lm(new_dist~DM[[metric]][1,])$coeff, col='grey')
(abs(DM$d0[1,] - new_dist) < 0.05) %>% table %>% prop.table

# stress
# how dist works 
dist(rbind(c(1,1),c(2,2)))
dist(rbind(c(1,1),c(2,2))) %>% as.matrix
# this needs more memory, working with a matrix
# we go the other way around and only use the relevant part of the dissimilarity matrix DM$d0 (or DM$d1,...), as a dist object

# stress 1
(as.dist(DM[[metric]]) - dist(cmdscale(DM[[metric]], eig = TRUE, k = 2)$points) )^2 %>% sum %>% (function(x) x/sum(as.dist(DM$d0)^2)) %>% sqrt
# obviously, stress of 0 is perfect, stress of 0.05 would be good, 0.2 would be poor
# compare to random perturbation
(as.dist(DM[[metric]]) - sample(as.dist(DM[[metric]])) )^2 %>% sum %>% (function(x) x/sum(as.dist(DM[[metric]])^2)) %>% sqrt
# random perturbation of dissimilarity results in more stress than the mds
# that is good

# required to get different values within the 0.1,0.2,0.3... range in DM$d0
no_duplicates <- function(x) {
  eps <- .Machine$double.eps
  # apply to each matrix entry individually
  apply(x,c(1,2),function(x1) x1 + abs(rnorm(1,0,0.001)))
}

dont_stress_me <- function(metric){
  stress1 <- (as.dist(DM[[metric]]) - dist(cmdscale(DM[[metric]], eig = TRUE, k = 2)$points) )^2 %>% sum %>% (function(x) x/sum(as.dist(DM[[metric]])^2)) %>% sqrt
  #stress2 <- (DM[[metric]]- as.matrix(dist(cmdscale(DM[[metric]], eig = TRUE, k = 2)$points)))^2 %>% sum %>% (function(x) x/sum(as.dist(DM[[metric]])^2)) %>% sqrt
  #assertthat::assert_that(stress1==stress2) %>% print # not equal, makes a difference if you work with the matrix or the dist object
  
  random_stress <- (as.dist(DM[[metric]]) - sample(as.dist(DM[[metric]])) )^2 %>% sum %>% (function(x) x/sum(as.dist(DM[[metric]])^2)) %>% sqrt
  
  # to do sammon dissimilarities may not be 0 and there may not be any dublicates (like two pairs with a dissimilarity of 0.1)
  if(metric=='d0'){
    d <- no_duplicates(DM[[metric]])
    diag(d) <- 0
  }else{
      d <- DM[[metric]]
    }
  sammon1 <- sammon(d=d, k=2)
  
  # explicitly working with a dist object => distance or dissimilarity, not the raw data
  t1 <-  tsne::tsne(as.dist(DM[[metric]]), k=2)
  stress3 <- (as.dist(DM[[metric]]) - dist(t1))^2 %>% sum %>% (function(x) x/sum(as.dist(DM[[metric]])^2)) %>% sqrt

  return(list(mds.GOF=cmdscale(DM[[metric]], eig = TRUE, k = 2)$GOF %>% unlist
              , mds_stress1=stress1
              , random_stress1=random_stress 
              , sammon_stress = sammon1$stress
              , tsne_stress = stress3)
         )
}

dont_stress_me('d0')

doc1 <- matrix(0, nrow=length(names(DM)), ncol=6)
ct <- 1
for(metric in names(DM)){
  doc1[ct,] <- dont_stress_me(metric) %>% unlist
  ct <- ct+1
}

df1 <- data.frame(doc1)
colnames(df1) <- c('mds GOF 1','mds GOF 2', 'mds stress','stress of perturbation','sammon stress','tsne stress')
rownames(df1) <- names(DM)
df1

# central trees : d0 9 , d1 179 , d2 159 , sb 21
central_trees <-  c(9,179,159,21)
for(i in 1:length(names(DM))){
  mds1 <- cmdscale(d=DM[[i]], eig=T, k=2) # with eig=T only the points are returned , no GOF returned
  plot(mds1$points , main=paste('mds', names(DM)[i]))
  xlim <- par('usr')[c(1,2)]
  ylim <-  par('usr')[c(3,4)]
  plot(mds1$points[c(1:3,central_trees[i]),]
       , col=c(1,1,1,2)
       , main=paste('mds', names(DM)[i])
       , xlim=xlim , ylim=ylim)
  text(mds1$points[c(1:3,central_trees[i]),] + 0.05*(xlim[2]-xlim[1])
       , labels = c('1','2','3','central') )
  
  plot(as.dist(DM[[i]]), dist(mds1$points)
       , main=paste('quality of mds, ', names(DM)[i])
       , xlab='dissimilarity'
       , ylab='plotted distance')
  
  plot(x= DM[[i]][central_trees[i],]
       , y=as.matrix(dist(mds1$points- mds1$points[central_trees[i],]))[central_trees[i],]
       , main=paste('quality of mds,', names(DM)[i])
       , xlab='dissimilarity to central tree'
       , ylab='plotted distance to central tree')
  
  plot(x=as.dist(DM[[i]][c(1,2,3,central_trees[i]), c(1,2,3,central_trees[i])])
       , y=dist(mds1$points[c(1,2,3,central_trees[i]),])
       , main=paste('quality of mds for selected trees', names(DM)[i])
       , xlab='dissimilarity'
       , ylab='plotted distance')
  
  ######
  t1 <- tsne(as.dist(DM[[i]]))
  plot(t1, main=paste('t-sne', names(DM)[i]))
  xlim <- par('usr')[c(1,2)]
  ylim <-  par('usr')[c(3,4)]
  plot(t1[c(1:3,central_trees[i]),]
       , col=c(1,1,1,2)
       , main=paste('t-sne', names(DM)[i])
       , xlim=xlim , ylim=ylim)
  text(t1[c(1:3,central_trees[i]),] + 2
       , labels = c('1','2','3','central') )
  
  plot(as.dist(DM[[i]]), dist(t1)
       , main=paste('quality of t-sne, ', names(DM)[i])
       , xlab='dissimilarity'
       , ylab='plotted distance')
  
  plot(as.dist(DM[[i]][c(1,2,3,central_trees[i]), c(1,2,3,central_trees[i])])
       , dist(t1[c(1,2,3,central_trees[i]),])
       , main=paste('quality of t-sne for selected trees' , names(DM)[i])
       , xlab='dissimilarity'
       , ylab='plotted distance')
  
  seriation::dissplot(DM[[i]][c(1,2,3,central_trees[i]),c(1,2,3,central_trees[i])])
}

install.packages('seriation')
library(seriation)

seriation::dissplot(DM$d0)

