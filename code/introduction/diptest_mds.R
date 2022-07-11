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

{
  metric <- 'sb'
  dm <- DM[[metric]]
  attr(dm,'metric') <- metric
}

dim(dm)

m <- cmdscale(d=dm,k=250,eig=T) # make k the maximum of dimensions you might consider
# cmdscale will always return all (=500, = nrow(dm)) eigenvalues
length(m$eig)
# goodness of fit is for the given k only
# its 2 components are for 2 different ways of calculating goodness of fit, see documentation
m$GOF 
# from the eigenvalues we plot the curve for goodness of fit over all eigenvalues , all possible dimensions k
m$eig %>% lapply(function(x) max(x,0)) %>% unlist %>% sum -> sum.pos.eig
m$eig %>% 
  lapply(function(x) max(x,0)) %>% 
  unlist %>% 
  cumsum %>% 
  (function(x) x/sum.pos.eig) %>% 
  plot(type='l', xlab='k , dimensions in mds', ylab='scaled sum of first k eigenvalues')
legend('topright', legend=metric, cex=0.8)
m50 <- cmdscale(d=dm,k=50,eig=T)
points(x=50, y=m50$GOF[2], col='grey', pch=0) # example, how GOF for k=50 lies on the curve
# in this kind of plot we look for a knee

##################################

m <- cmdscale(d=dm,k=250,eig=T)
doc.R <-  rep(NA, ncol(m$points))
for(k in 2:ncol(m$points)){
  if(k %% 50 == 0) print(k)
  mpk <-  m$points[,1:k]
  #dmp <- dim(mpk) ; dmp # DMP : Dim M$Points
  SSres <- (as.dist(dm) -dist(mpk))^2 %>% sum 
  # mean value for columns is basically 0
  # this is by construction in cmdscale: points are translation invariant !!
  mean_m_points <- apply(mpk,2,mean) ; mean_m_points
  SStot <- dist(mpk 
              - matrix(rep(mean_m_points,nrow(mpk))
                       ,byrow=T
                       ,ncol=k)
              )^2 %>% sum
  # so we get the same when removing it from the calculation
  # assertthat::assert_that(SStot==sum(dist(mpk)^2), msg=paste(k, SStot,sum(dist(mpk)^2)))
  doc.R[k] <- 1 - SSres/SStot
}
x <- 10:50
plot(x, doc.R[x], type='l',xlab='k, dimensions in mds', ylab='R squared')
legend('topright', legend=metric, cex=0.8)

##################################

m <- cmdscale(d=dm,k=250,eig=T)

# how many dimensions are enough? look for an elbow!
plot(m$eig[1:500], xlab='k , dimensions in mds' , ylab='eigenvalue', type='l')
abline(h=0, col='grey')
legend('topright', legend=metric, cex=0.8)
# negative eigen values are errors 

# focus
par(mar=c(4,4,1,1)+0.2)
plot(m$eig[1:60], xlab='dimension' , ylab='eigenvalue', type='l')
abline(h=0, col='grey')
legend('topright', legend=metric)

elbow <- 20 # 10 or 6 for d0 # for sb: 20 is too small (R square at 0.86), 200 too large (just too much), we go with 50 (R square a little below that of 200)
plot(m$eig[1:(elbow+5)], xlab='dimension' , ylab='eigenvalue', type='l')
abline(h=0, col='grey')
legend('topright', legend=metric, cex=0.8)


# R squared with 
# optimal dimension now fixed at knee
mp <- m$points[,1:elbow] # mp : m$points
dmp <- dim(mp) ; dmp # DMP : Dim M$Points
SSres <- (as.dist(dm) -dist(mp))^2 %>% sum 
# mean value for columns is basically 0
# this is by construction in cmdscale: points are translation invariant !!
mean_m_points <- apply(mp,2,mean) ; mean_m_points
SStot <- dist(mp 
              - matrix(rep(mean_m_points,dmp[1])
                                ,byrow=T
                                ,ncol=dmp[2])
              )^2 %>% sum
# so we get the same when removing it from the calculation
assertthat::assert_that(SStot==sum(dist(mp)^2))
1 - SSres/SStot

# we work with this spacial representation in 6 dimensions
mp %>% dim

# apply the dip test on data in R^6

my_diptest <-  function(mp){
  doc.dip <-  rep(NA,nrow(mp))
  dm2 <- as.matrix(dist(mp))
  for(i in 1:nrow(dm2)){
    doc.dip[i] <- diptest::dip.test(x=dm2[i,-i])$p.value 
  }
  plot(doc.dip
     , type='l'
     #, main=paste('dip test after dim reduction to ' , elbow, ',',metric)
     , ylab='p.value dip test distances to tree'
     , xlab='tree'
     )
  return(list('min p.value'=min(doc.dip)))
}

my_diptest(mp)
# d0 will always reject unimodality
# sb does not reject multimodality when transforming to R^499 and reducing to R^50
# smallest pvalue is 13.5%

# diptest for representing subforest
dim(mp)
stress(dm[R,R], list(points=mp)) # poor stress # WHICH PACKAGE ??

       