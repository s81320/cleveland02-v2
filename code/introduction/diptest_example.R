# checking out the Hartigans' dip test
# 5.5.2022
# example for the diptest applied to higher dimensional data
# in thesis : figure fig:intro:mds:example
# code / introduction / diptest example.R

# better: do mds with as many dimensions as needed to get stress below 0.1 or 0.05 or minimize stress.
# then test for multimodality in R^k

par(mar=c(4,4,3,1)+0.1)

#### 1 dimensional data ####

# construct as multimodal
set.seed(1)
data <- c(rnorm(6,-4,0.5) , rnorm(6,0,0.5) , rnorm(6,4,0.5)) # 3 modes at -4, 0,4
plot(density(data))
hist(data, breaks=20)
dip.test(data) # not recognised as multimodal, 0.14 p-value

#### data created as multimodal
#### unimodlity rejected for many observations per cluster, not rejected for few observations per cluster

set.seed(1)
n <- 20 # set to 20 or 2
data <- c(rnorm(n,-4,0.5) , rnorm(n,0,0.5) , rnorm(n,4,0.5))
plot(density(data)) # see modes
hist(data, breaks=40) # cannot see modes with default number of bars / breaks. Can see them with 20 breaks
plot(ecdf(data))
dip.test(data) # for n=50 reject unimodality on level 1% , for n=2 pvalue at 0.2599 cannot reject unimodality

#### other option : create data as multimodal. 
#### with centers close or far , fixed number of elements per cluster

#### high dimensional data : 2 dim as example

# install.packages('MASS')
library(MASS)
mu1 <- c(-4,0)
mu2 <-  c(0,4)
#sigma0 <- matrix(c(1,0.5,0.5,1),2,2)
sigma0 <- matrix(c(1,0,0,1),2,2)
set.seed(2)
data <- rbind(mvrnorm(10,mu1,sigma0),
              mvrnorm(10,mu2,sigma0))
par(mar=c(2,2,1,1)+0.1)
plot(data, type='n' , xlab='' , ylab='')
text(data, labels = 1:20, col=c(rep('green',10),rep('blue',10)))
points(rbind(mu1,mu2),col='red',pch=0)
legend('bottomright'
       , legend=c('mu1 , mu2','centered at mu1','centered at mu2') 
       , col=c('red','green','blue')
       , pch=c(0,1,1)
       , cex=0.8)

par(mar=c(2,2,1,1)+0.1)
pch <- rep(1,nrow(data))
pch[10] <- 8
pch[6] <- 8
plot(data , xlab='' , ylab='' 
     , col=c(rep('green',10),rep('blue',10))
     , pch=pch)
points(rbind(mu1,mu2),col='red',pch=0)
legend('bottomright'
       , legend=c('mu1 , mu2','centered at mu1','centered at mu2') 
       , col=c('red','green','blue')
       , pch=c(0,13,1)
       , cex=0.8)

# applying the dip test when working with a dissimilarity matrix

# creating the dissimilarity matrix based on points (R^k)
dm <- data %>% dist %>% as.matrix

doc.dip <-  rep(NA,nrow(dm))
for(i in 1:nrow(dm)){
  doc.dip[i] <- dip.test(dm[i,-i])$p.value 
} 
doc.dip
min(doc.dip) ; which.min(doc.dip) ; quantile(doc.dip,0.01)
plot(doc.dip[order(doc.dip)], type='l') 


hist(dm[6,-6],breaks=40, main='')
plot(density(dm[6,-6]),main='')
dip.test(dm[6,-6])

hist(dm[10,-10],breaks=40, main='')
plot(density(dm[10,-10]),main='')
dip.test(dm[10,-10])

# discreet variables are always multimodal !!
# we can get a density - but a density is not for discreet variables. It does not make sense!!

summary(as.dist(dm))
dip.test(c(rep(1,20),rep(2,20))) # multimodal
rep(c(1,2),times=10)
rep(c(1,2),each=10)
dip.test(rep(c(1,2),each=10)) # multimodal
dip.test(rep(c(1,2),each=10)/10000) # multimodal
# not exactly always. It depends on the number of occurences relative to the distance between points...
hist(rep(1:10,each=10)/100,breaks=20)
dip.test(rep(1:10,each=10)/100) # unimodal , pvalue 0.06

# same distance but more observations
hist(rep(1:10,each=100)/100, breaks=20)
dip.test(rep(1:10,each=100)/100) # multimodal

# should we standardize the dissimilarity matrix before clustering? 
# It seems to matter for the silhouette width

set.seed(1)
dm <-  matrix(abs(rnorm(100,0,1)), nrow=10) # not a dissimilarity matrix
dm1 <-  dm %>% as.dist %>% as.matrix # now symmetric and 0 on diagonal
# transformations of the dissimilarity matrix
# dm2 <-  dm*5 # does not change sil width (common factor, cancels in definition of silhouette)
# dm2 <- dm /(max(dm)-min(dm)) # 1/(max-min) is also a common factor
dm2 <-  10 + dm # not a dissimilarity matrix, 10 on the diagonal
dm2 <-  dm2 %>% as.dist %>% as.matrix # now symmetric and 0 on diagonal
# but dm2 does not represent the same dissimilarities, it is profoundly changed. Not allowed!
pam(x=dm1, k=2)$silinfo$avg.width
pam(x=dm2, k=2)$silinfo$avg.width

# what matters is the distribution of of distances
# min(dm) > 0 makes it harder to get a good clustering
# max(dm) close to min(dm) also makes it harder (especially / only if min(dm)>0)

load('data/nursery02/nursery02_01.rda')

DM <- doc[[1]]$DM

# ranges of dissimilarities, ratio of max to min in dissimilarities
for(i in 1:4){
  names(DM)[i] %>% print
  dm <-  DM[[i]]
  r <- range(as.dist(dm)) 
  r %>% print
  (r[2]/r[1]) %>% print
}

