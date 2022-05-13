# explain dissimilarities in the introduction chapter of thesis
# 21.3.2022

rm(list=ls())

source('code/source/plotTree.R')
source('code/source/distance-matrices.R')

?plotTreeFb

load('data/nursery/nursery01_01_50x500.rda')

rg <- doc[[1]]$rg

treeInfo(rg,1) %>% xtable

fbt1 <- plotTreeFb(rg,1)
fbt2 <-  plotTreeFb(rg,2)
fbt3 <- plotTreeFb(rg,3)

plotTree1(rg,1)
plotTree1(rg,2)
plotTree1(rg,3)

# how often are entries equal?
t1 <- fbt1$sb.info
t2 <- fbt2$sb.info
(t1 == t2 ) %>% table
# if they are equal, what value do they have? # 0!
(t1 == t2 ) %>% t1[.]
(t1==t2 ) %>% 
  (function(x) t1[x] ==0) %>%
  table
#if they are not equal, what value do they have?
(t1 != t2 ) %>% t1[.]
(t1 != t2 ) %>% t2[.]
(t1 != t2 ) %>% 
  (function(x) t1[x]) %>%
  table
(t1 != t2) %>% 
  (function(x) t2[x]) %>%
  table

# trees are not equal and one of them has a split node (one of them is 0)
(t1 != t2 ) %>% 
  (function(x) t2[x]*t1[x] == 0) %>%
  table

# trees 1 and 3
# how often are entries equal?
t3 <- fbt3$sb.info
(t1 == t3 ) %>% table
# if they are equal, what value do they have? # 0!
(t1 == t3 ) %>% t1[.]
(t1==t3 ) %>% 
  (function(x) t1[x] ==0) %>%
  table

# trees are not equal and one of them has a split node (one of them is 0)
(t1 != t3 ) %>% 
  (function(x) t3[x]*t1[x] == 0) %>%
  table

DM <-  doc[[1]]$DM
DM$d0[1:3,1:3]
DM$d1[1:3,1:3]
DM$d2[1:3,1:3]
DM$sb[1:3,1:3]

# look for indistinguishable trees
which(DM$d0[upper.tri(DM$d0, diag = F )]==0) %>% length/(500*499/2)
which(DM$d1[upper.tri(DM$d1, diag = F )]==0)
which(DM$d2[upper.tri(DM$d2, diag = F )]==0)
which(DM$sb[upper.tri(DM$sb, diag = F )]==0)

# indistinguishable trees
which(DM$d0[3,-3]==0)  %>% length #/499

lapply(1:500, function(i) which(DM$d0[i,-i]==0) %>% length) %>% unlist %>% table %>% plot(main='frequencies of number of indistinguishable trees')

# trees with maximum dissimilarity
which(DM$d0[2,]==1)  %>% length #/499
max(DM$d0[3,])

lapply(1:500, function(i) min(DM$d0[i,-i])) %>% unlist %>% plot
lapply(1:500, function(i) min(DM$d0[i,-i])) %>% unlist %>% plot
lapply(1:500, function(i) min(DM$d0[i,-i])>0) %>% unlist %>% table

DM$d0[upper.tri(DM$d0, diag=F)] %>% table %>% unlist %>% plot

#### d0 dissimilarity 
DM$d0[upper.tri(DM$d0, diag=F)] %>% summary %>% xtable %>% print
DM$d0[upper.tri(DM$d0, diag=F)] %>% (function(x) c(min(x),mean(x),max(x),sd(x)))

DM$d0[upper.tri(DM$d0, diag=F)] %>% table
DM$d0[upper.tri(DM$d0, diag=F)] %>% hist

#### d1 dissimilarity ####
DM$d1[1:3,1:3]
E3[upper.tri(E3, diag=F)] %>% (function(x) c(min(x),mean(x),max(x),sd(x)))
DM$d1[upper.tri(DM$d1, diag=F)] %>% summary %>% xtable %>% print
DM$d1[upper.tri(DM$d1, diag=F)] %>% (function(x) c(min(x),mean(x),max(x),sd(x)))
DM$d1[upper.tri(DM$d1, diag=F)] %>% hist
# single pair of most similar trees
(DM$d1[upper.tri(DM$d1, diag=F)]==min(DM$d1[upper.tri(DM$d1, diag=F)])) %>% table 
which.min(DM$d1[upper.tri(DM$d1, diag=F)])
# closest 2 trees , most similar trees
d1min <- min(DM$d1[upper.tri(DM$d1, diag=F)]) 
lapply(1:500, function(i) which(DM$d1[i,-i]<0.0457))  %>% unlist %>% table
DM$d1[35,102]
plotTree1(rg,35)
plotTree1(rg,102)


data.train <- doc[[1]][["bootstapped training data"]]
predict(rg, Cleve[unique(doc[[1]]$resample),] , type = 'terminalNodes') -> predd1
View(predd1$predictions[,c(35,102)])

#### d2 dissimilarity 
DM$d2[upper.tri(DM$d2, diag=F)] %>% summary %>% xtable %>% print
DM$d2[upper.tri(DM$d2, diag=F)] %>% (function(x) c(min(x),mean(x),max(x),sd(x)))

#### shannon banks
# whenever full binary trees agree, it is on a non-existent node. 
# They never agreee having a split node at a position and having the same split variable
fbt1$sb.info[which(fbt1$sb.info == fbt2$sb.info)] %>% table
fbt1$sb.info[which(fbt1$sb.info == fbt3$sb.info)] %>% table
fbt2$sb.info[which(fbt2$sb.info == fbt3$sb.info)] %>% table


DM$sb[upper.tri(DM$sb, diag=F)] %>% summary %>% xtable %>% print
DM$sb[upper.tri(DM$sb, diag=F)] %>% (function(x) c(min(x),mean(x),max(x),sd(x)))

DM$d0[upper.tri(DM$d0, diag=F)] %>% table
DM$d0[upper.tri(DM$d0, diag=F)] %>% hist

#### trees 1,2,3
predd1$predictions[,c(1,2,3)]

DM$d2[upper.tri(DM$d2, diag=F)] %>% summary
DM$d2[upper.tri(DM$d2, diag=F)] %>% hist

DM$sb[upper.tri(DM$sb, diag=F)] %>% summary
DM$sb[upper.tri(DM$sb, diag=F)] %>% hist
DM$sb[upper.tri(DM$sb, diag=F)] %>% table
