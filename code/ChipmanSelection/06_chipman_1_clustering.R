# trying to find the best clustering in Chipman 1
# and trying to find a test for uni / multimodality

rm(list=ls())

library(dplyr)
#library(caret)
library(ranger)
#library(effects)
library(cluster)
library(diptest)

# load doc from data / nursery / ..

# load Cleveland and test data
load('data/data_SupMat4.rda') # loads the data sets Cleve (for OOB), Hung, Swiss, VA

source('code/source/prep.R') # calcLogloss (on forests) , calcLogloss2 (on predicted probabilities)

# get dissimilarity matrix
i <- 1
metric <-  'd1'
dm <- doc[[i]][["DM"]][[metric]]

data.train <- doc[[i]]$`bootstapped training data`

OOB <-  base::setdiff(1:nrow(Cleve), unique(doc[[i]]$resample))
data.set.val <- Cleve[OOB,] # goes back to original Cleveland data

forest<- doc[[i]]$rg$forest

pp <- predict(forest 
              , data=data.set.val 
              , predict.all = T)$predictions[,2,]
#pp <- simplify2array(pp, higher=FALSE)

lapply(1:forest$num.trees
       , function(k){ 
         pp[,k] %>% 
           calcLogloss2( df=data.set.val ) %>% 
           unlist
       }) %>% 
  unlist -> LL

oLL <- order(LL) # tree indices , ordered by logloss on OOB

I <- 100
dm2 <-  dm[oLL,oLL][1:I,1:I] #; dm2

#m <- c( "average", "single", "complete", "ward")

#function to compute coefficient
#ac <- function(x) {
#  agnes(dm2[1:I,1:I], method = x, diss=T)$ac
#}
#hpo <-  lapply(m, ac)
#hpo

hc <- cluster::agnes(dm2, 
                     #method='complete',
                     #method='average', 
                     method="ward", # always optimal - whenever I did hpo
                     #method = m[which.max(hpo)] , 
                     diss=T)
#pltree(hc, cex = 0.6, hang = -1, main = paste('ward, ', metric))
#hc %>% summary
#hc.clus <- cutree(hc, k = sizeSF)

col1 <- list()
lapply(2:(I-1), 
       function(k){
         cutree(hc, k = k) %>% 
         silhouette(dmatrix=dm2) %>%
         .[,'sil_width'] %>% mean
         }
) -> col1
unlist(col1) %>% which.max +1

unlist(col1) %>% plot(x=2:(I-1)
                      , main=paste(metric, ', mean silhouette width'))


metric <-  'd2'
dm <- doc[[i]][["DM"]][[metric]]
dm2 <- dm[oLL,oLL][1:I,1:I]
dip.test(dm2)

m <- cmdscale(dm2[1:10,1:10], eig = TRUE, k = 2)
x <- m$points[, 1]
y <- m$points[, 2]
  
col <- heat.colors(10)[1:10]

plot(x,y
       , main=paste('Trees, dissimilarity in the', metric ,'metric')
       , col=col)
text(x,y,labels = 1:10)
#legend('bottomright'
#         , legend=c('small forest','default forest')
#         , col=c(1,2)
#         , pch=1
#         , cex=0.8)

MASS::sammon(dm2[1:10,1:10],k=2) -> sammon1 
sammon1$points %>% plot(type='n')
sammon1$points %>% text(labels=1:10)
sammon1$stress 
# how to interprete the Sammon stress ??
# Wikipedia https://de.wikipedia.org/wiki/Multidimensionale_Skalierung#STRESS-Maße
# but which stress measure is used im Sammon ? https://en.wikipedia.org/wiki/Sammon_mapping

# correlation of LL / ordered LL and dissimilarity to best tree
# I thought both should give the same result.
# not working :-( But both are small :-) and close to 0
cor(oLL[2:100], dm2[1,-1])
cor(LL[-oLL[1]][1:99], dm[oLL[1],-oLL[1]][1:99])

