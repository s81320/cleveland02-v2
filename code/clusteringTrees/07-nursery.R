# code/clusteringForests/07.R
# 09.02.2022

# train on bootstraps of Cleve, previously created in nursery
# loop over multiple nurseries using a collector

# test on Hung or VA

# cluster default forest, select trees from cluster into subforest

# using cluster methods 1,3,4
# always cluster proportionally
# 1) cluster into the desired number directly
# 2) build 2 clusters 
# 3) build 3 clusters

# measure of success : logloss

rm(list=ls())

library(dplyr)
library(caret)
library(ranger)
library(effects)
library(cluster)

source('code/source/prep.R') # calcLogloss
source('code/source/subforest.R') # subforest

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

# set test data by name : VA, Swiss or Hung
data.test.name <-  'Hung'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name


# load('data/nursery/nursery01_03_50x500.rda') # loads doc and info
# contains forests and their dissimilarity matrices
# info

t1<-function(tri){
  ti<-treeInfo(rg,tri)
  return(length(which(ti$terminal==TRUE))) # terminal nodes , -1 gives the number of splits
}


cf <- function(pam.obj, goodTrees, nClus, preselectMethod, targetSize){ # function needed for clustering strategies 2,3,4
  #' function for cluster strategies
  #' 
  #' 
  #' @param nClus From cluster nClus select a number of trees (howMany) if available
  #' choose only trees with a positive silhouette width
  #' @param preselectMethod According to preselectMethod select trees of a fixed size (howMany) or proportionally to the size of the cluster.
  #' @param targetSize the final size of the subforest aimed at, should be sizeSF[i]
  #' For a forest of 600 trees, a cluster of 60 trees (10% of the forest) and we want to build a subforest of size 50
  #' then collect 50*10% = 5 trees from this cluster
  #' If medoids are not in the selection, replace the first element by the medoid
  
  A1 <- which(pam.obj$clustering==nClus) 
  S1 <- base::intersect(A1,goodTrees) 
  
  # check the cluster the trees belong to
  #pam.obj$clustering[S1] 
  assertthat::validate_that(all(pam.obj$clustering[S1]==nClus))
  
  #check the sil width of the tree (should be positive!)
  #(pam.obj$silinfo$width[as.character(S1),'sil_width'] >0) %>% table
  
  # preselection 
  
  assertthat::assert_that(class(targetSize)=='numeric')
  assertthat::assert_that(targetSize>0)
  
  if(preselectMethod=='const'){
    howMany <- max(1,round(targetSize / pam.obj$call$k ,0))
  }else{
    if(preselectMethod=='prop'){
      howMany <- round(targetSize*pam.obj$clusinfo[nClus,1] / sum(pam.obj$clusinfo[,1]),0) 
      # pam.obj is created with forest and clusters forest$num.trees trees
      # how many <- target size (=how many you finally want) * cluster size / sum of trees in all clusters
    }
  }
    
  if(length(S1)>howMany){ # select from intersection if it has enough elements
    sample(S1,howMany) -> preselect
  }else{
    S1 -> preselect
  }
    
  # final selection , making sure the medoid is returned
  if(pam.obj$medoids[[nClus]] %in% preselect){ # has the medoid been selected?
      return(preselect) # if so, we're good
  }else{ 
    preselect[[1]] <- pam.obj$medoids[[nClus]] # else, we put the medoid as the first element in our selection
    return(preselect)
  }
}

calc_mDiss <-  function(dm , idcs){
  Dee <- dm[idcs,idcs]
  upper.tri(Dee, diag=F) %>% 
    Dee[.] %>%
    mean %>% 
    return
}

calc_LL_for_selection <- function(doc, sizeSF){

metrices<-c('d0','d1','d2','sb')
#metrices <- 'd0'

nT <- doc[[1]]$rg$num.trees
sizeDefault <- nT

nBs <- length(doc)

evalT <-matrix(NA,nrow=17*nBs,ncol=11) # table of evaluations

ct <- 1
for(i in 1:nBs){
  #print(paste(i/nBs , Sys.time()))
  if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
  # data frame , bootstrapped
  
  data.train <- doc[[i]]$`bootstapped training data`

  # grow small forest (size nT)
  rg <- doc[[i]]$rg
  forest<-rg$forest
  
  #### random subforest ####
  #### will be used later, when dissimilarity matrics are calculated
  
  #idcs.r <-  sample(1:forest$num.trees , sizeSF , replace=F)
  idcs.r <- 1:sizeSF # this will agree with the data generated in ranger_num_trees-v3-nusery
  LL.r <- calcLogloss(subforest(forest,idcs.r), data.test)
  
  #### default forest ####
  ########################
  
  idcs <- 1:sizeDefault
  LL.d <- calcLogloss(subforest(rg,idcs), data.test)

  c(i,nT
    , sizeDefault
    , 'default'
    , 'none'
    , NA
    , NA
    , NA
    , LL.d
    , NA
    , LL.d - LL.r
)  -> evalT[ct,]
  ct <- ct+1
  
  #### 1. strategy ###################################
  #### cluster, select medoids into the subforest ####
  ####################################################
  
  DM <- doc[[i]]$DM
  for(metric in metrices){
    #print(metric)
    #metric <- 'd2'
    dm <-  DM[[metric]]
    #print(dim(dm))
    
    pam.obj <- cluster::pam(dm
                          , k=sizeSF
                          , diss=TRUE
                          , medoids='random'
                          , nstart = 5
                          )
    #print(pam.obj)
  
    idcs <- pam.obj$medoids
    LL.s <- calcLogloss(subforest(forest,idcs), data.test)
  
    mDiss <- calc_mDiss(dm,idcs)

    c(i,nT, sizeSF, 'clustering1', metric
      , pam.obj$silinfo$avg.width # mSW mean silhouette width
      , pam.obj$silinfo$width[as.character(idcs),'sil_width'] %>% mean # mean silhouette width over the selected trees
      , mDiss
      , LL.s
      , LL.s - LL.d
      , LL.s - LL.r
      ) -> evalT[ct,]
    ct <- ct+1

    #### 3. strategy ##########################################################################################################
    #### cluster into 2 clusters, randomly select from clusters proportionally to their size until 50 trees are selected  ####
    ###########################################################################################################################
    
    pam.obj <- cluster::pam(dm
                            , k= 2 
                            , diss=TRUE
                            , medoids='random'
                            , nstart = 5)
    
    (pam.obj$silinfo$width[,'sil_width']>0) %>% 
      which %>%
      names %>% 
      as.integer -> goodTrees
    
    which(pam.obj$clusinfo[,'size'] >1) -> clustersSelected 
    
    cf3 <- function(nClus) cf(pam.obj,goodTrees, nClus,'prop',sizeSF)
    Vectorize(cf3, SIMPLIFY=F)(clustersSelected) %>% unlist -> idcs 
    
    LL.s <- calcLogloss(subforest(forest,idcs), data.test)
    
    assertthat::validate_that(sizeSF==length(idcs)
                              , msg =paste('Should have selected ',sizeSF,' into the subforest. Selected ', length(idcs)))
    
    mDiss <- calc_mDiss(dm,idcs)
  
    c(i,nT, length(idcs), 'clustering3', metric
      , pam.obj$silinfo$avg.width # mSW mean silhouette width 
      , pam.obj$silinfo$width[as.character(idcs),'sil_width'] %>% mean # mean silhouette width over the selected trees
      , mDiss
      , LL.s
      , LL.s - LL.d
      , LL.s - LL.r
    ) -> evalT[ct,]
    ct <- ct+1
    
    #### 4. strategy #####################################################################################################
    #### cluster into 3 clusters, randomly select from clusters according to their size until 50 trees are selected  ####
    ######################################################################################################################
    
    pam.obj <- cluster::pam(dm
                            , k= 3 # 50/5 = 10 # number needed for percentage of positive sil width!
                            , diss=TRUE
                            , medoids='random'
                            , nstart = 5)
    
    (pam.obj$silinfo$width[,'sil_width']>0) %>% 
      which %>%
      names %>% 
      as.integer -> goodTrees
    
    which(pam.obj$clusinfo[,'size'] >1) -> clustersSelected 
    
    cf4 <- function(nClus) cf(pam.obj,goodTrees,nClus,'prop',sizeSF)
    Vectorize(cf4, SIMPLIFY=F)(clustersSelected) %>% unlist -> idcs 
    
    LL.s <- calcLogloss(subforest(forest,idcs), data.test)
  
    assertthat::validate_that(sizeSF==length(idcs)
                              , msg =paste('Should have selected ',sizeSF,' into the subforest. Selected ', length(idcs)))
    
    mDiss <- calc_mDiss(dm , idcs)
    
    c(i,nT, length(idcs), 'clustering4', metric
      , pam.obj$silinfo$avg.width # mSW mean silhouette width, average over all trees
      , pam.obj$silinfo$width[as.character(idcs),'sil_width'] %>% mean # mean silhouette width over the selected trees
      , mDiss
      , LL.s
      , LL.s - LL.d
      , LL.s - LL.r
    ) -> evalT[ct,]
    ct <- ct+1
        
    #### random sub-forest ####
    ###########################
    
    mDiss <- calc_mDiss(dm , idcs.r)
    
    c(i,nT, sizeSF, 'random', metric
      , NA
      , NA
      , mDiss
      , LL.r
      , LL.r - LL.d
      , NA
    ) -> evalT[ct,]
    ct <- ct+1
  }
}

evalT %>% as.data.frame -> et
names(et)<-c('bootstrap','total','size', 'type','metric'
             , 'mSW.d', 'mSW.s'
             , 'mDiss'
             , 'logloss'
             , 'logloss.diff.d'
             , 'logloss.diff.r'
             )

for(A in names(et)[-c(4,5)]){
  et[,A] <- as.numeric(et[,A]) %>% round(5)
}

return(et)
}

# to base the result on more bootstraps
folder <- 'data/nursery'
files <- list.files(folder)[16:18]
# dir(folder)
collector <-  list()
ct <-  1 # counter for the above collector

sizeSF <- 5
#sizeDefault <- 500

for(file in files){
  # run loops over doc loaded from file
  load(paste(folder,file, sep='/'))
  LL <- calc_LL_for_selection(doc, sizeSF=sizeSF)
  # keep result from loop as list element
  collector[[ct]] <- LL
  ct <-  ct+1
}

# create single data.frame from list of data.
et <- bind_rows(collector)

View(et)

info<-paste('new preprocessing, tiny subforests, 3 cluster strategies, select positive sil width, document mean sil width for the whole forest and the selected subforest , clustering with all dissimilarities, separately. Cluster Quality. logloss, diff to default and random', sep='')
#moreInfo <-  list('nBs' = nBs, sizes=list('sizeSF'=sizeSF, 'nT'=nT , 'sizeDefault'=500), 'metrices'=metrices ,'rg' = rg , 'seed'=seed)
save(et, info, file=paste('data/cluster07-nursery/',sizeSF,'trees_',round(100*runif(1),0),'*.rda', sep=''))

#### what's in the data we created ?? ####
##########################################

et$type %>% table
et$metric %>% table

table(et$type ,et$metric)

table(et$size)

#### plots , visualize ####
###########################

et %>% filter(type %in% c('random','default')) -> et.dr # default and random 
# including the random forests and their mDiss in all 4 dissimilarity metrices
dim(et.dr)
et %>% filter(type %in% c('random',paste('clustering',1:4,sep=''))) -> et.s # selected : clustered and random

et %>% filter(type %in% paste('clustering',1:4,sep='')) -> et.c # clustered

et.s %>%
  group_by(metric , type) %>% 
  summarise(mean(logloss), sd(logloss), mean(mSW.s), mean(mSW.d))

et.s %>%
  group_by(metric , type) %>% 
  ggplot(aes(x=metric , y=logloss.diff.d, fill=type))+
  geom_boxplot()+
  labs(title='It\'s not easy to be better than random\n and quite impossible to be better than default')

# single plot for each metric
metrices <-  unique(et$metric)
for(metric in metrices){
p <- et.s[et.s$metric==metric,] %>% 
  group_by(type) %>% 
  ggplot(aes(x=type , y=logloss.diff.d))+
  geom_boxplot()+
  ggtitle(metric) 
plot(p)
}

et.s %>%
  group_by(metric , type) %>% 
  ggplot(aes(x=metric , y=logloss.diff.r, fill=type))+
  geom_boxplot()+
  labs(title='It\'s not easy to be better than random')

metrices <- c('d0','d1','d2','sb')
# Überblich, alle clustering Strategien zusammengefasst, nur nach metric unterschieden.
for(m in metrices){
  print(paste('************** ' , m , '************** '))
  et.s %>%
    filter(metric == m) %>%
    select(c("size","mDiss","mSW.d",'mSW.s','logloss')) -> da.ta
  lm(logloss~poly(mDiss,2), data=da.ta) -> lm1
  summary(lm1) %>% print
  lm(logloss~poly(mDiss,2)*mSW.s, data=da.ta) -> lm1
  summary(lm1) %>% print
  lm1 %>% 
    allEffects %>% 
    plot(multiline=T, main=paste(m , 'mean dissimilarity effect'))
  plot(da.ta$logloss~da.ta$mDiss , main=m)
}

for(m in metrices){
for(t in c('clustering1','clustering3','clustering4')){
  print(paste('************** ' , m ,' *** ', t , '************** '))
  et.s %>%
   filter(metric == m) %>%
   filter(type == t) %>%
   select(c("size","mDiss","mSW.d",'mSW.s','logloss')) -> da.ta
  lm(logloss~poly(mDiss,2), data=da.ta) %>% summary %>% print
  lm(logloss~poly(mDiss,2)*mSW.s, data=da.ta) -> lm1
  summary(lm1) %>% print
  lm1 %>% 
    allEffects %>% 
    plot(multiline=T, main=paste(m ,t , 'mean dissimilarity effect'))
  plot(da.ta$logloss~da.ta$mDiss , main=paste(m,t))
}

}
# not working
et[et$type=='default','logloss'] %>% mean

names(et)
et$type %>% unique

et %>% 
  filter(type=="default") %>%
  select(logloss)%>%
  apply(2,function(x)c(mean(x), sd(x)))

et %>% 
  filter(type=="random") %>%
  select(logloss)%>%
  apply(2,function(x)c(mean(x), sd(x)))

et %>% 
  filter(type=="clustering4") %>%
  filter(metric=='d0') %>%
  select(logloss)%>%
  apply(2,function(x)c(mean(x), sd(x)))
