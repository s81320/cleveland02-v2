# code/clusteringForests/08-simulation.R
# 18.02.2022

# simplified code, selecting from all trees in cluster, not only those with a positive silhouette width

# train on bootstraps of Cleve, previously created in nursery
# loop over multiple nurseries using a collector

# test on Hung or VA

# cluster default forest, select trees from cluster into subforest

# using cluster methods 1,3,4
# always cluster proportionally
# 1) cluster into the desired number directly
# 2) build 2 clusters 
# 3) build 3 clusters

# new: choose trees with lowest logloss on Clevland OOB

# measure of success : logloss

rm(list=ls())

library(dplyr)
library(caret)
library(ranger)
library(effects)
library(cluster)

library(xtable)

source('code/source/prep.R') # calcLogloss
source('code/source/subforest.R') # subforest

load('data/data_SupMat4.rda') # loads the data sets Cleve, Hung, Swiss, VA

# set test data by name : VA, Swiss or Hung
data.test.name <-  'Swiss'
data.test <-  get(data.test.name)
attr(data.test,'data.test.name') <- data.test.name


# load('data/nursery/nursery01_03_50x500.rda') # loads doc and info
# contains forests and their dissimilarity matrices
# info

t1<-function(tri){
  ti<-treeInfo(rg,tri)
  return(length(which(ti$terminal==TRUE))) # terminal nodes , -1 gives the number of splits
}


cf <- function(pam.obj, goodTrees, nClus, preselectMethod, targetSize, OOB, forest, dm){ # function needed for clustering strategies 2,3,4
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
  #S1 <- base::intersect(A1,goodTrees) 
  
  # check the cluster the trees belong to
  #pam.obj$clustering[S1] 
  #assertthat::validate_that(all(pam.obj$clustering[S1]==nClus)) %>% print
  
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
      # how many <- target size (=how many you finally want) * cluster size / sum of trees in all clusters (=num.trees of clustered object, but that would have to be passed as an argument)
    }
  }
    
  #if(length(S1)>howMany){ # select from intersection if it has enough elements
  #  sample(S1,howMany) -> preselect
  #}else{
  #  S1 -> preselect
  #}
    
  # final selection , making sure the medoid is returned
  #if(pam.obj$medoids[[nClus]] %in% preselect){ # has the medoid been selected?
  #    return(preselect) # if so, we're good
  #}else{ 
  #  preselect[[1]] <- pam.obj$medoids[[nClus]] # else, we put the medoid as the first element in our selection
  #  return(preselect)
  #}
  
  # final selection, no preselection
  #if(length(A1)>howMany){ # select from intersection if it has enough elements
  #  sample(A1,howMany) %>% return
  #}else{
  #  A1 %>% return
  #}
  
  
  # final selection by performance on validation set
  # get OOB observations as arguments
  if(length(A1) >howMany){
    data.set.val <- Cleve[OOB,] # goes back to original Cleveland data
    # pass full forest as argument? Select trees in cluster nClus
    # predictions of all trees in cluster nClus / indices in A1
    pp <- predict(subforest(forest,A1) 
                , data=data.set.val 
                , predict.all = T)$predictions[,2,]
    #pp <- simplify2array(pp, higher=FALSE)
    lapply(1:length(A1)
           , function(k){
             unlist(calcLogloss2(pp=pp[,k] , df=data.set.val ) ) 
             }) %>% 
             unlist -> LL
  
    #print('trees in cluster')
    #print(A1)
    #print(LL[order(LL)])
    
    # select close to medoid
    # look at: distances to cluster medoid and other trees of the cluster
    # order trees in A1 by this distance to mediod and take smallest
    # should include medoid (with dissimilarity 0)
    #print(paste('select ', howMany,' central trees from cluster ', nClus))
    A1[ order(dm[pam.obj$medoids[nClus],A1])[1:howMany] ] -> close2med 
    #close2med %>% print # indices in original forest , not in cluster
    #dm[pam.obj$medoids[nClus],close2med] %>% print
    #range(dm[pam.obj$medoids[nClus],A1]) %>% print
    #LL[order(dm[pam.obj$medoids[nClus],A1])[1:howMany]] %>% print # LL is cluster based , not suitable for indeces of forest
    return(close2med)
    
    # select trees with lowest rank(LL)
    #print(paste('select ', howMany,' high performers from cluster ', nClus))
    #which(rank(LL)<=howMany) %>% A1[.] %>% print
    #which(rank(LL)<=howMany) %>% LL[.] %>% print
    #which(rank(LL)<=howMany) %>% A1[.] %>% return
  }else{
    return(A1)
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
#nBs <- 5

evalT <-matrix(NA,nrow=(4*5+1)*nBs,ncol=11) # table of evaluations
# nrow = 4 :length(metrices)
#        5 : 4 strategies + 1 random
#        1 : default

ct <- 1
for(i in 1:nBs){
  #print(paste(i/nBs , Sys.time()))
  if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
  # data frame , bootstrapped
  
  data.train <- doc[[i]]$`bootstapped training data`

  # grow small forest (size nT)
  rg <- doc[[i]]$rg
  forest<-rg$forest
  
  OOB <-  base::setdiff(1:nrow(Cleve), unique(doc[[i]]$resample))
  #print('nr of OOB indices')
  #print(length(OOB))
  
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
  
  #### 1. strategy #######################################################
  #### cluster, select medoids (or best performer) into the subforest ####
  ########################################################################
  
  DM <- doc[[i]]$DM
  for(metric in metrices){
    #print(paste(i , ct , metric))
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

    #### 2. strategy ######################################################################################
    #### 50 : cluster into 25 clusters, take 2 best trees per cluster (best in terms of OOB of simulation)  ####
    #### 5 : cluster into 4 and select proportionally
    #######################################################################################################
    
    pam.obj <- cluster::pam(dm
                            , k=25
                          #  , k= 4 
                            , diss=TRUE
                            , medoids='random'
                            , nstart = 5)
    
    (pam.obj$silinfo$width[,'sil_width']>0) %>% 
      which %>%
      names %>% 
      as.integer -> goodTrees
    
    which(pam.obj$clusinfo[,'size'] >1) -> clustersSelected 
    
    cf3 <- function(nClus) cf(pam.obj
                              , goodTrees
                              , nClus
                              , 'const'
                              #,'prop'
                              , sizeSF
                              , OOB
                              , forest
                              , dm)
    Vectorize(cf3, SIMPLIFY=F)(clustersSelected) %>% unlist -> idcs 
    
    LL.s <- calcLogloss(subforest(forest,idcs), data.test)
    
    mDiss <- calc_mDiss(dm,idcs)
   
    c(i,nT, length(idcs), 'clustering2', metric
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
    
    #pam.obj <- cluster::pam(dm
    #                        , k= 2 
    #                        , diss=TRUE
    #                        , medoids='random'
    #                        , nstart = 5)
    
  #  (pam.obj$silinfo$width[,'sil_width']>0) %>% 
  #    which %>%
  #    names %>% 
  #    as.integer -> goodTrees
    
  #  which(pam.obj$clusinfo[,'size'] >1) -> clustersSelected 
    
   # cf3 <- function(nClus) cf(pam.obj,goodTrees, nClus,'prop',sizeSF,OOB,forest,dm)
  #  Vectorize(cf3, SIMPLIFY=F)(clustersSelected) %>% unlist -> idcs 
    
  #  LL.s <- calcLogloss(subforest(forest,idcs), data.test)
    
  #  assertthat::validate_that(sizeSF==length(idcs)
  #                            , msg =paste('Should have selected ',sizeSF,' into the subforest. Selected ', length(idcs)))
    
  #  mDiss <- calc_mDiss(dm,idcs)
    
  #  c(i,nT, length(idcs), 'clustering3', metric
  #    , pam.obj$silinfo$avg.width # mSW mean silhouette width 
  #    , pam.obj$silinfo$width[as.character(idcs),'sil_width'] %>% mean # mean silhouette width over the selected trees
  #    , mDiss
  #    , LL.s
  #    , LL.s - LL.d
  #    , LL.s - LL.r
  #  ) -> evalT[ct,]
  #  ct <- ct+1
    
    
    #### 4. strategy #####################################################################################################
    #### cluster into 3 clusters, randomly select from clusters according to their size until 50 trees are selected  ####
    ######################################################################################################################
    
    #pam.obj <- cluster::pam(dm
    #                        , k= 3 # 50/5 = 10 # number needed for percentage of positive sil width!
    #                        , diss=TRUE
    #                        , medoids='random'
    #                        , nstart = 5)
    
    #(pam.obj$silinfo$width[,'sil_width']>0) %>% 
    #  which %>%
    #  names %>% 
    #  as.integer -> goodTrees
    
    #which(pam.obj$clusinfo[,'size'] >1) -> clustersSelected 
    
    #cf4 <- function(nClus) cf(pam.obj,goodTrees,nClus,'prop',sizeSF,OOB, forest, dm)
    #Vectorize(cf4, SIMPLIFY=F)(clustersSelected) %>% unlist -> idcs 
    
    #LL.s <- calcLogloss(subforest(forest,idcs), data.test)
  
    #assertthat::validate_that(sizeSF==length(idcs)
    #                          , msg =paste('Should have selected ',sizeSF,' into the subforest. Selected ', length(idcs)))
    #
    #mDiss <- calc_mDiss(dm , idcs)
    
    #c(i,nT, length(idcs), 'clustering4', metric
    #  , pam.obj$silinfo$avg.width # mSW mean silhouette width, average over all trees
    #  , pam.obj$silinfo$width[as.character(idcs),'sil_width'] %>% mean # mean silhouette width over the selected trees
    #  , mDiss
    #  , LL.s
    #  , LL.s - LL.d
    #  , LL.s - LL.r
    #) -> evalT[ct,]
    #ct <- ct+1
        
    #### 5. strategy ######################################################################################
    #### 50 : cluster into 10 clusters, take 5 best trees per cluster (best in terms of OOB of simulation)  ####
    #### 5 : cluster into 4 and select proportionally
    #######################################################################################################
    
    pam.obj <- cluster::pam(dm
                            , k=10
                            #  , k= 4 
                            , diss=TRUE
                            , medoids='random'
                            , nstart = 5)
    
    (pam.obj$silinfo$width[,'sil_width']>0) %>% 
      which %>%
      names %>% 
      as.integer -> goodTrees
    
    which(pam.obj$clusinfo[,'size'] >1) -> clustersSelected 
    
    cf3 <- function(nClus) cf(pam.obj
                              , goodTrees
                              , nClus
                              , 'const'
                              #,'prop'
                              , sizeSF
                              , OOB
                              , forest
                              , dm)
    Vectorize(cf3, SIMPLIFY=F)(clustersSelected) %>% unlist -> idcs 
    
    LL.s <- calcLogloss(subforest(forest,idcs), data.test)
    
    mDiss <- calc_mDiss(dm,idcs)
    
    c(i,nT, length(idcs), 'clustering5', metric
      , pam.obj$silinfo$avg.width # mSW mean silhouette width 
      , pam.obj$silinfo$width[as.character(idcs),'sil_width'] %>% mean # mean silhouette width over the selected trees
      , mDiss
      , LL.s
      , LL.s - LL.d
      , LL.s - LL.r
    ) -> evalT[ct,]
    ct <- ct+1
    
    #### 6. strategy ######################################################################################
    #### 50 : cluster into 5 clusters, take 2 best trees per cluster (best in terms of OOB of simulation)  ####
    #######################################################################################################
    
    pam.obj <- cluster::pam(dm
                            , k=5
                            #  , k= 4 
                            , diss=TRUE
                            , medoids='random'
                            , nstart = 5)
    
    (pam.obj$silinfo$width[,'sil_width']>0) %>% 
      which %>%
      names %>% 
      as.integer -> goodTrees
    
    which(pam.obj$clusinfo[,'size'] >1) -> clustersSelected 
    
    cf3 <- function(nClus) cf(pam.obj
                              , goodTrees
                              , nClus
                              , 'const'
                              #,'prop'
                              , sizeSF
                              , OOB
                              , forest
                              , dm)
    Vectorize(cf3, SIMPLIFY=F)(clustersSelected) %>% unlist -> idcs 
    
    LL.s <- calcLogloss(subforest(forest,idcs), data.test)
    
    mDiss <- calc_mDiss(dm,idcs)
    
    c(i,nT, length(idcs), 'clustering6', metric
      , pam.obj$silinfo$avg.width # mSW mean silhouette width 
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
files <- list.files(folder)[1:5]
# dir(folder)
collector <-  list()
ct <-  1 # counter for the above collector

sizeSF <- 50
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

info <-  list(sizes=list('sizeSF'=sizeSF
                               , 'sizeDefault'=500)
                  , 'metrices'=c('d0','d1','d2','sb')
              , 'num.clusters' = c(50,25,10,5)
                  , 'cluster.strategies'=list('clustering1'='50 clusters, select 1 central tree: the medoid'
                                             , 'clustering2'='25 clusters, select 2 most central trees'
                                            , 'clustering5'='10 clusters, 5 most central trees'
                                            , 'clustering6'='5 clusters, 10 most trees')
                  , 'files'=files
                  , 'data.test'=attr(data.test,'data.test.name')
                  )
save(et, info, file=paste('data/simulation08/validation/selection_central_',sizeSF,'trees_',round(100*runif(1),0),'*.rda', sep=''))

#load('data/simulation08/validation/selection_hp_5trees_01.rda')
#### what's in the data we created ?? ####
##########################################

et$type %>% table
et$metric %>% table


table(et$type ,et$metric) %>%
  xtable -> et.xt
digits(et.xt) <- 3
et.xt

table(et$type ,et$metric) 

table(et$size) %>% xtable

#### plots , visualize ####
###########################

et %>% filter(type %in% c('random','default')) -> et.dr # default and random 
# including the random forests and their mDiss in all 4 dissimilarity metrices
dim(et.dr)
et %>% filter(type %in% c('random',paste('clustering',1:6,sep=''))) -> et.s # selected : clustered and random

et %>% filter(type %in% paste('clustering',1:6,sep='')) -> et.c # clustered

et.s %>%
  group_by(metric , type) %>% 
  # summarise(mean(logloss), sd(logloss), mean(mSW.s), mean(mSW.d)) -> xt
  summarise(mean(logloss), sd(logloss)) -> xt

xt%>%
  xtable -> xt
digits(xt) <- 4
xt


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
for(t in base::intersect(unique(et.s$type), paste('clustering',1:6,sep=''))){
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
  filter(type=="clustering2") %>%
  filter(metric=='d2') %>%
  select(logloss)%>%
  apply(2,function(x)c(mean(x), sd(x)))

et %>% 
  group_by(type) %>%
  summarise(m=mean(logloss), s=sd(logloss))

et %>% 
  group_by(metric) %>%
  summarise(m=mean(logloss), s=sd(logloss))

et %>% 
  group_by(type, metric) %>%
  summarise(m=mean(logloss), s=sd(logloss)) %>%
  as.data.frame -> res1

res1$mLL.rank <-  rank(res1$m)
res1$sdLL.rank <-  rank(res1$s)

res1  %>% 
  select(type, metric,mLL.rank , sdLL.rank) %>%
  arrange(sdLL.rank) %>% xtable -> res1.xt
digits(res1.xt) <- 0
res1.xt

names(res1)
unique(res1$type)
res2 <- res1 %>% 
  filter(type!='random') %>%
  filter(type!='default')

pch <- rep(0,nrow(res2))
t <- unique(res2$type)
pch.vec <- c(1,2,3,4) # colours
for(i in 1:length(pch.vec)){
  pch[which(res2$type==t[i])] <- pch.vec[i]
}

col <- rep(0,nrow(res2))
m <- unique(res2$metric)
col.vec <- c(2,3,4,7) # plot character
for(i in 1:length(col.vec)){
  col[which(res2$metric==m[i])] <- col.vec[i]
}

plot(res2$m
     , res2$s
     , col=col
     , pch=pch
     , xlab='mean logloss'
     , ylab='standard deviation of logloss'
     , xlim=range(res1$m)
     , ylim=range(res1$s)
     , main='mean and std dev of logloss for\ndifferent numbers of clusters (=cls) and dissimilarities')
# default forest
points(res1 %>% filter(type=='default') %>% select(c(m,s)), col=1,pch=8)
# text(res1 %>% filter(type=='default') %>% select(c(m,s)) + 0.01*c(1,1), 'default')
# text(res1 %>% filter(type=='default') %>% select(c(m,s)) - 0.0008*c(1,0), 'default')
text(res1 %>% filter(type=='default') %>% select(c(m,s)) + 0.001*c(0,1), 'default')

points(res1 %>% filter(type=='random') %>% select(c(m,s))%>% unique, col=1,pch=8)
# text(res1 %>% filter(type=='random') %>% select(c(m,s)) %>% unique + 0.01*c(1,0), 'random')
# text(res1 %>% filter(type=='random') %>% select(c(m,s)) %>% unique - 0.001*c(1,0), 'random')
text(res1 %>% filter(type=='random') %>% select(c(m,s)) %>% unique + 0.001*c(0,1), 'random')

legend('bottomright'
       #, legend=c('5 cls', '4 cls' , '2 cls' , '3 cls' , 'd0','d1','d2','sb')
       #, legend=c('50 cls', '25 cls' , '2 cls' , '3 cls' , 'd0','d1','d2','sb')
       , legend=c('50 cls', '25 cls' , '10 cls' , '5 cls' , 'd0','d1','d2','sb')
       , pch = c(pch.vec,rep(15,4))
       , col = c(rep(9,4),col.vec)
       , cex = 0.7)


