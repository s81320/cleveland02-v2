# code/clusteringForests/06.R
# 05.01.2022

# using cluster methods 1,3,4
# always cluster proportionally
# 1) cluster into the desired number directly
# 2) build 2 clusters 
# 3) build 3 clusters

# NEW 

# measure of success : 
# logloss # new # as discussed with A.Z.
# mae 
# 1-empirical distribution function of absolute error at 0.15 (how many get an absolute error more than 15%)
# all measures of success are designed to be better when smaller

# error analysis

rm(list=ls())

library(dplyr)
library(caret)
library(ranger)
library(effects)

source('code/source/sim-prep.R')
# reads Cleveland data set , returns enriched data frame df
# with probabilities generated and added as column prob

# 3) draw a bootstrap sample (or many...) from the enriched Cleveland data set
# dfb : data frame bootstrapped

source('code/source/distance-matrices.R') 
source('code/source/helper-functions.R')

t1<-function(tri){
  ti<-treeInfo(rg,tri)
  return(length(which(ti$terminal==TRUE))) # terminal nodes , -1 gives the number of splits
}

errors <-  function(dfb, train, forest, idcs){
  forest %>%
    subforest(idcs) %>%
    predict(data=dfb[-train,]) %>%
    .$predictions %>%
    .[,2] -> pp
  
  correctedpp <- ifelse(dfb[-train,'CAD']=='Yes',pp,1-pp) # problematic when this returns 0
  #which(correctedpp==0) %>% print
  
  # function copied from Calibration_Sup4_V01.R
  winsorize_probs <- function(x) {
     eps <- .Machine$double.eps
     pmax(pmin(x, 1 - eps), eps)
   }
  correctedpp <- winsorize_probs(correctedpp) # Vectorize?
  
  error <- pp-dfb[-train,'prob']
  # error by sex
  error.f <- error[dfb[-train,'Sex']=='Female'] 
  error.m <- error[dfb[-train,'Sex']=='Male'] 
  
    return(list(logloss= -mean(log(correctedpp))
                , mae=mean(abs(error))
                , unacceptable=1-ecdf(abs(error))(0.15)
                )
           )
}

cf <- function(nClus,preselectMethod,targetSize){ # function needed for clustering strategies 2,3,4
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
      howMany <- round(targetSize*pam.obj$clusinfo[nClus,1] / forest$num.trees,0) 
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

metrices<-c('d0','d1','d2','sb')
#metrices <- 'd0'
sizeSF <- 50

# used to be : nT <- c(rep(21,700),rep(41,700), rep(61,700), rep(81,700), rep(101,700))
nT <- 500
sizeDefault <- 500
nBs <-  50
# used to be length(sizeSF)

evalT <-matrix(NA,nrow=17*nBs,ncol=15) # table of evaluations

seed <- 2
set.seed(seed)
cr<-createResample(df$CAD
                   , times = nBs
                   , list = T)

set.seed(seed+1)

ct <- 1
for(i in 1:nBs){
  print(paste(i/nBs , Sys.time()))
  #if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
  # data frame, bootstrapped
  dfb<-df[cr[[i]],] 
  dim(dfb)
  
  dfb$CAD <- factor(rbinom(303,1,dfb$prob))
  levels(dfb$CAD)<-list("No"="0","Yes"="1")
  
  # split rows into training and test
  train <- createDataPartition(y=dfb$CAD , p=0.8 )[[1]]
  
  # grow small forest (size nT)
  rg <- ranger(CAD~.
               , data = dfb[train,-12] # exclude the probabilities
               , num.trees = nT
               , replace = F # neu nach Absprache mit AZ
               , mtry= 3 # default : 3
               , importance = 'impurity'
               , probability = T # this makes it a random forest of type 'Probability Estimation'
               , min.node.size = 13 # optimized in A.Z. paper: BiomJ 56, 2014, 
  )
  forest<-rg$forest
  
  #### random subforest ####
  #### will be used later, when dissimilarity matrics are calculated
  idcs.r <-  sample(1:forest$num.trees , sizeSF , replace=F)
  
  errors(dfb=dfb
         , train=train
         , forest=forest
         , idcs=idcs.r) -> e.r
  
  #### default forest ####
  ########################
  
  idcs <- 1:sizeDefault

  errors(dfb=dfb , train=train, forest=forest, idcs=idcs)-> e.d
  
  c(i,nT
    , sizeDefault
    , 'default'
    , 'none'
    , NA
    , NA
    , NA
    , e.d %>% unlist %>%  as.vector 
    , NA , NA
    , e.d$logloss - e.r$logloss
    , e.d$mae - e.r$mae)  -> evalT[ct,]
  ct <- ct+1
  

  #### 1. strategy ###################################
  #### cluster, select medoids into the subforest ####
  ####################################################
  
  for(metric in metrices){
    #print(metric)

    createDM(forest=forest , type=metric, dft=dfb[train,]) -> dm
  
    pam.obj <- cluster::pam(dm
                          , k=sizeSF
                          , diss=TRUE
                          , medoids='random'
                          , nstart = 5
                          )
  
    idcs <- pam.obj$medoids
  
    Dee <- dm[idcs,idcs]
    upper.tri(Dee, diag=F) %>% 
      Dee[.] %>%
      mean -> mDiss

    e.s <- errors(dfb=dfb, train=train, forest=forest, idcs=idcs) # error for selected sub-forest
    
    c(i,nT, sizeSF, 'clustering1', metric
      , pam.obj$silinfo$avg.width # mSW mean silhouette width
      , pam.obj$silinfo$width[as.character(idcs),'sil_width'] %>% mean # mean silhouette width over the selected trees
      , mDiss
      , e.s %>% unlist %>% as.vector
      , e.s$logloss - e.d$logloss
      , e.s$mae - e.d$mae
      , e.s$logloss - e.r$logloss
      , e.s$mae - e.r$mae
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
    
    cf3 <- function(nClus) cf(nClus,'prop',sizeSF)
    Vectorize(cf3, SIMPLIFY=F)(clustersSelected) %>% unlist -> idcs 
    
    assertthat::validate_that(sizeSF==length(idcs)
                              , msg =paste('Should have selected ',sizeSF,' into the subforest. Selected ', length(idcs)))
    
    Dee <- dm[idcs,idcs]
    upper.tri(Dee, diag=F) %>% 
      Dee[.] %>%
      mean -> mDiss
  
    e.s <- errors(dfb=dfb, train=train, forest=forest, idcs=idcs)
    
    c(i,nT, length(idcs), 'clustering3', metric
      , pam.obj$silinfo$avg.width # mSW mean silhouette width 
      , pam.obj$silinfo$width[as.character(idcs),'sil_width'] %>% mean # mean silhouette width over the selected trees
      , mDiss
      , e.s %>% unlist %>% as.vector
      , e.s$logloss - e.d$logloss
      , e.s$mae - e.d$mae
      , e.s$logloss - e.r$logloss
      , e.s$mae - e.r$mae
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
    
    cf4 <- function(nClus) cf(nClus,'prop',sizeSF)
    Vectorize(cf4, SIMPLIFY=F)(clustersSelected) %>% unlist -> idcs 
  
    assertthat::validate_that(sizeSF==length(idcs)
                              , msg =paste('Should have selected ',sizeSF,' into the subforest. Selected ', length(idcs)))
    
    Dee <- dm[idcs,idcs]
    upper.tri(Dee, diag=F) %>% 
      Dee[.] %>%
      mean -> mDiss
    
    e.s <- errors(dfb=dfb, train=train, forest=forest, idcs=idcs)
    
    c(i,nT, length(idcs), 'clustering4', metric
      , pam.obj$silinfo$avg.width # mSW mean silhouette width, average over all trees
      , pam.obj$silinfo$width[as.character(idcs),'sil_width'] %>% mean # mean silhouette width over the selected trees
      , mDiss
      , e.s %>% unlist %>% as.vector
      , e.s$logloss - e.d$logloss
      , e.s$mae - e.d$mae
      , e.s$logloss - e.r$logloss
      , e.s$mae - e.r$mae
    ) -> evalT[ct,]
    ct <- ct+1
        
    #### random sub-forest ####
    ###########################
    
    Dee <- dm[idcs.r,idcs.r]
    upper.tri(Dee, diag=F) %>% 
      Dee[.] %>%
      mean -> mDiss
    
    c(i,nT, sizeSF, 'random', metric
      , NA
      , NA
      , mDiss
      , e.r  %>% unlist %>% as.vector # was calculated when icds.r was set , calculated only once, used again for each dissim
      , e.r$logloss - e.d$logloss
      , e.r$mae - e.d$mae
      , NA
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
             , 'mae'
             , 'unacceptable'
             , 'logloss.diff.d' ,  'mae.diff.d'
             , 'logloss.diff.r' ,  'mae.diff.r'
             )

for(A in names(et)[-c(4,5)]){
  et[,A] <- as.numeric(et[,A]) %>% round(5)
}

View(et)

info<-paste('new preprocessing, tiny fubforests, 3 cluster strategies, select positive sil width, document mean sil width for the whole forest and the selected subforest. seed: ', seed, ' , clustering with all dissimilarities, separately. Cluster Quality. Error analysis. logloss, diff to default', sep='')
moreInfo <-  list('nBs' = nBs, sizes=list('sizeSF'=sizeSF, 'nT'=nT , 'sizeDefault'=sizeDefault ), 'metrices'=metrices ,  'pam.obj'=pam.obj , 'rg' = rg )
save(et, info, moreInfo , file=paste('data/cluster/06*',round(100*runif(1),0),'.rda', sep=''))

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
  ggplot(aes(x=metric , y=logloss.diff.r, fill=type))+
  geom_boxplot()+
  labs(title='It\'s not easy to be better than random')
# et.s[et.s$metric=='d0',] %>% group_by(type) %>% ggplot(aes(x=type , y=logloss.diff,))+geom_boxplot()

et %>%
  group_by(metric , type) %>% 
  ggplot(aes(x=metric , y=unacceptable, fill=type))+
  geom_boxplot()+
  labs(title='all strategies have a lower median than random')

et.s %>%
  group_by(metric , type) %>% 
  ggplot(aes(x=metric , y=logloss.diff.r, fill=type))+
  geom_boxplot()+
  labs(title='It\'s not easy to be better than random')

# works for d1 ?? only?
et.s %>%
  filter(metric =='d1') %>%
  select(c("size","mDiss","mSW.d",'mSW.s','logloss')) %>%
  lm(logloss~mDiss, data=.) -> lm1
summary(lm1)
lm1 %>% 
  allEffects %>% 
  plot(multiline=T)

# not working
et[et$type=='default','logloss'] %>% mean

names(et)
et$type %>% unique
et %>% 
  filter(type=="default") %>%
  select(logloss)%>%
  apply(2,mean)

et %>% 
  group_by(metric, size, type) %>%
  summarize(mean(logloss))%>% data.frame()

