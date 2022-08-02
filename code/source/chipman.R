print('sourcing c h i p m a n .R : functions for building and evaluating the subforests as descriped in Chipman, Georg, McCollouch (1998), based on diversity and representation.')

# functions to be called inside calc_LL_for selection
grow_chipForest_1 <-  function(dm, oLL, parameter, output=F){
  #' grows the Chipman 1 forest ('core technique' in the paper by Chipman, George, McCollouch) 
  #'
  #' needs data.test in parent environment
  #' @param parameter needs to have cutoff and selection
  #' @param output if TRUE the function outputs also the p-values of dip tests and the silhouette widths when clustering
  
  #print(dim(dm))
  #print(pa)
  
  # level of representation , remove the diagonal, it has only 0 values
  alpha <- quantile(dm[upper.tri(dm)],parameter$cutoff) # diagonal with all values 0 not included
  
  represented <-  function(I){
    #' testing for representation of best I trees (R= dense representing forest = {1,..,I}) at level alpha
    #' 
    #' uses N, dm2 and alpha from parent environment. i.e. function calc_chipForest_1
    #' for I >1 , 
    #' representation by one tree alone can happen only at level of largest dissimilarity
    #' no representation by almost all trees (I==N-1) only if last tree, 
    #' index N, is far away from all others, further than cutoff
    #assertthat::assert_that(class(I)=='numeric', 'argument I has to be numeric') %>% print
    #assertthat::assert_that( (I >0) , 'I must be at least 1 and at most one less than trees in forest') %>% print
    
    if(I==1){ # max
      dm2[1,2:N] %>% max %>%
        (function(x){x<=alpha}) 
    }else{
      if(I==N-1){ # min
        dm2[1:(N-1),N] %>% min %>%
          (function(x){x<=alpha})
      }else{ # min and max
        dm2[1:I,(I+1):N] %>%
          apply(2,min) %>% # not working for I = 1 or I=N-1
          max %>%
          (function(x){x<=alpha})
      }
    }
  }
  
  dm2 <-  dm[oLL,oLL]
  
  N <-  nrow(dm2) # used to be min(nrow(dm2),N) , N would stop early as soon as N trees were selected into the dense representing subforest.
  # But this function is for the original approach of the Chipmen, and they would not have done early stopping.
  
  if(alpha<max(dm2)){
    I <- 2
    while((I<N) && !represented(I)){
      I <- I+1
    }
  }else{I <- 1}
  # exits with smallest I with represented(I) TRUE
  
  # corrected error 10.7.22: Used to be R <-  1:I
  R <-  oLL[1:I] # (dense) representing subforest , dense as all trees inbetween (from best oLL[1] to I-th best oLL[I]) are included in the subforest
  kOpt <-  NA
  multimod <- NA
  
  # I will NOT be changed further down , R may be reduced by clustering
  # attempt clustering only if the representing forest is large (larger than 10)
  if(length(R)>10){ # corrected 10.7.22 , used to be I>10
    # attempting to cluster makes sense only for R relatively large. Larger than 10 trees, say. 
    # And number of cluster should be b/w 5 and 50? And not larger than I/4?
    # dissimilarity matrix needed only for R
    dm2 <-  dm[R,R] # same as dm[oLL[1:I],oLL[1:I]]
    # which 2 trees are meant in dm2[1,2] ?? the best and the second best , indexed oLL[1], oLL[2]
  
    {
      mds.dim <-  100
      mp <- cmdscale(d=dm,k=mds.dim) # returns points in R^k
      mp[R,] %>% dist %>% as.matrix -> dm3 # trees in representing subforest R as elements in R^50, distance matrix for them
    
      # dip test for all trees in R (dense representing forest)
      #lapply(R , function(i) dip.test(dm2[i,-i], simulate=T)$p.value) %>% # old code, replaced , 10.5.2022 , should run dip test on continuous data, on d0 , sb dip test will always reject hypothesis of unimodality because data is discreet
      # corrected 10.7.22 : used to be lapply(R, function...)
      lapply(1:length(R) , function(i) dip.test(dm3[i,-i], simulate=T)$p.value)  -> doc.pv # document intermediate results of p-values
      doc.pv %>% 
        unlist %>%
        min %>% 
        (function(x) x<0.05) -> multimod # logical , TRUE if pvalue of at least one dip test less than 0.05
    } # block to calculate multimod (should be a function! based on original dm)
  
    # cluster only if unimodality is rejected
    if(multimod){ # if hypothesis of unimodality is rejected : cluster
      #print('multimod')
      hc <- cluster::agnes(dm2 
                         , method="ward" # always optimal - whenever I did hpo
                         , diss=T
                         , keep.diss=F
                         , keep.data=F)
    
      col1 <- list()
      # to check all feasible numbers of clusters : lapply(2:I-1) but we only want sizes between 5 and 50 (arbitrary! faster!) or I-1 if that is smaller than 50
      # we started with I>10 in previous if condition
      lapply(5:min(50,I-1), 
           function(k){
             cutree(hc, k = k) %>% 
               silhouette(dmatrix=dm2)%>%
               .[,'sil_width'] %>%
               mean
           } ) %>% unlist -> col1
      
      # cluster only if cluster quality is somehow OK, silhouette width of 0.2
      if(max(col1)>0.2){
        col1 %>% 
          unlist %>% 
          which.max + 4 -> kOpt # which max returns a number in 1:(min(50,I-1)-4), we have to transform this back to the original range: 5 to 50 (or I-1)

        hc.clus <- cutree(hc, k = kOpt)
    
        # from each cluster select the tree with the smallest logloss (on OOB data)
        if(parameter$selection=='best'){
          lapply(1:kOpt,function(i) which(hc.clus==i)[1] %>% oLL[.]) %>% # we take the first, because indices are ordered, first is smallest
          # oLL to go back to original index
          unlist -> R # representing subset
          # if indices were not ordered, we'd do 
          #lapply(1:sizeSF,function(i) which(hc.clus==i) %>% min ) %>% unlist
          }
      
       "if(parameter$selection=='central'){
          lapply(1:kOpt,
             function(i){
               x <-  which(hc.clus==i) # trees of cluster i , based on dm2 , based on oLL
               lapply(1:length(x), 
                      # changed to dm2 in sum , used to be dm # 31.5.2022
                      function(j){sum(dm2[x[j],x])} %>% # sum of dissimilarities to all other trees in the same cluster (cluster i)
                        unlist ) %>% 
                 which.min %>% # smallest sum of dissimilarities
                 x[.] %>% # zugehöriges Element in x
                 oLL[.] # original index in default forest
               }
          ) %>% unlist -> R
        }" # old code , selects first (=best) element if trees are undistinguishable under given dissimilarity - needs random sampling!
        if(parameter$selection=='central'){
          lapply(1:kOpt, # i-th cluster
                 function(i){
                   tris.in.clus.i <-  which(hc.clus==i) # trees of cluster i , based on dm2 , based on oLL
                   lapply(1:length(tris.in.clus.i), 
                          # changed to dm2 in sum , used to be dm # 31.5.2022
                          function(j){
                            tri <- tris.in.clus.i[j]
                            sum(dm2[tri,tris.in.clus.i])} %>% # sum of dissimilarities to all other trees in the same cluster (cluster i)
                            unlist ) %>%
                     #(function(x){
                    #   print('dissimilarities')
                    #   print(unlist(x))
                    #   x}) %>%
                     which.min -> center.clus.i # smallest sum of dissimilarities
                 #  print(paste('center of cluster', i, ':', center.clus.i))
                   # is x.i.central unique?
                   center.indistinguishable <- which(dm2[tris.in.clus.i[center.clus.i],tris.in.clus.i]==0)
                  # print(paste('# indistinguishable from central',length(center.indistinguishable)))
                   # if not unique , select randomly from those undistinguishables , else we have the best performing by design but not intentionally
                   if( center.indistinguishable %>% length > 0){
                     center.clus.i <- sample(center.indistinguishable,1)
                   }
                   center.clus.i %>%
                     tris.in.clus.i[.] %>% # associated element in cluster
                     oLL[.] # original index in default forest
                 }
          ) %>% unlist -> R
        }
      }
    }
  }
  
  doc.ret <- list('forest'=R
                  , 'multimod'=multimod
                  , 'size.dense.representing.sf' = I
                  , 'alpha'= alpha
                  , 'call'=list(name='grow_chipForest_1' , 'parameter'=parameter) )
  
  if(output){
    doc.ret$mds.dim <- mds.dim
    doc.ret$pvalues <- unlist(doc.pv)
    if(multimod){
      doc.ret$sil <- unlist(col1) # col1 exists only in case of rejected unimodality
      }
    }
  
  return(doc.ret)
}

eval_chipForest_1 <- function(dm , forest, oLL, parameter){
  
  chip1 <- grow_chipForest_1(dm , oLL, parameter, output=F)
  
  return(list('logloss'=calcLogloss( subforest(forest, chip1$forest ), data.test)
              , 'size.dense.representing.sf'= chip1$size.dense.representing.sf # I
              , 'size.sf'=length(chip1$forest) # when multimodal, (opt.) number of clusters is the size
              #, kOpt # produces NA for unimodal R , allows to count unimodal vs multimodal cases
  )
  )
}

calc_chipForest_1 <- function(dm, forest , oLL, parameter){ 
  return(eval_chipForest_1(dm , forest, oLL, parameter)) 
}
 
grow_chipForest_1_simplified <-  function(dm, oLL, parameter, output=F){
  #' grows the Chipman 1 forest ('core technique' in the paper by Chipman, George, McCollouch) enforcing clustering with a fixed number of clusters
  #'
  #' needs data.test in parent environment
  #' @param dm full dissimilarity matrix to calculate the level of representation via quantils 
  #' @param parameter needs to have cutoff, selection, sizeSF
  #' @param output if TRUE the function outputs also ...
  
  #print(dim(dm))
  #print(pa)
  
  sizeSF <-  parameter$sizeSF
  # level of representation , remove the diagonal, it has only 0 values
  alpha <- quantile(dm[upper.tri(dm)],parameter$cutoff) # diagonal with all values 0 not included
  #oLL <-  order(LL)
  
  represented <-  function(I){
    #' testing for representation of best I trees (R= dense representing forest = {1,..,I}) at level alpha
    #' 
    #' uses N, dm2 and alpha from parent environment. i.e. function calc_chipForest_1
    #' for I >1 , 
    #' representation by one tree alone can happen only at level of largest dissimilarity
    #' no representation by almost all trees (I==N-1) only if last tree, 
    #' index N, is far away from all others, further than cutoff
    #assertthat::assert_that(class(I)=='numeric', 'argument I has to be numeric') %>% print
    #assertthat::assert_that( (I >0) , 'I must be at least 1 and at most one less than trees in forest') %>% print
    
    if(I==1){ # max
      dm2[1,2:N] %>% max %>%
        (function(x){x<=alpha}) 
    }else{
      if(I==N-1){ # min
        dm2[1:(N-1),N] %>% min %>%
          (function(x){x<=alpha})
      }else{ # min and max
        dm2[1:I,(I+1):N] %>%
          apply(2,min) %>% # not working for I = 1 or I=N-1
          max %>%
          (function(x){x<=alpha})
      }
    }
  }
  
  dm2 <-  dm[oLL,oLL]
  
  N <-nrow(dm2) # used to be min(nrow(dm2),N) , N would stop early as soon as N trees were selected into the dense representing subforest.
  # But this function is for the original approach of the Chipmen, and they would not have done early stopping.
  
  if(alpha<max(dm2)){
    I <- 2
    while((I<N) && !represented(I)){
      I <- I+1
    }
  }else{I <- 1}
  # exits with smallest I with represented(I) TRUE
  
  # corrected error 10.7.22: Used to be R <-  1:I
  R <-  oLL[1:I] # (dense) representing subforest , dense as all trees inbetween (from best oLL[1] to I-th best oLL[I]) are included in the subforest
  
  # I will NOT be changed further down , R may be reduced by clustering
  # attempt clustering only if the representing forest is large (larger than 10)
  if(length(R)>sizeSF){ # cluster only if you got more than you want!
    # attempting to cluster makes sense only for R relatively large. Larger than 10 trees, say. 
    # And number of cluster should be b/w 5 and 50? And not larger than I/4?
    # dissimilarity matrix needed only for R
    dm2 <-  dm[R,R] # using R # funny name in an R script ...
    # which 2 trees are meant in dm2[1,2] ?? the best and the second best , indexed oLL[1], oLL[2]
    
    # always cluster
    cluster::pam(x=dm2
                 , k=sizeSF
                 , diss=T
                 , medoids = 1:sizeSF # start with the best trees as medoids
                 , cluster.only=F # if T, are medoids returned? We need them!
    ) -> pam.clus
    
    # from each cluster select the tree with the smallest logloss (on OOB data)
    if(parameter$selection=='best'){
      lapply(1:sizeSF,function(i) which(pam.clus$clustering==i)[1] %>% oLL[.]) %>% # we take the first, because indices are ordered, first is smallest
        # oLL to go back to original index
        unlist -> cf # representing subset
      # if indices were not ordered, we'd do 
      #lapply(1:sizeSF,function(i) which(pam.clus$clustering==i) %>% min ) %>% unlist -> cf.alt
    }
    
    if(parameter$selection=='central'){
      pam.clus$medoids -> cf
    }
    
  }else{
    cf <-  R
  } # if more trees than needed for sizeSF
  
  
  doc.ret <- list('forest'=cf
                  , 'size.dense.representing.sf' = I
                  , 'alpha'= alpha
                  , 'call'=list(name='grow_chipForest_1b' , 'parameter'=parameter) )
  
  if(output){
    doc.ret$pam.sizes <- pam$clusinfo[,1] # cluster sizes
    doc.ret$quality <- pam$silinfo # silhouette width as quality measure
  }
  
  return(doc.ret)
}

eval_chipForest_1_simplified <-  function(dm , forest, oLL, parameter, output=F){
  cf1b <-  grow_chipForest_1_simplified(dm , oLL, parameter)
  
  chipForest <-  cf1b$forest
  
  #print('eval_chipForest_1_simplified will return')
  #print(list(calcLogloss( subforest(forest, chipForest ), data.test)
  #           , cf1b$size.dense.representing.sf # I
  #           , length(chipForest)))
  return(list(calcLogloss( subforest(forest, chipForest ), data.test)
              , cf1b$size.dense.representing.sf # I
              , length(chipForest))) 
}

calc_chipForest_2 <- function(dm , forest, oLL, parameter){
  #' old name for function to grow and evaluate a Chipman 2 forest.
  eval_chipForest_2(dm , forest, oLL, parameter, output=F)}

grow_chipForest_2 <- function(dm , oLL, parameter, output=F){ # 7.7.22 : removed forest as an argument ; it is not neccessary amd not used...‚
  #' grows the Chipman forest by adding diverse trees , Chipman 2
  #' 
  #' @param parameter is a list with items cutoff and sizeSF 
  #'
  #' returns the forest
  
  chipForest <- oLL[1] # start chipForest with the best tree
  # loDiv <- quantile(dm,parameter$cutoff)  # old
  loDiv<- quantile(dm[upper.tri(dm)],parameter$cutoff)
  
  # go through trees (ordered by performance) until sizeSF trees are collected / selected
  for(I in 2:length(oLL)){
    if(length(chipForest)==parameter$sizeSF){
      I<-I-1
      break
    }else{
      trindx <- oLL[I]
      # if new tree far from (current) chipForest , add it to chipForest
      # if cutoff =0 then each tree is added with certainty -> same as unimodal!
      if(min(dm[chipForest,trindx]) > loDiv){ # if adding trindex to R keeps it diverse: then add it
        chipForest <- c(trindx, chipForest)
      }
    }
  } # for exits with j the first index after all trees have been collected
  parameter$size.dense.representing.sf <- I # former, misleading name : parameter$represented : all trees are represented!
  return(list('forest'=chipForest, 'parameter'=parameter , 'type'='Chipman2'))
}

eval_chipForest_2 <-  function(dm , forest, oLL, parameter, output=F){
  
  cf2 <-  grow_chipForest_2(dm , oLL, parameter)
  
  chipForest <-  cf2$forest
  I <- cf2$parameter$size.dense.representing.sf
  
  return(list(calcLogloss( subforest(forest, chipForest ), data.test)
              , I
              , length(chipForest))) 
}

grow_meinForest<- function(dm , LL, parameter, output=T){
  #' calculates the Meiner forest
  #' 
  #' @param parameter is a list with items cutoff and sizeSF 
  #'
  #' returns its logloss on data.test (from higher environment)
  #' returns the rank of the last inspected tree (or the index +1 , like 6 when collecting the first 5 trees)
  #' returns the number of trees in the Meiner forest (enforced to be parameters$sizeSF or smaller)
  
  oLL <- order(LL)
  #print(parameter)
  meinForest <- oLL[1] # start chipForest with the best tree
  
  alpha <- quantile(dm[upper.tri(dm)], parameter$cutoff) # level of representation
  #print(cutoff)
  
  if(output==T){
    ct <- 1
    doc.build <- list()
  }
  # go through trees (ordered by performance) until sizeSF trees are collected / selected
  for(I in 2:length(oLL)){
    if(length(meinForest)==parameter$sizeSF){
      I <- I-1
      break
    }else{
      trindx <- oLL[I]
      # if new tree far from (current) meinForest , add the best tree (lowest logloss) to the meinForest that represents trindx.
      # if alpha==0 then each tree is added with certainty -> same as unimodal!
      
      if(min(dm[meinForest,trindx]) > alpha){ # if meinForest does not represent trindx
        allCandTrees <- base::setdiff(oLL[1:I], meinForest) # all candidate trees
        allCandTrees %>% # low logloss, trees not yet selected into meinForest
          dm[.,trindx] %>%
          (function(x) which(x<=alpha)) %>% # a mask , for representation
          allCandTrees[.] -> closeCandTrees # candidates close to trindex
        meinForest <- c(meinForest,closeCandTrees[[1]]) # order by LL is kept! index of best tree is first
        # if this is empty , add trindex!
        # will never be empty! trindex is always in closeCandTrees , it has dissim 0
        "if(length(closeCandTrees)==0){
          meinForest <- c(meinForest,trindx) # if nothing better can be found : add trindx
        }else{
          closeCandTrees %>%
            LL[.] %>%
            which.min %>%
            closeCandTrees[.] %>%
            c(meinForest) -> meinForest
        }" # removed code if(closeCandTrees is empty) because it never is! # 31.5.2022
        
        if(output){
          doc.build[[ct]] <-  list(LL=LL
                                   , oLL=order(LL)
                                   , I=I
                                   , trindx=oLL[I]
                                   , allCandTrees=allCandTrees
                                   , closeCandTrees=closeCandTrees
                                   , meinForest=meinForest)
          ct <- ct+1
        }
      }
    }
  } # for exits with I the first index after all trees have been collected
  # should be 500 for sizeSF=500 (unstopped Meiner) but may be smaller for sizeSF <500 (stopped Meiner)
    
  parameter$represented <- I
  doc.ret <- list(forest=meinForest
                #  , represented = I
                  , parameter=parameter
                  , call='grow_meinForest'
                  , type='Meiner')
  
  if(output==T){
    doc.ret$progress <-  doc.build
  }
 
  return(doc.ret)
  }

eval_meinForest <- function(dm , forest, LL, parameter){
  #' calls function to calculate the Meiner forest and evaluates it
  #' returns logloss for data.test from parent directory
  
  mf <-  grow_meinForest(dm , LL, parameter)
  
  meinForest <-  mf$forest
  I <- mf$parameter$represented
  
  return(list(calcLogloss( subforest(forest, meinForest ), data.test)
              , I
              , length(meinForest))) 
}

calc_meinForest <- function(dm , forest, LL, parameter){
  #' calls function to grow the Meiner forest and evaluates it
  return(eval_meinForest(dm , forest, LL, parameter)) 
  }

calc_LL <- function(forest , data){
  #' calculate logloss for trees in forest,
  #' logloss evaluated on predictions on data
  #' 
  #' naming is bad, should be LL_for_all_trees or something. 
  #' name calc_LL is too close to calc_LL_for_selection which is probably also not a good name ...
  
  pp <- predict(forest 
                , data=data 
                , predict.all = T)$predictions[,2,]
  
  lapply(1:forest$num.trees
         , function(t){ 
           pp[,t] %>% 
             calcLogloss2( df=data ) %>% 
             unlist
         }) %>% 
    unlist 
}

calc_oLL <- function(forest , data){
  #' calculate ordered logloss for trees in forest,
  #' logloss evaluated on predictions on on data
  
  pp <- predict(forest 
                , data=data 
                , predict.all = T)$predictions[,2,]
  
  lapply(1:forest$num.trees
         , function(k){ 
           pp[,k] %>% 
             calcLogloss2( df=data ) %>% 
             unlist
         }) %>% 
    unlist -> LL
  
  # LL can be calculated with vectorize :
  #((function(k){ 
  #  pp[,k] %>% 
  #    calcLogloss2( df=data.set.val ) %>% 
  #    unlist
  #} )%>%
  #Vectorize(.))(1:rg$num.trees) -> LL
  
  return(order(LL)) # tree indices , ordered by logloss on OOB
}

# name should change , it evaluates a packet of 50 simulations, it calcs the LL and sizes of forests of selected trees
calc_LL_for_selection <- function(method, doc, parameter, metrices=NULL){ # added metrices as optinal argument 14.7.2022
  #' gets arguments and calculates arguments needed for some inner function (funL) from doc
  #' then runs loops over the inner function (funL) for all metrices and all simulations (50) in doc
  #' 
  #' @param method specifies the inner function from a list of function (funL)
  #' @param parameter are parameters for function specified by method
  #' @param doc is the documentation of 50 simulations (forest, dissimilarity matrices, in-bag and OOB observations)
  #' @param metrices is optional, if NULL then we use all metrices from the (named!) dissimilarity matrix in doc
  
  nBs <- length(doc)
  #nBs <- 5
  
  DM <- doc[[1]]$DM
  if(is.null(metrices)){
    # no metrices given : use those / all for which we have a dissimilarity matrix
    metrices <- names(DM)
  }else{
    # metrices given : check that we have a dissimilarity matrix for each given metric
    assertthat::assert_that(all( metrices %in% names(DM) ) 
                            , msg=paste('passed/required metrices'
                              , paste(metrices,collapse=',')
                              ,'do not exist in doc[[1]]$DM: '
                              , paste(names(DM), collapse=',')
                              )
                            )
  }
  
  evalT <- matrix(NA,nrow=nBs,ncol=3*length(metrices)) # table of evaluations
  
  ct <- 1
  for(i in 1:nBs){
    #print(paste(i/nBs , Sys.time()))
    if(i%%10 ==0) print(paste(i/nBs , Sys.time()))
    
    DM <- doc[[i]]$DM
    assertthat::assert_that(all( metrices %in% names(DM) ), msg='not all passed metrices found in names of dissimilarity matrices, cf DM')
    
    forest <- doc[[i]]$rg$forest
    
    #data.train <- doc[[i]]$`bootstapped training data` # refers to nursery , we should work with nursery02 and the corrected typo!
    data.train <- doc[[i]]$`bootstrapped training data` # new 11.5.2022
    OOB <-  base::setdiff(1:nrow(Cleve), unique(doc[[i]]$resample))
    data.set.val <- Cleve[OOB,] # goes back to original Cleveland data
    
    LL <- calc_LL(forest, data.set.val)
    oLL <- order(LL)
    
    res <- list() # result from following loop
    
    # parameter will be used as a named list
    # with a sub-list for each metric , 
    # like this : parameter <- list('d0'=list('cutoff'=0.4, 'sizeSF'=500), 'd1'= list('cutoff'=0.3, 'sizeSF'=500) ...))
    # this allows for individual parameters for each metric
    
    # if parameter is passed as only a single list
    # parameter = list('cutoff'=0.4, 'sizeSF'=500)
    # then this list is repeated for each metric to become of the required format
    
    if(!is.null(parameter)){ 
    if(all(c('cutoff','sizeSF') %in% names(parameter))){
      P1 <-  list()
      for(metric in metrices){
        P1[[metric]] <- parameter
      }
      parameter <- P1
    }
    }
    
    # check for required format of named list with same names as the dissimilarity matrices
    assertthat::assert_that(all(names(parameter) %in% names(DM))
                            , msg = 'names for dissimilarities in parameter do not match names in dissimilarity matrices.')
    
    #print(parameter)
    
    funL <- list('meiner'=calc_meinForest
         ,'chip1'=calc_chipForest_1
         , 'chip2'=calc_chipForest_2
        #, 'chip2'=eval_chipForest_2
         , 'chip1_simplified'=eval_chipForest_1_simplified
        )
    
    for(metric in metrices){
      argsL <-  list('meiner'=list(dm=DM[[metric]], forest=forest, LL=LL, parameter=parameter[[metric]])
                     , 'chip1'=list(dm=DM[[metric]], forest=forest, oLL=oLL, parameter=parameter[[metric]])
                     , 'chip2'=list(dm=DM[[metric]], forest=forest, oLL=oLL, parameter=parameter[[metric]])
                     , 'chip1_simplified'=list(dm=DM[[metric]] , oLL=oLL, forest=forest, parameter=parameter[[metric]]))
      res[[metric]] <- do.call(funL[[method]]
                              , argsL[[method]]
                              )
    }
    #print('dim of evaluation table')
    #print(evalT %>% dim)
    #print('length of result')
    #print(unlist(res) %>% length)
    
    evalT[i,] <- unlist(res)
  }
  
  evalT <-  data.frame(evalT)
  names(evalT) <- c(paste(rep(c(paste('LL.test.', method,'.', sep=''),'I.','size.'),times=length(metrices))
                          #,rep(c('d0','d1','d2','sb'),each=3),sep='') # ok for metrices=NULL, and metrices==names(DM)==c('d0','d1','d2','sb')
                          ,rep(metrices,each=3),sep='')
                          )
  return(evalT)
}
# returns 12 columns
