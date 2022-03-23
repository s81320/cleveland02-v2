# dissimilarities between trees

# d0, d1, d2 : based on Banerjee 2012 paper
# identifying representative trees from ensembles
# https://deepblue.lib.umich.edu/bitstream/handle/2027.42/92082/sim4492.pdf;jsessionid=1BD56DC77B9AC9DA2B39632B5459528A?sequence=1

# Shannon-Banks-metric : introduced by Shannon and Banks, 1999
# https://pubmed.ncbi.nlm.nih.gov/10204200/
# cited in Chipman, George, McCulloch , 1998, available at
# http://www-stat.wharton.upenn.edu/~edgeorge/Research_papers/ , formula (5) on page 4

# d1, d2 metric use predictions and the terminal nodes observations fall into (d1) or the predicted class (d2). 
# Does this work with PE-forests?

# sb metric : when comparing 2 trees wrt their architecture
# we have to compare them as full binary trees
# The ranger forest has a shorter notation that does not allow for direct comparison
# trees will be of different depths since we removed the max.depth=5 argument.
# Don't know yet haw / where this is relevant...

# needed for d1
source("code/source/subforest.R") # loads function tNodes to get terminal nodes of a tree

#### split metric ####
createDMd0 <- function(forest){
  luc<-mapply(function(a,b) unique(a[which(b[[1]]!=0)]) # get split vars for inner nodes
              , forest$split.varIDs
              , forest$child.nodeIDs
              , SIMPLIFY=F # with the default SIMPLIFY=T we may sometimes get a matrix 
              # (if all lengths are equal or even a vector if all lengths are equally equal to 1)
              # this way we always get a list. May not really be efficient , sometimes ...
              # ... anytime we could have gotten a matrix :-)
  )
  d0<-function(a,b){length(base::union(a,b)) - length(base::intersect(a,b))}
  
  nVar<-length(forest$independent.variable.names)
  
  return(outer(luc,luc,Vectorize(d0))/nVar)
}

#### partition metric ####
createDMd1 <- function(forest, dft){
    #' partition metric d1
    #' 
    #' @description This function builds a dissimilarity matrix of the d1 partition metric 
    #' for a given forest and set of observations.
    #' 
    #' Each tree naturally partitions the observations by placing them in its terminal nodes. 
    #' Two trees are similar if they partition the observations similarly, putting pairs of observations in the same node or not.
    #' Two trees treating many pairs of observations differently 
    #' (one tree placing them in the same node the other tree in separate nodes) are less similar.
    #' 
    #' @param forest a ranger forest object
    #' @param dft the training data frame used to build the forest. 
    #' But it can be any set at all as long as the forest can handle it.
    #' 
    #' @return a matrix of d1 dissimilarities for the trees in the forest. matrix entries are in [0,1].
  
  nT<-forest$num.trees
  
  if (is.null(dft)){stop('data frame of training data missing. Please, use argument dft=')}
  
  nObs<-nrow(dft)
  
  pp<-predict(forest, # predict works with a ranger object or a ranger.forest object
              data=dft, # data.frame training
              #predict.all=TRUE, # 21.8.21
              type='terminalNodes'
  )$predictions
  
  A <- array(0, dim=c(nObs, nObs, nT)) # 0's are needed on the diagonal
  
  #docTN <- list() # documenting the terminal nodes and their observations
  for(tri in 1:nT){
    # loop through IDs of terminal nodes
    for(tn in tNodes(forest,tri)){
      sameTN <- which(pp[,tri]==tn) # observations mapped to this terminal node ID
      #docTN[[tn]] <- sameTN
      for(i in sameTN){
        triangleNodes <- base::intersect(sameTN , 1:i) # don't run through sameTN x sameTN , only through its triangle
        #print(paste(length(sameTN) , length(triangleNodes)))
        for(j in triangleNodes){
          A[i,j,tri]<-1  
          A[j,i,tri]<-1
        }
        # old loop , i,j looping over square of sameTN x sameTN
        #for(j in sameTN){
        #  A[i,j,tri]<-1  # making A[,,tri] a symmetric matrix , since looping over all combinations (i,j) in sameTN x sameTN
        #}
      }
    }
    #E1[[tri]] <<- docTN
  }
  #A[,,1]
  
  d1<-function(tri1,tri2){
    B<- A[,,tri1]-A[,,tri2]
    B[upper.tri(B)] %>% abs() %>% sum() # counting how many pairs of obs are treated differently by tri1 and tri2
  } # a function
  return(outer(1:nT,1:nT,Vectorize(d1)) * 2 / (nObs*(nObs-1))) # a matrix , dim  nT , nT
  # divide by n over 2 = divide by n*(n-1)/2 = multiply by 2 / (n*(n-1)) 
  # with n = number of observations nObs , number of rows in data set , here: the training data
}

#### fit or prediction metric ####
createDMd2 <- function(forest, dft){
  #' create a dissimilarity matrix for the prediction metric , d2
  #' 
  #' @param forest a ranger random forest of type probability estimation or classification
  #' @param dft training data (or could be any other data the trees can handle, ie that can fall down the tree, contains all variables needed at split nodes)
  
  pp<-predict(forest, # predict works with a ranger object or a ranger.forest object
              data=dft, # data.frame training
              predict.all = TRUE,
              # type='response' # it's the default
  )$predictions
  
  if(forest$treetype=='Classification'){
  #  assertthat::assert_that(length(dim(pp))==2, msg = paste('classification tree should return 2-dim predictions. got' , length(dim(pp))))
    # print('Classification tree, we got class predictions. do nothing')
    f1<-function(vec,vec2) mean(abs(vec-vec2))
  }else{
    if(forest$treetype=='Probability estimation'){
      pp <- pp[,2,] # only keep the probability for CAD==Yes
      # assertthat::assert_that(all(dim(pp)==c(nrow(dft),forest$num.trees)),msg='dimensions wrong in createDMd2')
      f1 <- function(vec,vec2) mean(abs(vec-vec2)) # mean abs error for the probabilities for CAD ==Yes
    }else{
        print(paste('got unknown tree type in forest : ' , forest$treetype))
      }
    }
  
  f2<-function(i,j) f1(pp[,i],pp[,j]) # trees i,j
  nT<-forest$num.trees
  return( outer(1:nT,1:nT,Vectorize(f2)) )
  
}

#### Shannon banks metric , tree metric ####
rt2fbt<-function(rgf=rg$forest, trindx=1){ 
  #' rt2fbt : ranger tree to full binary tree
  #' 
  #' @param rgf : forest object of a ranger object (a random forest)
  #' @param trindx : tree index , only 1 index, not a list of indices!
  
  t1<-rgf$child.nodeIDs[[trindx]] # a fixed tree
  ld<-t1[[1]]  # left children / daughters
  rd<-t1[[2]]  # right
  
  mf3<-function(a=0,b=1,k=0) {
    #' mf3 : my function 3 :-(
    #' 
    #' a ranger forest tree is coded with left and right daughters (2 lists) for each (non terminal) node.
    #' and an extra list for the split variables used at the inner nodes.
    #' a full binary tree will be coded as a sequence (single list) with 
    #' each list position referring to a fixed inner node in the layout of a full binary tree.
    #' the list value is the split variable applied at the inner node.
    #' 
    #' @param a a=0 is the root for a nodeID in a ranger forest-object
    #' @param b b=1 is the root for a full binary tree
    #' @param k distance to the root for current a,b
    #' 
    #' @return tree / subtree as a single list
    
    #### vector z created at last split nodes (no further split nodes, next are leaves) ####
    if(ld[a+1]==0){
      if(rd[a+1]==0){
        z<-c(rep(0,2**k-1)) # for this last split node (no further splitting, next is a leaf) the tree has (at least) 2**k-1 inner nodes. Local knowledge
        #print(paste(k, ' created z of length ' , length(z)))
        # z<-c(rep(0,2**(k+1)-1)) # unnecessary to create at the level of leaves / terminal nodes, there are no split variables to document! 
        # k is distance from the root , b in 2**k:(2**(k+1)-1)
        #print(paste('create z at node',a,b))
        return(z)}
      else{
        message(paste('error in mf3 at level ',k)) # error if only one child exists
        return(-1)
      }
    }
    #### at split nodes ####
    else{ #### recursive function calls for left and right sub-trees ####
      left.st<-mf3(ld[a+1], 2*b, k+1)
      right.st<-mf3(rd[a+1], 2*b+1, k+1)
      
      #print('merging')
      #print(z1)
      #print(z2)
      # adding z1 + z2 when they may be of unequal length
      l1<-length(left.st)
      l2<-length(right.st)
      if(l1!=l2){
        # append 0 to the shorter one, until of same length
        if(l1<l2){
          left.st <-c(left.st,c(rep(0,l2-l1))) }
        else{
          right.st<-c(right.st,c(rep(0,l1-l2))) }
      }
    }
    
    z<-left.st+right.st
    
    # split vars start at 0 so we have to +1 , z is initialized with 0 which is not in split.varIDs
    z[b]<-1+rgf$split.varIDs[[trindx]][a+1]
    #print(paste(b, z[b], rgf$split.varIDs[[trindx]][a+1]))
    
    #print('merged to')
    #print(z)
    return(z)
  }
  
  return(mf3())
}

createDMsb <- function(forest){
  
  lapply(1:forest$num.trees, function(i) rt2fbt(forest,i) ) ->
    A
  
  sb<-function(tri1,tri2){
    if(length(tri1)!=length(tri2)){
    # append 0 to the shorter one, until of same length
      if(length(tri1)<length(tri2)){
        tri1 <-c(tri1,c(rep(0,length(tri2)-length(tri1)))) 
      }else{
          tri2<-c(tri2,c(rep(0,length(tri1)-length(tri2)))) 
      }
    }
    # length(which(tri1!=tri2))/length(tri1) # since both trees now have the same length, we can divide by the length of any of them
    # the scaling factor (the division) may be different for each pair of trees, it is not uniform over the trees in the forest
    length(which(tri1!=tri2)) # do no scaling. If scaling has to be done it should be wrt the number of positions in the full binary tree that are non null when adding the 2 full binary trees 
    # removed the division 26.10.2021 (in Chipman paper there is no scaling factor)
    # possible scaling : length(which((tri1+tri2)!=0))
  }
  outer(A,A,Vectorize(sb)) %>% 
    as.matrix -> 
    result
  
  0 -> diag(result) 

  return(result)
}

#####################################################################

createDM <- function(forest, type, dft=NULL){
  # types d1 , d2 require the data.frame of the training set : dft
  if(type=='d0'){ createDMd0(forest) }
  else{ 
    if(type=='d1'){ 
      if(is.null(dft)){
        message('training data needed, please use argument dft=')}
      createDMd1(forest, dft) }
    else{
      if(type=='d2'){
        if(is.null(dft)){
          message('training data needed, please use argument dft=')}
        createDMd2(forest, dft) }
      else{
        if(type=='sb'){ createDMsb(forest) }
        else{ message(paste('undefined type ', as.character(type))) }
      } # end of else : when neither d0, d1 nor d2
      } # end of else : when neither d0 nor d1
    } # end of else : when not d0
  }

