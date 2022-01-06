# helper functions: 
## acc : accuracy
## apsf : accuracy for the predictions of a subforest (a p sf)
## pred.sf : generating a prediction by ensemble vote
## pred.sf.w : weighted version of pred.sf
## checkPartition

## test cases at end of file

acc <- function (vec1, vec2=df[test,'CAD']) {
  #' acc : accuracy
  #' 
  #' @param vec1 no default
  #' @param vec2 defaults to df[test,'CAD'] - what if df or test are not defined?
  if(length(vec1)==length(vec2))
  {return(length(which(vec1==vec2))/length(vec1))}
  else {return(-1)}
}

apsf<- function(forest=rg$forest, idx, newdata){
  #' apsf : accuracy predictions sub-forest
  #' 
  #' calculating accuracy for a subforest on given data
  #' 
  #' @param forest a ranger forest with treetype 'Classification' or 'Probability estimation'
  #' @param idx set of indices to select trees from a forest, thus selecting a sub-forest
  #' @param newdata with target variable 'CAD' (not generic!)
  #' @return in range 0 to 1 of all's well, -1 if something went wrong in acc (helper function to calculate accuracy)

  idx %>%
    subforest(forest,.) %>%
    predict(data=newdata) %>%
    .$predictions ->
    pp
  if(forest$treetype=='Probability estimation'){
    pp<-ifelse(pp[,1]>0.5,1,2) %>% as.factor
    levels(pp)<-forest$levels
  }
  acc(pp,newdata$CAD)
  
}


pred.sf <- function(mtrx,levels=rf$classes){
  # ensemble vote : return the class / name of the family most often returned
  # mtrx should be a matrix with row: observation to be voted on, column: trees that cast their vote. Matrix elements are the votes
  # when 2 classes get the maximum votes (tie) the vote goes always to the one first in line
  mf1<-function(vec) table(factor(vec, levels=levels))
  mf2<-function(vec) attr(which.max(vec),'names')
  
  x<-apply(mtrx,1, mf1) # turn each row into a table, counting how often each class is the classification result , counting votes for each class
  x<-apply(x,2,mf2) # get the maximum vote
  x<-factor(x, levels=rf$classes) # relevant when some families are not in the classification results, prevents the table to be reduced to a smaller size than 4x4.
  
  return(x)
}

pred.sf.w <- function(mtrx , wts=NULL , levels=rf$classes){
  # weighted ensemble vote : count the votes for each class, multiply the number of votes by the corresponding weight, return the class / name of the family with the highest value
  # mtrx should be a matrix with row: observation to be voted on, column: trees that cast their vote. Matrix elements are the votes
  # when 2 classes get the maximum votes (tie) the vote goes always to the one first in line
  
  if (is.null(wts)) return(pred.sf(mtrx,levels))
  else{
    if (length(wts)!=ncol(mtrx)) print('error in number of weights rel. to ncol for matrix mrtx')
    
    nC<- length(levels)
    nT<-ncol(mtrx)
    #print(paste('classes, trees',nC,nT))
    
    dimnames<-list()
    dimnames[1]<-list(1:nrow(mtrx))
    dimnames[2]<-list(levels)
    #dimnames<-list(dimnames)
    x<-matrix(0,nrow(mtrx),nC, dimnames=dimnames)
    # this x has a different layout than in predict.sf without weights
    
    # this should be done with apply
    for (r in 1:nrow(mtrx)){
      for (i in 1:nT){
        # x+=w
        x[r,as.character(mtrx[r,i])]<-x[r,as.character(mtrx[r,i])]+wts[i] 
      }
    }
    
    #print('for loop done')
    mf2<-function(vec) attr(which.max(vec),'names')
    x<-apply(x,1,mf2) # get the maximum vote
    # we apply mf2 row wise, in predict.sf we did it columns wise
    return(x) 
  }
}

checkPartition<-function(train_=train, val_=val, test_=test, df_=df){ 
  
  OK<-FALSE
  
  diff.in.strat <- (as.vector(prop.table(table(df_[train_,]$CAD))) +
                      as.vector(prop.table(table(df_[val_,]$CAD))) +
                      as.vector(prop.table(table(df_[test_,]$CAD))) - 
                      3*as.vector(prop.table(table(df_$CAD)))) %>%
    abs() %>%
    max()
  
  if(
    # check that the sizes of train, val and test set add to the size of the data set (number of rows)
    length(union(union(test_,val_),train_))==nrow(df_) &
    
    # check disjointness of train, val, test sets
    length(base::intersect(test_,val_))+length(base::intersect(train_,val_))+length(base::intersect(train_,test_))==0 & 
    
    # rough check that the class balance is about equal in all parts and in the overall data set , probably too rough ...
    diff.in.strat<0.05) { OK<-TRUE }
  
  return(OK)
}

# accuracy
acc.test<- function(){
print(acc(c(1,2,3,4),c(1,2,3,1)) == 3/4) # 3/4
print(acc(c(1,2,3,4))==-1) # -1 
}

checkPartition.test <- function(){
# checkPartition() # requires values for the defaults... and should return TRUE
# checkPartition(c(1,2,3)) # should return FALSE
# the following should return TRUE , any messing with the arguments should return FALSE
  checkPartition( c(1,2,3), c(4,5), c(6,7), as.data.frame(matrix(0,nrow=7,ncol=2)) )
}


