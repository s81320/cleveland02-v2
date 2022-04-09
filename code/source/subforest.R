
print(paste('sourcing s u b f o r e s t .R : functions subforest(forest , trindx) , forestHull(forest) and tNodes(forest , tri)'))

# constructor
subforest<-function(forest,trindx){
  sf<- forest # copy? point? # keep levels, class.values, treetype, indep.variable.names
  # all that is in names(forest) and not changes in the following
  # what is in names(forest) may depend on the treetype. Whether it is classification, regression or pe 
  sf$child.nodeIDs<-forest$child.nodeIDs[trindx]
  sf$split.varIDs<-forest$split.varIDs[trindx]
  sf$split.values<-forest$split.values[trindx]
  sf$num.trees<-length(trindx)
  if('terminal.class.counts' %in% names(forest)){
    sf$terminal.class.counts<-forest$terminal.class.counts[trindx]
  }
  class(sf)<-unique(c(class(forest),'subforest')) # one might do multiple sub-forests , and we do not want the class names to multiply
  return(sf)
}

# constructor
forestHull<- function(forest){
  # returns a ranger object with forest as its forest
  
  # sometimes we need a ranger object to pass as an argument though only forest is worked with
  # we create a minimal hull for a forest that is of class ranger
  
  hull<-list('forest'=list(),'num.trees'=0L)
  class(hull)<-c('ranger','forestHull')
  hull$forest<-forest
  hull$num.trees <- forest$num.trees
  hull$treetype="Classification" # added 21.8.21
  return(hull)
}

tNodes <- function(forest, tri){
  #' terminal nodes for th tri-th tree in a ranger forest
  #' 
  #' terminal nodes are those where child nodes (left and right) are 0
  #' 
  which(forest$child.nodeIDs[[tri]][[1]]==0)-1
}

