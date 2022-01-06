# 1) sb test
# 2) rt2fbt test
# 5.1.2022

# 1) sb test
# 2 full binary trees with 7 split nodes, depth 3. The terminal nodes / leaves are not encoded, they would never be used.
sb<-function(tri1,tri2){
  if(length(tri1)!=length(tri2)){
    # append 0 to the shorter one, until of same length
    if(length(tri1)<length(tri2)){
      tri1 <-c(tri1,c(rep(0,length(tri2)-length(tri1)))) 
    }else{
      tri2<-c(tri2,c(rep(0,length(tri1)-length(tri2)))) 
    }
  }
  length(which(tri1!=tri2))/length(tri1) # since both trees now have the same length, we can divide by the length of any of them
  # the scaling factor (the division) may be different for each pair of trees, it is not uniform over the trees in the forest
}

# compare results with hand drawn trees
# test case : 2 trees of same depth
tri1 <-  c(1
           ,2,3
           ,0,0,2,0)
tri2 <- c(2
          ,1,3
          ,2,0,0,0)

sb(tri1,tri2)
which(tri1!=tri2)
length(which(tri1!=tri2))

# test case : 2 trees of different depth
tri1 <-  c(1
           ,2,3
           ,0,0,2,0)
tri2 <- c(2
          ,1,3
          ,2,0,0,0
          ,1,2,0,0,0,0,0,0)
sb(tri1,tri2)

# count differences up to common depth
which(tri1!=tri2[1:7]) # 7 is length(tri1)
length(which(tri1!=tri2[1:7]))

# count differences for lower layer of deeper tree to all 0 entries
which(c(0,0,0,0,0,0,0,0)==tri2[8:15])
length(which(c(0,0,0,0,0,0,0,0)==tri2[8:15]))
# same as count non-zero entries in lower layers of deeper tree
length(which(tri2[8:15]!=0))
# add it up
length(which(tri1!=tri2[1:7])) + length(which(c(0,0,0,0,0,0,0,0)!=tri2[8:15]))
# or
length(which(tri1!=tri2[1:7])) + length(which(tri2[8:15]!=0))
6/length(tri2)

# test case : 2 trees of different depth, more than 1 layer difference
tri1 <-  c(1,2,3,0,0,2,0)
tri2 <- c(2 # 1 
          ,1,3 # 2
          ,2,rep(0,3) # 4
          ,1,2,rep(0,6) # 8
          ,2,3,rep(0,14)) # 16
sb(tri1,tri2)

# we should not divide by the number of nodes of the larger full binary tree
# but by the number of nodes a combined tree would have (whenever (at least) one tree has a split node, count this node in)