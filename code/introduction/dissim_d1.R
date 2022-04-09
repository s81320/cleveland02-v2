# needed in the introduction of d1


rg <- doc[[1]]$rg

rm(list='E1')
E1 <- list()
doc[[1]]$DM$d1 <- createDMd1(rg$forest,  Cleve[unique(rs),])
doc[[1]]$DM$d1 %>% dim
# E1 is now a list with for each tree a list of terminal nodes (key) and the observations mapped to it (values)

# if we want to check that observations fall to the correct terminal node
# we have to go back to the original data and the numbering there:
unique(rs)[E1[[1]][[10]]] # observations in original Cleveland data set

tNodes(rg$forest,1) %>% length
E1[[1]][[10]] # first terminal node
E1[[1]][[10]][[1]] # first observation in first terminal node
# find the terminal node in 2nd tree

# 
findMe <- function(obs,tri){
  #' find an observation (by number) in a tree (by number), returns the node number of the terminal node in the 0 based ranger tree (!= full binary tree!)
  tn <-  tNodes(rg$forest,tri)
  for(i in tn){
    if(length(which(E1[[tri]][[i]]==obs))>0) return(i)
  }
}

findMe(obs=1,tri=2)

f2 <- function(tri1, tri2, tnode){
  #' tnode should be a terminal node in tree tri1
  lapply(E1[[tri1]][[tnode]], FUN=findMe, tri=tri2) %>% 
    unlist %>% 
    table %>%
    unname %>% 
    lapply( function(y) y*(y-1)/2) %>% 
    unlist %>%
    sum
}

# needed in thesis, introduction of d1 dissim
f2(tri1=1,tri2=1,tnode=10) # test : look in tree 1 for obs that are in tnode 10 of tree 1
f2(tri1=1, tri2=2,tnode=10)  # tree 2 and observations from terminal node 10 : how many pairs fall together to the same terminal node in tree 2

# more
lapply(tNodes(rg$forest, 2) , FUN=f2, tri1=2, tri2=2) %>% unlist # how can there be a 0 ?? 0 pairs = 1 obs
lapply(tNodes(rg$forest, 1) , FUN=f2, tri1=1, tri2=1) %>% unlist # how can there be a 0 ??
lapply(tNodes(rg$forest, 1) , FUN=f2, tri1=1, tri2=2) %>% unlist
