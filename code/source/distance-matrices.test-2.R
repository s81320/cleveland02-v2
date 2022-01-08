# test for calcDMd0 , calcDMsb
# test calculation of dissimilarity matrices for d0 and sb against plotted trees
# compare visually

source('code/source/distance-matrices.R')

# grow a small forest
rg <- ranger(Species~. , num.trees=3 , data=iris , max.depth = 5)

# test that treeInfo and the plotted tree agree
# maybe draw tree by hand from treeInfo to compare

treeInfo(rg,1)
plotTree1(rg, 1, main='tree 1')

treeInfo(rg,2)
plotTree1(rg, 2, main='tree 2')

treeInfo(rg,2)
plotTree1(rg, 3, main='tree 3')

# use return values of plotTree1 of used split variables to calculate d0 by hand
createDMd0(rg$forest)

#########################################################

createDMsb(rg$forest)

# the following function is used in createDMsb
# full binary trees built from ranger trees
lapply(1:rg$forest$num.trees, function(i) rt2fbt(rg$forest,i) ) ->
  A
A

#################

# draw tree by hand from treeInfo , compare to plot
# plotTreeFb transforms the ranger tree to a full binary tree as used in calculation of the sb dissimilarity matrix.
# so if tree is drawn correctly, a step in the calculation of the DM is succussfully tested.
# the remaining part of calcDMsb is rather simple...

treeInfo(rg,1)
plotTreeFb(rg, 1, main='tree 1')

treeInfo(rg,2)
plotTreeFb(rg, 2, main='tree 2')

treeInfo(rg,3)
plotTreeFb(rg, 3, main='tree 3')

createDMd0(rg$forest)
createDMsb(rg$forest)
