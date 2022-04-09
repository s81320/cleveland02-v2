# new dissimilarities between trees

# modify d0: weight the split variables with weight 1/2^h for h the depth of the node
# or 5/2^h to make the difference larger for: never a split variable vs only a split variable at a low level (high depth)
# then we need a way to measure distances for the vector of weighted split variables
# mse sum((v1-v2)^2) ? mae sum(|v1-v2|)? KL divergence?
# mse or mae but with weights given by the variable importance of the full forest.

# modify d0: weight by the number of observations that go through a split node:
# how to calculate the number of observations that go through a split node:
# 1) find all terminal nodes below the split node.
# 2) do the prediction of type 'terminalNode' and count through terminal nodes found in 1)

# modify d0: Jaccard index , it is a similarity , not a dissimilarity .. ??!!
# does 1- similarity make it a dissimilarity?

# modify d1 : d1* of Banerjee : d1 * d0

# scale the dissimilarities to be in the interval 0,1 (put the min at 0, the max at 1) or to have mean 0, std 1
# Since we use quantiles to set the levels of representation or diversity this should not matter.

# modify d2 : use mse instead of mae as distance between predictions