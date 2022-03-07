# cluster only if hypothesis of unimodality is rejected

# select central tree per cluster
lapply(1:kOpt,
       function(i){
         x <-  which(hc.clus==i) # trees of cluster i
         lapply(1:length(x), 
                function(j){sum(dm[x[j],x])} %>% # Summe der Abstände zu allen Elementen des eigenen clusters
                  unlist ) %>% 
         which.min %>% # kleinste summe der Abstände
         x[.] # zugehöriges Element in x
         }
       ) %>% unlist -> R

# best performing tree per cluster
lapply(1:kOpt,function(i) which(hc.clus==i)[1]) %>% # we take the first, because indices are ordered, first is smallest
  unlist -> R