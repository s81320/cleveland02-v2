# load Cleveland
rm(list=ls())
load('data/data_SupMat4.rda')
source('code/source/distance-matrices.R')

read_folder <- 'data/nursery'
write_folder <- 'data/nursery02'
files <- list.files(read_folder)[2:3]

repair <- function(doc,info){
  for(i in 1:length(doc)){
    print(Sys.time())
    names(doc[[i]])[names(doc[[i]])=='bootstapped training data'] <- 'bootstrapped training data'
    rg <- doc[[i]]$rg
    rs <- doc[[i]]$resample
    doc[[i]]$DM$d1 <- createDMd1(rg$forest, Cleve[unique(rs),])
  }
  info$update <- paste('calculated DM$d1 anew on ', Sys.time() ,'. script code/nursery/update_d1.R')
  save(info, doc , file=paste(write_folder, file, sep='/'))
}

for(file in files){
  load(paste(read_folder, file, sep='/'))
  repair(doc,info)
}


# save doc , info to write_folder
