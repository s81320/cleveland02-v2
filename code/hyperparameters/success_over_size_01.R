# visualize data generated in hyperparameters / ranger_num_trees-v3-nursery.R
# and from code / ChipmanSelection / 03_hpo.R

library(dplyr)

load('data/LL_nested_regular_forests_size_2_to_500.rda') # loads LL

folder <- 'data/chipman'
files <- list.files(folder) ; files
file <- files[[3]] ; file
load(paste(folder,file,sep='/'))
et <- et03


# trade of for sucess and size for the regular random forest
apply(LL[1:100,],2,mean) -> LL.m # mean values only
par(mar=c(4,4,0,0)+0.3)
plot(colnames(LL)[1:9] %>% as.integer 
     , LL.m[1:9] 
     , type='l'
     , xlab='size'
     , ylab='logloss'
     #, main='trade off for regular random forests'
     )

plot(colnames(LL)[9:18] %>% as.integer 
     , LL.m[9:18] 
     , type='l'
     , xlab='size'
     , ylab='logloss'
     #, main='trade off for regular random forests'
     )

plot(colnames(LL)[18:58] %>% as.integer 
     , LL.m[18:58] 
     , type='l'
     , xlab='size'
     , ylab='logloss'
     #, main='trade off for regular random forests'
     )

# new methods in the context of regular random forests
matrix(nrow=0, ncol=5) %>%
  data.frame %>%
  setNames(c('logloss','I','size','dissim','cutoff')) -> et2

metrices <- c('d0','d1','d2','sb')
for(metric in metrices){
  et %>% 
    select(ends_with(metric) | 'cutoff') %>%
    cbind('dissim'=metric) %>%
    setNames(c('logloss','I','size','cutoff','dissim')) %>%
    rbind(et2) -> et2
}

par(mar=c(4,4,2,1)+0.3)
apply(LL,2,mean) -> LL.m # mean values only
plot(colnames(LL) %>% as.integer 
     , LL.m 
     , type='l'
     , xlab='size'
     , ylab='logloss'
    # , main='Chipman 2 forests in context of regular random forests'
     )
points(et2$size , et2$logloss, col=as.factor(et2$dissim))
legend('topright'
       # , legend=c('Chipman 1, d0','Chipman 1, d1','Chipman 1, d2','Chipman 1, sb','regular RF')
       , legend=c('Meiner, d0','Meiner, d1','Meiner, d2','Meiner, sb','regular RF')
       , pch=c(rep('o',4),'-')
       , col=c(1:4,1)
       , cex=0.8
)

par(mar=c(4,4,1,1)+0.01)
ylim <- range(0.54,0.60)
LL[,which((as.integer(colnames(LL)) > 9 ) & (as.integer(colnames(LL)) < 101 ) )] %>%
  apply(2,mean) %>%
  plot(x=names(.), y=. , 
       ylim=ylim
       , type='l' 
       , xlab='size'
       , ylab='logloss'
       #, main='success over size for Chipman 2 forests\nin context of regular random forests'
       )
abline(mean(LL[,'500']),0,col='grey')

text(et2[(et2$size > 9)&(et2$size<101),'size']
     , et2[(et2$size >9)&(et2$size<101),'logloss' ]
     , labels = et2[(et2$size >9)&(et2$size<101),'cutoff' ] 
     , cex=0.8
     , col=as.integer(as.factor(et2[(et2$size >9)&(et2$size<101),'dissim' ])))

legend('topright'
       #, legend=c('Chipman 2, d0','Chipman 2, d1','Chipman 2, d2','Chipman 2, sb','regular RF','default forest, size 500','cutoff parameter')
       , legend=c('Meiner, d0','Meiner, d1','Meiner, d2','Meiner, sb','regular RF','default forest, size 500','cutoff parameter')
       , pch=c(rep('o',4),'-','-','x')
       , col=c(1:4,1,'grey',1)
       , cex=0.8
)


LL[,which((as.integer(colnames(LL))>3) & (as.integer(colnames(LL))<11 ) )] %>%
  apply(2,mean) %>%
  plot(x=names(.), y=. , 
       ylim=range(0.5,1.2),  type='l' , xlab='size', ylab='logloss')
points(et2[(et2$size>3)&(et2$size<11),'size'] 
       , et2[(et2$size>3)&(et2$size<11),'logloss' ]
       , col=as.factor(et2[(et2$size>3)&(et2$size<11),'dissim' ]))
legend('topright'
      # , legend=c('regular RF, mean','Chipman, d0','Chipman, d1','Chipman, d2','Chipman, sb')
      , legend=c('regular RF, mean','Meiner, d0','Meiner, d1','Meiner, d2','Meiner, sb')
       , pch=c('-',rep('o',4))
       , col=c(1,1:4)
       , cex=0.8
)

# plot evolution of std dev
par(mar=c(2.5,2.5,0.5,0.5)+0.2) # axes labels : no || title / main : no
par(mar=c(2,2,3,1)+0.1) # axes labels : no ||  main / title : yes
par(mar=c(4,4,3,1)+0.1) # shows main / title , shows axes labels
plot(colnames(xt)[1:9]
     , xt['sd',1:9]
     , type='l'
     #, main='success over size, regular forests\nstandard deviation of logloss (mean over 1000 simulations)'
     #, xlab='size: number of trees in forest'
     , xlab='size'
     # , ylab='success: mean standard deviation'
     , ylab='standard deviation of logloss'
     #, cex.lab=2
     , cex.axis=2)
plot(colnames(xt)[9:18]
     , xt['sd',9:18]
     , type='l'
     #, main='success over size, regular forests\nstandard deviation of logloss (mean over 1000 simulations)'
     #, xlab='size: number of trees in forest'
     , xlab='size'
     #, ylab='success: mean standard deviation'
     , ylab='standard deviation of logloss'
     , cex.axis=2)
plot(colnames(xt)[18:ncol(xt)]
     , xt['sd',18:ncol(xt)]
     , type='l'
     #      , main='success over size, regular forests\nstandard deviation of logloss (mean over 1000 simulations)'
     # , xlab='size: number of trees in forest'
     , xlab='size'
     # , ylab='success: standard deviation of logloss'
     , ylab='standard deviation of logloss'
     , cex.axis=2)

