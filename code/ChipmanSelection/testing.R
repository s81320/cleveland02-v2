# chipman.R test

# load arguments (dm , forest, oLL, parameter) to feed the function
{load('data/nursery/nursery01_10_50x500.rda')
metric <- 'd1'
DM <- doc[[1]]$DM
dm <-  doc[[1]]$DM[[metric]]

data.train <- doc[[1]]$`bootstapped training data`
OOB <-  base::setdiff(1:nrow(Cleve), unique(doc[[1]]$resample))
data.set.val <- Cleve[OOB,] # goes back to original Cleveland data

forest<-doc[[1]]$rg$forest

pp <- predict(forest 
              , data=data.set.val 
              , predict.all = T)$predictions[,2,]

lapply(1:forest$num.trees
       , function(k){ 
         pp[,k] %>% 
           calcLogloss2( df=data.set.val ) %>% 
           unlist
       }) %>% 
  unlist -> LL

oLL <- order(LL) # tree indices , ordered by logloss on OOB

parameter <-  list('cutoff'=0.5 , 'sizeSF'=500)
pa <- parameter
}

## call function
calc_chipForest_1(dm , forest, oLL, pa)

mean(dm[1,])
summary(dm[1,])
quantile(dm[1,],0.5)
quantile(dm,0.5)

# 1)
dm.test <- matrix(0.4+rnorm(25,0,0.05),5,5)
# 2)
dm.test <- matrix(rnorm(25,0,1),5,5)
dm.test <- 3+dm.test
all(dm.test>0) # check
#
diag(dm.test) <- 0
dm.test[1,2] <- 0.1
dm.test[2,1] <- 0.1
dm.test[3,4] <- 0.1
dm.test[4,3] <- 0.1
dm2 <-  dm.test
LL.test <- c(0.5,0.7,0.45,0.8,0.6)
LL <-  LL.test

dip.test(c(-0.1,-0.2,0,2,1.5,1.6,1.7)) # very multimodal
dip.test(c(-4,-4.1,-4.21,-4.11,-0.1,-0.2,0,2,1.5,1.6,1.7)) # very multimodal
# construct as multimodal
data <- c(rnorm(5,-4,0.5) , rnorm(7,0,0.5) , rnorm(9,4,0.5))
plot(density(data))
plot(data)
dip.test(data) # recognise as multimodal, small p.value
#### overlapping , unimodlity rejected
data <- c(rnorm(50,-2,0.5) , rnorm(70,0,0.5) , rnorm(90,2,0.5))
plot(density(data))
plot(data)
dip.test(data)
#### more overlapping , unimodality not rejected
data <- c(rnorm(50,-1.7,0.5) , rnorm(70,0,0.5) , rnorm(90,1.7,0.5))
plot(density(data))
plot(data)
dip.test(data)

#### more overlapping and more (!) data (!), unimodality rejected
data <- c(rnorm(150,-1.7,0.5) , rnorm(150,0,0.5) , rnorm(150,1.7,0.5))
plot(density(data))
plot(data)
dip.test(data)
