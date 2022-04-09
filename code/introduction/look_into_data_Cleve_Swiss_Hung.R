library(ranger)

load('data/data_SupMat4.rda')

folder <- 'data/original'
files <- list.files(folder)
files %>% stringr::str_subset('.data') -> files # in alphabetical order!

# get original data
# run Cleve_preprocessing up to line 144 to get Cleve.data , Hung.data , ... before imputation

data.set.name <- 'Cleve.data'
ds <- get(data.set.name)
is.na(ds$STDepression) %>% table
(ds$RestingBP==0) %>% which

VA.data$RestingBP %>% hist

#### ST depression

VA.data$STDepression %>% table
(VA.data$STDepression==-0.5) %>% which
par(mar=c(2,4,2,1)+0.2)
VA.data$STDepression %>% hist(main='VA, STDepression',breaks=20, xlim=c(-3,6))
is.na(VA.data$STDepression) %>% table
(VA.data$STDepression==0) %>% table #%>% prop.table
(VA.data$STDepression>0) %>% table

# compare to other data sets
Cleve.data$STDepression %>% hist(main='Cleve, STDepression', breaks=20, xlim=c(-3,6))
(Cleve.data$STDepression==0) %>% table #%>% prop.table
(Cleve.data$STDepression>0) %>% table

Hung.data$STDepression %>% hist(main='Hung, STDepression', breaks=20, xlim=c(-3,6))
(Hung.data$STDepression==0) %>% table #%>% prop.table
(Hung.data$STDepression>0) %>% table

Swiss.data$STDepression %>% hist(main='Swiss, STDepression', breaks=20, xlim=c(-3,6))
(Swiss.data$STDepression==0) %>% table #%>% prop.table
(Swiss.data$STDepression>0) %>% table

# defects of VA
VA.data$STDepression %>% summary # probably due to imputing with normally distributed variables
VA.data$STDepression %>% (function(x) c(mean(x, na.rm=T), sd(x, na.rm=T)))
VA.data$STDepression %>% hist(breaks=20, main='VA STDepression, no imputation')

VA$STDepression %>% summary # probably due to imputing with normally distributed variables
VA$STDepression %>% (function(x) c(mean(x), sd(x)))
VA$STDepression %>% hist(breaks=20, main='VA STDepression, normal imputation')

# defects of Swiss : in Swiss STDepression there are negative values , in Cleve minimum is 0
# defects of Swiss: Cholesterol is always 0 (hidden missingness)

# move STDepression to have positive values only / undo centering
corSwiss1 <- Swiss[,1:11]
corSwiss1$STDepression %>% range
corSwiss1$STDepression <- Swiss$STDepression - min(Swiss$STDepression)

# VA has 49 observations with 0 Cholesterol. This has to be settled first.
# start with original VA, before any imputation

corVA <- VA
corVA$Cholesterol[corVA$Cholesterol==0] <- NA

corVA %>% mice(method='cart') %>%
  complete(1) -> corVA

# impute Cholesterol with help of VA
corSwiss2 <- corSwiss1
corSwiss2$Cholesterol <- NA
corSwiss2 %>% 
  rbind(VA[,1:11]) %>% # lengths
  as.data.frame %>%          
  mice(method='cart') %>% # used to be pmm , but I do not want to explain a new method. rf is available as method did not run well...
  complete(1) %>%
  .[1:nrow(corSwiss2),] -> corSwiss2

corSwiss <- corSwiss2

################################################################################
#### Verschiedenheit der Datensätze Cleveland , Swiss , Hungary
################################################################################

rg <-   rg <- ranger(CAD~.
                     , data = Cleve[,1:11]
                     , num.trees = 500
                     , replace = F 
                     , mtry= 3 # default : 3
                     , importance = 'impurity'
                     , probability = T 
                     , min.node.size = 13 # optimized in A.Z. paper: BiomJ 56, 2014, 
)
irg <- importance(rg) # 3 most important predictors : Chest pain type, max heart rate , ST depression
irg[order(irg)]

data_set <- c('Cleve','Hung','Hung.data.imp','Swiss' ,'Swiss.data.imp','VA','VA.data.imp')
# Cleve for imputation by AZ team , Cleve.data.imp my own imputation

# for qualitative predictor
g0 <- function(df,v) get(df)[,v] %>% table %>% prop.table %>% round(2)
# for quantitative predictor
f0 <- function(df,v) get(df)[,v] %>% (function(x) c(mean(x, na.rm=T),sd(x, na.rm=T)))
h0 <- function(df,v) get(df)[,v] %>% range(na.rm=T)

# compare 1 dim distributions of predictors
# most important : Chestpaintype
g1 <- function(df) g0(df,'Chestpaintype')
Vectorize(g1)(data_set) %>% t

# max heart rate
f1 <- function(df) f0(df,'MaxHeartRate') 
r1 <- Vectorize(f1)(data_set)%>% t
colnames(r1) <- c('mean','sd')
r1

# STDepression
f1 <- function(df) f0(df,'STDepression') 
r1 <- Vectorize(f1)(data_set)%>% t
colnames(r1) <- c('mean','sd')
r1

h1 <- function(df) h0(df,'STDepression')
r1 <- Vectorize(h1)(data_set) %>% t
colnames(r1) <-  c('min','max')
r1 #%>% xtable

((function(x) length(which(x==0))/length(x))%>% Vectorize)(data_set)

# Age
f1 <- function(df) f0(df,'Age') 
r1 <- Vectorize(f1)(data_set)%>% t
colnames(r1) <- c('mean','sd')
r1

# Cholesterol :  Missing values in Swiss
f1 <- function(df) f0(df,'Cholesterol') 
r1 <- Vectorize(f1)(data_set)%>% t
colnames(r1) <- c('mean','sd')
r1

# min values
data_set %>%
  ((function(x) min(get(x)$Cholesterol) ) %>%
  Vectorize)

# Sex # not important , ranked 8th out of 10
g1 <- function(df) g0(df,'Sex')
Vectorize(g1)(data_set) %>% t %>% xtable

# CAD
g1 <- function(df) g0(df,'CAD')
Vectorize(g1)(data_set) %>% t %>% xtable

###############################################################################
#### variable importance for simple rf (classification, simplest) on the different sets
###############################################################################
#### is a dissimilarity for the data sets #####################################
###############################################################################

f1 <-  function(df){
  rg <-   rg <- ranger(CAD~.
                     , data = df
                     , num.trees = 500
                     , replace = T
                     , mtry= 3 # default : 3
                     , importance = 'impurity'
                    , probability = T 
                     #, min.node.size = 13 # optimized in A.Z. paper: BiomJ 56, 2014, 
  )
  importance(rg) %>% (function(x) x/sum(abs(x))) # make the importance a measure on the predictor variables
  # importance(rg) # no scaling. Absolute value of a variable's importance is important!
}


doc.irg <-  list()
#data_set <-  c('Cleve','Swiss','corSwiss1', 'corSwiss2', 'VA', 'Hung')
for(ds in data_set){
doc.irg[[ds]] <-  f1(get(ds)[,1:11])
}

irg <- bind_rows(doc.irg) %>% data.frame
rownames(irg) <-  data_set
View(irg)

irg[, order(irg['Cleve',])]

d3L <-  function(i,j){
  I <- as.matrix(irg)[data_set[[i]],] 
  J <- as.matrix(irg)[data_set[[j]],] 
  mean((I-J)^2)
}

DMd3 <- outer(1:length(data_set), 1:length(data_set), Vectorize(d3L))
DMd3 <- data.frame(DMd3)
rownames(DMd3) <- data_set
colnames(DMd3) <- data_set
1000*DMd3 %>% round(5)

digits(xt) <- 3
DMd3

###############################################################################
#### prediction on unseen data ################################################
###############################################################################
#### is a dissimilarity for the data sets #####################################
###############################################################################

f1 <-  function(data.train, data.test){
  dtrain <-  get(data.train)
  if(ncol(dtrain)>11){
    dtrain <- dtrain[,1:11]
  }
  dtest <-  get(data.test)
  if(ncol(dtest)>11){
    dtest <- dtest[,1:11]
  }
  rg <- ranger(CAD~.
                       , data = dtrain
                       , num.trees = 500
                       , replace = T
                       , mtry= 3 # default : 3
                       , importance = 'impurity'
                       , probability = T 
                       , min.node.size = 13 # optimized in A.Z. paper: BiomJ 56, 2014, 
  )

  predict(rg, data=dtest)$predictions[,2] %>% calcLogloss2(dtest)
}

f1('Cleve', 'Swiss')
f2 <- function(dfname) f1('Cleve',dfname)
Vectorize(f2)(data_set)

