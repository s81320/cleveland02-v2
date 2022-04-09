
# Beispiel mice
library(mice)
library(dplyr)

df <-  data.frame('a'=c(1,2,3,4,2,3,4,3,4) 
                , 'b'=c(1,2,3,4,5,6,7,8,9)
                , 'c'=as.factor(c(rep('Y',4), rep('N',5)))) # factor is needed, does not work with strings

df[3,1]
df[3,1] <-  NA

df[5,1]
df[5,1] <-  NA

df[2,3]
df[2,3] <-  NA

md.pattern(df)
dfImp <- mice(data=df, method=c('cart','','sample'))
dfImp$where
data.frame(complete(dfImp,1))

##################

df <- cbind(rnorm(50), rnorm(50,1,1))
df <-  cbind(df,-2*df[,1]+2*df[,2]+rnorm(50,0,1))
df[2,1] <- NA
df[5,1] <-  NA
df[2,3] <- NA
md.pattern(df)
is.na(df) %>% table
dfImp <- mice(data=df, method=c('pmm','','pmm'))
dfImp$where
data.frame(complete(dfImp,1)) %>% is.na %>% table

###

load('data/data_SupMat4.rda')
Swiss$STDepression %>% summary
corSwiss <- Swiss[,1:11]
corSwiss$STDepression <- Swiss$STDepression - min(Swiss$STDepression)
corSwiss$STDepression %>% summary

corSwiss$Cholesterol <-  NA
df <- rbind(corSwiss , VA[,1:11])

md.pattern(df)
dfImp <- mice(data=df, method='pmm')
data.frame(complete(dfImp,1)) -> dfc

corSwiss %>%
  c(VA[,1:11]) %>%
  mice(method='pmm') %>%
  complete(1) %>%
  .[1:nrow(corSwiss),] -> corSwiss2

dfc %>% is.na %>% table

corSwiss <- dfc[1:nrow(corSwiss),]
md.pattern(corSwiss)

# does it improve performance for the forest?
# check look_into_data_Cleve_Swiss_Hung.R
