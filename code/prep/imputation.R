rm(list=ls())

library(mice)

# my own imputation for Cholesterol

# read data but do not impute:
# rum Cleve_preprocessing.R up to line 144 (about)

# this creates Cleve.data , Hung.data , Swiss.data and VA.data
# all should have 11 variables now , 10 of them predictive

# look at the descriptive statistics in Cleve_preprocessing.R

# look at explicitly missing data : NA resulting from ?
is.na(Cleve.data) %>% table
is.na(Hung.data) %>% table
is.na(Swiss.data) %>% table
is.na(VA.data) %>% table

# we add implicit missing data : Cholesterol = 0 , RestingBP = 0
(Hung.data$Cholesterol==0) %>% table
(Swiss.data$Cholesterol==0) %>% table
(VA.data$Cholesterol==0) %>% table
Swiss.data[!is.na(Swiss.data$Cholesterol) & Swiss.data$Cholesterol==0,'Cholesterol'] <- NA
VA.data[!is.na(VA.data$Cholesterol) & VA.data$Cholesterol==0,'Cholesterol'] <- NA

(VA.data$RestingBP==0) %>% table
VA.data[!is.na(VA.data$RestingBP) & VA.data$RestingBP==0,'RestingBP'] <- NA

VA.data.imp <- VA.data
VA.data.imp %>% mice(method='cart') %>%
  complete(1) -> VA.data.imp

# check : all should be FALSE
is.na(VA.data.imp) %>% table


Swiss.data.imp <- Swiss.data
Swiss.data.imp %>% mice(method='cart') %>% # oder 'rf' (with default ranger rf) oder 'pmm'
  complete(1) -> Swiss.data.imp

# check : all should be FALSE
is.na(Swiss.data.imp) %>% table 
# fewer is.na , but Cholesterol cannot be imputed if there are only NAs
# so we use VA to help

Swiss.data.imp %>% 
  rbind(VA.data.imp) %>% 
  as.data.frame %>%
  mice(method='cart') %>% 
  complete(1) %>%
  .[1:nrow(Swiss.data.imp),] -> Swiss.data.imp

# check : all should be FALSE
is.na(Swiss.data.imp) %>% table 

Hung.data.imp <-  Hung.data
# check : all should be FALSE
is.na(Hung.data.imp) %>% table

Hung.data.imp %>% mice(method='cart') %>%
  complete(1) -> Hung.data.imp
                       
# check : all should be FALSE
is.na(Hung.data.imp) %>% table

# no more missing values.
# hidden missing values for Cholesterol 0 are also removed

# run rename var in Cleve_preprocessing.R in lines 160 ff
# and do the same for my imputed guys:

Hung.data.imp<-renameVars(Hung.data.imp)
Swiss.data.imp<-renameVars(Swiss.data.imp)
VA.data.imp<-renameVars(VA.data.imp)

