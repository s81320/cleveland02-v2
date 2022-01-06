# Functions
impCon<-function (df, conVars){
  for (i in conVars){
    naPos <- which(is.na(df[,i]))
    if(length(naPos) != 0){
      mean_df <-mean(df[,i], na.rm=T)
      sd_df <-sd(df[,i], na.rm=T)
      n_df <-length(naPos)
      set.seed(50)
      df[naPos,i]<-rnorm(n=n_df, mean=mean_df, sd=sd_df)
    }
  }
  return (df)
}
impCat<-function (df, catVars){
  for (i in catVars){
    naPos <- which(is.na(df[,i]))
    if(length(naPos) != 0){
      temp<-df[,i]
      var=unique(temp[!is.na(temp)])
      pvar=NULL
      set.seed(50)
      for (j in 1:length(var)){
        pvar[j]<-length(temp[!is.na(temp)&temp==var[j]])/length(temp)
      }
      imp<-sample(var,size=length(naPos),replace=T, prob=pvar)
      df[naPos,i]<-imp
    }
  }
  return (df)
}
renameVars<-function(df){
  df$Sex<-ifelse(df$Sex=='Female',0,1)
  df$cp<-0
  df$cp[df$Chestpaintype=='Typical angina']=1
  df$cp[df$Chestpaintype=='Atypical angina']=2
  df$cp[df$Chestpaintype=='Non-anginal pain']=3
  df$cp[df$Chestpaintype=='Asymptomatic']=4
  df$HighFastBloodSugar<-ifelse(df$HighFastBloodSugar=='No',0,1)
  df$restECG<-0
  df$restECG[df$RestingECG=='Normal']=0
  df$restECG[df$RestingECG=='ST-T wave abnorm']=1
  df$restECG[df$RestingECG=='Leftvent hypertrophy']=2
  df$ExInducedAngina<-ifelse(df$ExInducedAngina=='No',0,1)
  df$CAD<-ifelse(df$CAD=='No',0,1)
  df<-df[,c(1,2,12,4:6,13,8:11)]
  names(df)[3]<-'Chestpaintype'
  names(df)[7]<-'RestingECG'
  return(df)
}

# Process the training and validation data
library(tidymodels)
library(dplyr)
library(sjmisc)
Cleve.data <- read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data", header = F, sep = ",") %>% na_if(., "?")
Cleve.data <- Cleve.data %>% rename(Age = V1, Sex = V2, Chestpaintype = V3, RestingBP = V4, Cholsterol = V5, HighFastBloodSugar = V6, 
                                    RestingECG = V7, MaxHeartRate = V8, ExInducedAngina = V9, STDepression = V10, SlopePeakST = V11, 
                                    NumbMajorVessels = V12, Thal = V13, CAD = V14) %>% 
  mutate(Chestpaintype = recode(as.factor(Chestpaintype), "1" = "Typical angina", "2" = "Atypical angina", "3" = "Non-anginal pain", 
                                "4" = "Asymptomatic") ) %>% 
  mutate(HighFastBloodSugar = recode(as.factor(HighFastBloodSugar), "1" = "Yes", "0" = "No") ) %>% 
  mutate(RestingECG = recode(as.factor(RestingECG), "0" = "Normal", "1" = "ST-T wave abnorm", "2" = "Leftvent hypertrophy") ) %>% 
  mutate(ExInducedAngina = recode(as.factor(ExInducedAngina), "1" = "Yes", "0" = "No") ) %>% 
  mutate(SlopePeakST = recode(as.factor(SlopePeakST), "1" = "Upsloping", "1" = "Flat", "2" = "Downsloping") ) %>% 
  mutate(Thal = recode(as.factor(Thal), "6.0" = "Fixed defect", "3.0" = "Normal", "7.0" = "Reversable defect") ) %>% 
  mutate(NumbMajorVessels = as.numeric(NumbMajorVessels) ) %>%
  mutate(Sex = recode(as.factor(Sex), "1" = "Male", "0" = "Female") ) %>%
  mutate( CAD = as.factor( (case_when( CAD > 0 ~ "Yes", CAD == 0 ~ "No" ) ) )) 
Hung.data <- read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.hungarian.data", header = F, sep = ",") %>% na_if(., "?")
Hung.data <- Hung.data %>% rename(Age = V1, Sex = V2, Chestpaintype = V3, RestingBP = V4, Cholsterol = V5, HighFastBloodSugar = V6, 
                                  RestingECG = V7, MaxHeartRate = V8, ExInducedAngina = V9, STDepression = V10, SlopePeakST = V11, 
                                  NumbMajorVessels = V12, Thal = V13, CAD = V14) %>% 
  mutate(Age = as.numeric(Age)) %>%
  mutate(Chestpaintype = recode(as.factor(Chestpaintype), "1" = "Typical angina", "2" = "Atypical angina", "3" = "Non-anginal pain", 
                                "4" = "Asymptomatic") ) %>% 
  mutate(RestingBP = as.numeric(RestingBP) ) %>%
  mutate(Cholsterol = as.numeric(Cholsterol) ) %>%
  mutate(HighFastBloodSugar = recode(as.factor(HighFastBloodSugar), "1" = "Yes", "0" = "No") ) %>% 
  mutate(RestingECG = recode(as.factor(RestingECG), "0" = "Normal", "1" = "ST-T wave abnorm", "2" = "Leftvent hypertrophy") ) %>% 
  mutate(MaxHeartRate = as.numeric(MaxHeartRate) ) %>%
  mutate(ExInducedAngina = recode(as.factor(ExInducedAngina), "1" = "Yes", "0" = "No") ) %>% 
  mutate(SlopePeakST = recode(as.factor(SlopePeakST), "1" = "Upsloping", "1" = "Flat", "2" = "Downsloping") ) %>% 
  mutate(Thal = recode(as.factor(Thal), "6.0" = "Fixed defect", "3.0" = "Normal", "7.0" = "Reversable defect") ) %>% 
  mutate(NumbMajorVessels = as.numeric(NumbMajorVessels) ) %>%
  mutate(Sex = recode(as.factor(Sex), "1" = "Male", "0" = "Female") ) %>%
  mutate( CAD = as.factor( (case_when( CAD > 0 ~ "Yes", CAD == 0 ~ "No" ) ) )) 
Swiss.data <- read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.switzerland.data", header = F, sep = ",") %>% na_if(., "?")
Swiss.data <- Swiss.data %>% rename(Age = V1, Sex = V2, Chestpaintype = V3, RestingBP = V4, Cholsterol = V5, HighFastBloodSugar = V6, 
                                    RestingECG = V7, MaxHeartRate = V8, ExInducedAngina = V9, STDepression = V10, SlopePeakST = V11, 
                                    NumbMajorVessels = V12, Thal = V13, CAD = V14) %>% 
  mutate(Age = as.numeric(Age)) %>%
  mutate(Chestpaintype = recode(as.factor(Chestpaintype), "1" = "Typical angina", "2" = "Atypical angina", "3" = "Non-anginal pain", 
                                "4" = "Asymptomatic") ) %>% 
  mutate(RestingBP = as.numeric(RestingBP) ) %>%
  mutate(Cholsterol = as.numeric(Cholsterol) ) %>%
  mutate(HighFastBloodSugar = recode(as.factor(HighFastBloodSugar), "1" = "Yes", "0" = "No") ) %>% 
  mutate(RestingECG = recode(as.factor(RestingECG), "0" = "Normal", "1" = "ST-T wave abnorm", "2" = "Leftvent hypertrophy") ) %>% 
  mutate(MaxHeartRate = as.numeric(MaxHeartRate) ) %>%
  mutate(ExInducedAngina = recode(as.factor(ExInducedAngina), "1" = "Yes", "0" = "No") ) %>% 
  mutate(STDepression = as.numeric(STDepression) ) %>%
  mutate(SlopePeakST = recode(as.factor(SlopePeakST), "1" = "Upsloping", "1" = "Flat", "2" = "Downsloping") ) %>% 
  mutate(Thal = recode(as.factor(Thal), "6.0" = "Fixed defect", "3.0" = "Normal", "7.0" = "Reversable defect") ) %>% 
  mutate(NumbMajorVessels = as.numeric(NumbMajorVessels) ) %>%
  mutate(Sex = recode(as.factor(Sex), "1" = "Male", "0" = "Female") ) %>%
  mutate( CAD = as.factor( (case_when( CAD > 0 ~ "Yes", CAD == 0 ~ "No" ) ) )) 
VA.data <- read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.va.data", header = F, sep = ",") %>% na_if(., "?")
VA.data <- VA.data %>% rename(Age = V1, Sex = V2, Chestpaintype = V3, RestingBP = V4, Cholsterol = V5, HighFastBloodSugar = V6, 
                              RestingECG = V7, MaxHeartRate = V8, ExInducedAngina = V9, STDepression = V10, SlopePeakST = V11, 
                              NumbMajorVessels = V12, Thal = V13, CAD = V14) %>% 
  mutate(Age = as.numeric(Age)) %>%
  mutate(Chestpaintype = recode(as.factor(Chestpaintype), "1" = "Typical angina", "2" = "Atypical angina", "3" = "Non-anginal pain", 
                                "4" = "Asymptomatic") ) %>% 
  mutate(RestingBP = as.numeric(RestingBP) ) %>%
  mutate(Cholsterol = as.numeric(Cholsterol) ) %>%
  mutate(HighFastBloodSugar = recode(as.factor(HighFastBloodSugar), "1" = "Yes", "0" = "No") ) %>% 
  mutate(RestingECG = recode(as.factor(RestingECG), "0" = "Normal", "1" = "ST-T wave abnorm", "2" = "Leftvent hypertrophy") ) %>% 
  mutate(MaxHeartRate = as.numeric(MaxHeartRate) ) %>%
  mutate(ExInducedAngina = recode(as.factor(ExInducedAngina), "1" = "Yes", "0" = "No") ) %>% 
  mutate(STDepression = as.numeric(STDepression) ) %>%
  mutate(SlopePeakST = recode(as.factor(SlopePeakST), "1" = "Upsloping", "1" = "Flat", "2" = "Downsloping") ) %>% 
  mutate(Thal = recode(as.factor(Thal), "6.0" = "Fixed defect", "3.0" = "Normal", "7.0" = "Reversable defect") ) %>% 
  mutate(NumbMajorVessels = as.numeric(NumbMajorVessels) ) %>%
  mutate(Sex = recode(as.factor(Sex), "1" = "Male", "0" = "Female") ) %>%
  mutate( CAD = as.factor( (case_when( CAD > 0 ~ "Yes", CAD == 0 ~ "No" ) ) ))  
## Exclude variables with too much missingness
library(jmv)
descriptives(Cleve.data, vars=vars('Sex','Chestpaintype','HighFastBloodSugar','RestingECG','ExInducedAngina','SlopePeakST','NumbMajorVessels','Thal'), freq=T)
descriptives(Cleve.data, vars=vars('Age','RestingBP','Cholsterol','MaxHeartRate','STDepression'))
descriptives(Hung.data, vars=vars('Sex','Chestpaintype','HighFastBloodSugar','RestingECG','ExInducedAngina','SlopePeakST','NumbMajorVessels','Thal'), freq=T)
descriptives(Hung.data, vars=vars('Age','RestingBP','Cholsterol','MaxHeartRate','STDepression'))
descriptives(Swiss.data, vars=vars('Sex','Chestpaintype','HighFastBloodSugar','RestingECG','ExInducedAngina','SlopePeakST','NumbMajorVessels','Thal'), freq=T)
descriptives(Swiss.data, vars=vars('Age','RestingBP','Cholsterol','MaxHeartRate','STDepression'))
descriptives(VA.data, vars=vars('Sex','Chestpaintype','HighFastBloodSugar','RestingECG','ExInducedAngina','SlopePeakST','NumbMajorVessels','Thal'), freq=T)
descriptives(VA.data, vars=vars('Age','RestingBP','Cholsterol','MaxHeartRate','STDepression'))
Cleve.data<-subset(Cleve.data, select=-c(SlopePeakST, Thal, NumbMajorVessels))
Hung.data<-subset(Hung.data, select=-c(SlopePeakST, Thal, NumbMajorVessels))
Swiss.data<-subset(Swiss.data, select=-c(SlopePeakST, Thal, NumbMajorVessels))
VA.data<-subset(VA.data, select=-c(SlopePeakST, Thal, NumbMajorVessels))
## Impute missing values
continuousVars <- c (1,4,5,8,10)
Cleve.data<-impCon(Cleve.data,continuousVars)
Hung.data<-impCon(Hung.data, continuousVars)
Swiss.data<-impCon(Swiss.data, continuousVars)
VA.data<-impCon(VA.data, continuousVars)
categoricalVars <- c(2,3,6,7,9)
Cleve.data<-impCat(Cleve.data, categoricalVars)
Hung.data<-impCat(Hung.data, categoricalVars)
Swiss.data<-impCat(Swiss.data, categoricalVars)
VA.data<-impCat(VA.data, categoricalVars)
## Rename variables
Cleve<-renameVars(Cleve.data)
Hung<-renameVars(Hung.data)
Swiss<-renameVars(Swiss.data)
VA<-renameVars(VA.data)
## No variables need to be removed due to high correlation
cor.cleve <- cor(Cleve, use = "complete.obs", method = "spearman")
which((cor.cleve> 0.7 | cor.cleve < -0.7) & cor.cleve != 1, arr.ind = TRUE)
cor.hung <- cor(Hung, use = "complete.obs", method = "spearman")
which((cor.hung> 0.7 | cor.hung < -0.7) & cor.hung != 1, arr.ind = TRUE)
cor.swiss <- cor(Swiss, use = "complete.obs", method = "spearman")
which((cor.swiss> 0.7 | cor.swiss < -0.7) & cor.swiss != 1, arr.ind = TRUE)
cor.va <- cor(VA, use = "complete.obs", method = "spearman")
