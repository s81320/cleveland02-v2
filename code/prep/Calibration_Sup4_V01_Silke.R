#Supplementary Material 4 for "A comprehensive comparison of approaches for the calibration of probability machines"
#by Francisco M. Ojeda *, Yao Hu *, Alexandre Thiery, Stefan Blankenberg,Christian Weimar, Matthias Schmid, and Andreas Ziegler

rm(list=ls())

# Attaching R packages ---------------------------------------------------------
library(dplyr)
library(tidymodels)

# Functions to perform imputation of missing values ----------------------------
# ==============================================================================
impCon <- function (df, conVars){
  for (i in conVars){
    naPos <- which(is.na(df[,i]))
    if(length(naPos) != 0){
      mean_df <- mean(df[,i], na.rm = TRUE)
      sd_df <- sd(df[,i], na.rm=T)
      n_df <- length(naPos)
      set.seed(50)
      df[naPos,i] <- rnorm(n=n_df, mean = mean_df, sd = sd_df)
    }
  }
  return (df)
}
impCat <- function (df, catVars){
  for (i in catVars){
    naPos <- which(is.na(df[,i]))
    if(length(naPos) != 0){
      temp <- df[,i]
      var <- unique(temp[!is.na(temp)])
      pvar <- NULL
      set.seed(50)
      for (j in 1:length(var)){
        pvar[j] <- length(temp[!is.na(temp) & temp == var[j]]) / length(temp)
      }
      imp <- sample(var, size=length(naPos), replace = TRUE, prob = pvar)
      df[naPos, i] <- imp
    }
  }
  return (df)
}

# Function that creates, recodes and reorders variables in Cleveland data
#===============================================================================
renameVars <- function(df){
  df$Sex <- ifelse(df$Sex == 'Female',0,1)
  df$cp <- 0
  df$cp[df$Chestpaintype == 'Typical angina'] <- 1
  df$cp[df$Chestpaintype == 'Atypical angina'] <- 2
  df$cp[df$Chestpaintype == 'Non-anginal pain'] <- 3
  df$cp[df$Chestpaintype == 'Asymptomatic'] <- 4
  df$HighFastBloodSugar<-ifelse(df$HighFastBloodSugar == 'No', 0, 1)
  df$restECG <- 0
  df$restECG[df$RestingECG == 'Normal'] <- 0
  df$restECG[df$RestingECG == 'ST-T wave abnorm'] <- 1
  df$restECG[df$RestingECG == 'Leftvent hypertrophy'] <- 2
  df$ExInducedAngina<-ifelse(df$ExInducedAngina == 'No', 0, 1)
  df$CAD <- ifelse(df$CAD == 'No', 0, 1)
  df <- df[,c(1, 2, 12, 4:6, 13, 8:11)]
  names(df)[3] <- 'Chestpaintype'
  names(df)[7] <- 'RestingECG'
  return(df)
}

# Function to read and process Cleveland data ---------------------
# ==============================================================================
read_process_df <- function(url) {
  col.names <- c("Age", "Sex", "Chestpaintype", "RestingBP", "Cholesterol",
                 "HighFastBloodSugar", "RestingECG", "MaxHeartRate", 
                 "ExInducedAngina", "STDepression", "SlopePeakST", 
                 "NumbMajorVessels", "Thal", "CAD")
  df <- read.table(url, header = FALSE, sep = ",", col.names = col.names, 
                   na.strings = "?") %>% 
    mutate(Chestpaintype = recode(as.factor(Chestpaintype), 
                                  "1" = "Typical angina", 
                                  "2" = "Atypical angina", 
                                  "3" = "Non-anginal pain", 
                                  "4" = "Asymptomatic") ) %>% 
    mutate(HighFastBloodSugar = recode(as.factor(HighFastBloodSugar), 
                                       "1" = "Yes", "0" = "No") ) %>% 
    mutate(RestingECG = recode(as.factor(RestingECG), "0" = "Normal", 
                               "1" = "ST-T wave abnorm", 
                               "2" = "Leftvent hypertrophy") ) %>% 
    mutate(ExInducedAngina = recode(as.factor(ExInducedAngina), 
                                    "1" = "Yes", "0" = "No") ) %>% 
    mutate(SlopePeakST = recode(as.factor(SlopePeakST), 
                                "1" = "Upsloping", "2" = "Flat", 
                                "3" = "Downsloping") ) %>% 
    mutate(Thal = recode(as.factor(Thal), "6" = "Fixed defect", "3" = "Normal", 
                         "7" = "Reversable defect") ) %>% 
    mutate(NumbMajorVessels = as.numeric(NumbMajorVessels) ) %>%
    mutate(Sex = recode(as.factor(Sex), "1" = "Male", "0" = "Female") ) %>%
    mutate( CAD = as.factor( (case_when( CAD > 0 ~ "Yes", CAD == 0 ~ "No" ) ) )) 
  # Variables with too many missing values
  df <- subset(df, select=-c(SlopePeakST, Thal, NumbMajorVessels))
  # Impute missing value
  continuousVars <- c (1, 4, 5, 8, 10)
  df<-impCon(df, continuousVars)
  categoricalVars <- c(2, 3, 6, 7, 9)
  df <- impCat(df, categoricalVars)
  # creating, recoding and reordering variables
  df <- renameVars(df)
  # Adding factor version of CAD
  df[["CAD_fac"]] <- factor(df[["CAD"]])
  df
}

# Reading and processing datasets ----------------------------------------
# ==============================================================================
url_root <- paste0("https://archive.ics.uci.edu/ml/machine-learning-databases/",
                   "heart-disease/processed.")
url_df <- data.frame(df_name = c("Cleve", "Hung", "Swiss", "VA"),
                     url = c(paste0(url_root, c("cleveland", "hungarian", 
                                                "switzerland", "va"), ".data")), 
                     stringsAsFactors = FALSE)
for(i in 1:nrow(url_df)) {
  assign(url_df[["df_name"]][i], read_process_df(url_df[["url"]][i]))
}

#### Silke : copied the mutate statement to mutate to factors
# removed previously removed features
mutate2factors <-  function(df){
  df %>% mutate(Chestpaintype = recode(as.factor(Chestpaintype), 
                                       "1" = "Typical angina", 
                                       "2" = "Atypical angina", 
                                       "3" = "Non-anginal pain", 
                                       "4" = "Asymptomatic") ) %>% 
    mutate(HighFastBloodSugar = recode(as.factor(HighFastBloodSugar), 
                                       "1" = "Yes", "0" = "No") ) %>% 
    mutate(RestingECG = recode(as.factor(RestingECG), "0" = "Normal", 
                               "1" = "ST-T wave abnorm", 
                               "2" = "Leftvent hypertrophy") ) %>% 
    mutate(ExInducedAngina = recode(as.factor(ExInducedAngina), 
                                    "1" = "Yes", "0" = "No") ) %>% 
    mutate(Sex = recode(as.factor(Sex), "1" = "Male", "0" = "Female") ) %>%
    mutate( CAD = as.factor( (case_when( CAD > 0 ~ "Yes", CAD == 0 ~ "No" ) ) )) 
}

Cleve <- mutate2factors(Cleve)
Swiss <- mutate2factors(Swiss)
Hung <-  mutate2factors(Hung)
VA <- mutate2factors(VA)

save(Cleve, Hung, Swiss,VA, file='data_SupMat4.rda')

# Fitting machines on the Cleveland data ---------------------------------------
# ==============================================================================
# Train and test data 
df_train <- Cleve
for(i in 2:nrow(url_df)) {
  study_name <- url_df[["df_name"]][i]
  df <- get(study_name)
  df[["study"]] <- study_name
  assign(study_name, df); rm(df, study_name)
}
df_test <- rbind(Hung, Swiss, VA)
# Creating formula
Formula <- paste0("CAD_fac ~ ", paste(setdiff(names(df_train), 
                                              c("CAD", "CAD_fac")), 
                                      collapse = "+")) %>% as.formula
# Creating folds for cross-validation
set.seed(50)
folds <- vfold_cv(df_train, v = 2, repeats = 2)

# Elastic net ------------------------------------------------------------------
# ==============================================================================
grid_size <- 50
set.seed(50)
# specifying engine and parameters to tune
spec_en <- logistic_reg(penalty = tune(),
                        mixture = tune()) %>% 
  set_engine("glmnet", family = "binomial") %>% 
  set_mode("classification")
# specifying grid for hyperparameter tuning
grid_en <- dials::grid_latin_hypercube(
  dials::penalty(), mixture(),
  size = grid_size)
# Creating workflow
wf_en <- workflows::workflow() %>%
  add_model(spec_en) %>% 
  add_formula(Formula)
# tuning model
en_tuned <- tune::tune_grid(
  object = wf_en,
  resamples = folds,
  grid = grid_en,
  metrics = metric_set(mn_log_loss),
  control = tune::control_grid(verbose = TRUE, 
                               allow_par = FALSE))
best_log_loss_en <- select_best(en_tuned, "mn_log_loss")
tuned_final_en <- finalize_workflow(wf_en, best_log_loss_en) 
fit_en <- tuned_final_en  %>% fit(data = df_train)
df_test[["prob_en"]] <- predict(fit_en, new_data = df_test, 
                                type = "prob")[[".pred_1"]]

# Random forest ----------------------------------------------------------------
# ==============================================================================
set.seed(50)
spec_rf <- rand_forest(mtry = tune(), trees = 500,
                       min_n = tune()) %>%
  set_engine("ranger") %>%
  set_mode("classification") %>% set_args(seed = NULL)
# random forest grid. Adding the parameters used in the original publication to the grid
vars <- all.vars(Formula)[-1] # these are the predictors
grid_rf <- dials::grid_latin_hypercube(
  finalize(dials::mtry(), 
           x = df_train[, vars, drop = FALSE]),
  dials::min_n(),
  size = grid_size
)
grid_rf_aux <- grid_rf[1, ]
grid_rf_aux[["mtry"]] <- floor(sqrt(ncol(df_train[, vars, 
                                                  drop = FALSE])))
grid_rf_aux[["min_n"]] <- floor(nrow(df_train) * 0.05)
grid_rf <- rbind(grid_rf, grid_rf_aux)
grid_rf <- grid_rf[!duplicated(grid_rf), ]
rm(grid_rf_aux)
# Creating workflow
wf_rf <- workflows::workflow() %>%
  add_model(spec_rf) %>% 
  add_formula(Formula)
# tuning model
rf_tuned <- tune::tune_grid(
  object = wf_rf,
  resamples = folds,
  grid = grid_rf,
  metrics = metric_set(mn_log_loss),
  control = tune::control_grid(verbose = TRUE, 
                               allow_par = FALSE))
best_log_loss_rf <- select_best(rf_tuned, "mn_log_loss")
tuned_final_rf <- finalize_workflow(wf_rf, best_log_loss_rf)
fit_rf <- tuned_final_rf  %>% fit(data = df_train)
df_test[["prob_rf"]] <- predict(fit_rf, new_data = df_test, 
                                type = "prob")[[".pred_1"]]

# Support vector machines ------------------------------------------------------
# ==============================================================================
# specifying engine and parameters to tune
set.seed(50)
spec_svm <- svm_rbf(cost = tune(),
                    rbf_sigma = tune()) %>% 
  set_engine('kernlab', type = "C-svc", kernel = "rbfdot",prob.model = TRUE) %>% 
  set_mode("classification")
# specifying grid for hyperparameter tuning
grid_svm <- dials::grid_latin_hypercube(
  dials::cost(), rbf_sigma(),
  size = grid_size)
# Creating workflow
wf_svm <- workflows::workflow() %>%
  add_model(spec_svm) %>% 
  add_formula(Formula)
# tuning model
svm_tuned <- tune::tune_grid(
  object = wf_svm,
  resamples = folds,
  grid = grid_svm,
  metrics = metric_set(mn_log_loss),
  control = tune::control_grid(verbose = TRUE, 
                               allow_par = FALSE))
best_log_loss_svm <- select_best(svm_tuned, "mn_log_loss")
tuned_final_svm <- finalize_workflow(wf_svm, best_log_loss_svm) 
fit_svm <- tuned_final_svm  %>% fit(data = df_train)
df_test[["prob_svm"]] <- predict(fit_svm, new_data = df_test, 
                                 type = "prob")[[".pred_1"]]

# Gradient boost ---------------------------------------------------------------
# ==============================================================================
grid_size <- 50
set.seed(50)
# specifying engine and parameters to tune
spec_gb <- boost_tree(learn_rate = tune(),
                      tree_depth = tune()) %>% 
  set_engine("xgboost") %>% 
  set_mode("classification")
# specifying grid for hyperparameter tuning
grid_gb <- dials::grid_latin_hypercube(
  dials::learn_rate(),tree_depth(),
  size = grid_size)
# Creating workflow
wf_gb <- workflows::workflow() %>%
  add_model(spec_gb) %>% 
  add_formula(Formula)
# tuning model
gb_tuned <- tune::tune_grid(
  object = wf_gb,
  resamples = folds,
  grid = grid_gb,
  metrics = metric_set(mn_log_loss),
  control = tune::control_grid(verbose = TRUE, 
                               allow_par = FALSE))
best_log_loss_gb <- select_best(gb_tuned, "mn_log_loss")
tuned_final_gb <- finalize_workflow(wf_gb, best_log_loss_gb) 
fit_gb <- tuned_final_gb  %>% fit(data = df_train)
df_test[["prob_gb"]] <- predict(fit_gb, new_data = df_test, 
                                type = "prob")[[".pred_1"]]

# Logistic regression ----------------------------------------------------------
# ==============================================================================
stepPvalue <- function(df, p.lim){
  fit <- glm(as.factor(CAD)~., data = df, family = binomial)
  coef <- summary(fit)$coefficients
  while(any(coef[-1,4] >= p.lim)){
    fit <- glm(as.factor(CAD)~., data = df, family = binomial)
    coef <- summary(fit)$coefficients
    removeVar <- which(max(coef[-1,4]) == coef[,4])
    df <- df[-which(names(df) == names(removeVar))]
  }
  returnList <- list(fit = fit, df = df, coef = coef)
  return(returnList)
}
library(mfp)
mfp(CAD ~ Age + RestingBP + Cholsterol + MaxHeartRate + STDepression, Cleve)
lr.cleve<-stepPvalue(Cleve[,-12],0.001)
df_test[["prob_lr"]] <- plogis(predict(lr.cleve[["fit"]], newdata=df_test))

# Functions for calibration ----------------------------------------------------
# ==============================================================================
# Functions to deal with 0 or 1 predictions
winsorize_probs <- function(x) {
  eps <- .Machine$double.eps
  pmax(pmin(x, 1 - eps), eps)
}
# Functions for calibration
cal_elkan <- function(response, traindata, caldata, newdata, prob) {
  br0 <- mean(traindata[, response])
  br1 <- mean(caldata[, response])
  prob.new <- newdata[, prob] %>% winsorize_probs()
  br1 * (prob.new - prob.new * br0) / (br0 - prob.new * br0 + br1 * prob.new - br0 * br1)
}
cal_beta <- function(response, traindata, caldata, newdata, prob) {
  prob1 <- caldata[, prob]%>% winsorize_probs()
  prob2 <- newdata[, prob]%>% winsorize_probs()
  caldata[["lp1"]] <- log(prob1)
  caldata[["lp2"]] <- -log(1 - prob1)
  newdata[["lp1"]] <- log(prob2)
  newdata[["lp2"]] <- -log(1 - prob2)
  Formula <- paste0(response, " ~ lp1 + lp2") %>% as.formula
  fit <- glm(Formula, data = caldata, family = binomial())
  predict(fit, newdata = newdata, type = "response")
} 
cal_iso <- function(response, traindata, caldata, newdata, prob) {
  fit <- isoreg(x = caldata[, prob], y = caldata[, response])
  stfun <- as.stepfun(fit)
  stfun(newdata[, prob])
}
cal_logistic <- function(response, traindata, caldata, newdata, prob) {
  prob1 <- caldata[, prob]%>% winsorize_probs()
  prob2 <- newdata[, prob]%>% winsorize_probs()
  caldata[["lp"]] <- log(prob1 / (1 - prob1))
  newdata[["lp"]] <- log(prob2 / (1 - prob2))
  Formula <- paste0(response, " ~ lp") %>% as.formula
  fit <- glm(Formula, data = caldata, family = binomial())
  predict(fit, newdata = newdata, type = "response")
}
cal_loginter <- function(response, traindata, caldata, newdata, prob) {
  prob1 <- caldata[, prob]%>% winsorize_probs()
  prob2 <- newdata[, prob]%>% winsorize_probs()
  caldata[["lp"]] <- log(prob1 / (1 - prob1))
  newdata[["lp"]] <- log(prob2 / (1 - prob2))
  Formula <- paste0(response, " ~ offset(lp)") %>% as.formula
  fit <- glm(Formula, data = caldata, family = binomial())
  predict(fit, newdata = newdata, type = "response")
}
cal_logslope <- function(response, traindata, caldata, newdata, prob) {
  prob1 <- caldata[, prob]%>% winsorize_probs()
  prob2 <- newdata[, prob]%>% winsorize_probs()
  caldata[["lp"]] <- log(prob1 / (1 - prob1))
  newdata[["lp"]] <- log(prob2 / (1 - prob2))
  Formula <- paste0(response, " ~ lp -1 ") %>% as.formula
  fit <- glm(Formula, data = caldata, family = binomial())
  predict(fit, newdata = newdata, type = "response")
}
cal_platt <- function(response, traindata, caldata, newdata, prob) {
  Formula <- paste0(response, " ~ ", prob) %>% as.formula
  fit <- glm(Formula, data = caldata, family = binomial())
  predict(fit, newdata = newdata, type = "response")
}
cal_chen <- function(response, traindata, caldata, newdata, prob) {
  br2 <- mean(caldata[, response])
  prob1 <- caldata[, prob]%>% winsorize_probs()
  prob2 <- newdata[, prob]%>% winsorize_probs()
  caldata[["lp"]] <- qlogis(prob1)
  newdata[["lp"]] <- qlogis(prob2)  
  for(current_value in 0:1) {
    v <- which(caldata[[response]] == current_value)
    if(length(v) == 0) stop("Response only takes one value")
    fit <- density(caldata[["lp"]][v], bw = "nrd0")
    assign(paste0("dens_fun_", current_value), approxfun(fit, rule = 2))
    rm(fit, v)
  }
  lp2 <- newdata[["lp"]]
  f1 <- dens_fun_1(lp2)
  f0 <- dens_fun_0(lp2)
  numerator <- f1 * br2 
  denominator <- (numerator + f0 * (1 - br2))
  numerator / denominator
}

# Creating folds for 10-fold cross-validation
# ==============================================================================
set.seed(50)
fold_hung<-vfold_cv(Hung, v=10, repeats = 1, strata='CAD')
Hung$cv<-0
for (i in 1:10){
  Hung$cv[-fold_hung$splits[[i]]$in_id]=i
}
fold_swiss<-vfold_cv(Swiss, v=10, repeats = 1, strata='CAD')
Swiss$cv<-0
for (i in 1:10){
  Swiss$cv[-fold_swiss$splits[[i]]$in_id]=i
}
fold_va<-vfold_cv(VA, v=10, repeats = 1, strata='CAD')
VA$cv<-0
for (i in 1:10){
  VA$cv[-fold_va$splits[[i]]$in_id]=i
}
df_test$cv<-c(Hung$cv, Swiss$cv, VA$cv)
df_test$id<-1:nrow(df_test)

# Calibration
# ==============================================================================
machine<-c('lr','rf','svm','en','gb')
cal.out<-df_test
set.seed(50)
for (j in machine){
  out<-NULL
  for (i in 1:10){
    cal<-cal.out[cal.out$cv!=i,]
    new<-cal.out[cal.out$cv==i,]
    new[[paste0('prob_',j,'_beta')]]<-cal_beta ('CAD', df_train, cal, new, paste0('prob_',j))
    new[[paste0('prob_',j,'_elkan')]]<-cal_elkan ('CAD', df_train, cal, new, paste0('prob_',j))
    new[[paste0('prob_',j,'_iso')]]<-cal_iso ('CAD', df_train, cal, new, paste0('prob_',j))
    new[[paste0('prob_',j,'_log')]]<-cal_logistic ('CAD', df_train, cal, new, paste0('prob_',j))
    new[[paste0('prob_',j,'_int')]]<-cal_loginter ('CAD', df_train, cal, new, paste0('prob_',j))
    new[[paste0('prob_',j,'_slope')]]<-cal_logslope ('CAD', df_train, cal, new, paste0('prob_',j))
    new[[paste0('prob_',j,'_platt')]]<-cal_platt ('CAD', df_train, cal, new, paste0('prob_',j))
    new[[paste0('prob_',j,'_chen')]]<-cal_chen('CAD', df_train, cal, new, paste0('prob_',j))
    out<-rbind(out,new)
  }
  out<-out[,20:28]
  df_test<-merge(df_test, out, by='id')
}
# Refit
refit<-df_test[,c(1:20)]
refit.out<-NULL
for (i in 1:10){
  cal<-refit[refit$cv!=i,]
  new<-refit[refit$cv==i,]
  set.seed(50)
  folds <- vfold_cv(cal, v = 2, repeats = 2)
  grid_size <- 50
  set.seed(50)
  spec_en <- logistic_reg(penalty = tune(),
                          mixture = tune()) %>% 
    set_engine("glmnet", family = "binomial") %>% 
    set_mode("classification")
  grid_en <- dials::grid_latin_hypercube(
    dials::penalty(), mixture(),
    size = grid_size)
  wf_en <- workflows::workflow() %>%
    add_model(spec_en) %>% 
    add_formula(Formula)
  en_tuned <- tune::tune_grid(
    object = wf_en,
    resamples = folds,
    grid = grid_en,
    metrics = metric_set(mn_log_loss),
    control = tune::control_grid(verbose = TRUE, 
                                 allow_par = FALSE))
  best_log_loss_en <- select_best(en_tuned, "mn_log_loss")
  tuned_final_en <- finalize_workflow(wf_en, best_log_loss_en) 
  fit_en <- tuned_final_en  %>% fit(data = cal)
  new[["prob_en_refit"]] <- predict(fit_en, new_data = new, 
                                    type = "prob")[[".pred_1"]]
  set.seed(50)
  spec_rf <- rand_forest(mtry = tune(), trees = 500,
                         min_n = tune()) %>%
    set_engine("ranger") %>%
    set_mode("classification") %>% set_args(seed = NULL)
  grid_rf <- dials::grid_latin_hypercube(
    finalize(dials::mtry(), 
             x = cal[, vars, drop = FALSE]),
    dials::min_n(),
    size = grid_size
  )
  grid_rf_aux <- grid_rf[1, ]
  grid_rf_aux[["mtry"]] <- floor(sqrt(ncol(cal[, vars, 
                                               drop = FALSE])))
  grid_rf_aux[["min_n"]] <- floor(nrow(cal) * 0.05)
  grid_rf <- rbind(grid_rf, grid_rf_aux)
  grid_rf <- grid_rf[!duplicated(grid_rf), ]
  rm(grid_rf_aux)
  wf_rf <- workflows::workflow() %>%
    add_model(spec_rf) %>% 
    add_formula(Formula)
  rf_tuned <- tune::tune_grid(
    object = wf_rf,
    resamples = folds,
    grid = grid_rf,
    metrics = metric_set(mn_log_loss),
    control = tune::control_grid(verbose = TRUE, 
                                 allow_par = FALSE))
  best_log_loss_rf <- select_best(rf_tuned, "mn_log_loss")
  tuned_final_rf <- finalize_workflow(wf_rf, best_log_loss_rf)
  fit_rf <- tuned_final_rf  %>% fit(data = cal)
  new[["prob_rf_refit"]] <- predict(fit_rf, new_data = new, 
                                    type = "prob")[[".pred_1"]]
  set.seed(50)
  spec_svm <- svm_rbf(cost = tune(),
                      rbf_sigma = tune()) %>% 
    set_engine('kernlab', type = "C-svc", kernel = "rbfdot",prob.model = TRUE) %>% 
    set_mode("classification")
  grid_svm <- dials::grid_latin_hypercube(
    dials::cost(), rbf_sigma(),
    size = grid_size)
  wf_svm <- workflows::workflow() %>%
    add_model(spec_svm) %>% 
    add_formula(Formula)
  svm_tuned <- tune::tune_grid(
    object = wf_svm,
    resamples = folds,
    grid = grid_svm,
    metrics = metric_set(mn_log_loss),
    control = tune::control_grid(verbose = TRUE, 
                                 allow_par = FALSE))
  best_log_loss_svm <- select_best(svm_tuned, "mn_log_loss")
  tuned_final_svm <- finalize_workflow(wf_svm, best_log_loss_svm) 
  fit_svm <- tuned_final_svm  %>% fit(data = cal)
  new[["prob_svm_refit"]] <- predict(fit_svm, new_data = new, 
                                     type = "prob")[[".pred_1"]]
  grid_size <- 50
  set.seed(50)
  spec_gb <- boost_tree(learn_rate = tune(),
                        tree_depth = tune()) %>% 
    set_engine("xgboost") %>% 
    set_mode("classification")
  grid_gb <- dials::grid_latin_hypercube(
    dials::learn_rate(),tree_depth(),
    size = grid_size)
  wf_gb <- workflows::workflow() %>%
    add_model(spec_gb) %>% 
    add_formula(Formula)
  gb_tuned <- tune::tune_grid(
    object = wf_gb,
    resamples = folds,
    grid = grid_gb,
    metrics = metric_set(mn_log_loss),
    control = tune::control_grid(verbose = TRUE, 
                                 allow_par = FALSE))
  best_log_loss_gb <- select_best(gb_tuned, "mn_log_loss")
  tuned_final_gb <- finalize_workflow(wf_gb, best_log_loss_gb) 
  fit_gb <- tuned_final_gb  %>% fit(data = cal)
  new[["prob_gb_refit"]] <- predict(fit_gb, new_data = new, 
                                    type = "prob")[[".pred_1"]]
  lr.refit<-stepPvalue(cal[,c(2:12)],0.001)
  new[["prob_lr_refit"]] <- plogis(predict(lr.refit[["fit"]], newdata=new))
  refit.out<-rbind(refit.out, new)
}
df_test<-merge(df_test,refit.out[,c(1, 21:25)], by='id')

# Bootstrap ----------------------------------------------------
# ==============================================================================
out_B<-NULL
for (k in 1:500) {
  set.seed(k)
  ind1 <- sample(x = seq_len(nrow(df_train)), size = nrow(df_train), replace = TRUE)
  ind2 <- sample(x = seq_len(nrow(df_test)), size = nrow(df_test), replace = TRUE)
  df.train <- df_train[ind1, ]
  df.test <- df_test[ind2,c(2:14)]
  folds <- vfold_cv(cal, v = 2, repeats = 2)
  grid_size <- 50
  set.seed(50)
  spec_en <- logistic_reg(penalty = tune(),
                          mixture = tune()) %>% 
    set_engine("glmnet", family = "binomial") %>% 
    set_mode("classification")
  grid_en <- dials::grid_latin_hypercube(
    dials::penalty(), mixture(),
    size = grid_size)
  wf_en <- workflows::workflow() %>%
    add_model(spec_en) %>% 
    add_formula(Formula)
  en_tuned <- tune::tune_grid(
    object = wf_en,
    resamples = folds,
    grid = grid_en,
    metrics = metric_set(mn_log_loss),
    control = tune::control_grid(verbose = TRUE, 
                                 allow_par = FALSE))
  best_log_loss_en <- select_best(en_tuned, "mn_log_loss")
  tuned_final_en <- finalize_workflow(wf_en, best_log_loss_en) 
  fit_en <- tuned_final_en  %>% fit(data = df.train)
  df.test[["prob_en"]] <- predict(fit_en, new_data = df.test, 
                                  type = "prob")[[".pred_1"]]
  set.seed(50)
  spec_rf <- rand_forest(mtry = tune(), trees = 500,
                         min_n = tune()) %>%
    set_engine("ranger") %>%
    set_mode("classification") %>% set_args(seed = NULL)
  vars <- all.vars(Formula)[-1]
  grid_rf <- dials::grid_latin_hypercube(
    finalize(dials::mtry(), 
             x = cal[, vars, drop = FALSE]),
    dials::min_n(),
    size = grid_size
  )
  grid_rf_aux <- grid_rf[1, ]
  grid_rf_aux[["mtry"]] <- floor(sqrt(ncol(df.train[, vars, 
                                                    drop = FALSE])))
  grid_rf_aux[["min_n"]] <- floor(nrow(df.train) * 0.05)
  grid_rf <- rbind(grid_rf, grid_rf_aux)
  grid_rf <- grid_rf[!duplicated(grid_rf), ]
  rm(grid_rf_aux)
  wf_rf <- workflows::workflow() %>%
    add_model(spec_rf) %>% 
    add_formula(Formula)
  rf_tuned <- tune::tune_grid(
    object = wf_rf,
    resamples = folds,
    grid = grid_rf,
    metrics = metric_set(mn_log_loss),
    control = tune::control_grid(verbose = TRUE, 
                                 allow_par = FALSE))
  best_log_loss_rf <- select_best(rf_tuned, "mn_log_loss")
  tuned_final_rf <- finalize_workflow(wf_rf, best_log_loss_rf)
  fit_rf <- tuned_final_rf  %>% fit(data = df.train)
  df.test[["prob_rf"]] <- predict(fit_rf, new_data = df.test, 
                                  type = "prob")[[".pred_1"]]
  set.seed(50)
  spec_svm <- svm_rbf(cost = tune(),
                      rbf_sigma = tune()) %>% 
    set_engine('kernlab', type = "C-svc", kernel = "rbfdot",prob.model = TRUE) %>% 
    set_mode("classification")
  grid_svm <- dials::grid_latin_hypercube(
    dials::cost(), rbf_sigma(),
    size = grid_size)
  wf_svm <- workflows::workflow() %>%
    add_model(spec_svm) %>% 
    add_formula(Formula)
  svm_tuned <- tune::tune_grid(
    object = wf_svm,
    resamples = folds,
    grid = grid_svm,
    metrics = metric_set(mn_log_loss),
    control = tune::control_grid(verbose = TRUE, 
                                 allow_par = FALSE))
  best_log_loss_svm <- select_best(svm_tuned, "mn_log_loss")
  tuned_final_svm <- finalize_workflow(wf_svm, best_log_loss_svm) 
  fit_svm <- tuned_final_svm  %>% fit(data = df.train)
  df.test[["prob_svm"]] <- predict(fit_svm, new_data = df.test, 
                                   type = "prob")[[".pred_1"]]
  grid_size <- 50
  set.seed(50)
  spec_gb <- boost_tree(learn_rate = tune(),
                        tree_depth = tune()) %>% 
    set_engine("xgboost") %>% 
    set_mode("classification")
  grid_gb <- dials::grid_latin_hypercube(
    dials::learn_rate(),tree_depth(),
    size = grid_size)
  wf_gb <- workflows::workflow() %>%
    add_model(spec_gb) %>% 
    add_formula(Formula)
  gb_tuned <- tune::tune_grid(
    object = wf_gb,
    resamples = folds,
    grid = grid_gb,
    metrics = metric_set(mn_log_loss),
    control = tune::control_grid(verbose = TRUE, 
                                 allow_par = FALSE))
  best_log_loss_gb <- select_best(gb_tuned, "mn_log_loss")
  tuned_final_gb <- finalize_workflow(wf_gb, best_log_loss_gb) 
  fit_gb <- tuned_final_gb  %>% fit(data = df.train)
  df.test[["prob_gb"]] <- predict(fit_gb, new_data = df.test, 
                                  type = "prob")[[".pred_1"]]
  lr.fit<-stepPvalue(df.train[,c(1:11)],0.001)
  df.test[["prob_lr"]] <- plogis(predict(lr.fit[["fit"]], newdata=df.test))
  set.seed(50)
  hung<-df.test[df.test$study=='Hung', c(1:13)]
  fold_hung<-vfold_cv(hung, v=10, repeats = 1, strata='CAD')
  hung$cv<-0
  for (i in 1:10){
    hung$cv[-fold_hung$splits[[i]]$in_id]=i
  }
  swiss<-df.test[df.test$study=='Swiss', c(1:13)]
  fold_swiss<-vfold_cv(swiss, v=10, repeats = 1, strata='CAD')
  swiss$cv<-0
  for (i in 1:10){
    swiss$cv[-fold_swiss$splits[[i]]$in_id]=i
  }
  va<-df.test[df.test$study=='VA', c(1:13)]
  fold_va<-vfold_cv(va, v=10, repeats = 1, strata='CAD')
  va$cv<-0
  for (i in 1:10){
    va$cv[-fold_va$splits[[i]]$in_id]=i
  }
  df.test$cv<-c(hung$cv, swiss$cv, va$cv)
  df.test$id<-1:nrow(df.test)
  machine<-c('lr','rf','svm','en','gb')
  cal.out<-df.test
  set.seed(50)
  for (j in machine){
    out<-NULL
    for (i in 1:10){
      cal<-cal.out[cal.out$cv!=i,]
      new<-cal.out[cal.out$cv==i,]
      new[[paste0('prob_',j,'_beta')]]<-cal_beta ('CAD', df.train, cal, new, paste0('prob_',j))
      new[[paste0('prob_',j,'_elkan')]]<-cal_elkan ('CAD', df.train, cal, new, paste0('prob_',j))
      new[[paste0('prob_',j,'_iso')]]<-cal_iso ('CAD', df.train, cal, new, paste0('prob_',j))
      new[[paste0('prob_',j,'_log')]]<-cal_logistic ('CAD', df.train, cal, new, paste0('prob_',j))
      new[[paste0('prob_',j,'_int')]]<-cal_loginter ('CAD', df.train, cal, new, paste0('prob_',j))
      new[[paste0('prob_',j,'_slope')]]<-cal_logslope ('CAD', df.train, cal, new, paste0('prob_',j))
      new[[paste0('prob_',j,'_platt')]]<-cal_platt ('CAD', df.train, cal, new, paste0('prob_',j))
      new[[paste0('prob_',j,'_chen')]]<-cal_chen('CAD', df.train, cal, new, paste0('prob_',j))
      out<-rbind(out,new)
    }
    out<-out[,c(20:28)]
    df.test<-merge(df.test, out, by='id')
  }
  refit<-df.test[,c(1:20)]
  refit.out<-NULL
  for (i in 1:10){
    cal<-refit[refit$cv!=i,]
    new<-refit[refit$cv==i,]
    set.seed(50)
    folds <- vfold_cv(cal, v = 2, repeats = 2)
    grid_size <- 50
    set.seed(50)
    spec_en <- logistic_reg(penalty = tune(),
                            mixture = tune()) %>% 
      set_engine("glmnet", family = "binomial") %>% 
      set_mode("classification")
    grid_en <- dials::grid_latin_hypercube(
      dials::penalty(), mixture(),
      size = grid_size)
    wf_en <- workflows::workflow() %>%
      add_model(spec_en) %>% 
      add_formula(Formula)
    en_tuned <- tune::tune_grid(
      object = wf_en,
      resamples = folds,
      grid = grid_en,
      metrics = metric_set(mn_log_loss),
      control = tune::control_grid(verbose = TRUE, 
                                   allow_par = FALSE))
    best_log_loss_en <- select_best(en_tuned, "mn_log_loss")
    tuned_final_en <- finalize_workflow(wf_en, best_log_loss_en) 
    fit_en <- tuned_final_en  %>% fit(data = cal)
    new[["prob_en_refit"]] <- predict(fit_en, new_data = new, 
                                      type = "prob")[[".pred_1"]]
    set.seed(50)
    spec_rf <- rand_forest(mtry = tune(), trees = 500,
                           min_n = tune()) %>%
      set_engine("ranger") %>%
      set_mode("classification") %>% set_args(seed = NULL)
    grid_rf <- dials::grid_latin_hypercube(
      finalize(dials::mtry(), 
               x = cal[, vars, drop = FALSE]),
      dials::min_n(),
      size = grid_size
    )
    grid_rf_aux <- grid_rf[1, ]
    grid_rf_aux[["mtry"]] <- floor(sqrt(ncol(cal[, vars, 
                                                 drop = FALSE])))
    grid_rf_aux[["min_n"]] <- floor(nrow(cal) * 0.05)
    grid_rf <- rbind(grid_rf, grid_rf_aux)
    grid_rf <- grid_rf[!duplicated(grid_rf), ]
    rm(grid_rf_aux)
    wf_rf <- workflows::workflow() %>%
      add_model(spec_rf) %>% 
      add_formula(Formula)
    rf_tuned <- tune::tune_grid(
      object = wf_rf,
      resamples = folds,
      grid = grid_rf,
      metrics = metric_set(mn_log_loss),
      control = tune::control_grid(verbose = TRUE, 
                                   allow_par = FALSE))
    best_log_loss_rf <- select_best(rf_tuned, "mn_log_loss")
    tuned_final_rf <- finalize_workflow(wf_rf, best_log_loss_rf)
    fit_rf <- tuned_final_rf  %>% fit(data = cal)
    new[["prob_rf_refit"]] <- predict(fit_rf, new_data = new, 
                                      type = "prob")[[".pred_1"]]
    set.seed(50)
    spec_svm <- svm_rbf(cost = tune(),
                        rbf_sigma = tune()) %>% 
      set_engine('kernlab', type = "C-svc", kernel = "rbfdot",prob.model = TRUE) %>% 
      set_mode("classification")
    grid_svm <- dials::grid_latin_hypercube(
      dials::cost(), rbf_sigma(),
      size = grid_size)
    wf_svm <- workflows::workflow() %>%
      add_model(spec_svm) %>% 
      add_formula(Formula)
    svm_tuned <- tune::tune_grid(
      object = wf_svm,
      resamples = folds,
      grid = grid_svm,
      metrics = metric_set(mn_log_loss),
      control = tune::control_grid(verbose = TRUE, 
                                   allow_par = FALSE))
    best_log_loss_svm <- select_best(svm_tuned, "mn_log_loss")
    tuned_final_svm <- finalize_workflow(wf_svm, best_log_loss_svm) 
    fit_svm <- tuned_final_svm  %>% fit(data = cal)
    new[["prob_svm_refit"]] <- predict(fit_svm, new_data = new, 
                                       type = "prob")[[".pred_1"]]
    grid_size <- 50
    set.seed(50)
    spec_gb <- boost_tree(learn_rate = tune(),
                          tree_depth = tune()) %>% 
      set_engine("xgboost") %>% 
      set_mode("classification")
    grid_gb <- dials::grid_latin_hypercube(
      dials::learn_rate(),tree_depth(),
      size = grid_size)
    wf_gb <- workflows::workflow() %>%
      add_model(spec_gb) %>% 
      add_formula(Formula)
    gb_tuned <- tune::tune_grid(
      object = wf_gb,
      resamples = folds,
      grid = grid_gb,
      metrics = metric_set(mn_log_loss),
      control = tune::control_grid(verbose = TRUE, 
                                   allow_par = FALSE))
    best_log_loss_gb <- select_best(gb_tuned, "mn_log_loss")
    tuned_final_gb <- finalize_workflow(wf_gb, best_log_loss_gb) 
    fit_gb <- tuned_final_gb  %>% fit(data = cal)
    new[["prob_gb_refit"]] <- predict(fit_gb, new_data = new, 
                                      type = "prob")[[".pred_1"]]
    lr.refit<-stepPvalue(cal[,c(2:12)],0.001)
    new[["prob_lr_refit"]] <- plogis(predict(lr.refit[["fit"]], newdata=new))
    refit.out<-rbind(refit.out, new)
  }
  df.test<-merge(df.test,refit.out[,c(1, 21:25)], by='id')
  df.test$bootstrap<-k
  out_B<-rbind(out_B, df.test)
  message(k, " out of 500 is done.")
}

# Scatter plots of predicted and observed probability
# ==============================================================================
library(ggplot2)
library(gridExtra)
library(ggpubr)
df.plot<-NULL
machine=c('lr','rf','svm','en','gb')
calib<-c('NoCal','Beta','Elkan','Isotonic','Logistic','Intercept','Slope','Platt','Chen','Refit')
for (j in machine){
  names<-c(paste0('prob_',j), paste0('prob_',j,'_beta'),paste0('prob_',j,'_elkan'),paste0('prob_',j,'_iso'),paste0('prob_',j,'_log'),paste0('prob_',j,'_int'),
           paste0('prob_',j,'_slope'),paste0('prob_',j,'_platt'),paste0('prob_',j,'_chen'),paste0('prob_',j,'_refit'))
  temp<-df_test[,names]
  temp<-data.frame(prob=c(as.matrix(temp)))
  temp$machine<-j
  temp$calibration=0
  for (i in 1:10){
    temp$calibration[c((1+(dim(df_test)[1])*(i-1)):((dim(df_test)[1])*i))]=calib[i]
  }
  df.plot<-rbind(df.plot,temp)
}
df.plot$CAD<-rep(df_test$CAD, times=50)
df.plot$machine[df.plot$machine=='lr']='Logistic regression'
df.plot$machine[df.plot$machine=='rf']='Random forest'
df.plot$machine[df.plot$machine=='svm']='Support vector machine'
df.plot$machine[df.plot$machine=='en']='Elastic net'
df.plot$machine[df.plot$machine=='gb']='Gradient boosting'
df.plot$calibration<-factor(df.plot$calibration, levels=c('NoCal','Logistic','Beta','Intercept','Slope','Refit','Platt','Isotonic','Elkan','Chen'))
df.plot$machine<-factor(df.plot$machine, levels=c('Logistic regression','Random forest','Support vector machine','Elastic net','Gradient boosting'))
pdf('scatter_lr_ci.pdf',width=18,height=6)
ggplot(df.plot[df.plot$machine=='Logistic regression',], aes(x = prob, y=CAD)) +
  geom_abline(slope = 1, intercept = 0, lty = 1, colour = "grey50", size=0.5)+
  geom_smooth (method=loess, se = TRUE, lty = 1, size=0.5, colour = "black", fill='grey60') +
  facet_wrap(vars(calibration),nrow = 2) + xlim(0,1)+ylim(0,1) +
  xlab("Predicted probability") + ylab("Observed proportion") + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))
dev.off()
pdf('scatter_rf_ci.pdf',width=18,height=6)
ggplot(df.plot[df.plot$machine=='Random forest',], aes(x = prob, y=CAD)) +
  geom_abline(slope = 1, intercept = 0, lty = 1, colour = "grey50", size=0.5)+
  geom_smooth (method=loess, se = TRUE, lty = 1, size=0.5, colour = "black", fill='grey60') +
  facet_wrap(vars(calibration),nrow = 2) + xlim(0,1)+ylim(0,1) + 
  xlab("Predicted probability") + ylab("Observed proportion") + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))
dev.off()
pdf('scatter_svm_ci.pdf',width=18,height=6)
ggplot(df.plot[df.plot$machine=='Support vector machines',], aes(x = prob, y=CAD)) +
  geom_abline(slope = 1, intercept = 0, lty = 1, colour = "grey50", size=0.5)+
  geom_smooth (method=loess, se = TRUE, lty = 1, size=0.5, colour = "black", fill='grey60') +
  facet_wrap(vars(calibration),nrow = 2) + xlim(0,1)+ylim(0,1) + 
  xlab("Predicted probability") + ylab("Observed proportion") + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))
dev.off()
pdf('scatter_en_ci.pdf',width=18,height=6)
ggplot(df.plot[df.plot$machine=='Elastic net',], aes(x = prob, y=CAD)) +
  geom_abline(slope = 1, intercept = 0, lty = 1, colour = "grey50", size=0.5)+
  geom_smooth (method=loess, se = TRUE, lty = 1, size=0.5, colour = "black", fill='grey60') +
  facet_wrap(vars(calibration),nrow = 2) + xlim(0,1)+ylim(0,1) + 
  xlab("Predicted probability") + ylab("Observed proportion") + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))
dev.off()
pdf('scatter_gb_ci.pdf',width=18,height=6)
ggplot(df.plot[df.plot$machine=='Gradient boost',], aes(x = prob, y=CAD)) +
  geom_abline(slope = 1, intercept = 0, lty = 1, colour = "grey50", size=0.5)+
  geom_smooth (method=loess, se = TRUE, lty = 1, size=0.5, colour = "black", fill='grey60') +
  facet_wrap(vars(calibration),nrow = 2) + xlim(0,1)+ylim(0,1) + 
  xlab("Predicted probability") + ylab("Observed proportion") + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))
dev.off()

# Boxplots for Brier score and logloss values after bootstrap
# ==============================================================================
BS<-NULL
logloss<-NULL
for (j in machine){
  output_BS<-matrix(NA, nrow=500, ncol=10)
  output_logloss<-matrix(NA, nrow=500, ncol=10)
  for (k in 1:500){
    boot1<-out_B[out_B$bootstrap==k,]
    name<-c('CAD', paste0('prob_',j), paste0('prob_',j,'_beta'), paste0('prob_',j,'_elkan'), paste0('prob_',j,'_iso'), paste0('prob_',j,'_log'),
            paste0('prob_',j,'_int'), paste0('prob_',j,'_slope'), paste0('prob_',j,'_platt'), paste0('prob_',j,'_chen'), paste0('prob_',j,'_refit'))
    boot2<-boot1[,name]
    boot2<-boot2[order(boot2[,2],decreasing=F),]
    for (i in 1:10){
      output_BS[k,i]<-sum((boot2$CAD - boot2[,(i+1)])^2)/length(boot2$CAD)
    }
    for (i in 1:10){
      boot2$CAD<-as.factor(boot2$CAD)
      output_logloss[k,i]<-mn_log_loss(boot2,CAD,name[i+1], event_level='second')$.estimate
    }
  }
  output_BS<-as.data.frame(output_BS)
  output_BS$machine<-j
  BS<-rbind(BS,output_BS)
  output_logloss<-as.data.frame(output_logloss)
  output_logloss$machine<-j
  logloss<-rbind(logloss,output_logloss)
}
bsboot<-NULL
calib<-c('NoCal','Beta','Elkan','Isotonic','Logistic','Intercept','Slope','Platt','Chen','Refit')
for (i in 1:10){
  temp<-BS[,c(i,11)]
  names(temp)[1]<-'BS'
  temp$calibration<-calib[i]
  bsboot<-rbind(bsboot, temp)
}
bsboot$calibration<-factor(bsboot$calibration, levels=rev(c('NoCal', 'Logistic', 'Beta', 'Intercept', 'Slope','Refit', 'Platt', 'Isotonic','Elkan','Chen')))
bsboot$machine[bsboot$machine=='lr']<-'Logistic regression'
bsboot$machine[bsboot$machine=='rf']<-'Random forest'
bsboot$machine[bsboot$machine=='svm']<-'Support vector machine'
bsboot$machine[bsboot$machine=='en']<-'Elastic net'
bsboot$machine[bsboot$machine=='gb']<-'Gradient boosting'
bsboot$machine<-factor(bsboot$machine, levels=c('Logistic regression','Random forest','Support vector machine','Elastic net','Gradient boosting'))
pdf('BS_boot_bw.pdf',width=12,height=8)
ggplot(bsboot, aes(x = calibration, y = BS)) + geom_boxplot() + 
  facet_wrap(vars(machine)) + xlab("") + ylab("") + coord_flip() +theme_bw() +
  theme(text = element_text(size = 15))
dev.off()
loglossboot<-NULL
for (i in 1:10){
  temp<-logloss[,c(i,11)]
  names(temp)[1]<-'logloss'
  temp$calibration<-calib[i]
  loglossboot<-rbind(loglossboot, temp)
}
loglossboot$calibration<-factor(loglossboot$calibration, levels=rev(c('NoCal', 'Logistic', 'Beta', 'Intercept', 'Slope','Refit', 'Platt', 'Isotonic','Elkan','Chen')))
loglossboot$machine[loglossboot$machine=='lr']<-'Logistic regression'
loglossboot$machine[loglossboot$machine=='rf']<-'Random forest'
loglossboot$machine[loglossboot$machine=='svm']<-'Support vector machine'
loglossboot$machine[loglossboot$machine=='en']<-'Elastic net'
loglossboot$machine[loglossboot$machine=='gb']<-'Gradient boosting'
loglossboot$machine<-factor(loglossboot$machine, levels=c('Logistic regression','Random forest','Support vector machine','Elastic net','Gradient boosting'))
pdf('logloss_boot_bw.pdf',width=12,height=8)
ggplot(loglossboot, aes(x = calibration, y = logloss)) + geom_boxplot() + 
  facet_wrap(vars(machine)) + xlab("") + ylab("") + coord_flip() +theme_bw() +
  theme(text = element_text(size = 15)) 
dev.off()

# Radar plots
# ==============================================================================
# Radar plots for effect estimates and 95%CI in logistic regression models
library(fmsb)
coef<-as.data.frame(lr.cleve$coef)
coef[4,]<-coef[4,]*10
coef$low<-coef$Estimate-1.96*coef$`Std. Error`
coef$high<-coef$Estimate+1.96*coef$`Std. Error`
coef$var<-rownames(coef)
coef<-coef[order(coef[,7],decreasing=F),]
var<-coef$var[-1]
coef<-coef[-1, c(1,5,6)]
coef<-as.data.frame(t(coef))
names(coef)<-var
names(coef)<-c('Chest pain type','Max heart rate/10','Sex','ST depression')
rownames(coef)<-c('Beta', 'Lower limit of the 95%CI', 'Upper limit of the 95%CI')
coef<-rbind(rep(3, dim(coef)[2]), rep(-1, dim(coef)[2]), coef)
val<-rbind(Swiss, Hung, VA)
fit.val <- glm(as.factor(CAD) ~ Sex+Chestpaintype+MaxHeartRate+STDepression, data = val[, c(1:11)], family = binomial())
coef.val<-as.data.frame(summary(fit.val)$coef)
coef.val[4,]<-coef.val[4,]*10
coef.val$low<-coef.val$Estimate-1.96*coef.val$`Std. Error`
coef.val$high<-coef.val$Estimate+1.96*coef.val$`Std. Error`
coef.val$var<-rownames(coef.val)
coef.val<-coef.val[order(coef.val[,7],decreasing=F),]
coef.val<-coef.val[-1, c(1,5,6)]
coef.val<-as.data.frame(t(coef.val))
names(coef.val)<-names(coef)
rownames(coef.val)<-c('Beta', 'Lower limit of the 95%CI', 'Upper limit of the 95%CI')
coef.val<-rbind(rep(3, dim(coef.val)[2]), rep(-1, dim(coef.val)[2]), coef.val)
pdf('LR_Cleve_radar.pdf',width=13,height=6)
par(mfrow=c(1,2))
line_col<-c('black', rgb(233,80,14, , maxColorValue = 255), rgb(64,133,198, , maxColorValue = 255))
line_type<-c(1, 3, 3)
radarchart(coef, pcol=line_col, plty=line_type, vlcex=0.8,cglcol="grey", cglty=1,axislabcol="grey", axistype=1, calcex=1, caxislabels=seq(-1,3,1),seg=4)
mtext('A',side=3,line=1,cex=1,adj=0)
radarchart(coef.val, pcol=line_col, plty=line_type, vlcex=0.8,cglcol="grey", cglty=1,axislabcol="grey", axistype=1, calcex=1, caxislabels=seq(-1,3,1),seg=4)
mtext('B',side=3,line=1,cex=1,adj=0)
dev.off()
# Radar plots for permutation importance and P values of each variable in random forest
rangerCleve<-ranger (as.factor(CAD) ~., data=Cleve[,-12], min.node.size=6, mtry=2, num.tree=500, seed=50, write.forest=T,importance='permutation', probability = T)
imp<-importance_pvalues(rangerCleve, method='altmann',num.permutations=10000,formula=as.factor(CAD) ~., data=Cleve[,-12])
imp$var<-rownames(imp)
imp<-imp[order(imp[,3],decreasing=F),]
imp1<-as.data.frame(t(imp[,1]))
imp2<-as.data.frame(t(imp[,2]))
imp1[1,]<-imp1[1,]*1000
imp2[1,]<--log(imp2[1,],base=10)
imp1<-rbind(rep(60,dim(imp1)[2]), rep(0,dim(imp1)[2]), imp1)
imp2<-rbind(rep(5,dim(imp2)[2]), rep(0,dim(imp2)[2]), imp2)
imp2[4,]<--log(0.05,base=10)
names(imp1)<-c('Age','Chest pain type','Cholesterol','Exercise induced angina','High blood glucose','Max heart rate','Resting blood pressure','Resting ECG','Sex','ST depression')
names(imp2)<-c('Age','Chest pain type','Cholesterol','Exercise induced angina','High blood glucose','Max heart rate','Resting blood pressure','Resting ECG','Sex','ST depression')
rangerVAL<-ranger (as.factor(CAD) ~., data=val[c(1:11)], min.node.size=6, mtry=2, num.tree=500, seed=50, write.forest=T,importance='permutation', probability = T)
impVAL<-importance_pvalues(rangerVAL, method='altmann',num.permutations=10000,formula=as.factor(CAD) ~., data=val[,c(1:11)])
impVAL$var<-rownames(impVAL)
impVAL<-impVAL[order(impVAL[,3],decreasing=F),]
val.imp<-as.data.frame(t(impVAL[,1]))
val.p<-as.data.frame(t(impVAL[,2]))
val.imp<-val.imp*1000
val.p<--log(val.p,base=10)
val.imp<-rbind(rep(60,dim(val.imp)[2]), rep(0,dim(val.imp)[2]), val.imp)
val.p<-rbind(rep(5,dim(val.p)[2]), rep(0,dim(val.p)[2]), val.p)
names(val.imp)<-c('Age','Chest pain type','Cholesterol','Exercise induced angina','High blood glucose','Max heart rate','Resting blood pressure','Resting ECG','Sex','ST depression')
names(val.p)<-c('Age','Chest pain type','Cholesterol','Exercise induced angina','High blood glucose','Max heart rate','Resting blood pressure','Resting ECG','Sex','ST depression')
val.p[4,]<--log(0.05,base=10)
pdf('RF_Cleve_radar.pdf',width=12,height=10)
par(mfrow=c(2,2))
radarchart(imp1, vlcex=1.2,cglcol="grey", cglty=1,axislabcol="grey", axistype=1, calcex=1, caxislabels=seq(0,60,10),seg=6)
mtext('A',side=3,line=1,cex=1,adj=0)
radarchart(val.imp, vlcex=1.2,cglcol="grey", cglty=1,axislabcol="grey", axistype=1, calcex=1, caxislabels=seq(0,60,10),seg=6)
mtext('B',side=3,line=1,cex=1,adj=0)
radarchart(imp2, vlcex=1.2,cglcol="grey", cglty=1,axislabcol="grey", axistype=1, calcex=1, caxislabels=seq(0,5,1),seg=5, pcol=c('black','red'),pty=c(16,32))
mtext('C',side=3,line=1,cex=1,adj=0)
radarchart(val.p, vlcex=1.2,cglcol="grey", cglty=1,axislabcol="grey", axistype=1, calcex=1, caxislabels=seq(0,5,1),seg=5, pcol=c('black','red'),pty=c(16,32))
mtext('D',side=3,line=1,cex=1,adj=0)
dev.off()

# CD plots
# ==============================================================================
## Functions
iman_davenport_test <- function(input_matrix, p.value = 0.05) {
  
  rankMatrix <- function(data) {
    f <- function(x) {
      rank(-x, ties.method = "average")
    }
    rankings <- t(apply(data, MARGIN = 1, FUN = f))
    colnames(rankings) <- colnames(data)
    rownames(rankings) <- rownames(data)
    rankings
  }
  
  nts <- nrow(input_matrix)
  nws <- ncol(input_matrix)
  mr <- colMeans(rankMatrix(input_matrix))
  
  chi <- 12 * nts/(nws * (nws + 1)) * (sum(mr^2) - (nws * (nws + 1)^2)/4)
  FF <- (nts - 1) * chi/(nts * (nws - 1) - chi)
  critVal <- df(1 - p.value, nws - 1, (nws - 1) * (nts - 1))
  rejNull <- FF > critVal
  pval_id <- 1 - pf(FF, df1 = nws - 1, df2 = (nws - 1) * (nts - 1))
  
  output <- list(F.test = list(chi = chi, FF = FF, 
                               critVal = critVal, rejNull = rejNull, 
                               pval_iman_davenport = pval_id),
                 Nemenyi.test = NULL,
                 avgRksWFs = mr
  )
  
  if (rejNull) {
    CD.n <- qtukey(1 - p.value, nws, df = Inf)/sqrt(2) * 
      sqrt(nws * (nws + 1)/(6 * nts))
    allRkDifs <- outer(mr, mr, function(x, y) abs(x - y))
    signifDifs <- allRkDifs >= CD.n
    output[["Nemenyi.test"]] <- list(critDif = CD.n, rkDifs = allRkDifs, 
                                     signifDifs = signifDifs)
  }
  
  list(logloss = output)
}
CD_diagram_Nemenyi <- function(r, metric = names(r)[1], coord_fixed_ratio = 0.5, 
                               bw = FALSE, return_ggplot2_object = FALSE) {
  
  if (is.null(r[[metric]]$F.test) | 
      is.null(r[[metric]]$Nemenyi.test)) stop(paste0("Results of both the F ",
                                                     "and Nemenyi tests are ",
                                                     "required for these ",
                                                     "diagrams."))
  
  o <- rank(r[[metric]]$avgRksWFs, ties.method = "first")
  mxl <- ceiling(length(r[[metric]]$avgRksWFs)/2)
  data <- data.frame(avgRk = r[[metric]]$avgRksWFs, 
                     invRk = length(r[[metric]]$avgRksWFs) + 
                       1 - r[[metric]]$avgRksWFs, 
                     sys = names(r[[metric]]$avgRksWFs), 
                     line = o%%mxl + ifelse(o%%mxl == 0, mxl, 0), 
                     side = ifelse(o <= mxl, -1, 1))
  data$line <- ifelse(data$side == 1, mxl + 1 - data$line, 
                      data$line)
  len <- length(r[[metric]]$avgRksWFs)
  cd <- r[[metric]]$Nemenyi.test$critDif
  
  col_arg <- "sys"; if(bw){col_arg <- NULL}
  
  g <- ggplot(data, aes_string(x = "invRk", y = "line")) + 
    geom_segment(aes_string(x = "invRk", y = 0, xend = "invRk", yend = "line", 
                            col = col_arg)) + 
    geom_text(data = data[data$side == -1, ], 
              aes_string(label = "sys", x = len, y = "line", col = col_arg), 
              hjust = 0, size = 4) + 
    geom_text(data = data[data$side == 1, ], 
              aes_string(label = "sys", x = 0, y = "line", 
                         col = col_arg), hjust = 1, size = 4) + 
    geom_segment(data = data[data$side == -1, ], 
                 aes_string(x = max(o),y = "line", xend = "invRk",yend = "line", 
                            col = col_arg)) + 
    geom_segment(data = data[data$side == 1, ], 
                 aes_string(x = "invRk", y = "line", 
                            xend = 0, yend = "line", col = col_arg)) + 
    # scale_x_continuous(limits = c(0, len), breaks = 0:(len + 1), 
    #                    labels = paste((len + 1):0)) + 
    # lines connecting a machine with ranking 1 (drawn at x = len), where
    # not being drawn and warnings like the following were being produced:
    # 'Removed 1 rows containing missing values (geom_segment).'
    # Increasing the width of the limits below solved the problem
    scale_x_continuous(limits = c(0 - 0.1, len + 0.1), breaks = 0:(len + 1), 
                       labels = paste((len + 1):0), expand=c(0.2,0.2)) + 
    scale_y_continuous(limits = c(0, mxl + 1)) + xlab("Average Rank") + 
    ylab("") + 
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
          legend.position = "none", panel.background = element_blank(), 
          panel.border = element_blank(), 
          axis.line = element_line(size = 1), 
          axis.line.y = element_blank(), panel.grid.major = element_blank(), 
          plot.background = element_blank(), 
          plot.margin = grid::unit(c(3, 6, 1, 5), "lines")) + 
    coord_fixed(ratio = coord_fixed_ratio) + 
    annotate("segment", x = 0, xend = cd, y = mxl + 1, yend = mxl + 1, 
             size = 1.5) + annotate("text", x = 0, y = mxl + 1, 
                                    label = paste("Critical Difference", 
                                                  round(cd, 1), sep = " = "), 
                                    vjust = -0.5, hjust = 0, size = 5)
  wfsOrd <- names(r[[metric]]$avgRksWFs[order(r[[metric]]$avgRksWFs)])
  ss <- r[[metric]]$Nemenyi.test$signifDifs[wfsOrd, wfsOrd]
  mx <- ncol(ss)
  pos <- rep(1:mxl, each = mxl) - seq(0, 1, by = 1/(mxl + 1))[-c(1, mxl + 2)]
  frees <- rep(1, mxl)
  from <- 1
  currTill <- 0
  for (i in 1:nrow(ss)) {
    #till <- which(ss[i, i:mx])
    #till <- if (length(till)) till - 1 else mx
    #if (till > currTill) {
    till <- which(ss[i,i:mx]) + (i-1) #add i-1 to get the correct which-values since the column selection removes i-1 values
    till <- if (length(till)>0) till[1]-1 else mx
    if (till > i && till > currTill) {
      theLine <- min(data[wfsOrd[i], "line"], data[wfsOrd[till], 
                                                   "line"])
      ypos <- pos[(theLine - 1) * mxl + frees[theLine]]
      frees[theLine] <- frees[theLine] + 1
      g <- g + annotate("segment", x = data[wfsOrd[till], "invRk"] - 0.1, 
                        xend = data[wfsOrd[i], "invRk"] + 0.1, y = ypos, 
                        yend = ypos, size = 1.2)
      currTill <- till
    }
    if (till == mx) break
  }
  
  if(return_ggplot2_object) return(g)
  
  #grid::grid.newpage()
  gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(g))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  gt
}
## Plots for each machine and averages across machines
library(yardstick)
library(performanceEstimation)
library(scmamp)
library(ggplot2)
library(ggpubr)
machine<-c('lr','rf','svm','en','gb')
calib<-c('NoCal','Beta','Elkan','Isotonic','Logistic','Intercept','Slope','Platt','Chen','Refit')
names(df)[1:10]<-calib
id_test_lr<-iman_davenport_test(logloss[logloss$machine=='lr',-11])
id_test_rf<-iman_davenport_test(logloss[logloss$machine=='rf',-11])
id_test_svm<-iman_davenport_test(logloss[logloss$machine=='svm',-11])
id_test_en<-iman_davenport_test(logloss[logloss$machine=='en',-11])
id_test_gb<-iman_davenport_test(logloss[logloss$machine=='gb',-11])
logloss.ave<-matrix(NA, ncol=10,nrow=500)
for (i in c(1:10)){
  for (j in c(1:500)){
    temp1<-logloss[logloss$machine=='lr',]
    temp2<-logloss[logloss$machine=='rf',]
    temp3<-logloss[logloss$machine=='svm',]
    temp4<-logloss[logloss$machine=='en',]
    temp5<-logloss[logloss$machine=='gb',]
    logloss.ave[j,i]=sum(temp1[j,i]+temp2[j,i]+temp3[j,i]+temp4[j,i]+temp5[j,i])/5
  }
}
logloss.ave<-as.data.frame(logloss.ave)
names(logloss.ave)<-calib
test.ave<-as.matrix(logloss.ave)
id_test_ave <- iman_davenport_test(test.ave)
pdf ('CD_Cleve.pdf', width=20,height=12)
p1<-CD_diagram_Nemenyi(id_test_lr, return_ggplot2_object = T, coord_fixed_ratio = 2)
p2<-CD_diagram_Nemenyi(id_test_rf, return_ggplot2_object = T, coord_fixed_ratio = 2)
p3<-CD_diagram_Nemenyi(id_test_svm, return_ggplot2_object = T, coord_fixed_ratio = 2)
p4<-CD_diagram_Nemenyi(id_test_en, return_ggplot2_object = T, coord_fixed_ratio = 2)
p5<-CD_diagram_Nemenyi(id_test_gb, return_ggplot2_object = T, coord_fixed_ratio = 2)
p6<-CD_diagram_Nemenyi(id_test_ave, return_ggplot2_object = T, coord_fixed_ratio = 2)
ggarrange (p1 + theme(plot.margin=unit(c(1,0,0,0),"cm")), p2 + theme(plot.margin=unit(c(1,0,0,0),"cm")), p3 + theme(plot.margin=unit(c(1,0,0,0),"cm")), 
           p4 + theme(plot.margin=unit(c(1,0,0,0),"cm")), p5 + theme(plot.margin=unit(c(1,0,0,0),"cm")), p6 + theme(plot.margin=unit(c(1,0,0,0),"cm")), 
           labels = 'AUTO', ncol=3, nrow=2, font.label = list(size = 20))
dev.off()
