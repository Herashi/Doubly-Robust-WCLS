setwd("~/Documents/DML-WCLS/Case Study")
# library(geepack)
library(dplyr)
library(ggplot2)
library(splines)
library(tidyr)
# library(foreach)
# library(parallel)
# library(doParallel)
library(randomForest)
library(stringr)
load("Data_missing.RData")
# load("imputation_list_daily_separated_20.RData")
# IHSdata_2018 <- read_csv("IHSdata_2018.csv")
# cur_list = impute_list[[9]]
# 
# all_baseline = cur_list$all_baseline
# full_data = cur_list$full_data
# colnames(full_data)[1] = "UserID"
# 
# 
# baseline_step = apply(all_baseline[,c(10,13,16)], MARGIN = 1, FUN = mean)
# baseline_sleep = apply(all_baseline[,c(11,14,17)], MARGIN = 1, FUN = mean)
# baseline_mood = apply(all_baseline[,c(12,15,18)], MARGIN = 1, FUN = mean)
# baseline_average = cbind(all_baseline$UserID, baseline_step, baseline_sleep, baseline_mood)
# baseline_average = data.frame(baseline_average)
# names(baseline_average)[1] = 'UserID'
# baseline_average$study_week = 1
# names(baseline_average)[2:4] = c('STEP_COUNTprev', 'SLEEP_COUNTprev', 'MOODprev')
# baseline_average$week_categoryprev = NA
# 
# aggregate_weekly = aggregate(full_data[, 4:6], by = full_data[,c(1,3,7)], FUN = mean,na.rm=TRUE)
# aggregate_weekly2 = aggregate_weekly
# names(aggregate_weekly2)[3:6] = paste(names(aggregate_weekly2)[3:6], "prev", sep = '')
# aggregate_weekly2$study_week = aggregate_weekly2$study_week + 1
# aggregate_weekly2 = rbind(aggregate_weekly2, baseline_average)
# aggregate_weekly_new = merge(aggregate_weekly, aggregate_weekly2, by = c('UserID', 'study_week'), all.x = TRUE)
# 
# aggregate_weekly_new1 = merge(aggregate_weekly_new, all_baseline[,1:9], by = 'UserID', all.x = TRUE)
# 
# analysis_dat = aggregate_weekly_new1[, -7]
# analysis_dat = analysis_dat[analysis_dat$week_category != 'unsure', ]
# analysis_dat$week_category_new = as.numeric(analysis_dat$week_category != 'None')
# analysis_dat$week_category = relevel(analysis_dat$week_category, ref = 'None')
# analysis_dat_gee = analysis_dat[order(analysis_dat$UserID, analysis_dat$study_week), ]
# analysis_dat_gee$week_category= droplevels(analysis_dat_gee$week_category)
# 
# data = analysis_dat_gee
# 
# data[sapply(data, is.nan)] <- NA
# 
# save(data, file = "Data_missing.RData")
# unequal fold-size

CV_id = data.frame(UserID = unique(data$UserID),
                   cv_id = sample(1:5, length(unique(data$UserID)), replace = TRUE))


# equal fold size

# CV_id = data.frame(UserID = unique(data$UserID),
#   cv_id = rep(1:3, length.out =length(unique(data$UserID)) ))

# table(CV_id$cv_id)
data = left_join(data, CV_id,by = "UserID")

data = data %>%
  mutate(R = if_else(is.na(MOOD),0,1))%>%
  arrange(UserID,study_week)%>%
  group_by(UserID)%>%
  mutate(R_prev =lag(R, n = 1, default = 1),
         R_cum_avg = cummean(R),
         R_factor = as.factor(R))

data$week_category_new = ifelse(data$week_category == "mood",1,0)
data$pn = mean(data$week_category_new)

# data = data %>% group_by(UserID) %>% mutate(week_category_new_lag1 = lag(week_category_new))
# data$week_category_new_lag1[is.na(data$week_category_new_lag1)] <- -1

data = data %>% mutate(week_category_new_c = week_category_new - pn,
                       prob = 1/4,
                       w = if_else(week_category_new==1, pn/prob, (1-pn)/(1-prob)) )

# r_train = randomForest(as.formula(paste0("R_factor ~ R_cum_avg + R_prev + PHQtot0 + Neu0 + depr0 + EFE0 + pre_intern_mood + pre_intern_sleep +pre_intern_sqrt_step")),
#                        data = data, importance = TRUE, ntree=300,mtry = 3)

Cov_set = colnames(data)[c(7:17)]
Cov_set = c(Cov_set,"study_week","R_cum_avg","R_prev")
Control_var = paste0(Cov_set, collapse = " + ",sep = "")

r_train = randomForest(as.formula(paste0("R_factor ~", Control_var)),
                                   data = data%>% drop_na(all_of(Cov_set)), importance = TRUE, ntree=300,mtry = 3)
#
data$prob_r_full <- unname(predict(r_train,data, type = "prob")[,"1"])
data$r_w = data$R / data$prob_r_full

data[sapply(data, is.nan)] <- NA

# ################ 
# # R-WCLS
# ################
# 

# 
data$MOOD_minus = rep(NA, nrow(data))
data$MOOD_pred = rep(NA, nrow(data))

for (i in 1:5){
  
  train = data[data$cv_id !=i,]
  test = data[data$cv_id ==i,]
  
  test_id = CV_id$UserID[which(CV_id$cv_id== i)]
  
  fita_train = glm(week_category_new~ 1, family = binomial(), data = train)
  p_tilde = unique(fita_train[["fitted.values"]])
  
  rf_1 <- randomForest(as.formula(paste0("MOOD ~ ", Control_var)),
                       data = train[train$week_category_new==1,] %>% drop_na(MOOD,all_of(Cov_set)), 
                       importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_1$mse) > 300){
    print("increase number of trees")
    break
  }
  
  rf_0 <- randomForest(as.formula(paste0("MOOD ~ ", Control_var)),
                       data = train[train$week_category_new==0,] %>% drop_na(MOOD,all_of(Cov_set)), 
                       importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_0$mse) > 300){
    print("increase number of trees")
    break
  }
  
  test.pred_1 <- predict(rf_1,test)
  test.pred_0 <- predict(rf_0,test)
  test.pred <- test.pred_1*p_tilde + test.pred_0*(1-p_tilde)
  
  
  data$MOOD_minus[data$UserID %in% test_id] = test$MOOD - test.pred
  data$MOOD_pred[data$UserID %in% test_id] = test.pred
  data$pn[data$UserID %in%test_id] = p_tilde
}

data = data %>% mutate(weights = if_else(week_category_new==1, pn/prob, (1-pn)/(1-prob)) )

data[sapply(data, is.nan)] <- NA

###########
# DR-WCLS
##########

data$MOOD_DR = rep(NA, nrow(data))
data$weights_DR = rep(NA, nrow(data))
data$prob_DR = rep(NA, nrow(data))
data$prob_r = rep(NA, nrow(data))
data$r_weight = rep(NA, nrow(data))

for (i in 1:5){
  
  train = data[data$cv_id !=i,]
  test = data[data$cv_id ==i,]
  
  test_id = CV_id$UserID[which(CV_id$cv_id== i)]
  
  fita_train = glm(week_category_new~ 1, family = binomial(), data = train)
  test$p_tilde = unique(fita_train[["fitted.values"]])
  denom = test$p_tilde * (1-test$p_tilde)
  
  ## random forest on the missing prob
  
  r_train = randomForest(as.formula(paste0("R_factor ~", Control_var)),
                         data = train%>% drop_na(all_of(Cov_set)), importance = TRUE, ntree=300,mtry = 3)
  #
  test$prob_r <- unname(predict(r_train,test, type = "prob")[,"1"])
  
  ## glm on the missing prob
  
  # r_train = glm(R ~ R_cum_avg + R_prev + PHQtot0 + Neu0 + depr0 + EFE0 + pre_intern_mood + pre_intern_sleep +pre_intern_sqrt_step,
  #               family = binomial(), data = train)
  # test$prob_r= predict(r_train,test, type = "response")
  
  test$r_weight = test$R /test$prob_r
  
  rf_1 <- randomForest(as.formula(paste0("MOOD ~ ", Control_var)),
                       data = train[train$week_category_new==1,] %>% drop_na(MOOD,all_of(Cov_set)), 
                       importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_1$mse) > 300){
    print("increase number of trees")
    break
  }
  
  rf_0 <- randomForest(as.formula(paste0("MOOD ~ ", Control_var)),
                       data = train[train$week_category_new==0,] %>% drop_na(MOOD,all_of(Cov_set)), 
                       importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_0$mse) > 300){
    print("increase number of trees")
    break
  }
  
  test.pred_1 <- predict(rf_1,test)
  test.pred_0 <- predict(rf_0,test)
  test.pred <- test.pred_1*test$week_category_new + test.pred_0*(1-test$week_category_new)
  
  test$pseudo = test$weights * test$r_weight* (test$week_category_new - test$p_tilde) * (test$MOOD - test.pred)/denom + test.pred_1 -test.pred_0
  test$pseudo_beta = test.pred_1 -test.pred_0
  
  data$MOOD_DR[data$UserID %in% test_id] = if_else(is.na(test$MOOD),test$pseudo_beta,test$pseudo)
  data$weights_DR[data$UserID %in%test_id] = denom 
  data$prob_DR[data$UserID %in%test_id] = test$prob
  data$prob_r[data$UserID %in%test_id] = test$prob_r
  data$r_weight[data$UserID %in%test_id]  = test$r_weight
  
}
data[sapply(data, is.nan)] <- NA

save(data, file = "data_rfm_mood.RData")

########### 
#Step 
##########

load("Data_missing.RData")

# unequal fold-size

CV_id = data.frame(UserID = unique(data$UserID),
                   cv_id = sample(1:5, length(unique(data$UserID)), replace = TRUE))

data = left_join(data, CV_id,by = "UserID")

data = data %>%
  mutate(R = if_else(is.na(STEP_COUNT),0,1))%>%
  arrange(UserID,study_week)%>%
  group_by(UserID)%>%
  mutate(R_prev =lag(R, n = 1, default = 1),
         R_cum_avg = cummean(R),
         R_factor = as.factor(R))

data$week_category_new = ifelse(data$week_category == "activity",1,0)
data$pn = mean(data$week_category_new)

# data = data %>% group_by(UserID) %>% mutate(week_category_new_lag1 = lag(week_category_new))
# data$week_category_new_lag1[is.na(data$week_category_new_lag1)] <- -1

data = data %>% mutate(week_category_new_c = week_category_new - pn,
                       prob = 1/4,
                       w = if_else(week_category_new==1, pn/prob, (1-pn)/(1-prob)) )

Cov_set = colnames(data)[c(7:17)]
Cov_set = c(Cov_set,"study_week","R_cum_avg","R_prev")
Control_var = paste0(Cov_set, collapse = " + ",sep = "")

r_train = randomForest(as.formula(paste0("R_factor ~", Control_var)),
                       data = data%>% drop_na(all_of(Cov_set)), importance = TRUE, ntree=300,mtry = 3)
#
data$prob_r_full <- unname(predict(r_train,data, type = "prob")[,"1"])
data$r_w = data$R / data$prob_r_full

data[sapply(data, is.nan)] <- NA

################ 
# R-WCLS
################

data$STEP_COUNT_minus = rep(NA, nrow(data))
data$STEP_COUNT_pred = rep(NA, nrow(data))

for (i in 1:5){
  
  train = data[data$cv_id !=i,]
  test = data[data$cv_id ==i,]
  
  test_id = CV_id$UserID[which(CV_id$cv_id== i)]
  
  fita_train = glm(week_category_new~ 1, family = binomial(), data = train)
  p_tilde = unique(fita_train[["fitted.values"]])
  
  #### tuning parameters
  
  # bestmtry <- tuneRF(train[,Cov_set],train[,"MOOD"], stepFactor=1.5, improve=1e-5, ntree=300)
  # print(bestmtry)
  
  #### 
  
  rf_1 <- randomForest(as.formula(paste0("STEP_COUNT ~ ", Control_var)),
                       data = train[train$week_category_new==1,] %>% drop_na(STEP_COUNT,all_of(Cov_set)), 
                       importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_1$mse) > 300){
    print("increase number of trees")
    break
  }
  
  rf_0 <- randomForest(as.formula(paste0("STEP_COUNT ~ ", Control_var)),
                       data = train[train$week_category_new==0,] %>% drop_na(STEP_COUNT,all_of(Cov_set)), 
                       importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_0$mse) > 300){
    print("increase number of trees")
    break
  }
  
  test.pred_1 <- predict(rf_1,test)
  test.pred_0 <- predict(rf_0,test)
  test.pred <- test.pred_1*p_tilde + test.pred_0*(1-p_tilde)
  
  
  data$STEP_COUNT_minus[data$UserID %in% test_id] = test$STEP_COUNT - test.pred
  data$STEP_COUNT_pred[data$UserID %in% test_id] = test.pred
  data$pn[data$UserID %in%test_id] = p_tilde
}

data = data %>% mutate(weights = if_else(week_category_new==1, pn/prob, (1-pn)/(1-prob)) )

###########
# DR-WCLS
##########

data$STEP_COUNT_DR = rep(NA, nrow(data))
data$weights_DR = rep(NA, nrow(data))
data$prob_DR = rep(NA, nrow(data))
data$prob_r = rep(NA, nrow(data))
data$r_weight = rep(NA, nrow(data))

for (i in 1:5){
  
  train = data[data$cv_id !=i,]
  test = data[data$cv_id ==i,]
  
  test_id = CV_id$UserID[which(CV_id$cv_id== i)]
  
  fita_train = glm(week_category_new~ 1, family = binomial(), data = train)
  test$p_tilde = unique(fita_train[["fitted.values"]])
  denom = test$p_tilde * (1-test$p_tilde)
  
  r_train = randomForest(as.formula(Cov_set = colnames(data)[c(7:17)]),
                         data = train%>% drop_na(all_of(Cov_set)), importance = TRUE, ntree=300,mtry = 3)
  #
  test$prob_r <- unname(predict(r_train,test, type = "prob")[,"1"])
  
  # r_train = glm(R ~ R_cum_avg + R_prev + PHQtot0 + Neu0 + depr0 + EFE0 + pre_intern_mood + pre_intern_sleep +pre_intern_sqrt_step,
  #               family = binomial(), data = train)
  # test$prob_r= predict(r_train,test, type = "response")
  # 
  
  test$r_weight = test$R /test$prob_r
  
  rf_1 <- randomForest(as.formula(paste0("STEP_COUNT ~ ", Control_var)),
                       data = train[train$week_category_new==1,] %>% drop_na(STEP_COUNT,all_of(Cov_set)), 
                       importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_1$mse) > 300){
    print("increase number of trees")
    break
  }
  
  rf_0 <- randomForest(as.formula(paste0("STEP_COUNT ~ ", Control_var)),
                       data = train[train$week_category_new==0,] %>% drop_na(STEP_COUNT,all_of(Cov_set)), 
                       importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_0$mse) > 300){
    print("increase number of trees")
    break
  }
  
  test.pred_1 <- predict(rf_1,test)
  test.pred_0 <- predict(rf_0,test)
  test.pred <- test.pred_1*test$week_category_new + test.pred_0*(1-test$week_category_new)
  
  test$pseudo = test$weights * test$r_weight* (test$week_category_new - test$p_tilde) * (test$STEP_COUNT - test.pred)/denom + test.pred_1 -test.pred_0
  test$pseudo_beta = test.pred_1 -test.pred_0
  
  
  data$STEP_COUNT_DR[data$UserID %in% test_id] = if_else(is.na(test$STEP_COUNT),test$pseudo_beta,test$pseudo)
  data$weights_DR[data$UserID %in%test_id] = denom 
  data$prob_DR[data$UserID %in%test_id] = test$prob
  data$prob_r[data$UserID %in%test_id] = test$prob_r
  data$r_weight[data$UserID %in%test_id]  = test$r_weight
  
}

# data$r_weight[is.infinite(data$r_weight)] <- 0

data[sapply(data, is.nan)] <- NA
save(data, file = "data_rfm_step.RData")