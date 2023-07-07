# library(geepack)
library(dplyr)
library(ggplot2)
library(splines)
# library(foreach)
# library(parallel)
# library(doParallel)
library(randomForest)
setwd("~/Documents/DML-WCLS/Case Study")
load("IHS_MRT.RData")
library(stringr)
data = IHS_MRT[[9]]

# unequal fold-size

CV_id = data.frame(UserID = unique(data$UserID),
                   cv_id = sample(1:5, length(unique(data$UserID)), replace = TRUE))


# equal fold size

# CV_id = data.frame(UserID = unique(data$UserID),
#   cv_id = rep(1:3, length.out =length(unique(data$UserID)) ))

# table(CV_id$cv_id)
data = left_join(data, CV_id,by = "UserID")

data$week_category_new = ifelse(data$week_category == "mood",1,0)
data$pn = mean(data$week_category_new)

# data = data %>% group_by(UserID) %>% mutate(week_category_new_lag1 = lag(week_category_new))
# data$week_category_new_lag1[is.na(data$week_category_new_lag1)] <- -1

data = data %>% mutate(week_category_new_c = week_category_new - pn,
                       prob = 1/4,
                       w = if_else(week_category_new==1, pn/prob, (1-pn)/(1-prob)) )

Cov_set = colnames(data)[c(7:17)]

# Cov_set = c("MOODprev","pre_intern_mood","Neu0","depr0")
# Control_var = paste0(Cov_set, collapse = " + ",sep = "")

################ 
# R-WCLS
################

Cov_set = c(Cov_set,"study_week")
Control_var = paste0(Cov_set, collapse = " + ",sep = "")

data$MOOD_minus = rep(NA, nrow(data))
data$MOOD_pred = rep(NA, nrow(data))

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
  
  rf_1 <- randomForest(as.formula(paste0("MOOD ~ ", Control_var)),
                       data = train[train$week_category_new==1,], importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_1$mse) > 300){
    print("increase number of trees")
    break
  }
  
  rf_0 <- randomForest(as.formula(paste0("MOOD ~ ", Control_var)),
                       data = train[train$week_category_new==0,], importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_0$mse) > 300){
    print("increase number of trees")
    break
  }
  
  test.pred_1 <- predict(rf_1,test)
  test.pred_0 <- predict(rf_0,test)
  test.pred <- test.pred_1*p_tilde + test.pred_0*(1-p_tilde)
  
  
  data$MOOD_minus[data$UserID %in% test_id] = test$MOOD - test.pred
  data$MOOD_pred[data$UserID %in% test_id] =  test.pred
  data$pn[data$UserID %in%test_id] = p_tilde
}

data = data %>% mutate(weights = if_else(week_category_new==1, pn/prob, (1-pn)/(1-prob)) )

###########
# DR-WCLS
##########

data$MOOD_DR = rep(NA, nrow(data))
data$weights_DR = rep(NA, nrow(data))
data$prob_DR = rep(NA, nrow(data))

for (i in 1:5){
  
  train = data[data$cv_id !=i,]
  test = data[data$cv_id ==i,]
  
  test_id = CV_id$UserID[which(CV_id$cv_id== i)]
  
  
  fita_train = glm(week_category_new~ 1, family = binomial(), data = train)
  test$p_tilde = unique(fita_train[["fitted.values"]])
  denom = test$p_tilde * (1-test$p_tilde)
  
  # train$week_category_new = as.factor(train$week_category_new)
  # fita_train = randomForest(as.formula(paste0("week_category_new ~ week_category_new_lag1 + ", Control_var)),
  #                           data = train, importance = TRUE, ntree=300,mtry = 3)
  # 
  # test$prob <- unname(predict(fita_train,test, type = "prob")[,"1"])
  # 
  # 
  # test = test %>% 
  #   mutate(weights = if_else(week_category_new==1, p_tilde/prob, (1-p_tilde)/prob))%>%
  #   mutate(weights = if_else(is.finite(weights), weights, 0))
  
  
  #### 
  rf_1 <- randomForest(as.formula(paste0("MOOD ~ ", Control_var)),
                       data = train[train$week_category_new==1,], importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_1$mse) > 300){
    print("increase number of trees")
    break
  }
  
  rf_0 <- randomForest(as.formula(paste0("MOOD ~ ", Control_var)),
                       data = train[train$week_category_new==0,], importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_0$mse) > 300){
    print("increase number of trees")
    break
  }
  
  test.pred_1 <- predict(rf_1,test)
  test.pred_0 <- predict(rf_0,test)
  test.pred <- test.pred_1*test$week_category_new + test.pred_0*(1-test$week_category_new)
  
  
  data$MOOD_DR[data$UserID %in% test_id] = test$weights * (test$week_category_new - test$p_tilde) * (test$MOOD - test.pred)/denom + test.pred_1 -test.pred_0
  data$weights_DR[data$UserID %in%test_id] = denom 
  data$prob_DR[data$UserID %in%test_id] = test$prob
  
}

# data$MOOD_DR[data$MOOD_DR == "NaN"] <- NA

save(data, file = "data_rfa_mood.RData")


########### 
#Step 
##########
data = IHS_MRT[[9]]

# unequal fold-size

CV_id = data.frame(UserID = unique(data$UserID),
                   cv_id = sample(1:5, length(unique(data$UserID)), replace = TRUE))

data = left_join(data, CV_id,by = "UserID")

data$week_category_new = ifelse(data$week_category == "activity",1,0)
data$pn = mean(data$week_category_new)

# data = data %>% group_by(UserID) %>% mutate(week_category_new_lag1 = lag(week_category_new))
# data$week_category_new_lag1[is.na(data$week_category_new_lag1)] <- -1

data = data %>% mutate(week_category_new_c = week_category_new - pn,
                       prob = 1/4,
                       w = if_else(week_category_new==1, pn/prob, (1-pn)/(1-prob)) )

Cov_set = colnames(data)[c(7:17)]

# Cov_set = c("MOODprev","pre_intern_mood","Neu0","depr0")
# Control_var = paste0(Cov_set, collapse = " + ",sep = "")

################ 
# R-WCLS
################

Cov_set = c(Cov_set,"study_week")
Control_var = paste0(Cov_set, collapse = " + ",sep = "")

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
                       data = train[train$week_category_new==1,], importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_1$mse) > 300){
    print("increase number of trees")
    break
  }
  
  rf_0 <- randomForest(as.formula(paste0("STEP_COUNT ~ ", Control_var)),
                       data = train[train$week_category_new==0,], importance = TRUE, ntree=300,mtry = 3)
  
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


for (i in 1:5){
  
  train = data[data$cv_id !=i,]
  test = data[data$cv_id ==i,]
  
  test_id = CV_id$UserID[which(CV_id$cv_id== i)]
  
  fita_train = glm(week_category_new~ 1, family = binomial(), data = train)
  test$p_tilde = ifelse(test$prob == 0, 0, unique(fita_train[["fitted.values"]]))
  denom = test$p_tilde * (1-test$p_tilde)
  
  # train$week_category_new = as.factor(train$week_category_new)
  # fita_train = randomForest(as.formula(paste0("week_category_new ~ week_category_new_lag1 + ", Control_var)),
  #                           data = train, importance = TRUE, ntree=300,mtry = 3)
  # 
  # test$prob <- unname(predict(fita_train,test, type = "prob")[,"1"])
  # 
  # test = test %>% 
  #   mutate(weights = if_else(week_category_new==1, p_tilde/prob, (1-p_tilde)/prob))%>%
  #   mutate(weights = if_else(is.finite(weights), weights, 0))
  
  
  
  rf_1 <- randomForest(as.formula(paste0("STEP_COUNT ~ ", Control_var)),
                       data = train[train$week_category_new==1,], importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_1$mse) > 300){
    print("increase number of trees")
    break
  }
  
  rf_0 <- randomForest(as.formula(paste0("STEP_COUNT ~ ", Control_var)),
                       data = train[train$week_category_new==0,], importance = TRUE, ntree=300,mtry = 3)
  
  if(which.min(rf_0$mse) > 300){
    print("increase number of trees")
    break
  }
  
  test.pred_1 <- predict(rf_1,test)
  test.pred_0 <- predict(rf_0,test)
  test.pred <- test.pred_1*test$week_category_new + test.pred_0*(1-test$week_category_new)
  
  
  data$STEP_COUNT_DR[data$UserID %in% test_id] = test$weights * (test$week_category_new - test$p_tilde) * (test$STEP_COUNT - test.pred)/denom + test.pred_1 -test.pred_0
  data$weights_DR[data$UserID %in%test_id] = denom 
  data$prob_DR[data$UserID %in%test_id] = test$prob
  
}

# data$STEP_COUNT_DR[data$STEP_COUNT_DR == "NaN"] <- NA

save(data, file = "data_rf_step.RData")

