library(randomForest)
library(glmnet)

ML-fitter_cv = function(d,group_ls, ml_method = "random forest"){
  
  runin <- runin.fity
  r <- which(d$time >= runin)
  d = d[r,]
  
  group = group_str(group_ls)
  beta = rep(0,group[["groups"]])
  
  for (k in 1:group[["groups"]]){
    train_id = which(group[["group_id"]]!= k)
    train = d[d$id %in%train_id, ]
    
    test_id = which(group[["group_id"]]== k)
    test = d[d$id %in%test_id, ]
    
    if(ml_method == "random forest"){
      rf_weight =  ifelse(train[, "a"] == 1, train[, "pn"]/ train[, "prob"],
                          (1 - train[, "pn"]) / (1 - train[ ,"prob"]))
      
      rf <- randomForest(as.formula(paste0("y ~ state + ", Control_var)), 
                         data = train, importance = TRUE, ntree=500, weights = rf_weight)
      
      if(which.min(rf$mse) > 500){
        print("increase number of trees")
        break
      }
      
      test.pred <- predict(rf,test)
    }else{
      train_x = model.matrix(as.formula(paste0("y ~ state + ", Control_var)), data = train)
      test_x = model.matrix(as.formula(paste0("y ~ state + ", Control_var)), data = test)
      
      glm_net = glmnet(x = train_x[,-1], y = train[, "y"])
      test.pred = predict(glm_net, newx = test_x[,-1], s = glm_net[["lambda"]][glm_net[["dim"]][2]])
    }
    
    test$y_minus = test$y - test.pred
    
    if (!is.null(a.formula)){
      l = list(x = model.matrix(as.formula(a ~ 1), data = test), y = test[, "a"])
      fita <- do.call(glm.fit,l )
    }
    
    l <- list(x = model.matrix(as.formula(y_minus ~ I(a - pn)), data = test), y = test[, "y_minus"])
    
    test[ , args[["w"]][["wn"]]] = fita[["fitted.values"]]
    l$w <- ifelse(test[, "a"] == 1, test[, "pn"]/ test[, "prob"],
                  (1 - test[, "pn"]) / (1 - test[ ,"prob"]))
    
    fun <- "lm.wfit"
    fit <- do.call(fun, l)
    
    beta[k] = fit[["coefficients"]][["I(a - pn)"]]
  }
  beta_mean = mean(beta)
  
}


ML-fitter = function(d,group_ls, ml_method = "random forest"){
  
  runin <- runin.fity
  r <- which(d$time >= runin)
  d = d[r,]
  d$y_minus = rep(NA, nrow(d))
  
  group = group_str(group_ls)
  
  for (k in 1:group[["groups"]]){
    train_id = which(group[["group_id"]]!= k)
    train = d[d$id %in%train_id, ]
    
    test_id = which(group[["group_id"]]== k)
    test = d[d$id %in%test_id, ]
    
    if(ml_method == "random forest"){
      rf_weight =  ifelse(train[, "a"] == 1, train[, "pn"]/ train[, "prob"],
                          (1 - train[, "pn"]) / (1 - train[ ,"prob"]))
      
      rf <- randomForest(as.formula(paste0("y ~ state + ", Control_var)), 
                         data = train, importance = TRUE, ntree=500, weights = rf_weight)
      
      if(which.min(rf$mse) > 500){
        print("increase number of trees")
        break
      }
      
      test.pred <- predict(rf,test)
    }else{
      train_x = model.matrix(as.formula(paste0("y ~ state + ", Control_var)), data = train)
      test_x = model.matrix(as.formula(paste0("y ~ state + ", Control_var)), data = test)
      
      glm_net = glmnet(x = train_x[,-1], y = train[, "y"])
      test.pred = predict(glm_net, newx = test_x[,-1], s = glm_net[["lambda"]][glm_net[["dim"]][2]])
    }

    d$y_minus[d$id %in%test_id] = test$y - test.pred
    
  }
  
  if (!is.null(a.formula)){
    l = list(x = model.matrix(as.formula(a ~ 1), data = d), y = d[, "a"])
    fita <- do.call(glm.fit,l )
  }
  
  l <- list(x = model.matrix(as.formula(y_minus ~ I(a - pn)), data = d), y = d[, "y_minus"])
  
  d[ , "pn"] = fita[["fitted.values"]]
  l$w <- ifelse(d[, "a"] == 1, d[, "pn"]/ d[, "prob"],
                (1 - d[, "pn"]) / (1 - d[ ,"prob"]))
  
  fun <- "lm.wfit"
  fit <- do.call(fun, l)
  
  beta = fit[["coefficients"]][["I(a - pn)"]]
  
}





  # Create a random forest with 1000 trees


# How many trees are needed to reach the minimum error estimate? 
# This is a simple problem; it appears that about 100 trees would be enough. 


#   # Plot rf to see the estimated error as a function of the number of trees
#   # (not running it here)
# plot(rf) 
#       
#         # Using the importance()  function to calculate the importance of each variable
# imp <- as.data.frame(sort(importance(rf)[,1],decreasing = TRUE),optional = T)
# names(imp) <- "% Inc MSE"
# imp
# 
# 
#   # As usual, predict and evaluate on the test set
# test.pred.forest <- predict(rf,test)
# RMSE.forest <- sqrt(mean((test.pred.forest-test$y)^2))
# RMSE.forest

# RMSE.lasso <- sqrt(mean((test.pred-test$y)^2))
# RMSE.lasso

#  
# 
# MAE.forest <- mean(abs(test.pred.forest-test$y))
# MAE.forest
