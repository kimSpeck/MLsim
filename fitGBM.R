# gbm fitting

# original gbm package (https://github.com/cran/gbm)
#   is retired/no longer under active development
#   but functions from the gbm-package are used inside of caret
#   see line 68: https://github.com/topepo/caret/blob/master/models/files/gbm.R

# xgboost
#   approx. 10x faster than gbm (implemented in C++) + still under development
#   (+ options for regularised gbms to prevent overfitting)
#   https://github.com/dmlc/xgboost/tree/master/R-package

# parameter information
# https://xgboost.readthedocs.io/en/latest/parameter.html#learning-task-parameters
# install.packages('xgboost')
library(xgboost)

# only works for matrices that contain all numeric variables
#   dummy coding for categorical/binary variables beforehand!
# scaling and dummy coding is not necessary (so far) because we only have continuous
#     predictors on identical scales

# train model
# input data with xgboosts own class (recommended)
# https://xgboost.readthedocs.io/en/stable/R-package/xgboostPresentation.html
trainData <- xgboost::xgb.DMatrix(data = Xtrain, label = ytrain, missing = NA) 
# xgboost::setinfo(x, "label", y)

# (tuning) parameter grid
# number of trees, depth of trees, learning rate, subsampling?
tuneGrid <- expand.grid(eta = seq(.001, .201, .02),              # shrinkage/learning rate
                        max_depth = c(1,2,3,4,5),                # tree depth; interaction depth hängt von simulation ab
                        min_child_weight = c(5,10,20,50),        # min node size
                        # subsample? # percent of training data to sample fpr each tree?/bag.fraction in gbm?
                        # colsample_bytree? # percent of columns to sample from for each tree
                        # nTrees = c(50,100,150,300,500,1000),     # max number of trees to run
                        # does max number of trees become irrelevant to tune as soon as we introduce early stopping?
                        optimalTrees = NA,                        # fill in cv with nTrees for each grid line    
                        minRMSE = NA)                             # fill in cv with RMSE for each grid line

# run tuning grid search 
for (iGline in seq_len(nrow(tuneGrid))) {

  # reproducibility
  set.seed(123)
  
  # train model
  fit_cv <- xgb.cv(
    data = trainData,                 # X
    # label = response_train,         # y
    nrounds = 1000,                   # max nTrees
    nfold = setParam$fit$nfolds,      # k-fold cross validation
    eta = tuneGrid$eta[iGline],
    max_depth = tuneGrid$max_depth[iGline],
    min_child_weight = tuneGrid$min_child_weight[iGline],
    # subsample = 1,                  # percent of training data to sample fpr each tree?/bag.fraction in gbm?
    # colsample_bytrees: percent of columns to sample from for each tree
    objective = "reg:squarederror",   # for regression models
    # eval_metric = "rmse",           # default for regression models 
    verbose = FALSE,                  # silent
    early_stopping_rounds = 10)       # stop running if the cross validated error 
                                      # does not improve for n continuous trees
  
  # identify the optimal number of trees and the corresponding minimum RMSE for the cv 
  # add number of Trees and corresponding training error to grid
  tuneGrid$optimalTrees[iGline] <- which.min(fit_cv$evaluation_log$test_rmse_mean)
  tuneGrid$minRMSE[iGline] <- min(fit_cv$evaluation_log$test_rmse_mean)
  
}

# cross-validated tuning parameters
idxOptim <- which.min(tuneGrid$minRMSE)

# to do: make tuned parameters ready to use in model
# to do: save tuned parameters in output

# fit model
fit <- xgboost(
  data = trainData,                 # X
  # label = response_train,         # y
  eta = tuneGrid$eta[idxOptim],
  max_depth = tuneGrid$max_depth[idxOptim],
  min_child_weight = tuneGrid$min_child_weight[idxOptim],
  nrounds = tuneGrid$optimalTrees[idxOptim],
  objective = "reg:squarederror",   # for regression models
  # eval_metric = "rmse",           # default for regression models 
  verbose = FALSE)

# performance test for train data
predTrain <- predict(fit, Xtrain)
performTrain <- evalPerformance(predTrain, ytrain)
performTrain <- matrix(performTrain, ncol = length(performTrain), nrow = 1)

# test data
# performance test for test data
# get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model testing
performTest <- lapply(paste0("test", seq_len(setParam$dgp$nTest)), function(iTest) {
  Xtest <- as.matrix(data[[iTest]][["X_int"]])
  ytest <- data[[iTest]][["yMat"]][,iCond]
  pred <- predict(fit, Xtest)  
  
  evalPerformance(pred, ytest)
})
performTest <- do.call(rbind, performTest)
