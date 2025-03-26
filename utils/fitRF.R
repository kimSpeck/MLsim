fitRF <- function(Xtrain, ytrain, Xtest, ytest, setParam) {
  # # test it
  # Xtrain <- data[["1"]][["X_int"]]
  # Xtest <- data[["test1"]][["X_int"]]
  # ytrain <- data[[1]][["yMat"]][,"R20.8lin_inter0.5_0.5"]
  # ytest <- data[["test1"]][["yMat"]][,"R20.8lin_inter0.5_0.5"]
  # 
  # # train data
  # idx_rmInter.train <- stringr::str_detect(colnames(Xtrain), ":")
  # Xtrain <- Xtrain[,colnames(Xtrain)[!idx_rmInter.train]]
  # # test data
  # idx_rmInter.test <- stringr::str_detect(colnames(Xtest), ":")
  # Xtest <- Xtest[,colnames(Xtest)[!idx_rmInter.test]]
  
  # combine predictors with outcome
  trainData <- cbind(Xtrain, ytrain)
  
  # cross validation for training model
  model <- caret::train(ytrain ~ .,
                        trainData,
                        method = "ranger", # random forest
                        metric = "RMSE",
                        num.trees = setParam$fit$numTreesRF,
                        #model training using ten-fold cross-validation
                        trControl = trainControl(method="cv",
                                                 number=setParam$fit$nfolds,
                                                 selectionFunction = "oneSE"),
                        tuneGrid = setParam$fit$tuneGrid_RF, #setting the tuning grid
                        verbose = FALSE)
  
  # save tuning parameters! (results from cross validation)
  # cross-validated tuning parameters
  tunedMtry <- model$finalModel$mtry
  tunedSplitrule <- model$finalModel$splitrule
  tunedMinNodeSize <- model$finalModel$min.node.size
  
  # R² across all test sets in the CV procedure (for the final hyperparameters)
  performCVtest <- model$results[rownames(model$results) == rownames(model$bestTune), 
                                 c("RMSE", "Rsquared", "MAE")]

  # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model training
  # performance test for train data
  # each observation is predicted using all trees, rather than only out-of-bag trees
  predTrain <- predict(model, Xtrain) 
  performTrain <- evalPerformance(predTrain, ytrain)
  # performTrain <- matrix(performTrain, ncol = length(performTrain), nrow = 1)
  # colnames(performTrain) <- c("RMSE", "Rsquared", "MAE")
  
  
  # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model testing
  # test data
  # performance test for test data
  pred <- predict(model, Xtest)  
  performTest <- evalPerformance(pred, ytest)
  
  # save dependent variables for each sample (in a list)
  return(list(
    performTrain = performTrain, # train performance (full train sample)
    performTest = performTest, # test performance
    performCVtest = performCVtest, # within CV test performance
    # oob predictions and R2
    oobPredictions = model[["finalModel"]][["predictions"]], # the OOB predictions
    oobR2 = model$finalModel$r.squared, # oob R²
    # tuning parameters (RF)
    tunedMtry = tunedMtry, 
    tunedSplitrule = tunedSplitrule,
    tunedMinNodeSize = tunedMinNodeSize))
}
