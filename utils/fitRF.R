fitRF <- function(Xtrain, ytrain, Xtest, ytest, setParam, iSample) {
  # # test it
  # iSim =  7 # inter, RF, 100, 10, 0.6, 451935
  # iSample = sample(size = 1, 1:1000)
  # iCond = sample(size = 1, 1:9)
  # 
  # fileName <- paste0(condGrid[iSim, "data"], 
  #                    "/simDataN", condGrid[iSim, "N"],
  #                    "_pTrash", condGrid[iSim, "pTrash"],
  #                    "_rel", condGrid[iSim, "reliability"])
  # 
  # load(paste0(dataFolder, "/", fileName, "/sample_", iSample, ".rda"))
  # 
  # # get predictor variable matrix (identical across all R2 x lin_inter conditions)
  # Xtrain <- as.matrix(dataList[["X_int"]]) 
  # # get criterion/target variable (column in yMat depending on R2 x lin_inter condition)
  # ytrain <- dataList[["yMat"]][,iCond]
  # 
  # # load test data
  # # by doing this in advane condGrid and full data not an additional argument for model functions
  # # ... but, only if number of test samples is 1!
  # # setParam$dgp$nTest == 1
  # load(paste0(dataFolder, "/", fileName, "/sample_test1.rda"))
  # 
  # Xtest <- as.matrix(dataList[["X_int"]]) 
  # ytest <- dataList[["yMat"]][,iCond] 
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
                                                 selectionFunction = "best"),
                        tuneGrid = setParam$fit$tuneGrid_RF, #setting the tuning grid
                        verbose = FALSE)

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
    # cross-validated tuning parameters
    tunedMtry = model$finalModel$mtry, 
    tunedSplitrule = model$finalModel$splitrule,
    tunedMinNodeSize = model$finalModel$min.node.size, 
    iSample = iSample))
}
