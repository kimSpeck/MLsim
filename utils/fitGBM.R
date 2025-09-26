fitGBM <- function(Xtrain, ytrain, Xtest, ytest, setParam, iSample, explanation = FALSE) {
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
                        method = "gbm", #gradient boosting machines
                        metric = "RMSE",
                        #model training using ten-fold cross-validation
                        trControl = trainControl(method="cv", 
                                                 number=setParam$fit$nfolds,
                                                 selectionFunction = "best"), 
                        tuneGrid = setParam$fit$tuneGrid_GBM, #setting the tuning grid
                        verbose = FALSE)
  
  # R² across all test sets in the CV procedure (for the final hyperparameters)
  performCVtest <- model$results[rownames(model$results) == rownames(model$bestTune), 
                                 c("RMSE", "Rsquared", "MAE")]
  
  # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model training
  # performance test for train data
  predTrain <- predict(model, Xtrain)
  performTrain <- evalPerformance(predTrain, ytrain)
  # performTrain <- matrix(performTrain, ncol = length(performTrain), nrow = 1)
  # colnames(performTrain) <- c("RMSE", "Rsquared", "MAE")
  
  # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model testing
  # test data
  # performance test for test data
  pred <- predict(model, Xtest)  
  performTest <- evalPerformance(pred, ytest)
  
  if (setParam$fit$explanation) {
    # multiple dependent measures from the iml package
    predictor.gbm <- Predictor$new(
      model = model, 
      data = data.frame(Xtrain), 
      y = ytrain)
    
    # permutation variable importance
    #       increase n.repetitions? (How often should the shuffling of the feature 
    #       be repeated? The higher the number of repetitions the more stable 
    #       and accurate the results become.)
    imp.gbm <- FeatureImp$new(predictor.gbm, loss = "mse")
    pviRank <- imp.gbm$results$feature
    pviValue <- imp.gbm$results$importance
  } else {
    pviRank <- NA
    pviValue <- NA
  }
  
  # evaluate interaction strength with H-statistic
  if (setParam$fit$InterStrength) {
    interStrength <- lapply(paste0("Var", seq_len(setParam$dgp$p)), function(iVar) {
      tmp <- iml::Interaction$new(predictor.gbm, feature = iVar,
                                  grid.size = setParam$fit$nInterStrength)
      idxTmp <- sort(tmp$results$.interaction, decreasing = T)
      cbind(var = rep(iVar, length(idxTmp)),
            feature = tmp$results$.feature[match(idxTmp, tmp$results$.interaction)],
            interaction = tmp$results$.interaction[match(idxTmp, tmp$results$.interaction)])
      
    })
    interStrength <- do.call(rbind, interStrength)
  }
  
  
  # save dependent variables for each sample (in a list)
  return(list(
    performTrain = performTrain, # train performance (full train sample)
    performTest = performTest, # (holdout) test performance
    performCVtest = performCVtest, # within CV test performance
    # model agnostic measures
    pvi = cbind(pviRank, pviValue), # permutation variable importance
    # h-statistic two-way interactions
    interStrength = if (setParam$fit$InterStrength) interStrength else NA,
    # cross-validated tuning parameters (GBM)
    tunedShrinkage = model$finalModel$shrinkage,
    tunedMax_depth = model$finalModel$interaction.depth,
    tunedMin_child_weight =  model$finalModel$n.minobsinnode,
    tunedNrounds = model$finalModel$n.trees,
    iSample = iSample))
}