fitGBM <- function(Xtrain, ytrain, Xtest, ytest, setParam) {
  
  # combine predictors with outcome
  trainData <- cbind(Xtrain, ytrain)
  
  # cross validation for training model
  model <- caret::train(ytrain ~ .,
                        trainData,
                        method = "gbm", #gradient boosting machines
                        metric = "RMSE",
                        #model training using ten-fold cross-validation
                        trControl = trainControl(method="cv", 
                                                 number=setParam$fit$nfolds), 
                        tuneGrid = setParam$fit$tuneGrid, #setting the tuning grid
                        verbose = FALSE)
  
  # save tuning parameters! (results from cross validation)
  # cross-validated tuning parameters
  tunedShrinkage <- model$finalModel$shrinkage
  tunedMax_depth <- model$finalModel$interaction.depth
  tunedMin_child_weight <- model$finalModel$n.minobsinnode
  tunedNrounds <- model$finalModel$n.trees
  
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
    performTrain = performTrain, # train performance
    performTest = performTest, # test performance
    # model agnostic measures
    pvi = cbind(pviRank, pviValue), # permutation variable importance
    # h-statistic two-way interactions
    if (setParam$fit$InterStrength) {interStrength = interStrength 
    } else {interStrength = NA},
    # tuning parameters (GBM)
    tunedShrinkage = tunedShrinkage,
    tunedMax_depth = tunedMax_depth,
    tunedMin_child_weight = tunedMin_child_weight,
    tunedNrounds = tunedNrounds))
}