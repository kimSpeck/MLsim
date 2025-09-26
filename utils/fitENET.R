fitENET <- function(Xtrain, ytrain, Xtest, ytest, setParam, iSample, explanation = FALSE) {
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
  
  # cross-validation on full sample ("training")
  #     ... alpha (elastic net) & lambda-grid (overall penalty strength) for fitting Enet
  
  # chose lambda grid based on glmnet for alpha = 0.5 (i.e., warm start option in caret)
  # https://github.com/topepo/caret/blob/5f4bd2069bf486ae92240979f9d65b5c138ca8d4/models/files/glmnet.R#L15C23-L23C58
  if (setParam$fit$warmStart) {
    init <- glmnet::glmnet(x = Xtrain, 
                           y = ytrain,
                           family = "gaussian",
                           nlambda = length(setParam$fit$lambda)+2, # + 2 lambdas removed 
                           alpha = .5)
    
    lambdaWS <- unique(init$lambda)
    lambdaWS <- lambdaWS[-c(1, length(lambdaWS))] # remove first and last lambda
    lambdaWS <- lambdaWS[1:min(length(lambdaWS), length(setParam$fit$lambda))]
  }
  
  # lambda values come from warm start procedure if setParam$fit$warmStart == TRUE; 
  #   otherwise lambda vector from parameter set are used
  lambdaValues <- if (setParam$fit$warmStart) lambdaWS else setParam$fit$lambda
  
  # tune alpha by iterating trough alphas
  set.seed(89101) # seed for foldid
  foldid <- sample(rep(seq_len(setParam$fit$nfolds), length.out = length(ytrain)), replace = FALSE)
  
  fit_cv <- lapply(setParam$fit$alpha, function(iAlpha) {
    # within iteration through different alpha values tune lambda via cv.glmnet function
    tmp_fit <- cv.glmnet(x = Xtrain, 
                         y = ytrain, 
                         foldid = foldid, # as suggested for alpha tuning via cv
                         alpha = iAlpha, # iteration through different alpha values
                         lambda = lambdaValues, # tune lambda within cv.glmnet
                         family = "gaussian", 
                         standardize = TRUE, 
                         nfolds = setParam$fit$nfolds, # 10 fold cross validation  
                         type.measure = "mse") 
    
    # cv results that we need for parameter tuning
    #     extract mse and lambda for lambda criterions according to parameter list
    #     setParam$fit$lambdaCrit = {1se, min} or any subset
    tmpLambdaCV <- sapply(seq_along(setParam$fit$lambdaCrit), function(iCrit) {
      idxLambdaCrit <- match(setParam$fit$lambdaCrit[iCrit], rownames(tmp_fit$index)) # idx according to parameter (min or se)
      idxLambda <- tmp_fit$index[idxLambdaCrit] # optimal lambda index for lambda criterion (min or se)
      tmp_mse <- tmp_fit$cvm[idxLambda] # choose MSE value for lambda criterion (min or se)
      tmp_lambda <- tmp_fit[[paste0("lambda.", setParam$fit$lambdaCrit[iCrit])]]   
      c(cvm = tmp_mse, 
        lambda = tmp_lambda)
    })
    
    colnames(tmpLambdaCV) <- setParam$fit$lambdaCrit
    tmpLambdaCV
  })
  
  # choose alpha & lambda based on cross validation tuning of lambda given 
  #     specific alpha values results
  #   -> mse from fit_cvs to choose alpha since lambda in cv is chosen based on mse as well
  tuneParam <- lapply(seq_along(fit_cv), function(iAlpha) {
    t(sapply(setParam$fit$lambdaCrit, function(iCrit) {
      tmp_mse <- fit_cv[[iAlpha]]["cvm", iCrit] # MSE vector
      tmp_lambda <- fit_cv[[iAlpha]]["lambda", iCrit] # rather conservative
      c(alpha = setParam$fit$alpha[iAlpha],
        MSE = tmp_mse,
        lambda = tmp_lambda)
    }))
  })
  
  tuneParam <- do.call(rbind, tuneParam)
  
  # for every lambda criterion applied find the respective optimum, i.e., the 
  #     minimum MSE and return alpha & lambda for this optimum
  tunedParams <- sapply(setParam$fit$lambdaCrit, function(iCrit) {
    
    idxCrit <- row.names(tuneParam) == iCrit
    # ! if there are multiple optima: which optimum should we choose?
    # in the case of a tie the tied variable with lowest index is selected.
    idxOptim <- which(tuneParam[idxCrit,"MSE"] == min(tuneParam[idxCrit,"MSE"]))[1]
    
    # Best parameters
    c(tunedAlpha = tuneParam[idxCrit,][idxOptim, "alpha"],
      tunedLambda = tuneParam[idxCrit,][idxOptim, "lambda"])
  })
  
  # ridge: alpha = 0; lasso: alpha = 1
  # here!!! avoid refitting
  fit <- glmnet(x = Xtrain,
                y = ytrain, 
                alpha = tunedParams["tunedAlpha", "min"],
                lambda = tunedParams["tunedLambda", "min"],
                family = "gaussian", 
                standardize = TRUE)
  
  estBeta <- as.matrix(fit$beta)
  selectedVars <- estBeta != 0
  
  # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model training
  predTrain <- predict(fit, Xtrain)
  performTrain <- evalPerformance(predTrain, ytrain)
  performTrain <- matrix(performTrain, ncol = length(performTrain), nrow = 1)
  
  # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model testing
  pred <- predict(fit, Xtest)  
  performTest <- evalPerformance(pred, ytest)
  
  if (setParam$fit$explanation) {
    # permutation variable importance
    predictor.enet <- Predictor$new(
      model = fit, 
      data = data.frame(Xtrain), 
      y = ytrain, 
      predict.fun = predEnet)
    
    imp.enet <- FeatureImp$new(predictor.enet, loss = "mse")
    pviRank <- imp.enet$results$feature
    pviValue <- imp.enet$results$importance
    
  } else {
    pviRank <- NA
    pviValue <- NA
  }
  
  return(list(estB = estBeta, 
              selectedVars = selectedVars,
              performTrain = performTrain,  # train performance (full train sample)
              performTest = performTest, # (holdout) test performance
              # model agnostic measures
              pvi = cbind(pviRank, pviValue), # permutation variable importance
              tunedAlpha = tunedParams["tunedAlpha", "min"], 
              tunedLambda = tunedParams["tunedLambda", "min"], 
              iSample = iSample))
}
