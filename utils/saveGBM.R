saveGBM <- function(estRes, iCond, setParam, resList){
  # in general: concatenate outcome from different conditions!
  ## H-Statistic
  if (setParam$fit$explanation & setParam$fit$InterStrength) {
    resList[["interStrength"]] <- do.call(rbind, lapply(seq_len(length(estRes)), function(iSample) {
      cbind(idxCondLabel = rep(iCond, dim(estRes[[iSample]][["interStrength"]])[1]),
            sample = iSample,
            estRes[[iSample]][["interStrength"]])
    }))
  } else {
    resList[["interStrength"]] <- NA
  }
  
  # R² across all test sets in the CV procedure (for the final hyperparameters)
  performCVtestMat <- do.call(rbind, lapply(estRes, function(X) X[["performCVtest"]])) # each sample
  colnames(performCVtestMat) <- c("RMSE_CVtest", "Rsq_CVtest", "MAE_CVtest")
  # aggregated performance measures (M, SD, SE)
  resList[["performCVtestStats"]] <- cbind(getStats(performCVtestMat, 2, setParam$dgp$nSamples),
                                         idxCondLabel = iCond)
  
  # tuning parameters (GBM)
  tunedShrinkage <- unname(do.call(c, lapply(estRes, 
                                       function(X) X[["tunedShrinkage"]])))
  tunedMax_depth <- unname(do.call(c, lapply(estRes, 
                                             function(X) X[["tunedMax_depth"]])))
  tunedMin_child_weight <- unname(do.call(c, lapply(estRes, 
                                                    function(X) X[["tunedMin_child_weight"]]))) 
  tunedNrounds <- unname(do.call(c, lapply(estRes, 
                                           function(X) X[["tunedNrounds"]])))
  
  resList[["selectionPerSample"]] <- cbind(shrinkage = tunedShrinkage, max_depth = tunedMax_depth,
                              min_child_weight = tunedMin_child_weight, 
                              Nrounds = tunedNrounds,
                              idxCondLabel = iCond)
  
  return(resList)
}