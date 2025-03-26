saveRF <- function(estRes, iCond, setParam, resList) {
  # in general: concatenate outcome from different conditions!
  # save oob predictions (oob based predicted criterion values) and oob R^2
  resList[["oobPredictions"]] <- do.call(cbind, lapply(estRes, function(X) X[["oobPredictions"]]))
  
  resList[["oobR2"]] <- do.call(c, lapply(estRes, function(X) X[["oobR2"]]))
  
  # R² across all test sets in the CV procedure (for the final hyperparameters)
  performCVtestMat <- do.call(rbind, lapply(estRes, function(X) X[["performCVtest"]])) # each sample
  colnames(performCVtestMat) <- c("RMSE_CVtest", "Rsq_CVtest", "MAE_CVtest")
  # aggregated performance measures (M, SD, SE)
  resList[["performCVtestStats"]] <- cbind(getStats(performCVtestMat, 2, setParam$dgp$nSamples),
                                           idxCondLabel = iCond) 
  
  # tuning parameter choice
  tunedMtryVec <- unname(do.call(c, lapply(estRes, function(X) X[["tunedMtry"]])))
  tunedSplitruleVec <- unname(do.call(c, lapply(estRes, function(X) X[["tunedSplitrule"]])))
  tunedMinNodeVec <- unname(do.call(c, lapply(estRes, function(X) X[["tunedMinNodeSize"]])))
  
  # variable selection stats per sample
  resList[["selectionPerSample"]] <- cbind(mtry = tunedMtryVec, 
                                           splitRule = tunedSplitruleVec,
                                           minNode = tunedMinNodeVec,
                                           idxCondLabel = iCond)
  
  return(resList)
}
