saveRF <- function(estRes, iCond, setParam, resList) {
  # here!
  # what else do we need to extract specifically for the random forest?
   
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
