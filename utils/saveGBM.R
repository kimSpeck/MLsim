saveGBM <- function(estRes, iCond, setParam, resList){
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