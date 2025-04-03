saveENET <- function(estRes, iCond, setParam, resList) {
  # in general: concatenate outcome from different conditions!
  # coefficients
  estBeta <- do.call(cbind, lapply(estRes, function(X) X[["estB"]]))
  # remove variables which are not selected by Enet from mean calculation:
  #   coefficients that are exactly 0 -> NA
  resList[["estBetaFull"]] <- ifelse(estBeta == 0, NA, estBeta)
  
  resList[["estBeta"]] <- cbind(getStats(estBeta, 1, setParam$dgp$nSamples),
                                     idxCondLabel = iCond)
  
  # selected variables (how often are predictors selected in model)
  selectedVars <- do.call(cbind, lapply(estRes, function(X) X[["selectedVars"]]))
  
  nSelection <- apply(selectedVars, MARGIN = 1, sum) # frequency of variable selection 
  percSelection <- nSelection / ncol(selectedVars) * 100 # relative frequency 
  resList[["varSelection"]] <- cbind(nSelection = nSelection, 
                        percSelection = percSelection,
                        idxCondLabel = iCond) # stats per predictor

  # only true predictors in model (8 predictors with simulated effects and every other variable == 0)
  # how many of the linear effects are recovered?
  nSelectedLin <- sapply(seq_len(ncol(selectedVars)), function(iCol) { 
    sum(setParam$dgp$linEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
  # how many of the interaction effects are recovered?
  nSelectedInter <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
    sum(setParam$dgp$interEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
  
  # all simulated effects selected in model?
  selectedAll <- nSelectedLin + nSelectedInter == length(c(setParam$dgp$linEffects, setParam$dgp$linEffects))
  
  # only simulated effects selected (i.e., every other predictor is not selected!)
  nSelectedOthers <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
    sum(!(rownames(selectedVars)[selectedVars[,iCol]] %in% c(setParam$dgp$interEffects, setParam$dgp$linEffects)))})
  
  # final Enet-Parameters for each sample
  tunedAlphaVec <- unname(do.call(c, lapply(estRes, function(X) X[["tunedAlpha"]])))
  tunedLambdaVec <- unname(do.call(c, lapply(estRes, function(X) X[["tunedLambda"]])))
  
  # variable selection stats per sample
  # for all: 1 = TRUE, 0 = FALSE
  resList[["selectionPerSample"]] <- cbind(nLin = nSelectedLin, nInter = nSelectedInter, 
                                           all.T1F0 = selectedAll, nOthers = nSelectedOthers, 
                                           alpha = tunedAlphaVec, lambda = tunedLambdaVec,
                                           idxCondLabel = iCond)
  
  return(resList)
}
