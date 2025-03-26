# (actual) function to simulate linear vs. non-linear data

sampleNonlinearData <- function() {
  
  createFolder(paste0(dataFolder,"/nonlinear"))
  
  if (setParam$dgp$singleSamples) {
    sampleFolder <- paste0("/simDataN", N, "_pTrash", pTrash, "_rel", reliability)
    createFolder(paste0(dataFolder,"/nonlinear", sampleFolder))
  }
  
  # sample data in parallel; 
  # generate samples in parallel as samples are drawn as random & independent
  data <- parLapply(cl, seq_len(setParam$dgp$nSamples), function(iSample) {
    
    P <- setParam$dgp$p + pTrash # total number of variables
    # generate matrix of (almost) uncorrelated predictors
    if (iSample > setParam$dgp$nTrain) {
      # test sample with fixed sample size across all simulated conditions (independent of N)
      N <- setParam$dgp$testN
    }
    
    # get predictor values without (!) measurement error
    X <- createPredictors(N = N, P = P, 
                          corMat = setParam$dgp$predictorCorMat[seq_len(P), seq_len(P)])
    
    # add names to variables
    colnames(X) <- paste0("Var", seq_len(P))
    
  })
  
  if (!setParam$dgp$singleSamples) {
    # save to one big rda file with nSamples list entries! 
    #   all training and test sample in one big rda file  
    names(data) <- c(seq_len(setParam$dgp$nTrain), 
                     paste0("test", seq_len(setParam$dgp$nTest)))
    
    fileName <- paste0("simDataN", N, "_pTrash", pTrash, "_rel", reliability, ".rda")
    save(data, file = paste0(dataFolder, "/nonlinear/", fileName))
  }
}