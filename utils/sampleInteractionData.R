# (actual) function to simulate linear vs. interaction data

sampleInteractionData <- function() {
  
  # # generate one big test data set
  # N <- 1000000
  # pTrash <- 10
  # dgpFolder <- paste0(dataFolder, "/bigTestSamples")
  # reliability <- 1
  # data <- "inter"
  
  dgpFolder <- paste0(dataFolder, "/", data)
  createFolder(dgpFolder)
  
  if (setParam$dgp$singleSamples) {
    sampleFolder <- paste0("/simDataN", N, "_pTrash", pTrash, "_rel", reliability)
    createFolder(paste0(dgpFolder, sampleFolder))
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
    
    # create model formula (allows polynomial and interaction effects of any degree/depth)
    popModel <- genModel(colnames(X), setParam$dgp$interDepth, setParam$dgp$poly)
    
    # predictor matrix that allows for polynomials and interactions
    X_int <- model.matrix(as.formula(popModel), data.frame(X))
    
    # remove first degree polynomials from data (they are duplicates!)
    #   only if poly in model matrix, else error 
    if (setParam$dgp$poly > 0) {
      X_int <- rmDuplicatePoly(X_int)
    }
    
    # generate matrix of regression coefficients (matrix includes all conditions)
    # rows represent predictors (thus, number of rows depends on pTrash which varies 
    #     between simulated conditions)
    # columns represent conditions (= combination of R2 and lin/inter effect balance)
    bMatrix <- genBmat(X_int, data, setParam)
    
    # calculate R^2 for every combination of R2 and lin/inter effect balance
    # print R^2 as a quick sanity check (removed for speed sake) 
    R2 <- sapply(seq_len(ncol(bMatrix)), function(x) getR2(X_int, bMatrix[,x], setParam$dgp$sigmaE))
    
    # calculate dependent variable for every combination of R2 and lin/inter effect balance
    # dependent variable is simulated from predictors without measurement error!
    yMatrix <- sapply(seq_len(ncol(bMatrix)), function(x) {
      calcDV(X = X_int, b = bMatrix[,x],
             sigmaE = setParam$dgp$sigmaE, N = N)
    })
    colnames(yMatrix) <- setParam$dgp$condLabels
    
    # add measurement error/reliability manipulatio to data
    #   add measurement error only to X (poly & interactions are calculated based on X)
    #   measurement error ...
    #     ... independent for each predictor 
    #     ... normally distributed with M = 0 & SD according to reliability 
    
    # error variance according to reliability
    covMatError <- diag(P) * (1 - reliability)/reliability
    measureError <- rmvnorm(n = N, mean = rep(0, P), sigma = covMatError)
    
    # add measurement error to predictors
    X_wME <- X + measureError
    
    # predictor matrix that allows for polynomials and interactions
    X_final <- model.matrix(as.formula(popModel), data.frame(X_wME))
    
    # remove first degree polynomials from data (they are duplicates!)
    if (setParam$dgp$poly > 0) {
      X_final <- rmDuplicatePoly(X_final)
    }
    
    # recalculate R2 for predictors with measurement error
    R2_wME <- sapply(seq_len(ncol(bMatrix)), function(x) {
      var(X_final %*% bMatrix[,x]) / (var(X_int %*% bMatrix[,x]) + setParam$dgp$sigmaE^2)
    })
    
    # save ...
    #     ... yMat with dependent variable for all R2 - lin/inter effect conditions in columns
    #     ... X_int predictor matrix (identical for all R2 - lin/inter effect conditions)
    #         including measurement error in the data
    #     ... R2 based on simulated regression coefficients
    #     ... R2_wME based on data with measurement error
    dataList <- list(yMat = yMatrix, # criterion (DV) with R2 x lin_inter in columns
                     X_int = X_final, # these are the predictors (IV) with measurement error
                     R2 = R2, # without measurement error
                     R2_wME = R2_wME) # with measurement error
    
    # testFileName <- paste0("simDataN", N, "_pTrash", pTrash, "_rel", reliability, "_", data, ".rda")
    # save(dataList, file = paste0(dgpFolder, "/", testFileName))
    
    if (setParam$dgp$singleSamples){
      sampleNames <- c(seq_len(setParam$dgp$nTrain),
                       paste0("test", seq_len(setParam$dgp$nTest)))
      fileName <- paste0("sample_", sampleNames[iSample], ".rda")
      save(dataList, file = paste0(dgpFolder, sampleFolder, "/", fileName))
    } else {
      dataList
    }
  })
  
  if (!setParam$dgp$singleSamples) {
    # save to one big rda file with nSamples list entries! 
    #   all training and test sample in one big rda file  
    names(data) <- c(seq_len(setParam$dgp$nTrain), 
                     paste0("test", seq_len(setParam$dgp$nTest)))
    
    fileName <- paste0("simDataN", N, "_pTrash", pTrash, "_rel", reliability, ".rda")
    save(data, file = paste0(dgpFolder, "/", fileName))
  }
}