# (actual) function to simulate linear vs. piecewise-linear data

samplePiecewiseLinearData <- function() {
  
  # # generate one big test data set
  # N <- 1000000
  # pTrash <- 10
  # dgpFolder <- paste0(dataFolder, "/bigTestSamples")
  # reliability <- 1
  # data <- "pwlinear"

  dgpFolder <- paste0(dataFolder, "/", data)
  createFolder(dgpFolder)
  
  # create folder to save data in single sample files
  if (setParam$dgp$singleSamples) {
    sampleFolder <- paste0("/simDataN", N, "_pTrash", pTrash, "_rel", reliability)
    createFolder(paste0(dgpFolder, sampleFolder))
  }
  
  # sample data in parallel; 
  # generate samples in parallel as samples are drawn as random & independent
  data <- parLapply(cl, seq_len(setParam$dgp$nSamples), function(iSample) {
    
    P <- setParam$dgp$p + setParam$dgp$pPWL + pTrash # total number of variables
    # generate matrix of (almost) uncorrelated predictors
    if (iSample > setParam$dgp$nTrain) {
      # test sample with fixed sample size across all simulated conditions (independent of N)
      N <- setParam$dgp$testN
    }
    
    # get predictor values without (!) measurement error
    # additional variable for stepwise effect should be uncorrelated with other 
    #   variables as interactions are with linear effects (Bohrnstedt and Goldberger, 1969)
    # here: add additional variable for nonlinear effect to predictor matrix or
    #       use one of the pTrash variables to generate dummy? 
    X <- createPredictors(N = N, P = P, 
                          corMat = setParam$dgp$predictorCorMat_pwl[seq_len(P), seq_len(P)])
    
    # add names to variables
    colnames(X) <- paste0("Var", seq_len(P))
    
    # create model formula (allows polynomial and interaction effects of any degree/depth)
    #   only create model formula to put into model.matrix
    popModel <- genModel(colnames(X), setParam$dgp$interDepth, setParam$dgp$poly)
    
    # predictor matrix that allows for polynomials and interactions
    #   here we need the numeric non-linear variable to create all possible interactions
    X_int <- model.matrix(as.formula(popModel), data.frame(X))
    
    # remove first degree polynomials from data (they are duplicates!)
    #   only if poly in model matrix, else error 
    if (setParam$dgp$poly > 0) {
      X_int <- rmDuplicatePoly(X_int)
    }
    
    # # check data generating mechanism    
    # apply(X_int, 2, mean)
    # apply(X_int, 2, var)
    # cor(X_int)
    
    # here: we need to switch to the piecewise linear variables
    # add ...
    #   ... the dummy coded version of the threshold and
    #   ... the resulting x value for the second segment to data matrix
    # median split of the original variable means splitting at 0
    #   ... therefore beta0 and beta0.2nd has to be 0! 
    # median split at 0 essentially means we are applying a ReLu function to the original variable
    # -> how to calculate (expected mean and) expected variance of a ReLu transformed variable
    #     ... there are closed form solutions!
    # we do not need to overwrite the linear variable because we need it for the 
    #     linear regression as predictor later and for the piecewise linear regression model, too
    # here: do not use effect coding on the dummy as the dummy is only intermediate step
    #       but dummy will not be kept in predictor matrix to model y; therefore only vector
    
    # dummy to calculate 2nd segment x
    dumVar5 = createDummy(X_int[, "Var5"], q = 0.5, effectCoding = F) 
    dumVar6 = createDummy(X_int[, "Var6"], q = 0.5, effectCoding = F) 
    dumVar7 = createDummy(X_int[, "Var7"], q = 0.5, effectCoding = F) 
    X_int <- cbind(X_int, 
                   Var5.2nd = (X_int[, "Var5"] - quantile(X_int[, "Var5"], 0.5))*dumVar5,
                   Var6.2nd = (X_int[, "Var6"] - quantile(X_int[, "Var6"], 0.5))*dumVar6,
                   Var7.2nd = (X_int[, "Var7"] - quantile(X_int[, "Var7"], 0.5))*dumVar7)
    
    # apply(X_int[,c("Var5", "Var6", "Var5.2nd", "Var6.2nd")], 2, mean)
    # # Var5          Var6      Var5.2nd      Var6.2nd 
    # # 0.0006299155 -0.0002065031  0.3991402179  0.3982295602 
    # apply(X_int[,c("Var5", "Var6", "Var5.2nd", "Var6.2nd")], 2, var)
    # # Var5      Var6  Var5.2nd  Var6.2nd 
    # # 1.0002533 0.9974187 0.3412924 0.3394212 
    # cor(X_int[, c(setParam$dgp$linEffects, setParam$dgp$pwlinEffects)])
    # # nonlineare Variablen sind unkorreliert miteinander; 
    # #   first and second Segment ignorieren, weil Koeffizient für Var5, etc. = 0
    
    # generate matrix of regression coefficients (matrix includes all conditions)
    # rows represent predictors (thus, number of rows depends on pTrash which varies 
    #     between simulated conditions)
    # columns represent conditions (= combination of R2 and lin/inter effect balance)
    
    # here: use data with the original variable and the second segment variable
    # here: add weights (= regression-coefficients) for the second segment variable
    bMatrix <- genBmat(X_int, data, setParam)
    
    # # quick check
    # bMatrix[rownames(bMatrix) %in% c(setParam$dgp$linEffects, setParam$dgp$pwlinEffects),]
    
    # calculate R^2 for every combination of R2 and lin/inter effect balance
    # print R^2 as a quick sanity check (removed for speed sake) 
    # here: use data without the original variables, but with dummies in the data
    R2 <- sapply(seq_len(ncol(bMatrix)), function(x) getR2(X_int, bMatrix[,x], setParam$dgp$sigmaE))
    
    # calculate dependent variable for every combination of R2 and lin/inter effect balance
    # dependent variable is simulated based on ...
    #   ... predictors without measurement error!
    #   ... dummy variables instead of the original variable
    
    # here: use data without the original variables, but with dummies in the data
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
    
    # save ...
    #     ... yMat with dependent variable for all R2 - lin/inter effect conditions in columns
    #     ... X_int predictor matrix (identical for all R2 - lin/inter effect conditions)
    #         including measurement error in the data
    #     ... R2 based on simulated regression coefficients
    #     ... R2_wME based on data with measurement error
    dataList <- list(yMat = yMatrix, # criterion (DV) with R2 x lin_inter in columns
                     X_int = X_final, # these are the predictors (IV) with measurement error
                     R2 = R2) #, # without measurement error
    # R2_wME = R2_wME) # with measurement error
    
    # # to save big test sample
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
  
  # save data in one big rda file (with all samples)
  if (!setParam$dgp$singleSamples) {
    # save to one big rda file with nSamples list entries! 
    #   all training and test sample in one big rda file  
    names(data) <- c(seq_len(setParam$dgp$nTrain), 
                     paste0("test", seq_len(setParam$dgp$nTest)))
    
    fileName <- paste0("simDataN", N, "_pTrash", pTrash, "_rel", reliability, ".rda")
    save(data, file = paste0(dgpFolder, "/", fileName))
  }
}