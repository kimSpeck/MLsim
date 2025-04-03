# (actual) function to simulate linear vs. non-linear data

# simulate two additional variables 
#   - nonlinear "original" variables are uncorrelated with ...
#     ... each other and with the linear predictor variables
#     - nonlinear effects should suit the random forest which typically works better
#       for uncorrelated variables 
#     - technical reason: engineering the correlation of the dummies for which the 
#       effects are simulated is rather difficult as we sample data for a continuous 
#       variable first
#   - additional to p and pTrash; but! keep total number of predictors in the ENET 
#     constant to linear:interaction condition! 
#     thus, remove two randomly chosen noise variables and their interactions with 
#       every other variable from the predictor matrix of the ENET
# separate both additional variables into one dummy (based on quantiles)
#   one dummy (= 2 steps): 50% quantile; effect-coded {-1, 1} 
#     due to effect coding each dummy variable has unit variance
#     their interaction has smaller variance
# effect coding of the dummy variable has the additional benefit that the main effect
#   of the dummy variable is reflected in the regression coefficient of the dummy
#   therefore, we do not need to add an intercept!
# simulate mean differences for the dummies as regression coefficients of the effect
#   coded dummy variable
# interaction between these two dummies

# Vergleich ist linear:Interaktions und linear:non-linear; nicht interactions vs.non-linear; 
#   therefore reliability of non-linear effects shoulde be the same as for linear effects
#   which is more than for interactions (thus, non-linear effects should be easier to find)
# Wie halte ich die Varianz der Dummies identisch mit der Varianz in den linearen Prädiktoren? 
# Wie kontrolliere ich das R² wenn es Mittelwertsdifferenzen und Reg.koeffizienten gibt?
#   - Nein: Mittelwertsunterschiede bestimmen die Varianz der Dummies; Balanciertheit bestimmt die Varianz der Dummies
#   - Reg.koeffizient bestimmt das R² und den Mittelwertsunterschied!
# Was nehme ich dann als Prädiktor in den Modellen? immer die "normale" Variable
#   - ENET muss Gerade durch Gruppen ziehen
#   - Bäume können splitten
# Jetzt durch 2 Stufen, durch die ENET einfach Regressionsgerade ziehen kann?! 

#####
# NO: separate additional variable into two dummies (based on quantiles)
#   two dummies (= 3 steps): 33.3 quantiles
# simulate mean differences for the dummies
# interaction between these two dummies is meaningless as no observation will 
#   fall into both dummies (D1 = 1 and D2 = 1 never occurs by design)
#####

sampleNonlinearData <- function() {
  
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
    
    P <- setParam$dgp$p + setParam$dgp$pNL + pTrash # total number of variables
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
                          corMat = setParam$dgp$predictorCorMat_nl[seq_len(P), seq_len(P)])
    
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
    
    # here: we need to switch to the dummy coded variables
    # add the dummy coded version of the non-linear effect variables to the data
    #   overwriting var does not work because we need the numeric variable later on
    # here: set effect coding to get dummies with variance = 1 
    X_int <- cbind(X_int, 
                   dumVar5.1 = createDummy(X_int[, "Var5"], q = 0.5, effectCoding = T),
                   dumVar6.1 = createDummy(X_int[, "Var6"], q = 0.5, effectCoding = T))
    # add the interaction between the dummies
    X_int <- cbind(X_int, 
                   `dumVar5.1:dumVar6.1` = X_int[, "dumVar5.1"] * X_int[, "dumVar6.1"])
    # # var(dumVar5.1) = 1; var(dumVar6.1) = 1; var(dumVar5.1:dumVar6.1) = 1
    # # dummies are uncorrelated with everything
    # var(X_int[, "dumVar5.1"]); var(X_int[, "dumVar6.1"]); var(X_int[, "dumVar5.1:dumVar6.1"])
    # cor(X_int[, c(setParam$dgp$linEffects, setParam$dgp$nonlinEffects)])
    
    # detect dummies to remove the corresponding original variables temporarily
    detectDum <- colnames(X_int)[stringr::str_detect(colnames(X_int), "^dum") ] 
    varNamesDummies <- paste0("Var", unique(stringr::str_extract(detectDum, "(?<=dumVar)\\d+")))
    X_intDummy <- X_int[, !colnames(X_int) %in% c(varNamesDummies)] # remove original variables
    
    # generate matrix of regression coefficients (matrix includes all conditions)
    # rows represent predictors (thus, number of rows depends on pTrash which varies 
    #     between simulated conditions)
    # columns represent conditions (= combination of R2 and lin/inter effect balance)
    
    # here: use data without the original variables, but with dummies in the data
    # here: add weights (= regression-coefficients) for the dummies and their interaction
    bMatrix <- genBmat(X_intDummy, data, setParam)
    
    # # quick check
    # bMatrix[rownames(bMatrix) %in% c(setParam$dgp$linEffects, setParam$dgp$nonlinEffects),]
    
    # calculate R^2 for every combination of R2 and lin/inter effect balance
    # print R^2 as a quick sanity check (removed for speed sake) 
    # here: use data without the original variables, but with dummies in the data
    R2 <- sapply(seq_len(ncol(bMatrix)), function(x) getR2(X_intDummy, bMatrix[,x], setParam$dgp$sigmaE))
    
    # calculate dependent variable for every combination of R2 and lin/inter effect balance
    # dependent variable is simulated based on ...
    #   ... predictors without measurement error!
    #   ... dummy variables instead of the original variable
    
    # here: use data without the original variables, but with dummies in the data
    yMatrix <- sapply(seq_len(ncol(bMatrix)), function(x) {
      calcDV(X = X_intDummy, b = bMatrix[,x],
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
    # we cannot calculate R2 with measurement error directly: coefficients only apply to
    #   predictor matrix with dummies; we could however recreate the dummies based on 
    #   the variable with measurement error (?)
    # R2_wME <- sapply(seq_len(ncol(bMatrix)), function(x) {
    #   var(X_final %*% bMatrix[,x]) / (var(X_int %*% bMatrix[,x]) + setParam$dgp$sigmaE^2)
    # })
    
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