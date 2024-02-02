# simulate data

# ToDo: simulate test data sample for each dgp (50% size in terms of observations)
# ! interactions are specified between variables that do not show main effects
#     unrealistic assumption but due to correlations easier to find coefficients 
#     for simulating effects given fixed R2 target
# ! correlation matrices to sample predictor variables from multivariate normal are
#     randomly drawn in each sample 
#     ->  in real simulation: one randomly drawn but fixed correlation matrix. 
#         calculate regression coefficients for this random but specific correlation matrix
#         and check that R^2 is met exactly


# fixed: depending in seed we receive non PSM matrices sigma as covariance matrices in multivariate normal
#     Warnmeldungen: 
#       1: In rmvnorm(n = N, mean = mP, sigma = rX) :
#       sigma is numerically not positive semidefinite
# however, rmvnorm somehow deals with this issue (see: https://stats.stackexchange.com/questions/267908/clarification-regarding-rmvnorm-in-r)

# ToDo: write big matrices directly to their final place of storage (memory intensive and time consuming)

# load packages
library(mvtnorm)
library(truncnorm)
library(parallel)

# load parameters & custom functions 
source("setParameters.R") # parameter values
source("simTools.R") # functions for data simulation

timeStampFolder <- format(Sys.time(), "%d%m%y_%H%M%S")

# generate folder for data files
dataFolder <- "data"
if (!file.exists(dataFolder)){
  dir.create(dataFolder)
}

# generate folder for log files
logFolder = "log"
if (!file.exists(logFolder)){
  dir.create(logFolder)
}

# nCoresSampling <- detectCores() - 1 # simulate data
nCoresSampling <- 5 # simulate data

# grid to simulate data with mapply later
# simulation via mapply to easily simulate subsets of parameter combinations
# create factor structure for predictor variables to increase reliability
gridFull <- expand.grid(N = setParam$dgp$N,
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability,
                        factors = c(TRUE, FALSE))
str(gridFull) # ! factors are interpreted as level numbers!; only character variables are interpreted by their name!
# create seed number for parallel cluster (reproducibility of generated data)
set.seed(8967369)
seedNum <- sample(1:999999, dim(gridFull)[1], replace = FALSE) 
gridFull$sampleSeed <- seedNum[1:dim(gridFull)[1]]

# sample data in parallel
createData <- function(N, pTrash, reliability, factors, sampleSeed){
  
  environment(sampleData) <- environment()
  
  # Initiate cluster; type = "FORK" only on Linux/MacOS: contains all environment variables automatically
  cl <- makeCluster(nCoresSampling, type = "FORK",
                    outfile = paste0(logFolder, "/", "simulateDataStatus",
                                     timeStampFolder, ".txt"))
  
  # set seed that works for parallel processing
  set.seed(sampleSeed)
  s <- .Random.seed
  clusterSetRNGStream(cl = cl, iseed = s)
  
  sampleData() # run function to actually create dataset
  
  # close cluster to return resources (memory) back to OS
  stopCluster(cl)
}

# actual function to simulate data
sampleData <- function() {
  
  # sample data in parallel
  data <- parLapply(cl, seq_len(setParam$dgp$nSamples), function(iSample) {
    
    # # test it
    # pTrash <- 10
    # N <- 10000
    # reliability <- 0.6
    # iSample <- 1
    # factors <- TRUE
    
    P <- setParam$dgp$p + pTrash # total number of variables
    # generate matrix of (almost) uncorrelated predictors
    if (iSample > setParam$dgp$nTrain) {
      N <- setParam$dgp$testNpc * N
    }
    
    # either create single indicators or factors as predictor
    if (factors) {
      # get predictor values without (!) measurement error
      X_full <- createFactorPredictors(N = N, p = setParam$dgp$p, 
                                       nIndicator = setParam$dgp$nIndicator, 
                                       pTrash = pTrash, reliability = 1,
                                       corMat = setParam$dgp$predictorCorMat[seq_len(P), seq_len(P)])
      
      # add names to variables
      varGrid <- expand.grid(indicator = seq_len(setParam$dgp$nIndicator),
                             factor = seq_len(setParam$dgp$p))
      colnames(X_full) <- c(paste0("F", varGrid$factor, "_", varGrid$indicator),
                            paste0("Var", setParam$dgp$p+seq_len(pTrash)))
      
      # calculate sum scores (i.e., average from all indicators; 
      #     otherwise range of factors as predictors would deviate from other predictors)
      X <- matrix(NA, ncol = P, nrow = N)
      X[,seq_len(setParam$dgp$p)] <- sapply(seq_len(setParam$dgp$p), function(iF) {
        tmp_indicators <- stringr::str_which(colnames(X_full), paste0("F", iF))
        rowSums(X_full[,tmp_indicators], na.rm = T) / length(tmp_indicators)
      })
      X[,c(setParam$dgp$p + seq_len(pTrash))] <- X_full[,c((setParam$dgp$p * setParam$dgp$nIndicator) + seq_len(pTrash))]
      
      # X_full: all variables (i.e., indicators of every factor)
      # X: only sum scores
      
    } else { # create predictors without factor structure!
      # get predictor values without (!) measurement error
      X <- createPredictors(N = N, P = P, 
                            corMat = setParam$dgp$predictorCorMat[seq_len(P), seq_len(P)])
    }
    
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
    
    # to do: check if code runs without these lines
    # # generate vector regression weights
    # b <- rep(0, ncol(X_int)) # initiate all b-values with value of zero 
    # names(b) <- colnames(X_int) # assign names
    #
    # linEffects <- setParam$dgp$linEffects
    # interEffects <- setParam$dgp$interEffects
    # nEffects <- length(c(setParam$dgp$linEffects, setParam$dgp$interEffects))
    
    # generate matrix of regression coefficients (matrix includes all conditions)
    # rows represent predictors (thus, number of rows depends on pTrash which varies 
    #     between simulated conditions)
    # columns represent conditions (= combination of R2 and lin/inter effect balance)
    bMatrix <- genBmat(X_int, setParam)
    
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
    if (factors) {
      # here!
      covMatError <- diag(setParam$dgp$p * setParam$dgp$nIndicator + pTrash) * (1 - reliability)/reliability
      measureError <- rmvnorm(n = N, 
                              mean = rep(0, setParam$dgp$p * setParam$dgp$nIndicator + pTrash), 
                              sigma = covMatError)
      
      # add measurement error to predictors
      X_full_wME <- X_full + measureError
      
      colnames(X_full_wME) <- c(paste0("F", varGrid$factor, "_", varGrid$indicator),
                                paste0("Var", setParam$dgp$p+seq_len(pTrash)))
      
      # calculate sum scores (i.e., average from all indicators; 
      #     otherwise range of factors as predictors would deviate from other predictors)
      X_wME <- matrix(NA, ncol = P, nrow = N)
      X_wME[,seq_len(setParam$dgp$p)] <- sapply(seq_len(setParam$dgp$p), function(iF) {
        tmp_indicators <- stringr::str_which(colnames(X_full_wME), paste0("F", iF))
        rowSums(X_full_wME[,tmp_indicators], na.rm = T) / length(tmp_indicators)
      })
      X_wME[,c(setParam$dgp$p + seq_len(pTrash))] <- X_full_wME[,c((setParam$dgp$p * setParam$dgp$nIndicator) + seq_len(pTrash))]
      
      # add names to variables
      colnames(X_wME) <- paste0("Var", seq_len(P))
      
    } else {
      covMatError <- diag(P) * (1 - reliability)/reliability
      measureError <- rmvnorm(n = N, mean = rep(0, P), sigma = covMatError)
      
      # add measurement error to predictors
      X_wME <- X + measureError
    }
    
    # predictor matrix that allows for polynomials and interactions
    X_final <- model.matrix(as.formula(popModel), data.frame(X_wME))
    
    # remove first degree polynomials from data (they are duplicates!)
    X_final <- rmDuplicatePoly(X_final)
    
    # recalculate R2 for predictors with measurement error
    R2_wME <- sapply(seq_len(ncol(bMatrix)), function(x) {
      var(X_final %*% bMatrix[,x]) / (var(X_int %*% bMatrix[,x]) + setParam$dgp$sigmaE^2)
    })
    
    # # run single indicator SEM to check if reliabilities are simulated correctly
    # SImodel <- genSingleIndicatorModel(P, reliability)
    # checkSimParam <- lapply(seq_len(dim(yMatrix)[2]), function(iR2_LI) {
    #   iR2_LI <- 6
    #   X_check <- cbind(X_final[,1:P], y = yMatrix[,iR2_LI])
    # 
    #   fit <- lavaan::sem(SImodel, data=X_check)
    #   # lavaan::summary(fit) # check lavaan output 
    #   
    #   # save path coefficients for predictors with simulated effects
    #   estBeta <- fit@Model@GLIST[["beta"]][(P+1),seq_along(setParam$dgp$linEffects)]
    #   
    #   # save correlations between latent variables of correlated predictors
    #   estPsi <- fit@Model@GLIST[["psi"]][seq_along(setParam$dgp$linEffects), seq_along(setParam$dgp$linEffects)]
    #   estPsi <- estPsi[upper.tri(estPsi)] # F1F2, F1F3, F2F3, F1F4, F2F4, F3F4 
    #   
    #   list(estBeta = estBeta, 
    #        estPsi = estPsi)
    # })
    # estPsi <- do.call(rbind, lapply(seq_along(checkSimParam), function(subList) {
    #   rbind(checkSimParam[[subList]][["estPsi"]])
    # }))
    # colnames(estPsi) <- c("F1F2", "F1F3", "F2F3", "F1F4", "F2F4", "F3F4")
    # 
    # estBeta <- do.call(rbind, lapply(seq_along(checkSimParam), function(subList) {
    #   rbind(checkSimParam[[subList]][["estBeta"]])
    # }))
    # estBeta <- cbind(estBeta, matrix(setParam$dgp$trueEffects$lin, ncol = 1))
    # colnames(estBeta) <- c(setParam$dgp$linEffects, "trueBeta")
    
    # save ...
    #     ... yMat with dependent variable for all R2 - lin/inter effect conditions in columns
    #     ... X_int predictor matrix (identical for all R2 - lin/inter effect conditions)
    #     ... trueB simulated regression coefficients
    #     ... R2 based on simulated regression coefficients
    list(yMat = yMatrix,
         X_int = X_final, # save predictors with measurement error (with sum score for factors)
         # all indicators as predictors but NA matrix if factors = FALSE
         X_full_wME = if(factors) X_full_wME else matrix(NA, ncol = (setParam$dgp$p * setParam$dgp$nIndicator) + pTrash, nrow = N),
         # X_int = X_int, # data without measurement error
         # X_full = X_full, # with single indicators for each factor but without measurement error
         R2 = R2, # without measurement error
         R2_wME = R2_wME) # with measurement error
  })
  
  names(data) <- c(seq_len(setParam$dgp$nTrain), 
                   paste0("test", seq_len(setParam$dgp$nTest)))
  
  #save to rda file
  fileName <- paste0("simDataN", N, "_pTrash", pTrash, "_rel", reliability, "_f", ifelse(factors, 1, 0), ".rda")
  save(data, file = paste0(dataFolder, "/", fileName))
}

# out <- do.call(mapply, c(FUN = createData, gridFull[1,])) # test
out <- do.call(mapply, c(FUN = createData, gridFull))
