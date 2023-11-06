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

# ToDo: use simulated effects from setParam!

# for reproducibility purposes
# set.seed(42) # not positive semidefinite covariance matrix ):
# set.seed(7382) 

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
gridFull <- expand.grid(N = setParam$dgp$N,
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability)
str(gridFull) # ! factors are interpreted as level numbers!; only character variables are interpreted by their name!
# create seed number for parallel cluster (reproducibility of generated data)
set.seed(8967369)
seedNum <- sample(1:999999, dim(gridFull)[1], replace = FALSE) 
gridFull$sampleSeed <- seedNum[1:dim(gridFull)[1]]

# sample data in parallel
createData <- function(N, pTrash, reliability, sampleSeed){
  
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
    # reliability <- 0.7
    # iSample <- 1
    
    P <- setParam$dgp$p + pTrash # total number of variables
    # generate matrix of (almost) uncorrelated predictors
    if (iSample > setParam$dgp$nTrain) {
      N <- setParam$dgp$testNpc * N
    }
    # to do: remove function and use rmvnorm directly now?!
    X <- createPredictors(N = N, P = P, 
                          corMat = setParam$dgp$predictorCorMat[seq_len(P), seq_len(P)])
    
    # add names to variables
    colnames(X) <- paste0("Var", seq_len(P))
    
    # create model formula (allows polynomial and interaction effects of any degree/depth)
    popModel <- genModel(colnames(X), setParam$dgp$interDepth, setParam$dgp$poly)
    
    # predictor matrix that allows for polynomials and interactions
    X_int <- model.matrix(as.formula(popModel),data.frame(X))
    
    # remove first degree polynomials from data (they are duplicates!)
    X_int <- rmDuplicatePoly(X_int)
    
    # generate vector regression weights
    b <- rep(0, ncol(X_int)) # initiate all b-values with value of zero 
    names(b) <- colnames(X_int) # assign names
    
    # polyPreds <- names(b)[stringr::str_which(names(b), "^poly")]
    # interPreds <- names(b)[stringr::str_which(names(b), ":")]
    # linPreds <- names(b)[!(names(b) %in% polyEffects) & !(names(b) %in% interEffects)]
    
    linEffects <- sapply(seq_len(setParam$dgp$p), function(x) paste0("Var", x))
    # choose variables for interaction that have no linear effects (R2 budget)
    # interEffects <- c("Var1:Var2", "Var1:Var4", "Var2:Var3", "Var3:Var4")
    interEffects <- c("Var5:Var6", "Var5:Var8", "Var6:Var7", "Var7:Var8")
    nEffects <- length(c(linEffects, interEffects))
    
    # generate matrix of regression coefficients (matrix includes all conditions)
    # rows represent predictors
    # columns represent conditions (= combination of R2 and lin/inter effect balance)
    bMatrix <- matrix(0, 
                      ncol = length(setParam$dgp$Rsquared) * length(setParam$dgp$percentLinear),
                      nrow = length(b))
    
    rownames(bMatrix) <- colnames(X_int) 
    colnames(bMatrix) <- setParam$dgp$condLabels 
    
    # 50 - 50 (all effects with identical regresssion coefficients)
    # R2 <- var(X %*% b) / (var(X %*% b) + sigmaE^2)
    # (R2 <- sapply(seq_len(ncol(bMatrix)), function(x) getR2(X_int, bMatrix[,x], setParam$dgp$sigmaE)))
    
    bMatrix[c(linEffects, interEffects), 1] <- rep(0.096, times = nEffects) # R2 = 0.1
    bMatrix[c(linEffects, interEffects), 2] <- rep(0.187, times = nEffects) # R2 = 0.3
    bMatrix[c(linEffects, interEffects), 3] <- rep(0.286, times = nEffects) # R2 = 0.5
    bMatrix[c(linEffects, interEffects), 4] <- rep(0.572, times = nEffects) # R2 = 0.8
    
    # 80 - 20 (most effects in linEffects)
    # setParam$dgp$Rsquared/100*80
    # setParam$dgp$Rsquared/100*20
    
    # R^2 = 0.1
    bMatrix[c(linEffects), 5] <- rep(0.101, times = length(linEffects)) # R2 = 0.08
    bMatrix[c(interEffects), 5] <- rep(0.074, times = length(interEffects)) # R2 = 0.02
    # R^2 = 0.3
    bMatrix[c(linEffects), 6] <- rep(0.192, times = length(linEffects)) # R2 = 0.24
    bMatrix[c(interEffects), 6] <- rep(0.128, times = length(interEffects)) # R2 = 0.06
    # R^2 = 0.5
    bMatrix[c(linEffects), 7] <- rep(0.32, times = length(linEffects)) # R2 = 0.4
    bMatrix[c(interEffects), 7] <- rep(0.2, times = length(interEffects)) # R2 = 0.1
    # R^2 = 0.8
    bMatrix[c(linEffects), 8] <- rep(0.655, times = length(linEffects)) # R2 = 0.64
    bMatrix[c(interEffects), 8] <- rep(0.27, times = length(interEffects)) # R2 = 0.16
    
    # 20 - 80 (most effects in interEffects)
    # R^2 = 0.1
    bMatrix[c(interEffects), 9] <- rep(0.151, times = length(interEffects)) # R2 = 0.08
    bMatrix[c(linEffects), 9] <- rep(0.049, times = length(linEffects)) # R2 = 0.02
    # R^2 = 0.3
    bMatrix[c(interEffects), 10] <- rep(0.305, times = length(interEffects)) # R2 = 0.24
    bMatrix[c(linEffects), 10] <- rep(0.095, times = length(linEffects)) # R2 = 0.06
    # R^2 = 0.5
    bMatrix[c(interEffects), 11] <- rep(0.48, times = length(interEffects)) # R2 = 0.4
    bMatrix[c(linEffects), 11] <- rep(0.128, times = length(linEffects)) # R2 = 0.1
    # R^2 = 0.8
    bMatrix[c(interEffects), 12] <- rep(1.0, times = length(interEffects)) # R2 = 0.64
    bMatrix[c(linEffects), 12] <- rep(0.17, times = length(linEffects)) # R2 = 0.16
    
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
    covMatError <- diag(P)                              # identity matrix
    diag(covMatError) <- (1 - reliability)/reliability  # error variance according to reliability
    measureError <- rmvnorm(n = N, mean = rep(0, P), sigma = covMatError)
    
    # add measurement error to predictors
    XwME <- X + measureError
    
    # predictor matrix that allows for polynomials and interactions
    X_final <- model.matrix(as.formula(popModel),data.frame(XwME))
    
    # remove first degree polynomials from data (they are duplicates!)
    X_final <- rmDuplicatePoly(X_final)
    
    # # recalculate R2 for predictors with measurement error
    # R2_wME <- sapply(seq_len(ncol(bMatrix)), function(x) getR2(X_final, bMatrix[,x], setParam$dgp$sigmaE))
    
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
         # X_int = X_int, # without measurement error
         X_int = X_final, # save predictors with measurement error
         trueB = bMatrix, 
         R2 = R2)
  })
  
  names(data) <- c(seq_len(setParam$dgp$nTrain), 
                   paste0("test", seq_len(setParam$dgp$nTest)))
  
  #save to rda file
  fileName <- paste0("simDataN", N, "_pTrash", pTrash, "_rel", reliability, ".rda")
  save(data, file = paste0(dataFolder, "/", fileName))
}

# out <- do.call(mapply, c(FUN = createData, gridFull[1,])) # test
out <- do.call(mapply, c(FUN = createData, gridFull))
