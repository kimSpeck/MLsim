# simulate data
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
                        pTrash = setParam$dgp$pTrash)
str(gridFull) # ! factors are interpreted as level numbers!; only character variables are interpreted by their name!
# create seed number for parallel cluster (reproducibility of generated data)
set.seed(8967369)
seedNum <- sample(1:999999, dim(gridFull)[1], replace = FALSE) 
gridFull$sampleSeed <- seedNum[1:dim(gridFull)[1]]

# ToDo: write big matrices directly to their final place of storage (memory intensive and time consuming)
# ToDo: depending in seed we receive non PSM matrices sigma as covariance matrices in multivariate normal
#     Warnmeldungen: 
#       1: In rmvnorm(n = N, mean = mP, sigma = rX) :
#       sigma is numerically not positive semidefinite
# however, rmvnorm somehow deals with this issue (see: https://stats.stackexchange.com/questions/267908/clarification-regarding-rmvnorm-in-r)

# sample data in parallel
createData <- function(N, pTrash, sampleSeed){
  
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
    
    P <- setParam$dgp$p + pTrash # total number of variables
    # generate matrix of (almost) uncorrelated predictors
    X <- createPredictors(N = N, P = P, 
                          mR = setParam$dgp$meanR, sdR = setParam$dgp$sdR)
    
    # add names to variables
    colnames(X) <- paste0("Var", seq_len(P))
    
    # create model formula (allows polynomial and interaction effects of any degree/depth)
    popModel <- genModel(colnames(X), setParam$dgp$interDepth, setParam$dgp$poly)
    
    # predictor matrix that allows for polynomials and interactions
    X_int <- model.matrix(as.formula(popModel),data.frame(X))
    
    # remove first degree polynomials from data (they are duplicates!)
    # colnames(X_int)[stringr::str_detect(colnames(X_int), pattern = "^(poly\\().+(\\)1)$")]
    rmColsIdx <- which(stringr::str_detect(colnames(X_int), pattern = "^(poly\\().+(\\)1)$"))
    X_int <- X_int[,-rmColsIdx] 
    
    # ToDo: generate regression coefficients, check that R^2 matches fixed R^2
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
    bMatrix[c(linEffects, interEffects), 1] <- rep(0.116, times = nEffects) # R2 = 0.1
    bMatrix[c(linEffects, interEffects), 2] <- rep(0.227, times = nEffects) # R2 = 0.3
    bMatrix[c(linEffects, interEffects), 3] <- rep(0.346, times = nEffects) # R2 = 0.5
    bMatrix[c(linEffects, interEffects), 4] <- rep(0.693, times = nEffects) # R2 = 0.8
    
    # 80 - 20 (most effects in linEffects)
    # setParam$dgp$Rsquared/100*80
    # setParam$dgp$Rsquared/100*20
    
    # R^2 = 0.1
    bMatrix[c(linEffects), 5] <- rep(0.144, times = length(linEffects)) # R2 = 0.08
    bMatrix[c(interEffects), 5] <- rep(0.079, times = length(interEffects)) # R2 = 0.02
    # R^2 = 0.3
    bMatrix[c(linEffects), 6] <- rep(0.28, times = length(linEffects)) # R2 = 0.24
    bMatrix[c(interEffects), 6] <- rep(0.15, times = length(interEffects)) # R2 = 0.06
    # R^2 = 0.5
    bMatrix[c(linEffects), 7] <- rep(0.43, times = length(linEffects)) # R2 = 0.4
    bMatrix[c(interEffects), 7] <- rep(0.235, times = length(interEffects)) # R2 = 0.1
    # R^2 = 0.8
    bMatrix[c(linEffects), 8] <- rep(0.85, times = length(linEffects)) # R2 = 0.64
    bMatrix[c(interEffects), 8] <- rep(0.45, times = length(interEffects)) # R2 = 0.16
    
    # 20 - 80 (most effects in interEffects)
    # R^2 = 0.1
    bMatrix[c(interEffects), 9] <- rep(0.144, times = length(interEffects)) # R2 = 0.08
    bMatrix[c(linEffects), 9] <- rep(0.079, times = length(linEffects)) # R2 = 0.02
    # R^2 = 0.3
    bMatrix[c(interEffects), 10] <- rep(0.28, times = length(interEffects)) # R2 = 0.24
    bMatrix[c(linEffects), 10] <- rep(0.15, times = length(linEffects)) # R2 = 0.06
    # R^2 = 0.5
    bMatrix[c(interEffects), 11] <- rep(0.43, times = length(interEffects)) # R2 = 0.4
    bMatrix[c(linEffects), 11] <- rep(0.235, times = length(linEffects)) # R2 = 0.1
    # R^2 = 0.8
    bMatrix[c(interEffects), 12] <- rep(0.85, times = length(interEffects)) # R2 = 0.64
    bMatrix[c(linEffects), 12] <- rep(0.45, times = length(linEffects)) # R2 = 0.16
    
    # calculate R^2 for every combination of R2 and lin/inter effect balance
    # print R^2 as a quick sanity check (removed for speed sake) 
    R2 <- sapply(seq_len(ncol(bMatrix)), function(x) getR2(X_int, bMatrix[,x], setParam$dgp$sigmaE))
    
    # calculate dependent variable for every combination of R2 and lin/inter effect balance
    yMatrix <- sapply(seq_len(ncol(bMatrix)), function(x) {
      calcDV(X = X_int, b = bMatrix[,x],
             sigmaE = setParam$dgp$sigmaE, N = N)
    })
    colnames(yMatrix) <- setParam$dgp$condLabels
    
    # save ...
    #     ... yMat with dependent variable for all R2 - lin/inter effect conditions in columns
    #     ... X_int predictor matrix (identical for all R2 - lin/inter effect conditions)
    #     ... trueB simulated regression coefficients
    #     ... R2 based on simulated regression coefficients
    list(yMat = yMatrix, 
         X_int = X_int, 
         trueB = bMatrix, 
         R2 = R2)
  })
  
  #save to rda file
  fileName <- paste0("simDataN", N, "_pTrash", pTrash, ".rda")
  save(data, file = paste0(dataFolder, "/", fileName))
}

# out <- do.call(mapply, c(FUN = createData, gridFull[1,])) # test
out <- do.call(mapply, c(FUN = createData, gridFull))
