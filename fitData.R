# Model fitting
library(parallel) # fit models in parallel across samples
library(glmnet)   # ENET
library(caret)    # GBM
library(gbm)      # GBM
library(iml)      # permutation variable importance, h-statistic

# load parameters & custom functions 
source("setParameters.R") # parameter values
source("simTools.R")

source("fitENET.R")
source("fitGBM.R")

# build folder structure
dataFolder <- "data"

resFolder <- "results"
createFolder(resFolder)

logFolder <- "log"
createFolder(logFolder)

timeStamp <- format(Sys.time(), "%d%m%y_%H%M%S")
# nCoresSampling <- detectCores() - 1 
nCoresSampling <- 10

# get condGrid from parameter set 
condGrid <- setParam$fit$condGrid
# remove lines in condGrid for which there already are results
# check which files are already in the results folder
resFileList <- list.files(resFolder)
resFileList <- sub('results', "", resFileList); resFileList <- sub('.rds', "", resFileList)

allDataFiles <- paste0("simDataN", condGrid[, "N"], 
                       "_pTrash", condGrid[, "pTrash"], 
                       "_rel", condGrid[,"reliability"], ".rda")
allDataFiles <- sub("simData", "", allDataFiles); allDataFiles <- sub(".rda", "", allDataFiles)

fitIdx <- which(!(allDataFiles %in% resFileList)) # index of conditions that need to be fitted 
condGrid <- condGrid[fitIdx, ] # remove conditions that are already done

# Initiate cluster; type = "FORK" only on Linux/MacOS: contains all environment variables automatically
cl <- makeCluster(nCoresSampling, type = "FORK",
                  outfile = paste0(logFolder, "/", "fitDataStatus",
                                   timeStamp, ".txt"))

results <- lapply(seq_len(nrow(condGrid)), function(iSim) {
  
  # test it!
  # iSim = 37 # ENETwo x 100 x 10 x 1.0
  # iSim = 38 # ENETw x 100 x 10 x 1.0
  # iSim = 39 # GBM x 100 x 10 x 1.0
  
  fileName <- paste0("simDataN", condGrid[iSim, "N"],
                     "_pTrash", condGrid[iSim, "pTrash"],
                     "_rel", condGrid[iSim, "reliability"], ".rda")
  
  load(paste0(dataFolder, "/", fileName))
  
  # set seed that works for parallel processing
  set.seed(condGrid[iSim, "sampleSeed"])
  s <- .Random.seed
  clusterSetRNGStream(cl = cl, iseed = s)
  
  # fitte eine regularisierte Regression
  tmp_estRes <- lapply(seq_along(setParam$dgp$condLabels), function(iCond) { # R2 x lin_inter combinations
    
    # estRes <- lapply(seq_len(setParam$dgp$nTrain), function(iSample) {
    estRes <- parLapply(cl, seq_len(setParam$dgp$nTrain), function(iSample) {
      
      # # test it; set lapply variables to any possible value 
      # iCond <- 1
      # iSample <- 1
      
      # monitor progress of the simulation study (irrespective of model)
      # write progress into log file to keep track of simulation progress
      if (iSample == 1) {
        print(paste0(fileName,
                     "; model: ", condGrid[iSim, "model"], 
                     "; iCond: ", setParam$dgp$condLabels[iCond], 
                     "; first Sample at ", format(Sys.time(), "%d%m%y_%H%M%S")))  
      }
      
      # irrespective of model
      # get predictor variable matrix (identical across all R2 x lin_inter conditions)
      Xtrain <- as.matrix(data[[iSample]][["X_int"]]) 
      # get criterion/target variable (column in yMat depending on R2 x lin_inter condition)
      ytrain <- data[[iSample]][["yMat"]][,iCond]
      
      # extract test data 
      # by doing this in advane condGrid and full data not an additional argument for model functions
      # ... but, only if number of test samples is 1!
      # setParam$dgp$nTest == 1
      Xtest <- as.matrix(data[["test1"]][["X_int"]]) 
      ytest <- data[["test1"]][["yMat"]][,iCond] 
      
      # if GBM or ENETwo remove interactions from predictor matrix data (train & test) 
      if (condGrid[iSim, "model"] != "ENETw") {
        # train data
        idx_rmInter.train <- stringr::str_detect(colnames(Xtrain), ":")
        Xtrain <- Xtrain[,colnames(Xtrain)[!idx_rmInter.train]]
        # test data
        idx_rmInter.test <- stringr::str_detect(colnames(Xtest), ":")
        Xtest <- Xtest[,colnames(Xtest)[!idx_rmInter.test]]
      }
      
      # depending on condGrid[iSim, "model"] run ENET or GBM code
      if (condGrid[iSim, "model"] != "GBM") {
        fitENET(Xtrain, ytrain, Xtest, ytest, setParam)
      } else if (condGrid[iSim, "model"] == "GBM") {
        fitGBM(Xtrain, ytrain, Xtest, ytest, setParam)
      }
    
    }) # end of parallel fitting for samples
    
    # coefficients
    estBeta <- do.call(cbind, lapply(estRes, function(X) X[["estB"]]))
    # remove variables which are not selected by Enet from mean calculation:
    #   coefficients that are exactly 0 -> NA
    estBeta <- ifelse(estBeta == 0, NA, estBeta)
    
    estBetaStats <- getStats(estBeta, 1, setParam$dgp$nSamples)
    estBetaStats <- cbind(estBetaStats, idxCondLabel = iCond)
    
    # selected variables (how often are predictors selected in model)
    selectedVars <- do.call(cbind, lapply(estRes, function(X) X[["selectedVars"]]))
    
    # AV: Wie oft bleiben Terme (lineare Praediktoren/Interaktionen) im Modell?
    nSelection <- apply(selectedVars, MARGIN = 1, sum) # frequency of variable selection 
    percSelection <- nSelection / ncol(selectedVars) * 100 # relative frequency 
    varSelection <- cbind(nSelection = nSelection, 
                          percSelection = percSelection,
                          idxCondLabel = iCond) # stats per predictor
    
    # here! to do! only linear and interaction effects extracted; not indicators!
    
    # only true predictors in model (8 predictors with simulated effects and every other variable == 0)
    # how many of the linear effects are recovered?
    nSelectedLin <- sapply(seq_len(ncol(selectedVars)), function(iCol) { 
      sum(setParam$dgp$linEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
    # how many of the interaction effects are recovered?
    nSelectedInter <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
      sum(setParam$dgp$interEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
    
    # how many of the single indicators were recovered?
    nSelectedInd <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
      sum(setParam$dgp$indEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
    
    # how many of the single indicator interactions were recovered?
    nSelectedIndInter <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
      sum(setParam$dgp$indInterEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
    
    # all simulated effects selected in model?
    selectedAll <- nSelectedLin + nSelectedInter == length(c(setParam$dgp$linEffects, setParam$dgp$linEffects))
    selectedAllInd <- nSelectedInd + nSelectedIndInter == length(c(setParam$dgp$indEffects, setParam$dgp$indInterEffects))
    
    # only simulated effects selected (i.e., every other predictor is not selected!)
    nSelectedOthers <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
      sum(!(rownames(selectedVars)[selectedVars[,iCol]] %in% c(setParam$dgp$interEffects, setParam$dgp$linEffects)))})
    
    nSelectedOthersInd <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
      sum(!(rownames(selectedVars)[selectedVars[,iCol]] %in% c(setParam$dgp$indInterEffects, setParam$dgp$indEffects)))})
    
    
    # final Enet-Parameters for each sample
    tunedAlphaVec <- unname(do.call(c, lapply(estRes, function(X) X[["tunedAlpha"]])))
    tunedLambdaVec <- unname(do.call(c, lapply(estRes, function(X) X[["tunedLambda"]])))
    
    # variable selection stats per sample
    # for all: 1 = TRUE, 0 = FALSE
    selectionPerSample <- cbind(nLin = nSelectedLin, nInter = nSelectedInter, 
                                nInd = nSelectedInd, nIndInter = nSelectedIndInter,
                                all.T1F0 = selectedAll, all.Ind = selectedAllInd,
                                nOthers = nSelectedOthers, nOthersInd = nSelectedOthersInd,
                                alpha = tunedAlphaVec, lambda = tunedLambdaVec,
                                idxCondLabel = iCond)
    
    # training performance (RMSE, Rsquared, MAE)
    performTrainMat <- do.call(rbind, lapply(estRes, function(X) X[["performTrain"]])) # each sample
    colnames(performTrainMat) <- c("RMSE_train", "Rsq_train", "MAE_train")
    performTrainStats <- getStats(performTrainMat, 2, setParam$dgp$nSamples) # M, SD, SE
    rownames(performTrainStats) <- c("RMSE", "Rsquared", "MAE")
    performTrainStats <- cbind(performTrainStats, idxCondLabel = iCond)
    
    # test performance (RMSE, Rsquared, MAE)
    performTestMat <- do.call(rbind, lapply(estRes, function(X) X[["performTest"]])) # each sample
    colnames(performTestMat) <- c("RMSE_test", "Rsq_test", "MAE_test")
    performTestStats <- getStats(performTestMat, 2, setParam$dgp$nSamples) # M, SD, SE
    performTestStats <- cbind(performTestStats, idxCondLabel = iCond)
    
    perfromPerSample <- cbind(performTrainMat, performTestMat, idxCondLabel = iCond)
    
    # join results 
    list(estBeta = estBetaStats, # M, SD, SE for coefficients across samples
         estBetaFull = estBeta, # coefficients for every predictor in every sample 
         varSelection = varSelection, 
         selectSample = selectionPerSample, 
         perfromPerSample = perfromPerSample, 
         performTrainStats = performTrainStats, 
         performTestStats = performTestStats)
  }) # nested lapplys take ~ 15 mins
  
  # manage data 
  names(tmp_estRes) <- setParam$dgp$condLabels
  iCondRes <- do.call(Map, c(f = rbind, tmp_estRes)) # das?
  
  # clean up working memory
  rm(data)
  gc()
  
  # save in nested loop across all N x pTrash conditions 
  #   -> does this fill up wm leading the algorithm to slow down or does it only 
  #       become slower due to larger number of predictors in later conditions?
  #       sieht eigentlich nicht danach aus, wenn man Zeiten in log file und wm load betrachtet
  # iCondRes 
  
  # save separate each N x pTrash condition in rds files  
  #   -> pro: skip fitted condGrid conditions if algorithm does not run till the end... 
  #   -> con: writing data to rds file is slow and data needs to be united later on
  resFileName <- paste0(resFolder, "/", "resultsModel", condGrid[iSim, "model"], 
                        "_N", condGrid[iSim, "N"], 
                        "_pTrash", condGrid[iSim, "pTrash"], 
                        "_rel", condGrid[iSim,"reliability"], ".rds")
  saveRDS(iCondRes, file = resFileName)
})

# close cluster to return resources (memory) back to OS
stopCluster(cl)



