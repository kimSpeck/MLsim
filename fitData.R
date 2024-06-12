# Model fitting
# to do:
# fitting of only unfitted files does not work (modelGBM in filename now!)
# with 10 cores at simDataN300_pTrash50_rel0.6.rda; model: GBM; iCond: R20.5lin_inter0.2_0.8
# Error in serialize(data, node$con, xdr = FALSE) : error writing to connection
library(parallel) # fit models in parallel across samples
library(glmnet)   # ENET
library(caret)    # GBM
library(gbm)      # GBM
library(iml)      # permutation variable importance, h-statistic

# load parameters & custom functions 
source("setParameters.R") # parameter values
source("simTools.R")

# functions to fit different models
source("fitENET.R")
source("fitGBM.R")

# functions to calculate and save outcome measures from models
source("saveENET.R")
source("saveGBM.R")

# build folder structure
dataFolder <- "data"

resFolder <- "results"
createFolder(resFolder)

logFolder <- "log"
createFolder(logFolder)

timeStamp <- format(Sys.time(), "%d%m%y_%H%M%S")
# nCoresSampling <- detectCores() - 1 
nCoresSampling <- 10
# 10 cores seem to be max for ENET memory usage wise

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

results <- lapply(seq_len(nrow(condGrid)), function(iSim) {
  
  # test it!
  # iSim = 37 # ENETwo x 100 x 10 x 1.0
  # iSim = 38 # ENETw x 100 x 10 x 1.0
  # iSim = 39 # GBM x 100 x 10 x 1.0
  
  fileName <- paste0("simDataN", condGrid[iSim, "N"],
                     "_pTrash", condGrid[iSim, "pTrash"],
                     "_rel", condGrid[iSim, "reliability"], ".rda")
  
  load(paste0(dataFolder, "/", fileName))
  
  # fitte eine regularisierte Regression
  tmp_estRes <- lapply(seq_along(setParam$dgp$condLabels), function(iCond) { # R2 x lin_inter combinations
  
    # test it: 
    # iCond <- 1
    
    # Initiate cluster; type = "FORK" only on Linux/MacOS: contains all environment variables automatically
    cl <- makeCluster(nCoresSampling, type = "FORK",
                      outfile = paste0(logFolder, "/", "fitDataStatus",
                                       timeStamp, ".txt"))
    
    # set seed that works for parallel processing
    set.seed(condGrid[iSim, "sampleSeed"])
    s <- .Random.seed
    clusterSetRNGStream(cl = cl, iseed = s)
      
    # estRes <- lapply(seq_len(setParam$dgp$nTrain), function(iSample) {
    # estRes <- parLapply(cl, seq_len(setParam$dgp$nTrain), function(iSample) {
    tStart <- Sys.time()
    estRes <- parLapply(cl, seq_len(10), function(iSample) {
      
      # # test it; set lapply variables to any possible value 
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
    tEnd <- Sys.time()
    difftime(tEnd, tStart)
    
    # close cluster to return resources (memory) back to OS
    stopCluster(cl)
    
    ##### save measures #####
    # model agnostic measures
    #   - training and test performance (RMSE, Rsquared, MAE)
    #   - permutation variable importance (pvi)
    
    # generate list for all results; write results in list immediatly after generating! 
    if (condGrid[iSim, "model"] != "GBM") {
      resList <- vector(mode = "list", 
                        length = setParam$fit$out + setParam$fit$outENET)
      names(resList) <- c(setParam$fit$outLabels, setParam$fit$enetLabels)
    } else if (condGrid[iSim, "model"] == "GBM") {
      resList <- vector(mode = "list", 
                        length = setParam$fit$out + setParam$fit$outGBM)
      names(resList) <- c(setParam$fit$outLabels, setParam$fit$gbmLabels)
    }
    
    # training performance (RMSE, Rsquared, MAE)
    performTrainMat <- do.call(rbind, lapply(estRes, function(X) X[["performTrain"]])) # each sample
    colnames(performTrainMat) <- c("RMSE_train", "Rsq_train", "MAE_train")
    # aggregated performance measures (M, SD, SE)
    resList[["performTrainStats"]] <- cbind(getStats(performTrainMat, 2, setParam$dgp$nSamples),
                          idxCondLabel = iCond)
    rownames(resList[["performTrainStats"]]) <- c("RMSE", "Rsquared", "MAE")
    
    
    # test performance (RMSE, Rsquared, MAE)
    performTestMat <- do.call(rbind, lapply(estRes, function(X) X[["performTest"]])) # each sample
    colnames(performTestMat) <- c("RMSE_test", "Rsq_test", "MAE_test")
    # aggregated performance measures (M, SD, SE)
    resList[["performTestStats"]] <- cbind(getStats(performTestMat, 2, setParam$dgp$nSamples),
                          idxCondLabel = iCond)
  
    
    # test & train performance for each sample
    resList[["performPerSample"]] <- cbind(performTrainMat, performTestMat, idxCondLabel = iCond)
    
    # PVI - permutation variable importance from every sample
    resList[["pvi"]] <- do.call(rbind, lapply(seq_len(length(estRes)), function(iSample) {
      cbind(idxCondLabel = rep(iCond, dim(estRes[[iSample]][["pvi"]])[1]),
            sample = iSample, 
            estRes[[iSample]][["pvi"]])
    }))
    
    resList <- if (condGrid[iSim, "model"] != "GBM") {
      saveENET(estRes, iCond, setParam, resList)
    } else if (condGrid[iSim, "model"] == "GBM") {
      saveGBM(estRes, iCond, setParam, resList)
    }
    
    resList
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
  resFileName <- paste0(resFolder, "/", 
                        "resultsModel", condGrid[iSim, "model"], 
                        "_N", condGrid[iSim, "N"], 
                        "_pTrash", condGrid[iSim, "pTrash"], 
                        "_rel", condGrid[iSim,"reliability"], ".rds")
  saveRDS(iCondRes, file = resFileName)
})





