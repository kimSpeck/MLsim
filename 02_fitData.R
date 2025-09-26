# Model fitting
library(parallel) # fit models in parallel across samples
library(glmnet)   # ENET
library(caret)    # GBM
library(gbm)      # GBM
library(ranger)   # random forests
library(iml)      # permutation variable importance, h-statistic

# load parameters & custom functions 
source("utils/setParameters.R") # parameter values
source("utils/simTools.R")

# functions to fit different models
source("utils/fitENET.R")
source("utils/fitGBM.R")
source("utils/fitRF.R")

# functions to calculate and save outcome measures from models
source("utils/saveENET.R")
source("utils/saveGBM.R")
source("utils/saveRF.R")

# build folder structure
dataFolder <- "data"

resFolder <- "results"
createFolder(resFolder)

createFolder(paste0(resFolder, "/inter"))
createFolder(paste0(resFolder, "/nonlinear3"))
createFolder(paste0(resFolder, "/pwlinear"))

logFolder <- "log"
createFolder(logFolder)

timeStamp <- format(Sys.time(), "%d%m%y_%H%M%S")
# nCoresSampling <- detectCores() - 1 
nCoresSampling <- 30

# get condGrid from parameter set 
condGrid <- setParam$fit$condGrid

setParam$dgp$nTrain <- 10
setParam$dgp$nSamples <- setParam$dgp$nTrain + setParam$dgp$nTest
# choose one single condition to test code/timing
condGrid <- condGrid[condGrid$data == "inter" &
                       condGrid$model == "ENETwo" &
                       condGrid$N == 100 &
                       condGrid$pTrash == 10 &
                       condGrid$reliability == 0.6,]


# flag to save single conditions (TRUE) or only full set of fitted results (FALSE)
setParam$fit$saveConds <- FALSE
# setParam$dgp$condLabels # check 

testFlag <- F

options(warn = -1)

results <- lapply(seq_len(nrow(condGrid)), function(iSim) {
  
  # Initiate cluster; type = "FORK" only on Linux/MacOS: contains all environment variables automatically
  cl <- makeCluster(nCoresSampling, type = "FORK",
                    outfile = paste0(logFolder, "/", "fitDataStatus",
                                     timeStamp, ".txt"))
  
  # set seed that works for parallel processing
  set.seed(condGrid[iSim, "sampleSeed"])
  s <- .Random.seed
  clusterSetRNGStream(cl = cl, iseed = s)
  
  # file name for all samples in one big rda file or folder name for single sample rda files
  fileName <- paste0(condGrid[iSim, "data"], 
                     "/simDataN", condGrid[iSim, "N"],
                     "_pTrash", condGrid[iSim, "pTrash"],
                     "_rel", condGrid[iSim, "reliability"])
  
  if (!setParam$dgp$singleSamples) { # load big rda file data if technical argument applies
    load(paste0(dataFolder, "/", fileName, ".rda"))
  }
  
  if (condGrid[iSim, "model"] == "RF") {
    nPredictors <- condGrid[iSim, "pTrash"] + setParam$dgp$p
    setParam$fit$tuneGrid_RF <- setParam$fit$setTuningGrid_RF(nPredictors)
  }
  
  # fit data (ENET, GBM, RF depending on fitting function; see below)
  tmp_estRes <- lapply(seq_along(setParam$dgp$condLabels), function(iCond) { # R2 x lin_inter combinations
      
    tStart <- Sys.time() # time monitoring
    # iterate over different training samples in parallel
    estRes <- parLapply(cl, seq_len(setParam$dgp$nTrain), function(iSample) {
      
      # monitor progress of the simulation study (irrespective of model)
      # write progress into log file to keep track of simulation progress
      if (iSample == 1) {
        print(paste0(fileName,
                     "; model: ", condGrid[iSim, "model"], 
                     "; iCond: ", setParam$dgp$condLabels[iCond], 
                     "; first Sample at ", format(Sys.time(), "%d%m%y_%H%M%S")))  
      }
      
      # load data for single sample rda files or extract relevant data from big rda file 
      if (setParam$dgp$singleSamples) { # data saved in single files
        # load single sample
        load(paste0(dataFolder, "/", fileName, "/sample_", iSample, ".rda"))
        
        # irrespective of model
        # get predictor variable matrix (identical across all R2 x lin_inter conditions)
        Xtrain <- as.matrix(dataList[["X_int"]]) 
        # get criterion/target variable (column in yMat depending on R2 x lin_inter condition)
        ytrain <- dataList[["yMat"]][,setParam$dgp$condLabels[iCond]]
        
        # load test data
        # by doing this in advane condGrid and full data not an additional argument for model functions
        # ... but, only if number of test samples is 1!
        # setParam$dgp$nTest == 1
        load(paste0(dataFolder, "/", fileName, "/sample_test1.rda"))
        
        Xtest <- as.matrix(dataList[["X_int"]]) 
        ytest <- dataList[["yMat"]][,setParam$dgp$condLabels[iCond]]
        
      } else { # data in one big rda file
        # irrespective of model
        # get predictor variable matrix (identical across all R2 x lin_inter conditions)
        Xtrain <- as.matrix(data[[iSample]][["X_int"]]) 
        # get criterion/target variable (column in yMat depending on R2 x lin_inter condition)
        ytrain <- data[[iSample]][["yMat"]][,setParam$dgp$condLabels[iCond]]
        
        # extract test data 
        # by doing this in advane condGrid and full data not an additional argument for model functions
        # ... but, only if number of test samples is 1!
        # setParam$dgp$nTest == 1
        Xtest <- as.matrix(data[["test1"]][["X_int"]]) 
        ytest <- data[["test1"]][["yMat"]][,setParam$dgp$condLabels[iCond]]
      }
      
      # if GBM, ENETwo or RF remove interactions from predictor matrix data (train & test) 
      if (condGrid[iSim, "model"] != "ENETw") {
        # train data
        idx_rmInter.train <- stringr::str_detect(colnames(Xtrain), ":")
        Xtrain <- Xtrain[,colnames(Xtrain)[!idx_rmInter.train]]
        # test data
        idx_rmInter.test <- stringr::str_detect(colnames(Xtest), ":")
        Xtest <- Xtest[,colnames(Xtest)[!idx_rmInter.test]]
      }
      
      # depending on condGrid[iSim, "model"] run GBM, RF or ENET code
      if (condGrid[iSim, "model"] == "GBM") {
        fitGBM(Xtrain, ytrain, Xtest, ytest, setParam, iSample = iSample, setParam$fit$explanation)
      } else if (condGrid[iSim, "model"] == "RF") {
        fitRF(Xtrain, ytrain, Xtest, ytest, setParam, iSample = iSample)
      } else if (condGrid[iSim, "model"] != "GBM" & condGrid[iSim, "model"] != "RF") {
        fitENET(Xtrain, ytrain, Xtest, ytest, setParam, iSample = iSample, setParam$fit$explanation)
      }
    
    }) # end of parallel fitting for samples
    tEnd <- Sys.time() # time monitoring
    difftime(tEnd, tStart)
    
    ##### save measures #####
    # model agnostic measures
    #   - training and test performance (RMSE, Rsquared, MAE)
    #   - permutation variable importance (pvi)
    
    # generate list for all results; write results in list immediately after generating! 
    if (condGrid[iSim, "model"] == "GBM") {
      resList <- vector(mode = "list", 
                        length = setParam$fit$out + setParam$fit$outGBM)
      names(resList) <- c(setParam$fit$outLabels, setParam$fit$gbmLabels)
    } else if (condGrid[iSim, "model"] == "RF") {
      resList <- vector(mode = "list", 
                        length = setParam$fit$out + setParam$fit$outRF)
      names(resList) <- c(setParam$fit$outLabels, setParam$fit$rfLabels)
    } else if (condGrid[iSim, "model"] != "GBM" & condGrid[iSim, "model"] != "RF") {
      resList <- vector(mode = "list", 
                        length = setParam$fit$out + setParam$fit$outENET)
      names(resList) <- c(setParam$fit$outLabels, setParam$fit$enetLabels)
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
    sampleVec <- do.call(rbind, lapply(estRes, function(X) X[["iSample"]])) # each sample
    resList[["performPerSample"]] <- cbind(performTrainMat, performTestMat, 
                                           idxCondLabel = iCond, sample = sampleVec)
    
    # PVI - permutation variable importance from every sample
    if (setParam$fit$explanation) {
      resList[["pvi"]] <- do.call(rbind, lapply(seq_len(length(estRes)), function(iSample) {
        cbind(idxCondLabel = rep(iCond, dim(estRes[[iSample]][["pvi"]])[1]),
              sample = iSample, 
              estRes[[iSample]][["pvi"]])
      }))
    } else {resList[["pvi"]] <- NA} 

    # add model specific measures to resList and return full list
    resList <- if (condGrid[iSim, "model"] == "GBM") {
      saveGBM(estRes, iCond, setParam, resList)
    } else if (condGrid[iSim, "model"] == "RF") {
      saveRF(estRes, iCond, setParam, resList)
    } else if (condGrid[iSim, "model"] != "GBM" & condGrid[iSim, "model"] != "RF") {
      saveENET(estRes, iCond, setParam, resList)
    }
    
    # if flag TRUE: save every single condition
    if (setParam$fit$saveConds) {
      resCondFileName <- paste0(resFolder, "/", condGrid[iSim, "data"], "/",
                            "resultsModel", condGrid[iSim, "model"], 
                            "_N", condGrid[iSim, "N"], 
                            "_pTrash", condGrid[iSim, "pTrash"], 
                            "_rel", condGrid[iSim,"reliability"], 
                            "_cond", iCond, ".rds")  
      saveRDS(resList, file = resCondFileName)
    }
    
    resList
  }) 
  
  # close cluster to return resources (memory) back to OS
  stopCluster(cl)
  
  # manage data 
  names(tmp_estRes) <- setParam$dgp$condLabels
  iCondRes <- do.call(Map, c(f = rbind, tmp_estRes)) 

  # save in nested loop across all N x pTrash conditions 
  
  # save separate each N x pTrash condition in rds files  
  #   -> pro: skip fitted condGrid conditions if algorithm does not run till the end... 
  #   -> con: writing data to rds file is slow and data needs to be united later on
  resFileName <- paste0(resFolder, "/", condGrid[iSim, "data"], "/",
                        "resultsModel", condGrid[iSim, "model"], 
                        "_N", condGrid[iSim, "N"], 
                        "_pTrash", condGrid[iSim, "pTrash"], 
                        "_rel", condGrid[iSim,"reliability"], ".rds")
  saveRDS(iCondRes, file = resFileName)
})

options(warn = 0)
