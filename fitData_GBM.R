# Model fitting
# only for GBM (quick and dirty)

# GBM packages
library(xgboost)
library(iml)
# library(EIX)

# elastic net packages
# library(glmnet)

# parallelisation
library(parallel)

# load parameters & custom functions 
source("setParameters.R") # parameter values
source("simTools.R")

dataFolder <- "data"

resFolder <- "results"
createFolder(resFolder)

logFolder <- "log"
createFolder(logFolder)

timeStamp <- format(Sys.time(), "%d%m%y_%H%M%S")
# nCoresSampling <- detectCores() - 1 
nCoresSampling <- 10

# test it!
# iSim = 1

# # different modes for Enet
# ## variables included in Enet (TRUE = include poly/inter; FALSE = without poly/inter)
# includePoly <- FALSE
# includeInter <- TRUE # run Enet with and without interactions

# run GBM without interactions!
includePoly <- FALSE
includeInter <- FALSE

## warm start to arrive at lambda parameters for cross validation
warmStart <- setParam$fit$warmStart

# iterate through these combinations of data conditions
# only use data without simulated factors
condGrid <- expand.grid(N = setParam$dgp$N, 
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability,
                        factors = FALSE)

# remove lines in condGrid for which there already are results
# check which files are already in the results folder
resFileList <- list.files(resFolder)
resFileList <- sub('results', "", resFileList); resFileList <- sub('.rds', "", resFileList)

allDataFiles <- paste0("simDataN", condGrid[, "N"], 
                       "_pTrash", condGrid[, "pTrash"], 
                       "_rel", condGrid[,"reliability"], 
                       "_f", ifelse(condGrid[,"factors"], 1, 0), ".rda")
allDataFiles <- sub("simData", "", allDataFiles); allDataFiles <- sub(".rda", "", allDataFiles)

fitIdx <- which(!(allDataFiles %in% resFileList)) # index of conditions that need to be fitted 
condGrid <- condGrid[fitIdx, ] # remove conditions that are already done

# Initiate cluster; type = "FORK" only on Linux/MacOS: contains all environment variables automatically
cl <- makeCluster(nCoresSampling, type = "FORK",
                  outfile = paste0(logFolder, "/", "fitDataStatus",
                                   timeStamp, ".txt"))

# set seed that works for parallel processing
set.seed(342890)
s <- .Random.seed
clusterSetRNGStream(cl = cl, iseed = s)

results <- lapply(seq_len(nrow(condGrid)), function(iSim) {
  
  fileName <- paste0("simDataN", condGrid[iSim, "N"],
                     "_pTrash", condGrid[iSim, "pTrash"],
                     "_rel", condGrid[iSim, "reliability"], 
                     "_f", ifelse(condGrid[iSim, "factors"], 1, 0), ".rda")
  # test it
  fileName <- "simDataN1000_pTrash10_rel1_f0.rda"
  load(paste0(dataFolder, "/", fileName))
  
  # fitte eine regularisierte Regression
  # R2 x lin_inter combinations
  tmp_estRes <- lapply(seq_along(setParam$dgp$condLabels), function(iCond) { 
    
    # test it
    iCond <- 3
    estRes <- lapply(seq_len(2), function(iSample) {
    # estRes <- lapply(seq_len(setParam$dgp$nTrain), function(iSample) {
    # estRes <- parLapply(cl, seq_len(setParam$dgp$nTrain), function(iSample) {
      
      if (iSample == 1) {
        print(paste0(fileName,
                     "; iCond: ", setParam$dgp$condLabels[iCond], 
                     "; first Sample at ", format(Sys.time(), "%d%m%y_%H%M%S")))  
      }
      
      # # test it
      # iCond <- 1
      # iSample <- 1
      
      tstart <- Sys.time()
      
      # get predictors
      Xtrain <- as.matrix(data[[iSample]][["X_int"]])
      
      # remove polynomials for Enet and GBM 
      if (!includePoly) {
        idx_rmPoly <- stringr::str_detect(colnames(Xtrain), "^poly")
        Xtrain <- Xtrain[,colnames(Xtrain)[!idx_rmPoly]]
      }
      # remove interactions from predictor matrix in GBM (and depending on condition in Enet)
      if (!includeInter) {
        idx_rmInter <- stringr::str_detect(colnames(Xtrain), ":")
        Xtrain <- Xtrain[,colnames(Xtrain)[!idx_rmInter]]
      }
      
      # this changes the simulated R2 x lin_inter condition we are analysing
      ytrain <- data[[iSample]][["yMat"]][,iCond]
      
      # input data with xgboosts own class (recommended)
      trainData <- xgboost::xgb.DMatrix(data = Xtrain, label = ytrain, missing = NA)
      
      # copy tuning grid from setParam
      tuneGrid <- setParam$fit$tuneGrid
      
      # cross validation loop for training model
      for (iGline in seq_len(nrow(tuneGrid))) {
        
        # reproducibility
        set.seed(3829)
        
        # train model
        fit_cv <- xgb.cv(
          data = trainData,                 # X & y in xgboost data format
          nrounds = tuneGrid$nTrees[iGline],# max nTrees
          nfold = setParam$fit$nfolds,      # k-fold cross validation
          eta = tuneGrid$eta[iGline],
          max_depth = tuneGrid$max_depth[iGline],
          min_child_weight = tuneGrid$min_child_weight[iGline],
          objective = "reg:squarederror",   # for regression models
          verbose = FALSE,                  # silent
          early_stopping_rounds = 10)       # stop running if the cross validated error 
                                            # does not improve for n continuous trees
        
        # identify the optimal number of trees and the corresponding minimum RMSE for the cv 
        # add number of Trees and corresponding training error to grid
        tuneGrid$optimalTrees[iGline] <- which.min(fit_cv$evaluation_log$test_rmse_mean)
        tuneGrid$minRMSE[iGline] <- min(fit_cv$evaluation_log$test_rmse_mean)
        
      }
      
      # save tuning parameters! (results from cross validation)
      # cross-validated tuning parameters
      idxOptim <- which.min(tuneGrid$minRMSE)
      
      # to do: save tuned parameters in output
      # Best parameters
      tunedEta <- tuneGrid$eta[idxOptim]
      tunedMax_depth <- tuneGrid$max_depth[idxOptim]
      tunedMin_child_weight <- tuneGrid$min_child_weight[idxOptim]
      tunedNrounds <- tuneGrid$optimalTrees[idxOptim]
      
      # fit model on training data with tuned hyperparameters
      fit <- xgboost(
        data = trainData,                 # X & y in xgboost data format
        eta = tuneGrid$eta[idxOptim],
        max_depth = tuneGrid$max_depth[idxOptim],
        min_child_weight = tuneGrid$min_child_weight[idxOptim],
        nrounds = tuneGrid$optimalTrees[idxOptim],
        objective = "reg:squarederror",   # for regression models
        verbose = FALSE)
      
      # calculate dependent measures
      
      # in Enet:
      #     - beta coefficients
      # estBeta <- as.matrix(fit$beta)
      #
      #     - selected predictors
      # selectedVars <- estBeta != 0
      
      # # model agnostic measures:
      #     - train performance 
      # to do: write this into a function to reuse it for gbm code as well (?)
      # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model training
      # performance test for train data
      predTrain <- predict(fit, Xtrain)
      performTrain <- evalPerformance(predTrain, ytrain)
      performTrain <- matrix(performTrain, ncol = length(performTrain), nrow = 1)
      colnames(performTrain) <- c("RMSE", "Rsquared", "MAE")
      
      #     - test performance
      # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model testing
      # test data
      # performance test for test data
      # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model testing
      performTest <- lapply(paste0("test", seq_len(setParam$dgp$nTest)), function(iTest) {
        Xtest <- as.matrix(data[[iTest]][["X_int"]])
        if (!includePoly) {
          idx_rmPoly <- stringr::str_detect(colnames(Xtest), "^poly")
          Xtest <- Xtest[,colnames(Xtest)[!idx_rmPoly]]
        }
        if (!includeInter) {
          idx_rmInter <- stringr::str_detect(colnames(Xtest), ":")
          Xtest <- Xtest[,colnames(Xtest)[!idx_rmInter]]
        }
        
        ytest <- data[[iTest]][["yMat"]][,iCond]
        pred <- predict(fit, Xtest)  
        
        evalPerformance(pred, ytest)
      })
      performTest <- do.call(rbind, performTest)
      
      # multiple dependent measures from the iml package
      #   prepare data
      Xtest <- as.matrix(data[["test1"]][["X_int"]])
      if (!includePoly) {
        idx_rmPoly <- stringr::str_detect(colnames(Xtest), "^poly")
        Xtest <- Xtest[,colnames(Xtest)[!idx_rmPoly]]
      }
      if (!includeInter) {
        idx_rmInter <- stringr::str_detect(colnames(Xtest), ":")
        Xtest <- Xtest[,colnames(Xtest)[!idx_rmInter]]
      }
      
      ytest <- data[["test1"]][["yMat"]][,iCond]
      
      predictor.gbm <- Predictor$new(
        model = fit, 
        data = data.frame(Xtrain), 
        y = ytrain, 
        predict.fun = predGBM,
        type = "prob")
      
      #     - permutation variable importance
      #       increase n.repetitions? (How often should the shuffling of the feature 
      #       be repeated? The higher the number of repetitions the more stable 
      #       and accurate the results become.)
      imp.gbm <- FeatureImp$new(predictor.gbm, loss = "mse")
      pviRank <- imp.gbm$results$feature
      pviValue <- imp.gbm$results$importance
      
      #     - evaluate interaction strength with H-statistic
      # Molnar, C. (2020). Interpretable machine learning. Lulu. com.
      # Friedman & Popescu (2005). 
      #   only for variables with simulated main effects otherwise interaction strength is overestimated
      #   see also Greenwell et al. (2018) and Henninger et al. (2023) 
      # interStrength <- lapply(paste0("Var", seq_len(setParam$dgp$p)), function(iVar) {
      #   tmp <- iml::Interaction$new(predictor.gbm, feature = iVar, 
      #                          grid.size = setParam$fit$nInterStrength)
      #   idxTmp <- sort(tmp$results$.interaction, decreasing = T)
      #   cbind(var = rep(iVar, length(idxTmp)),
      #         feature = tmp$results$.feature[match(idxTmp, tmp$results$.interaction)],
      #         interaction = tmp$results$.interaction[match(idxTmp, tmp$results$.interaction)])
      #   
      # })
      # interStrength <- do.call(rbind, interStrength)
      
      tend <- Sys.time()
      difftime(tend, tstart)
      # timing: 
      #   ! lower estimate data has only 10 trash variables (simDataN1000_pTrash10_rel1_f0.rda)
      #   1 sample with interStrength (setParam$fit$nInterStrength = 100) = 2.098334 mins
      #     -> 2 min x 100 samples = 200 Minuten = 3h 20min pro Bedingung 
      #     -> 3h 20min x 18 Bedingungen = 60h
      #   1 sample without interStrength = 23.39594 secs
      #     -> 24 sec x 100 samples = 40 Minuten 
      #     -> 40 Minuten x 18 Bedingungen = 12h
      
      # save dependent variables for each sample (in a list)
      list(#estB = estBeta, 
           #selectedVars = selectedVars,
           performTrain = performTrain, # train performance
           performTest = performTest, # test performance
           # model agnostic measures
           pvi = cbind(pviRank, pviValue), # permutation variable importance
           interStrength = interStrength, # h-statistic two-way interactions
           # tuning parameters (GBM)
           tunedEta = tunedEta,
           tunedMax_depth = tunedMax_depth,
           tunedMin_child_weight = tunedMin_child_weight,
           tunedNrounds = tunedNrounds)
      
      # # remove fitted sample to reduce working memory load
      # data[[iSample]] <- NULL
    }) # end of loop across samples
    
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
    
    # model agnostic measures
    ## PVI
    pvi <- do.call(rbind, lapply(seq_len(length(estRes)), function(iSample) {
      cbind(idxCondLabel = rep(iCond, dim(estRes[[iSample]][["pvi"]])[1]),
            sample = iSample, 
            estRes[[iSample]][["pvi"]])
    }))
    
    ## H-Statistic
    interStrength <- do.call(rbind, lapply(seq_len(length(estRes)), function(iSample) {
      cbind(idxCondLabel = rep(iCond, dim(estRes[[iSample]][["interStrength"]])[1]),
            sample = iSample,
            estRes[[iSample]][["interStrength"]])
    }))
    
    # tuning parameters (GBM)
    tunedEta <- unname(do.call(c, lapply(estRes, 
                                         function(X) X[["tunedEta"]])))
    tunedMax_depth <- unname(do.call(c, lapply(estRes, 
                                               function(X) X[["tunedMax_depth"]])))
    tunedMin_child_weight <- unname(do.call(c, lapply(estRes, 
                                                      function(X) X[["tunedMin_child_weight"]]))) 
    tunedNrounds <- unname(do.call(c, lapply(estRes, 
                                             function(X) X[["tunedNrounds"]])))
    
    selectionPerSample <- cbind(eta = tunedEta, max_depth = tunedMax_depth,
                                min_child_weight = tunedMin_child_weight, 
                                Nrounds = tunedNrounds,
                                idxCondLabel = iCond)
    # # coefficients
    # estBeta <- do.call(cbind, lapply(estRes, function(X) X[["estB"]]))
    # # remove variables which are not selected by Enet from mean calculation:
    # #   coefficients that are exactly 0 -> NA
    # estBeta <- ifelse(estBeta == 0, NA, estBeta)
    # 
    # estBetaStats <- getStats(estBeta, 1, setParam$dgp$nSamples)
    # estBetaStats <- cbind(estBetaStats, idxCondLabel = iCond)
    # 
    # # selected variables (how often are predictors selected in model)
    # selectedVars <- do.call(cbind, lapply(estRes, function(X) X[["selectedVars"]]))
    # 
    # # AV: Wie oft bleiben Terme (lineare Praediktoren/Interaktionen) im Modell?
    # nSelection <- apply(selectedVars, MARGIN = 1, sum) # frequency of variable selection 
    # percSelection <- nSelection / ncol(selectedVars) * 100 # relative frequency 
    # varSelection <- cbind(nSelection = nSelection, 
    #                       percSelection = percSelection,
    #                       idxCondLabel = iCond) # stats per predictor
    # 
    # # here! to do! only linear and interaction effects extracted; not indicators!
    # 
    # # only true predictors in model (8 predictors with simulated effects and every other variable == 0)
    # # how many of the linear effects are recovered?
    # nSelectedLin <- sapply(seq_len(ncol(selectedVars)), function(iCol) { 
    #   sum(setParam$dgp$linEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
    # # how many of the interaction effects are recovered?
    # nSelectedInter <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
    #   sum(setParam$dgp$interEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
    # 
    # # how many of the single indicators were recovered?
    # nSelectedInd <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
    #   sum(setParam$dgp$indEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
    # 
    # # how many of the single indicator interactions were recovered?
    # nSelectedIndInter <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
    #   sum(setParam$dgp$indInterEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
    # 
    # # all simulated effects selected in model?
    # selectedAll <- nSelectedLin + nSelectedInter == length(c(setParam$dgp$linEffects, setParam$dgp$linEffects))
    # selectedAllInd <- nSelectedInd + nSelectedIndInter == length(c(setParam$dgp$indEffects, setParam$dgp$indInterEffects))
    # 
    # # only simulated effects selected (i.e., every other predictor is not selected!)
    # nSelectedOthers <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
    #   sum(!(rownames(selectedVars)[selectedVars[,iCol]] %in% c(setParam$dgp$interEffects, setParam$dgp$linEffects)))})
    # 
    # nSelectedOthersInd <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
    #   sum(!(rownames(selectedVars)[selectedVars[,iCol]] %in% c(setParam$dgp$indInterEffects, setParam$dgp$indEffects)))})
    # 
    # 
    # # variable selection stats per sample
    # # for all: 1 = TRUE, 0 = FALSE
    # selectionPerSample <- cbind(nLin = nSelectedLin, nInter = nSelectedInter, 
    #                             nInd = nSelectedInd, nIndInter = nSelectedIndInter,
    #                             all.T1F0 = selectedAll, all.Ind = selectedAllInd,
    #                             nOthers = nSelectedOthers, nOthersInd = nSelectedOthersInd,
    #                             alpha = tunedAlphaVec, lambda = tunedLambdaVec,
    #                             idxCondLabel = iCond)
    
    
    
    # join results 
    list(# estBeta = estBetaStats, # M, SD, SE for coefficients across samples
         # estBetaFull = estBeta, # coefficients for every predictor in every sample 
         # varSelection = varSelection, 
         # selectSample = selectionPerSample, 
         perfromPerSample = perfromPerSample, 
         performTrainStats = performTrainStats, 
         performTestStats = performTestStats,
         pvi = pvi, 
         interStrength = interStrength,
         selectionPerSample = selectionPerSample)
  }) # end of loop across conditions # nested lapplys take ~ 15 mins
  
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
  resFileName <- paste0(resFolder, "/", "resultsN", condGrid[iSim, "N"], 
                        "_pTrash", condGrid[iSim, "pTrash"], 
                        "_rel", condGrid[iSim,"reliability"],
                        "_f", ifelse(condGrid[iSim,"factors"], 1, 0), ".rds")
  saveRDS(iCondRes, file = resFileName)
})

# close cluster to return resources (memory) back to OS
stopCluster(cl)

## on my machine
# used to take ~ 3.5 hours (without N = 10.000) (old version) - without parallelisation

## on dalek server
# 10 cores parallel, on server: 2.5 hours (only N = {100, 300, 1000} and pTrash = {10, 50, 100})
#     + no GBM estimation
# really seems twice as fast on 10 compared to 5 cores 
# in N1000xpTrash100 max. 14 cores?! or less ~ 10?
# more cores in early conditions possible (memory is not exhausted)



