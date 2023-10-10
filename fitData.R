# Model fitting
# to do: 
#     -> Rsquared often NA; this happens if 0 variables are chosen in Enet!  
library(glmnet)
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
nCoresSampling <- 2

# test it!
# iSim = 9

condGrid <- expand.grid(N = setParam$dgp$N, 
                        pTrash = setParam$dgp$pTrash)

# remove lines in condGrid for which there already are results
# check which files are already in the results folder
resFileList <- list.files(resFolder)
resFileList <- sub('results', "", resFileList); resFileList <- sub('.rds', "", resFileList)

allDataFiles <- paste0("simDataN", condGrid[, "N"], "_pTrash", condGrid[, "pTrash"], ".rda")
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
  
  fileName <- paste0("simDataN", condGrid[iSim, "N"], "_pTrash", condGrid[iSim, "pTrash"], ".rda")
  load(paste0(dataFolder, "/", fileName))
  
  # fitte eine regularisierte Regression
  tmp_estRes <- lapply(seq_along(setParam$dgp$condLabels), function(iCond) { # R2 x lin_inter combinations
    
    # estRes <- lapply(seq_len(setParam$dgp$nTrain), function(iSample) {
    estRes <- parLapply(cl, seq_len(setParam$dgp$nTrain), function(iSample) {
      
      if (iSample == 1) {
        print(paste0(fileName,
                     "; iCond: ", setParam$dgp$condLabels[iCond], 
                     "; first Sample at ", format(Sys.time(), "%d%m%y_%H%M%S")))  
      }
      
      # # test it
      # iCond <- 1
      # iSample <- 1
      
      Xtrain <- as.matrix(data[[iSample]][["X_int"]])
      ytrain <- data[[iSample]][["yMat"]][,iCond]
      
      # fit data on full sample ("training")
      #     ... alpha (elastic net) & lambda-grid (overall penalty strength) for fitting Enet
      # glmnet package only allows to tune lambda via cross-validation for fixed alpha
      # see source code of glmnet: (https://github.com/cran/glmnet/blob/master/R/cv.glmnet.R) 
      #   Note that \code{cv.glmnet} does NOT search for values for \code{alpha}. 
      #   A specific value should be supplied, else \code{alpha=1} is assumed by default. 
      #   If users would like to cross-validate \code{alpha} as well, they should call 
      #   \code{cv.glmnet} with a pre-computed vector \code{foldid}, and then use this 
      #   same fold vector in separate calls to \code{cv.glmnet} with different values of \code{alpha}.
      # see solution in this source code: https://github.com/hongooi73/glmnetUtils/blob/master/R/cvaGlmnetFormula.r
      # see example code here: https://stats.stackexchange.com/questions/268885/tune-alpha-and-lambda-parameters-of-elastic-nets-in-an-optimal-way
      

      # to do: think about "grid" size (are 100 lambda values too much?)
      #     -> alphaVec ~20 values?
      #     -> choose lambda in caret style? lambda vector based on enet with alpha = 0.5?
      # # provide lambda values by "warm start" for alpha = .5 as done in caret? 
      # #   -> kein warm start wie in caret (weil Kontrolle)
      # # https://stackoverflow.com/questions/48280074/r-how-to-let-glmnet-select-lambda-while-providing-an-alpha-range-in-caret
      # init <- glmnet::glmnet(Matrix::as.matrix(X), y,
      #                        family = "gaussian",
      #                        nlambda = len+2,
      #                        alpha = .5)
      # 
      # lambda <- unique(init$lambda)
      # lambda <- lambda[-c(1, length(lambda))]
      # lambda <- lambda[1:min(length(lambda), len)]
      
      # lambdaVec = 10^seq(-1, 1, length = 100) # tune lambda within cv.glmnet
      # alphaVec = seq(0, 1, .1) # tune alpha by iterating trough alphas
      
      set.seed(89101)
      
      # tune alpha by iterating trough alphas
      foldid <- sample(1:setParam$fit$nfolds, size = length(ytrain), replace = TRUE)
      
      fit_cv <- lapply(setParam$fit$alpha, function(iAlpha) {
        tmp_fit <- cv.glmnet(x = Xtrain, 
                             y = ytrain, 
                             foldid = foldid, # as suggested for alpha tuning via cv
                             alpha = iAlpha, 
                             lambda = setParam$fit$lambda, # tune lambda within cv.glmnet
                             family = "gaussian", 
                             standardize = TRUE, 
                             nfolds = setParam$fit$nfolds, # 10 fold cross validation  
                             type.measure = "mse") 
        
        # cv results that we need for parameter tuning
        idxLambdaCrit <- match(setParam$fit$lambdaCrit, rownames(tmp_fit$index))
        idxLambda <- tmp_fit$index[idxLambdaCrit] # optimal lambda index for lambda criterion (min or se)
        tmp_mse <- tmp_fit$cvm[idxLambda] # choose MSE value for lambda criterion (min or se)
        tmp_lambda <- tmp_fit[[paste0("lambda.", setParam$fit$lambdaCrit)]] 
        
        list(cvm = tmp_mse,
             lambda.1se = tmp_lambda)
      })
    
      # Groesse der Objekte mitschreiben lassen?
      
      # choose alpha & lambda based on cross validation tuning of lambda given 
      #     specific alpha values results
      #   -> mse from fit_cvs to choose alpha since lambda in cv is chosen based on mse as well
      
      tuneParam <- lapply(seq_along(fit_cv), function(iAlpha) {
        tmp_mse <- fit_cv[[iAlpha]]$cvm # MSE vector
        tmp_l1se <- fit_cv[[iAlpha]]$lambda.1se # rather conservative
        c(alpha = setParam$fit$alpha[iAlpha],
             MSE = tmp_mse, 
             lambda1SE = tmp_l1se)
      })
      
      tuneParam <- do.call(rbind, tuneParam)
      
      # ! to do: if there are multiple optima: which optimum should we choose?
      # in the case of a tie the tied variable with lowest index is selected.
      idxOptim <- which(tuneParam[,"MSE"] == min(tuneParam[,"MSE"]))[1]
      
      # Best parameters
      tunedAlpha <- tuneParam[idxOptim, "alpha"]
      tunedLambda <- tuneParam[idxOptim, "lambda1SE"]
      
      # "lambda.min" = the lambda at which the smallest MSE is achieved.
      # "lambda.1se" = the largest λ at which the MSE is within one SE of the smallest MSE (default).
      # here: "one-standard-error" rule for choosing lambda (Hastie et al. 2009)
      #   Friedman et al. 2010. Regularization Paths for Generalized Linear Models via Coordinate Descent.
      fit <- glmnet(x = Xtrain,
                    y = ytrain, 
                    alpha = tunedAlpha,
                    lambda = tunedLambda,
                    family = "gaussian", 
                    standardize = TRUE)
      
      # variable selection strongly depends on alpha ... 
      estBeta <- as.matrix(fit$beta)
      selectedVars <- estBeta != 0
      
      # to do: write this into a function to reuse it for gbm code as well (?)
      # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model training
      predTrain <- predict(fit, Xtrain)
      performTrain <- evalPerformance(predTrain, ytrain)
      performTrain <- matrix(performTrain, ncol = length(performTrain), nrow = 1)
        
      # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model testing
      performTest <- lapply(paste0("test", seq_len(setParam$dgp$nTest)), function(iTest) {
        Xtest <- as.matrix(data[[iTest]][["X_int"]])
        ytest <- data[[iTest]][["yMat"]][,iCond]
        pred <- predict(fit, Xtest)  
        
        evalPerformance(pred, ytest)
      })
      performTest <- do.call(rbind, performTest)
      
      list(estB = estBeta, 
           selectedVars = selectedVars,
           performTrain = performTrain,
           performTest = performTest, 
           tunedAlpha = tunedAlpha, 
           tunedLambda = tunedLambda)
      
      # # remove fitted sample to reduce working memory load
      # data[[iSample]] <- NULL
    })
    
    # coefficients
    estBeta <- do.call(cbind, lapply(estRes, function(X) X[["estB"]]))
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

    # only true predictors in model (8 predictors with simulated effects and every other variable == 0)
    # how many of the linear effects are recovered?
    nSelectedLin <- sapply(seq_len(ncol(selectedVars)), function(iCol) { 
      sum(setParam$dgp$linEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
    # how many of the interaction effects are recovered?
    nSelectedInter <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
      sum(setParam$dgp$interEffects %in% rownames(selectedVars)[selectedVars[,iCol]])})
    
    # all simulated effects selected in model?
    selectedAll <- nSelectedLin + nSelectedInter == length(c(setParam$dgp$linEffects, setParam$dgp$linEffects))
    
    # only simulated effects selected (i.e., every other predictor is not selected!)
    nSelectedOthers <- sapply(seq_len(ncol(selectedVars)), function(iCol) {
      sum(!(rownames(selectedVars)[selectedVars[,iCol]] %in% c(setParam$dgp$interEffects, setParam$dgp$linEffects)))})
    
    # final Enet-Parameters for each sample
    tunedAlphaVec <- unname(do.call(c, lapply(estRes, function(X) X[["tunedAlpha"]])))
    tunedLambdaVec <- unname(do.call(c, lapply(estRes, function(X) X[["tunedLambda"]])))
    
    # variable selection stats per sample
    # for all: 1 = TRUE, 0 = FALSE
    selectionPerSample <- cbind(nLin = nSelectedLin, nInter = nSelectedInter, 
                                all.T1F0 = selectedAll, nOthers = nSelectedOthers,
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
    list(estBeta = estBetaStats, 
         varSelection = varSelection, selectSample = selectionPerSample, 
         perfromPerSample = perfromPerSample, 
         performTrainStats = performTrainStats, performTestStats = performTestStats)
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
  resFileName <- paste0(resFolder, "/", "resultsN", condGrid[iSim, "N"], "_pTrash", condGrid[iSim, "pTrash"], ".rds")
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

# memory leakage due to outer loop (condGrid-loop?) 
#   -> write data to file instead of saving in growing list?

# fileName = "initialFullResults.rda"
# save(results, file = paste0(resFolder, "/", fileName))


