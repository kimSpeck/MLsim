# Model fitting

library(glmnet)

# load parameters & custom functions 
source("setParameters.R") # parameter values
source("simTools.R")
 
# N <- setParam$dgp$N[1]
# pTrash <- setParam$dgp$pTrash[1]
dataFolder <- "data"
resFolder <- "results"
if (!file.exists(resFolder)){
  dir.create(resFolder)
}
# test it!
# iSim = 3

condGrid <- expand.grid(N = setParam$dgp$N, 
                        pTrash = setParam$dgp$pTrash)

results <- lapply(seq_len(nrow(condGrid)), function(iSim) {
  
  fileName <- paste0("simDataN", condGrid[iSim, "N"], "_pTrash", condGrid[iSim, "pTrash"], ".rda")
  load(paste0(dataFolder, "/", fileName))
  
  # fitte eine regularisierte Regression
  tmp_estRes <- lapply(seq_along(setParam$dgp$condLabels), function(iCond) {
    
    estRes <- lapply(seq_len(setParam$dgp$nTrain), function(iSample) {
      
      # # test it
      # iCond <- 1
      # iSample <- 1
      
      X <- as.matrix(data[[iSample]][["X_int"]])
      y <- data[[iSample]][["yMat"]][,iCond]
      
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
      foldid <- sample(1:setParam$fit$nfolds, size = length(y), replace = TRUE)
      
      fit_cv <- lapply(setParam$fit$alpha, function(iAlpha) {
        cv.glmnet(x = X, 
                  y = y, 
                  foldid = foldid, # as suggested for alpha tuning via cv
                  alpha = iAlpha, 
                  lambda = setParam$fit$lambda, # tune lambda within cv.glmnet
                  family = "gaussian", 
                  standardize = TRUE, 
                  nfolds = setParam$fit$nfolds, # 10 fold cross validation  
                  type.measure = "mse") 
      })
    
      # choose alpha & lambda based on cross validation tuning of lambda given 
      #     specific alpha values results
      #   -> mse from fit_cvs to choose alpha since lambda in cv is chosen based on mse as well
      
      tuneParam <- lapply(seq_along(fit_cv), function(iAlpha) {
        idxLambda1se <- fit_cv[[iAlpha]]$index[2] # optimal lambda index for lambda.1se
        tmp_mse <- fit_cv[[iAlpha]]$cvm[idxLambda1se] # MSE vector
        tmp_l1se <- fit_cv[[iAlpha]]$lambda.1se # rather conservative
        c(alpha = setParam$fit$alpha[iAlpha],
             MSE = tmp_mse, 
             lambda1SE = tmp_l1se)
      })
      
      tuneParam <- do.call(rbind, tuneParam)
      idxOptim <- which(tuneParam[,"MSE"] == min(tuneParam[,"MSE"]))
      
      # Best parameters
      tunedAlpha <- tuneParam[idxOptim, "alpha"]
      tunedLambda <- tuneParam[idxOptim, "lambda1SE"]
      
      # "lambda.min" = the lambda at which the smallest MSE is achieved.
      # "lambda.1se" = the largest λ at which the MSE is within one SE of the smallest MSE (default).
      # here: "one-standard-error" rule for choosing lambda (Hastie et al. 2009)
      #   Friedman et al. 2010. Regularization Paths for Generalized Linear Models via Coordinate Descent.
      fit <- glmnet(x = as.matrix(data[[iSample]][["X_int"]]),
                    y = data[[iSample]][["yMat"]][,iCond], 
                    alpha = tunedAlpha, 
                    lambda = tunedLambda,
                    family = "gaussian", 
                    standardize = TRUE)
      
      estBeta <- as.matrix(fit$beta)
      
      # evaluate performance (R2, RMSE, MAE) for new sample from same dgp (size is 50% of training sample)
      
      # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model training
      predTrain <- predict(fit, X)  
      performTrain <- evalPerformance(predTrain, y)
        
      # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model testing
      performTest <- lapply(paste0("test", seq_len(setParam$dgp$nTest)), function(iTest) {
        Xtest <- as.matrix(data[[iTest]][["X_int"]])
        obs <- data[[iTest]][["yMat"]][,iCond]
        pred <- predict(fit, Xtest)  
        
        evalPerformance(pred, obs)
      })
      performTest <- do.call(rbind, performTest)
      
      list(estB = estBeta, 
           performTrain = performTrain,
           performTest = performTest)
    })
    
    # coefficients
    estBeta <- do.call(cbind, lapply(estRes, function(X) X[["estB"]]))
    estBetaStats <- getStats(estBeta, 1, setParam$dgp$nSamples)
    
    # training performance (RMSE, Rsquared, MAE)
    performTrainMat <- do.call(rbind, lapply(estRes, function(X) X[["performTrain"]])) # each sample
    performTrainStats <- getStats(performTrainMat, 2, setParam$dgp$nSamples) # M, SD, SE
    
    # test performance (RMSE, Rsquared, MAE)
    performTestMat <- do.call(rbind, lapply(estRes, function(X) X[["performTest"]])) # each sample
    performTestStats <- getStats(performTestMat, 2, setParam$dgp$nSamples) # M, SD, SE
    
    # here!
    # list(estB_M = estB_M, estB_SD = estB_SD, estB_SE = estB_SE)
  })
  
  names(tmp_estRes) <- setParam$dgp$condLabels
  tmp_estRes <- do.call(Map, c(f = cbind, tmp_estRes)) # das?
  rm(data)
  gc()
  tmp_estB
})
# takes ~ 3.5 hours (without N = 10.000)

fileName = "estB_initialResults.rda"
save(res_estB, file = paste0(resFolder, "/", fileName))


