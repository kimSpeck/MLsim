createFolder <- function(folderPath) {
  if (!file.exists(folderPath)){
    dir.create(folderPath)
  }
}

# generate matrix of coefficients to simulate data depending on predictor matrix
genBmat <- function(X_int, dgp, setParam) {
  ###
  # create matrix with beta coefficients for each simulation condition R2 x lin_inter
  ## input:
  # X_int     - [matrix] predictor value matrix
  # dgp       - [char] indicating current dgp {inter, nonlinear} 
  # setParam  - [vector] list with parameter values
  ## output:
  # bMatrix   - [matrix] matrix with beta coefficients for each simulation condition R2 x lin_inter
  #                       with predictors in rows and simulated conditions in columns
  ###
  
  if (!(dgp %in% c("inter", "nonlinear", "pwlinear", "nonlinear3"))) {
    stop("We can only simulate inter, nonlinear or piecewise linear data!")
  }
  
  # generate empty matrix of beta coefficients with ...
  #   ... number of rows = total number of terms (linear effects, poly, interactions)
  #   ... number of columns = R2 {0.1, 0.3, 0.5, 0.8} x effect size splitting {0.5_0.5, 0.8_0.2, 0.2_0.8}
  bMatrix <- matrix(0, 
                    ncol = length(setParam$dgp$Rsquared) * length(setParam$dgp$percentLinear),
                    nrow = ncol(X_int))
  
  rownames(bMatrix) <- colnames(X_int) # Var1, ..., poly(...), Var1:Var2 
  colnames(bMatrix) <- setParam$dgp$condLabels # e.g., R20.1lin_inter0.5_0.5
  
  # match structure of bMatrix with arrangement of values in matrix of true effects 
  #   (i.e., population beta coefficients)
  #   R2 indicated by rows in matrices of population beta coefficients in setParam
  bMatrixR2 <- stringr::str_sub(colnames(bMatrix), start = 3L, end = 5L)
  
  # # organization of R² in rows is identical for lin vs. inter and for lin vs. nonlin 
  # idxRowInter <- match(bMatrixR2, rownames(setParam$dgp$trueB$inter$inter))
  # idxRowInter <- match(bMatrixR2, rownames(setParam$dgp$trueB$nonlinear$nonlinear))
  # identical(idxRow, idxRowInter)
  
  # how much R2 gets the linear effect: get lin from lin_inter in R20.2lin_inter0.5_0.5
  bMatrixLin <- stringr::str_sub(colnames(bMatrix), start = -7L, end = -5L)
  # how much R2 gets the interaction or nonlinear effect: get inter from lin_inter in R20.2lin_inter0.5_0.5
  bMatrixOther <- stringr::str_sub(colnames(bMatrix), start = -3L, end = -1L)
  
  if (dgp == "inter") {
    idxRow <- match(bMatrixR2, rownames(setParam$dgp$trueB$inter$lin)) 
    idxColLin <- match(bMatrixLin, colnames(setParam$dgp$trueB$inter$lin))
    idxColOther <- match(bMatrixOther, colnames(setParam$dgp$trueB$inter$inter))
    
    # fill up bMatrix with different beta coefficients for linear and interaction
    #   effects in each combination of ...
    #   R2 {0.1, 0.3, 0.5, 0.8} and 
    #   effect size splitting between linear & interaction effects {0.5_0.5, 0.8_0.2, 0.2_0.8}
    for (iCond in seq_len(dim(bMatrix)[2])) {
      # linear effects
      bMatrix[c(setParam$dgp$linEffects), iCond] <- rep(
        setParam$dgp$trueB$inter$lin[idxRow[iCond],idxColLin[iCond]], 
        times = length(setParam$dgp$linEffects)) 
      # interaction effects
      bMatrix[c(setParam$dgp$interEffects), iCond] <- rep(
        setParam$dgp$trueB$inter$inter[idxRow[iCond],idxColOther[iCond]], 
        times = length(setParam$dgp$interEffects)) # R2 = 0.02
    }
  } else if (dgp == "nonlinear"){
    idxRow <- match(bMatrixR2, rownames(setParam$dgp$trueB$nonlinear$lin))
    idxColLin <- match(bMatrixLin, colnames(setParam$dgp$trueB$nonlinear$lin))
    idxColOther <- match(bMatrixOther, colnames(setParam$dgp$trueB$nonlinear$nonlinear))
    
    # fill up bMatrix with different beta coefficients for linear and interaction
    #   effects in each combination of ...
    #   R2 {0.1, 0.3, 0.5, 0.8} and 
    #   effect size splitting between linear & interaction effects {0.5_0.5, 0.8_0.2, 0.2_0.8}
    for (iCond in seq_len(dim(bMatrix)[2])) {
      # linear effects
      bMatrix[c(setParam$dgp$linEffects), iCond] <- rep(
        setParam$dgp$trueB$nonlinear$lin[idxRow[iCond],idxColLin[iCond]], 
        times = length(setParam$dgp$linEffects)) 
      # interaction effects
      bMatrix[c(setParam$dgp$nonlinEffects), iCond] <- rep(
        setParam$dgp$trueB$nonlinear$nonlinear[idxRow[iCond],idxColOther[iCond]], 
        times = length(setParam$dgp$nonlinEffects)) # R2 = 0.02
    }
  } else if (dgp == "nonlinear3"){
    idxRow <- match(bMatrixR2, rownames(setParam$dgp$trueB$nonlinear3$lin))
    idxColLin <- match(bMatrixLin, colnames(setParam$dgp$trueB$nonlinear3$lin))
    idxColOther <- match(bMatrixOther, colnames(setParam$dgp$trueB$nonlinear3$nonlinear))
    
    # fill up bMatrix with different beta coefficients for linear and interaction
    #   effects in each combination of ...
    #   R2 {0.1, 0.3, 0.5, 0.8} and 
    #   effect size splitting between linear & interaction effects {0.5_0.5, 0.8_0.2, 0.2_0.8}
    for (iCond in seq_len(dim(bMatrix)[2])) {
      # linear effects
      bMatrix[c(setParam$dgp$linEffects), iCond] <- rep(
        setParam$dgp$trueB$nonlinear3$lin[idxRow[iCond],idxColLin[iCond]], 
        times = length(setParam$dgp$linEffects)) 
      # interaction effects
      bMatrix[c(setParam$dgp$nonlinEffects3), iCond] <- rep(
        setParam$dgp$trueB$nonlinear3$nonlinear[idxRow[iCond],idxColOther[iCond]], 
        times = length(setParam$dgp$nonlinEffects3)) # R2 = 0.02
    }
  } else if (dgp == "pwlinear"){
    idxRow <- match(bMatrixR2, rownames(setParam$dgp$trueB$pwlinear$lin))
    idxColLin <- match(bMatrixLin, colnames(setParam$dgp$trueB$pwlinear$lin))
    idxColOther <- match(bMatrixOther, colnames(setParam$dgp$trueB$pwlinear$nonlinear))
    
    # fill up bMatrix with different beta coefficients for linear and interaction
    #   effects in each combination of ...
    #   R2 {0.1, 0.3, 0.5, 0.8} and 
    #   effect size splitting between linear & interaction effects {0.5_0.5, 0.8_0.2, 0.2_0.8}
    for (iCond in seq_len(dim(bMatrix)[2])) {
      # linear effects
      bMatrix[c(setParam$dgp$linEffects), iCond] <- rep(
        setParam$dgp$trueB$pwlinear$lin[idxRow[iCond],idxColLin[iCond]], 
        times = length(setParam$dgp$linEffects)) 
      # interaction effects
      bMatrix[c(setParam$dgp$pwlinEffects), iCond] <- rep(
        setParam$dgp$trueB$pwlinear$nonlinear[idxRow[iCond],idxColOther[iCond]], 
        times = length(setParam$dgp$pwlinEffects)) # R2 = 0.02
    }
  }
  
  return(bMatrix)
}

createPredictors <- function(N, P, mP = rep(0, P), corMat) {
  ###
  # create matrix with predictor values
  ## input:
  # N       - [scalar] number of observations
  # P       - [scalar] number of predictors
  # M       - [vector] mean values for predictors (default: centered predictors)
  # corMat  - [matrix] correlation between predictors
  ## output:
  # X       - [matrix] matrix of predictor values [N x P] 
  ###
  
  # check if mean vector has enough values
  if (P != length(mP)) stop("mismatch in number of predictors and provided mean values")
  if (P != dim(corMat)[1] && P != dim(corMat)[2]) stop("mismatch in number of predictors and provided correlation matrix")
  
  # sample predictors from multivariate normal distribution
  X <- rmvnorm(n = N, mean = mP, sigma = corMat)
  return(X)
}

# 
genModel <- function(varVec, interDepth, polyDegree) {
  ###
  # create model formula
  ## input:
  # varVec      - [vector] vector with all predictor names
  # interDepth  - [scalar] depth of interaction (applies to all predictors)
  # polyDegress - [scalar] which degree should polynomial take? (applies to all predictors)
  ## output:
  # model       - [char] formula as character string
  ###
  # basic model (only predictors)
  model <- paste0("~ (0 +", paste(varVec, sep = "", collapse = "+"), ")") # all variables 
  
  # add interactions to model
  if (interDepth > 0) model <- paste0(model, "^", interDepth) # interactions
  
  # add polynomial terms to model
  # if polyDegree == 1, add polynomial of degree 1 although its the same as predictors 
  #   to remove it later (as it is done for higher polynomials as well)
  if (polyDegree > 0) model <- paste0(model, " +", 
                                      paste("poly(", varVec, ", degree = ", polyDegree, ", raw = 2)", 
                                            sep = "", collapse = "+")) # polynomials)
  return(model)
}

#   1. compute median/quantile in general, 2. create dummy based on median/quantile thresh
createDummy <- function(var, q = 0.5, effectCoding = T) {
  #####
  # takes a numeric vector, computes a specified quantile, and returns a dummy (0/1) 
  #   variable indicating which observations exceed that cutoff.
  ## Input:
  # var - numeric vector.
  # q   - numeric value in [0, 1] specifying the quantile cutoff (default 0.5)
  ## Output
  #     - numeric vector of 0s and 1s (same length as var)
  #####
  # Safety checks for quantile
  if (q > 1 | q < 0) {
    stop("Please set quantile to a value between 0 and 1.")
  } 
  
  # compute cutoff
  thresh <- quantile(var, probs = q)
  # create dummy: -1 for values <= median, 1 for values > median
  if (effectCoding) {
    ifelse(var > thresh, 1, -1)
  } else { 
    ifelse(var > thresh, 1, 0)
  }
}

getR2 <- function(X, b, sigmaE) {
  ###
  # calculate proportion of explained variance relative to total variance (R squared)
  ## input:
  # X       - [matrix] matrix with all predictor values [predictors, interactions and polynomials]
  # b       - [vector] regression coefficients [for all predictor, interactions and polynomials]
  # sigmaE  - [scalar] error standard deviation
  ## output:
  # R2      - [scalar] R squared
  ###
  R2 <- var(X %*% b) / (var(X %*% b) + sigmaE^2)
  return(R2)
}

calcDV <- function(X, b, sigmaE, N) {
  ###
  # calculate dependent variable from predictors and regression coefficients and add noise
  #   noise: we assume normally distributed noise with an expected value of 0
  ## input:
  # X       - [matrix] matrix with all predictor values [predictors, interactions and polynomials]
  # b       - [vector] regression coefficients [for all predictor, interactions and polynomials]
  # sigmaE  - [scalar] error standard deviation
  # N       - [scalar] number of observations to sample error
  ## output:
  # y       - [vector] 
  ###
  
  y <- X %*% b + rnorm(n = N, mean = 0, sd = sigmaE)  
  return(y)
}

evalPerformance <- function(pred, obs) {
  ###
  # evaluate performance (R2, RMSE, MAE) for train or test data
  # see postResample function of the caret package
  # https://github.com/topepo/caret/blob/5f4bd2069bf486ae92240979f9d65b5c138ca8d4/pkg/caret/R/postResample.R#L126C3-L144C49
  # adjusted caret package function to return 0 instead of NA if there are no 
  #     predictors chosen in elastic net or lasso regression
  ## input:
  # pred    - [vector] outcome as predicted based on model 
  # obs     - [vector] observed outcome 
  ## output:
  #         - [vector] c(RMSE, Rsquared, MAE)
  ###
  
  # # this is the original caret code; it produces NA if no predictor is chosen
  # #     with (silent) warning message standard deviation is zero
  # #   -> NAs bias the estimation of Rsquared in simulated conditions with low R2
  # resamplCor <- try(cor(pred, obs, use = "pairwise.complete.obs"), silent = TRUE)
  
  if (!is.factor(obs) && is.numeric(obs)) {
    # empirical R2
    # use tryCatch to set correlation and therefore R2 to zero if no predictor is 
    #     chosen in elastic net or lasso regression
    resamplCor <- tryCatch(cor(pred, obs, use = "pairwise.complete.obs"),
                           warning = function(w) 0, silent = TRUE) 
    mse <- mean((pred - obs)^2)
    mae <- mean(abs(pred - obs))
  } else {
    resamplCor <- NA
    mse <- NA
    mae <- NA
  }
  
  return(c(RMSE = sqrt(mse), Rsquared = resamplCor^2, MAE = mae))
}

getStats <- function(data, aggrDim, nSamples, na.rm = T) {
  ###
  # get mean, sd and standard error of the mean
  ## input:
  # data      - [matrix] data from different samples 
  # aggrDim   - [scalar] dimension to aggregate over (1 = rows, 2 = columns)
  # nSamples  - [scalar] number of samples that contribute to the mean
  # na.rm     - [boolean] remove NAs to avoid NAs as outcome
  ## output:
  #         - [vector] c(M, SD, SE)
  ###
  
  M <- apply(data, MARGIN = aggrDim, mean, na.rm = na.rm)
  SD <- apply(data, MARGIN = aggrDim, sd, na.rm = na.rm)
  SE <- SD / sqrt(nSamples) # standard error of the mean
  return(cbind(M = M, SD = SD, SE = SE))
}

rmDuplicatePoly <- function(X) {
  ###
  # remove first degree polynomials from data (they are duplicates!)
  # for testing purposes:
  # colnames(X_int)[stringr::str_detect(colnames(X_int), pattern = "^(poly\\().+(\\)1)$")]
  ## input:
  # X       - [matrix] predictor variables from model.matrix
  ## output:
  # X       - [matrix] predictor variables without "duplicate" due to first polynomial
  ###
  
  rmColsIdx <- which(stringr::str_detect(colnames(X), pattern = "^(poly\\().+(\\)1)$"))
  X <- X[,-rmColsIdx] 
  return(X)
}

# check simulated reliabilities with single-indicator SEM
# check correlations between predictors (i.e., latent variables) and compare simulated 
#     regression coefficients with path-coefficients of the structure model
# see checkDataSimulation.R
genSingleIndicatorModel <- function(P, reliability) { 
  # define latent variables with single indicator each
  lV <- paste0(sapply(seq_len(P), function(iP) {
    paste0("F", iP, " =~ 1 * Var", iP, "\n")}, simplify = "array"), collapse = " ")
  
  # fix error variance in the observed variable according to reliability 
  varError <- (1-reliability)/reliability 
  varOV <- paste0(sapply(seq_len(P), function(iP) {
    paste0("Var", iP, " ~~ ", varError," * Var", iP, "\n")}, simplify = "array"), collapse = " ")
  
  # fix variance in the latent variable to 1 (as modeled in the predictors)
  # apply(X, 2, var)
  varLV <- paste0(sapply(seq_len(P), function(iP) {
    paste0("F", iP, " ~~ ", 1 ," * F", iP, "\n")}, simplify = "array"), collapse = " ")
  
  # # model dependent variable (all latent formulation)
  # # this version is equivalent to simply using y as an observed variable in the structure model
  # yLV <- "Fy =~ 1 * y\n y ~~ 0*y"
  
  # structure model (regression part to recover true simulated regression coefficients
  # use latent factors to get simulated regression coefficients (without measurement error
  #     and hence without attenuated estimated regression coefficients)
  # if manifest variables in structure model we can observe the attenuation of the
  #     regression coefficients due to measurement error as: estBeta = trueBeta * sqrt(reliability)
  structureModel <- paste0(sapply(seq_len(P-1), function(iP) {
    paste0("F", iP, " + ")}, simplify = "array"), collapse = "")
  # structureModel <- paste0("Fy ~ ", structureModel, paste0("Var", P)) # all latent
  structureModel <- paste0("y ~ ", structureModel, paste0("F", P)) # y as OV in structure model
  
  # SImodel <- paste0(lV, "\n", varOV, "\n", varLV, "\n", yLV, "\n", structureModel)
  SImodel <- paste0(lV, "\n", varOV, "\n", varLV, "\n", structureModel)
  return(SImodel)
}


# evaluate ENET (for xgboost fit models)
predEnet <- function(model, newdata)  {
  results <- predict(model, as.matrix(newdata))
  return(results)
}
