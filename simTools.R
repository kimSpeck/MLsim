# ToDo: implement function to vary correlations and include correlation function to process
createFolder <- function(folderPath) {
  if (!file.exists(folderPath)){
    dir.create(folderPath)
  }
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
  # see plorResample function of the caret package
  # https://github.com/topepo/caret/blob/master/pkg/caret/R/postResample.R
  ## input:
  # pred    - [vector] outcome as predicted based on model 
  # obs     - [vector] observed outcome 
  ## output:
  #         - [vector] c(RMSE, Rsquared, MAE)
  ###
  
  resamplCor <- try(cor(pred, obs, use = "pairwise.complete.obs"), silent = TRUE)
  mse <- mean((pred - obs)^2)
  mae <- mean(abs(pred - obs))
  
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
  structureModel <- paste0(sapply(seq_len(P-1), function(iP) {
    paste0("Var", iP, " + ")}, simplify = "array"), collapse = "")
  # structureModel <- paste0("Fy ~ ", structureModel, paste0("Var", P)) # all latent
  structureModel <- paste0("y ~ ", structureModel, paste0("Var", P)) # y as OV in structure model
  
  # SImodel <- paste0(lV, "\n", varOV, "\n", varLV, "\n", yLV, "\n", structureModel)
  SImodel <- paste0(lV, "\n", varOV, "\n", varLV, "\n", structureModel)
  return(SImodel)
}

