# ToDo: implement function to vary correlations and include correlation function to process

createPredictors <- function(N, P, mP = rep(0, P), mR = 0, sdR) {
  ###
  # create matrix with predictor values
  ## input:
  # N       - [scalar] number of observations
  # P       - [scalar] number of predictors
  # M       - [vector] mean values for predictors (default: centered predictors)
  # mR      - [scalar] mean of predictor correlations
  # sdR     - [scalar] sd predictor correlations
  ## output:
  # X       - [matrix] matrix of predictor values [N x P] 
  ###
  
  # check if mean vector has enough values
  if (P != length(mP)) stop("mismatch in number of predictors and provided mean values")
  
  # randomly sample correlations between predictors 
  # correlations are sampled from truncated normal and are limited to [-0.5, 0.5]
  corN <- P*(P-1)/2
  corVec <- rtruncnorm(corN, mean = mR, sd = sdR, a = -0.5, b = 0.5)
  rX <- lavaan::lav_matrix_upper2full(corVec, diagonal = F) # fill up upper triangle
  diag(rX) <- 1
  
  # sample predictors from multivariate normal distribution
  X <- rmvnorm(n = N, mean = mP, sigma = rX)
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
