setParam <- list()

# data generating process
setParam$dgp$nTrain <- 100
setParam$dgp$nTest <- 1
setParam$dgp$nSamples <- setParam$dgp$nTrain + setParam$dgp$nTest

# setParam$dgp$N <- c(100, 300, 1000, 10000) # number of observations
setParam$dgp$N <- c(100, 300, 1000) # number of observations

# !if 50% of training sample size is used as test sample, test sample sizes vary across N observation conditions
#   as a result, the empirical SE of Rsquared changes according to N 
#   -> fix test sample size 
setParam$dgp$testNpc <- 0.5 
setParam$dgp$p <- 4
setParam$dgp$pTrash <- c(10, 50, 100) # number of "trash" predictors

setParam$dgp$interDepth <- c(2) # depth of interactions (so far: only two-way interaction)
setParam$dgp$poly <- c(2) # degree of polynomials (so far: only quadratic effects)

# predictors and their polynomials + all interactions of depth
P <- (setParam$dgp$p + setParam$dgp$pTrash) # all predictors
setParam$dgp$nModelPredictors <- P * setParam$dgp$poly + choose(P, setParam$dgp$interDepth)
# P * setParam$dgp$poly + (P * (P-1) / 2) # this only works for interactionDepth = 2
rm(P)

################################################################################
# simulated effects
################################################################################
# linear effects
setParam$dgp$linEffects <- sapply(seq_len(setParam$dgp$p), function(x) paste0("Var", x))
# choose variables for interaction that have no linear effects (R2 budget)
# interEffects <- c("Var1:Var2", "Var1:Var4", "Var2:Var3", "Var3:Var4")
setParam$dgp$interEffects <- c("Var5:Var6", "Var5:Var8", "Var6:Var7", "Var7:Var8")

# proportion of effect explained by linear effects vs. interaction
setParam$dgp$percentLinear <- c(0.5, 0.8, 0.2) 
setParam$dgp$percentInter <- c(0.5, 0.2, 0.8)
setParam$dgp$percentPoly <- c(0, 0, 0)

# check
if (!all(setParam$dgp$percentLinear+
         setParam$dgp$percentInter +
         setParam$dgp$percentPoly == 1)) stop("Proportion of different kinds of effects do not sum up to 1! Check values!")

################################################################################
# R squared
################################################################################
setParam$dgp$Rsquared <- c(.10, .30, .50, .80)

setParam$dgp$trueEffects <- cbind(p0.5 = c(0.116, 0.227, 0.346, 0.693),
                                  p0.8 = c(0.144, 0.28, 0.43, 0.85),
                                  p0.2 = c(0.079, 0.15, 0.235, 0.45))  
rownames(setParam$dgp$trueEffects) <- setParam$dgp$Rsquared

comboGrid <- expand.grid(setParam$dgp$Rsquared, 
                         paste(setParam$dgp$percentLinear, setParam$dgp$percentInter, sep = "_"))
setParam$dgp$condLabels <- sapply(seq_len(length(setParam$dgp$Rsquared) * length(setParam$dgp$percentLinear)), 
                                  function(x) paste0("R2", comboGrid$Var1[x], "lin_inter", comboGrid$Var2[x]))
rm(comboGrid)

################################################################################
# predictor correlations
# ToDo: manipulate correlations between different variables (no separate simulated conditions)
################################################################################
# do not randomly sample correlation matrix in each sample!
# instead: randomly choose one correlation matrix and use the same correlation matrix in each sample
setParam$dgp$Reffects <- 0.4 # correlation between predictors with simulated effects

# average correlations of trash-predictors and their SD 
setParam$dgp$meanR <- 0
setParam$dgp$sdR <- 0.05 # 0.1 leads to non-PSM correlation matrix
# setParam$dgp$sdR <- 0.2 # original value

# use this correlation matrix in data simulation! 
# set seed
set.seed(42)

# randomly choose correlation matrix for biggest set of predictors (max(pTrash))
# (use subsets for pTrash < max(pTrash))
P <- max(setParam$dgp$p + setParam$dgp$pTrash)
corN <- P*(P-1)/2
corNp <- setParam$dgp$p*(setParam$dgp$p-1)/2

repeat{
  # correlations around 0 for trash variables (# alternatively exactly 0 correlations)
  corVec <- rtruncnorm(corN, mean = setParam$dgp$meanR, sd = setParam$dgp$sdR, 
                       a = -0.5, b = 0.5) 
  rX <- lavaan::lav_matrix_upper2full(corVec, diagonal = F) # fill up upper triangle
  # rX <- lavaan::lav_matrix_upper2full(rep(0, corN), diagonal = F) # fill up upper triangle
  
  # correlations between predictor variables (all identical atm)
  rP <- lavaan::lav_matrix_upper2full(rep(setParam$dgp$Reffects, corNp), diagonal = F)
  rX[1:dim(rP)[1], 1:dim(rP)[1]] <- rP
  diag(rX) <- 1
  
  # check if correlation matrix is positive semi-definite
  # thereby avoid warning and correction procedure in rmvnorm
  if(all(eigen(rX)$values >= 0)) {
    break
  }
}
setParam$dgp$predictorCorMat <- rX # save PSM correlation matrix
rm(P, corN, corNp, corVec, rP, rX) # remove temporary variables

# error "variance" (= standard deviation)
setParam$dgp$sigmaE <- 1

# measurement error in predictors
setParam$dgp$reliability <- c(0.5, 0.7, 1)

# parameter for gbm or model fitting more general
setParam$fit$lambda <- 10^seq(-1, 1, length = 100)
setParam$fit$alpha <- seq(0, 1, length.out = 20) # alpha with at least 20 steps
setParam$fit$nfolds <- 10
setParam$fit$lambdaCrit <- "1se" # min or 1se
if (!(setParam$fit$lambdaCrit %in% c("min", "1se"))) {
  stop("value for 'setParam$fit$lambdaCrit' not available! Choose 'min' or '1se'.")
}