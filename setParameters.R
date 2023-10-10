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

# 
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

# r squared
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
# average correlations of predictors and their SD 
# ToDo: manipulate correlations between different variables (no separate simulated conditions)
setParam$dgp$meanR <- 0
setParam$dgp$sdR <- 0.1
# setParam$dgp$sdR <- 0.2 # original value

# error "variance" (= standard deviation)
setParam$dgp$sigmaE <- 1

# parameter for gbm or model fitting more general
setParam$fit$lambda <- 10^seq(-1, 1, length = 100)
setParam$fit$alpha <- seq(0, 1, .1) 
setParam$fit$nfolds <- 10
setParam$fit$lambdaCrit <- "1se" # min or 1se
if (!(setParam$fit$lambdaCrit %in% c("min", "1se"))) {
  stop("value for 'setParam$fit$lambdaCrit' not available! Choose 'min' or '1se'.")
}