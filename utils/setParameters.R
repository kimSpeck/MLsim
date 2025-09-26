setParam <- list()

################################################################################
# data generating process
################################################################################
setParam$dgp$nTrain <- 1000
setParam$dgp$nTest <- 1
setParam$dgp$nSamples <- setParam$dgp$nTrain + setParam$dgp$nTest

# this is only a technical argument which determines if data is saved in ...
#   ... either one big rda file which heavily stresses RAM in parallelisation
#   ... or in rda files for every individual sample (singleSamples = T) 
setParam$dgp$singleSamples <- TRUE

setParam$dgp$N <- c(100, 300, 1000) # number of observations

# !if 50% of training sample size is used as test sample, test sample sizes vary across N observation conditions
#   as a result, the empirical SE of Rsquared changes according to N 
#   -> fix test sample size 
setParam$dgp$testN <- 1000            # fixed test sample size across all N conditions 

setParam$dgp$p <- 4               # number of latent variables
setParam$dgp$interDepth <- c(2) # depth of interactions (so far: only two-way interaction)
setParam$dgp$poly <- c(0) # degree of polynomials (so far: no polynomials)

# setParam$dgp$pNL <- 2             # number of original variables for nonlinear effects   
setParam$dgp$pNL3 <- 3            # number of original variables for nonlinear effects        
setParam$dgp$pPWL <- 3            # number of original variables for piecewise-linear effects
setParam$dgp$pTrash <- c(10, 50)  # number of noise variables

# predictors and their polynomials + all interactions of depth
P <- (setParam$dgp$p + setParam$dgp$pTrash) # all predictors
setParam$dgp$nModelPredictors <- P * setParam$dgp$poly + choose(P, setParam$dgp$interDepth)
# P * setParam$dgp$poly + (P * (P-1) / 2) # this only works for interactionDepth = 2
rm(P)

##### simulated effects #####
# linear effects
setParam$dgp$linEffects <- sapply(seq_len(setParam$dgp$p), function(x) paste0("Var", x))
# choose variables for interaction that have no linear effects (R2 budget)
setParam$dgp$interEffects <- c("Var1:Var2", "Var1:Var4", "Var2:Var3", "Var3:Var4")

# stepwise DGP

# setParam$dgp$nonlinEffects <- sapply((setParam$dgp$p+1):(setParam$dgp$p+setParam$dgp$pNL), 
#                                      function(x) paste0("dumVar", x, ".1"))
# setParam$dgp$nonlinEffects <- c(setParam$dgp$nonlinEffects, "dumVar5.1:dumVar6.1")

setParam$dgp$nonlinEffects3 <- sapply((setParam$dgp$p+1):(setParam$dgp$p+setParam$dgp$pNL3), 
                                     function(x) paste0("dumVar", x, ".1"))

# piecewise linear effect
setParam$dgp$pwlinEffects <- c("Var5.2nd", "Var6.2nd", "Var7.2nd")

# proportion of effect explained by linear effects vs. interaction
setParam$dgp$percentLinear <- c(0.5, 0.8, 0.2) 
setParam$dgp$percentInter <- c(0.5, 0.2, 0.8)
setParam$dgp$percentPoly <- c(0, 0, 0)

# # simulate conditions in which interactions, nonlinear and piecewise linear effects get 100% R² 
# setParam$dgp$percentLinear <- c(0.5, 0.8, 0.2, 0.0) 
# setParam$dgp$percentInter <- c(0.5, 0.2, 0.8, 1.0)
# setParam$dgp$percentPoly <- c(0, 0, 0, 0)

# check
if (!all(setParam$dgp$percentLinear+
         setParam$dgp$percentInter +
         setParam$dgp$percentPoly == 1)) stop("Proportion of different kinds of effects do not sum up to 1! Check values!")

##### R squared #####
setParam$dgp$Rsquared <- c(.20, .50, .80) 

# beta coefficients brute forced based on correlation matrix of predictors
#     via gibbs sampling procedure using optim to derive beta coefficients
### these are the parameters that are used in bruteForceB.R and bruteForceNonLinearB.R
setParam$bruteForceB$pTrash <- 0
setParam$bruteForceB$N <- 100000
setParam$bruteForceB$reliability <- 1
setParam$bruteForceB$poly <- setParam$dgp$poly

# parameters for beta coefficient estimation
setParam$fit$optimLowerLimit <- 0
setParam$fit$optimUpperLimit <- 2
setParam$fit$optimTol <- 1e-5
setParam$fit$optimBetaTol <- 1e-5

# read in coefficients as results from bruteForceB.R and bruteForceNonLinearB.R
bruteForceB_inter <- read.table("utils/bruteForceBcoeff_inter.csv", header = T, sep = ",")
# bruteForceB_nl <- read.table("utils/bruteForceBcoeff_nonlinear.csv", header = T, sep = ",")
bruteForceB_nl3 <- read.table("utils/bruteForceBcoeff_nonlinear3.csv", header = T, sep = ",")
# bruteForceB_nl <- read.table("utils/bruteForceBcoeff_nonlinear_plus.csv", header = T, sep = ",")
bruteForceB_pwl <- read.table("utils/bruteForceBcoeff_piecewise.csv", header = T, sep = ",")

# reorganize coefficients
# # lin vs. interaction
setParam$dgp$trueB$inter$lin <- reshape2::dcast(bruteForceB_inter, R2 ~ lin, value.var = "betaLin")
setParam$dgp$trueB$inter$inter <- reshape2::dcast(bruteForceB_inter, R2 ~ inter, value.var = "betaInter")

row.names(setParam$dgp$trueB$inter$lin) <- setParam$dgp$trueB$inter$lin[["R2"]]
setParam$dgp$trueB$inter$lin[["R2"]]    <- NULL
row.names(setParam$dgp$trueB$inter$inter) <- setParam$dgp$trueB$inter$inter[["R2"]]
setParam$dgp$trueB$inter$inter[["R2"]]    <- NULL

# # # lin vs. nonlinear
# setParam$dgp$trueB$nonlinear$lin <- reshape2::dcast(bruteForceB_nl, R2 ~ lin, value.var = "betaLin")
# setParam$dgp$trueB$nonlinear$nonlinear <- reshape2::dcast(bruteForceB_nl, R2 ~ inter, value.var = "betaInter")
# 
# row.names(setParam$dgp$trueB$nonlinear$lin) <- setParam$dgp$trueB$nonlinear$lin[["R2"]]
# setParam$dgp$trueB$nonlinear$lin[["R2"]]    <- NULL
# row.names(setParam$dgp$trueB$nonlinear$nonlinear) <- setParam$dgp$trueB$nonlinear$nonlinear[["R2"]]
# setParam$dgp$trueB$nonlinear$nonlinear[["R2"]]    <- NULL

# # naming conventions for 0.0 vs. 1.0
# colnames(setParam$dgp$trueB$nonlinear$lin) <- formatC(sort(setParam$dgp$percentLinear), format = "f", digits = 1)
# colnames(setParam$dgp$trueB$nonlinear$nonlinear) <- formatC(sort(setParam$dgp$percentInter), format = "f", digits = 1)

# # lin vs. nonlinear (with 3 dummy variables)
setParam$dgp$trueB$nonlinear3$lin <- reshape2::dcast(bruteForceB_nl3, R2 ~ lin, value.var = "betaLin")
setParam$dgp$trueB$nonlinear3$nonlinear <- reshape2::dcast(bruteForceB_nl3, R2 ~ inter, value.var = "betaInter")

row.names(setParam$dgp$trueB$nonlinear3$lin) <- setParam$dgp$trueB$nonlinear3$lin[["R2"]]
setParam$dgp$trueB$nonlinear3$lin[["R2"]]    <- NULL
row.names(setParam$dgp$trueB$nonlinear3$nonlinear) <- setParam$dgp$trueB$nonlinear3$nonlinear[["R2"]]
setParam$dgp$trueB$nonlinear3$nonlinear[["R2"]]    <- NULL

# # lin vs. piecewise linear
setParam$dgp$trueB$pwlinear$lin <- reshape2::dcast(bruteForceB_pwl, R2 ~ lin, value.var = "betaLin")
setParam$dgp$trueB$pwlinear$nonlinear <- reshape2::dcast(bruteForceB_pwl, R2 ~ inter, value.var = "betaInter")

row.names(setParam$dgp$trueB$pwlinear$lin) <- setParam$dgp$trueB$pwlinear$lin[["R2"]]
setParam$dgp$trueB$pwlinear$lin[["R2"]]    <- NULL
row.names(setParam$dgp$trueB$pwlinear$nonlinear) <- setParam$dgp$trueB$pwlinear$nonlinear[["R2"]]
setParam$dgp$trueB$pwlinear$nonlinear[["R2"]]    <- NULL

rm(bruteForceB_inter, bruteForceB_nl3, bruteForceB_pwl) # rm temporary matrices 

comboGrid <- expand.grid(setParam$dgp$Rsquared, 
                         paste(formatC(setParam$dgp$percentLinear, format = "f", digits = 1), 
                               formatC(setParam$dgp$percentInter, format = "f", digits = 1), sep = "_"))
setParam$dgp$condLabels <- sapply(seq_len(length(setParam$dgp$Rsquared) * length(setParam$dgp$percentLinear)), 
                                  function(x) paste0("R2", comboGrid$Var1[x], "lin_inter", comboGrid$Var2[x]))
rm(comboGrid)

##### predictor correlations #####
# do not randomly sample correlation matrix in each sample!
# instead: randomly choose one correlation matrix and use the same correlation matrix in each sample
setParam$dgp$rLinEffects <- 0.4 # correlation for linear predictors with simulated effects
setParam$dgp$rNonLinEffects <- 0 # correlation for nonlinear predictors with simulated effects

# average correlations of trash-predictors and their SD 
setParam$dgp$meanR <- 0
setParam$dgp$sdR <- 0.05 # 0.1 leads to non-PSM correlation matrix

##### linear vs. interaction effects #####
# use this correlation matrix in data simulation! 
set.seed(42)

# randomly choose correlation matrix for biggest set of predictors (max(pTrash))
# (use subsets for pTrash < max(pTrash))
P <- max(setParam$dgp$p + setParam$dgp$pTrash)
corN <- P*(P-1)/2
corNp <- setParam$dgp$p*(setParam$dgp$p-1)/2

repeat{
  # correlations around 0 for trash variables (# alternatively exactly 0 correlations)
  corVec <- truncnorm::rtruncnorm(corN, mean = setParam$dgp$meanR, sd = setParam$dgp$sdR, 
                       a = -0.5, b = 0.5) 
  rX <- lavaan::lav_matrix_upper2full(corVec, diagonal = F) # fill up upper triangle
  # rX <- lavaan::lav_matrix_upper2full(rep(0, corN), diagonal = F) # fill up upper triangle
  
  # correlations between predictor variables (all identical atm)
  rP <- lavaan::lav_matrix_upper2full(rep(setParam$dgp$rLinEffects, corNp), diagonal = F)
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

##### linear vs. nonlinear effects #####
set.seed(420)

# # randomly choose correlation matrix for biggest set of predictors (max(pTrash))
# # (use subsets for pTrash < max(pTrash))
# P_nl <- max(setParam$dgp$p + setParam$dgp$pNL + setParam$dgp$pTrash)
# corN_nl <- P_nl*(P_nl-1)/2
# corNp <- setParam$dgp$p*(setParam$dgp$p-1)/2
# corNpNL <- setParam$dgp$pNL*(setParam$dgp$pNL-1)/2
# 
# repeat{
#   # correlations around 0 for all variables (but finally for trash variables) 
#   #   alternatively exactly 0 correlations
#   corVec <- truncnorm::rtruncnorm(corN_nl, mean = setParam$dgp$meanR, sd = setParam$dgp$sdR, 
#                                   a = -0.5, b = 0.5) 
#   rX <- lavaan::lav_matrix_upper2full(corVec, diagonal = F) # fill up upper triangle
#   # rX <- lavaan::lav_matrix_upper2full(rep(0, corN), diagonal = F) # fill up upper triangle
#   
#   # correlations between linear predictor variables (all identical atm)
#   #   replace around 0 correlations at the corresponding positions in the matrix
#   rP <- lavaan::lav_matrix_upper2full(rep(setParam$dgp$rLinEffects, corNp), diagonal = F)
#   rX[1:dim(rP)[1], 1:dim(rP)[1]] <- rP
#   
#   # correlations between nonlinear predictor variables (all identical atm)
#   #   replace around 0 correlations at the corresponding positions in the matrix
#   rPnl <- lavaan::lav_matrix_upper2full(rep(setParam$dgp$rNonLinEffects, corNpNL), diagonal = F)
#   rX[(setParam$dgp$p+1):(setParam$dgp$p+dim(rPnl)[1]), 
#      (setParam$dgp$p+1):(setParam$dgp$p+dim(rPnl)[1])] <- rPnl
#   diag(rX) <- 1
#   
#   # check if correlation matrix is positive semi-definite
#   # thereby avoid warning and correction procedure in rmvnorm
#   if(all(eigen(rX)$values >= 0)) {
#     break
#   }
# }
# setParam$dgp$predictorCorMat_nl <- rX # save PSM correlation matrix
# rm(P_nl, corN_nl, corNp, corNpNL, corVec, rP, rPnl, rX) # remove temporary variables

##### linear vs. piecewise linear effects #####
set.seed(42420)

# randomly choose correlation matrix for biggest set of predictors (max(pTrash))
# (use subsets for pTrash < max(pTrash))
P_pwl <- max(setParam$dgp$p + setParam$dgp$pPWL + setParam$dgp$pTrash)
corN_pwl <- P_pwl*(P_pwl-1)/2
corNp <- setParam$dgp$p*(setParam$dgp$p-1)/2
corNpPWL <- setParam$dgp$pPWL*(setParam$dgp$pPWL-1)/2

repeat{
  # correlations around 0 for all variables (but finally for trash variables) 
  #   alternatively exactly 0 correlations
  corVec <- truncnorm::rtruncnorm(corN_pwl, mean = setParam$dgp$meanR, sd = setParam$dgp$sdR, 
                                  a = -0.5, b = 0.5) 
  rX <- lavaan::lav_matrix_upper2full(corVec, diagonal = F) # fill up upper triangle
  # rX <- lavaan::lav_matrix_upper2full(rep(0, corN), diagonal = F) # fill up upper triangle
  
  # correlations between linear predictor variables (all identical atm)
  #   replace around 0 correlations at the corresponding positions in the matrix
  rP <- lavaan::lav_matrix_upper2full(rep(setParam$dgp$rLinEffects, corNp), diagonal = F)
  rX[1:dim(rP)[1], 1:dim(rP)[1]] <- rP
  
  # correlations between nonlinear predictor variables (all identical atm)
  #   replace around 0 correlations at the corresponding positions in the matrix
  rPpwl <- lavaan::lav_matrix_upper2full(rep(setParam$dgp$rNonLinEffects, corNpPWL), diagonal = F)
  rX[(setParam$dgp$p+1):(setParam$dgp$p+dim(rPpwl)[1]), 
     (setParam$dgp$p+1):(setParam$dgp$p+dim(rPpwl)[1])] <- rPpwl
  diag(rX) <- 1
  
  # check if correlation matrix is positive semi-definite
  # thereby avoid warning and correction procedure in rmvnorm
  if(all(eigen(rX)$values >= 0)) {
    break
  }
}
setParam$dgp$predictorCorMat_pwl <- rX # save PSM correlation matrix
rm(P_pwl, corN_pwl, corNp, corNpPWL, corVec, rP, rPpwl, rX) # remove temporary variables

# error "variance" (= standard deviation when sampling via rnorm due to fun. argument)
setParam$dgp$sigmaE <- 1

# measurement error in predictors 
setParam$dgp$reliability <- c(0.6, 0.8, 1)

################################################################################
# model fitting
################################################################################
# always create the full set of conditions including sample seeds and remove conditions 
#   -> ensure reproducibility
# iterate through these combinations of data conditions
# setParam$fit$condGrid <- expand.grid(data = c("inter", "nonlinear", "pwlinear", "nonlinear3"),
setParam$fit$condGrid <- expand.grid(data = c("inter", "pwlinear", "nonlinear3"),
                                     model = c("ENETwo", "ENETw", "GBM", "RF"),
                                     N = setParam$dgp$N, 
                                     pTrash = setParam$dgp$pTrash,
                                     reliability = setParam$dgp$reliability)

# create seed number for parallel cluster (reproducibility of results)
set.seed(7849380)
seedNum <- sample(1:999999, dim(setParam$fit$condGrid)[1], replace = FALSE) 
setParam$fit$condGrid$sampleSeed <- seedNum[1:dim(setParam$fit$condGrid)[1]]


##### Hyperparameter Tuning #####
setParam$fit$nfolds <- 10
setParam$fit$explanation <- FALSE
setParam$fit$saveConds <- TRUE

##### ENET #####
# parameter for gbm or model fitting more general
# fixed set of lambda parameters (only used if warmStart = FALSE)
setParam$fit$lambda <- 10^seq(-1, 1, length = 100) 
# use warmStart?
#  pro: this is how caret operates (i.e., get lambda values for alpha = 0.5)
#         -> thus, lots of users use warm start for lambda tuning
#       with warmStart exact model is chosen more frequently
#  con: less control over lambda grid to be searched      
setParam$fit$warmStart <- TRUE  # thus, lambda parameter are not used
setParam$fit$alpha <- seq(0, 1, length.out = 20) # alpha with at least 20 steps

# the code can handle setParam$fit$lambdaCrit = {"1se", "min"} or every subset
# "lambda.min" = the lambda at which the smallest MSE is achieved.
# "lambda.1se" = the largest lambda at which the MSE is within one SE of the smallest MSE (default).
# here: "one-standard-error" rule for choosing lambda (Hastie et al. 2009)
#   Friedman et al. 2010. Regularization Paths for Generalized Linear Models via Coordinate Descent.
setParam$fit$lambdaCrit <- c("min") # not 1se anymore, but min instead!; but only for fitENET function!
# -> despite the more conservative criterion (lambda.1se instead of lambda.min) all
#     predictors are found often but too much trash is extracted to find the exact model
for (iCrit in setParam$fit$lambdaCrit) {
  if (!(iCrit %in% c("min", "1se"))) {
    stop("value for 'setParam$fit$lambdaCrit' not available! Choose 'min' or '1se'.")
  }
}

##### GBM #####
# hyperparameter tuning grid for gbm 
setParam$fit$tuneGrid_GBM <- expand.grid(
  interaction.depth = c(1,2,3), # tree depth (= max_depth in xgboost)
  n.minobsinnode = c(5, 10, 20),    # end node size (= min_child_weight in xgboost)
  # n.trees = c(seq(20, 500, 30)),      # max number of trees (= nTrees in xgboost)
  n.trees = c(5, seq(10, 40, 10), seq(50, 500, 30)),      # max number of trees (= nTrees in xgboost)
  shrinkage = c(0.001, .011, 0.031,seq(.051, .401, .05))) # shrinkage/learning rate (= eta in xgboost)
  # shrinkage = c(.011, 0.031,seq(.051, .201, .05))) # shrinkage/learning rate (= eta in xgboost)

# number of samples observations to calculate partial dependencies and h-statistic based on pds
setParam$fit$nInterStrength <- 50
setParam$fit$InterStrength <- FALSE
if (setParam$fit$InterStrength & !setParam$fit$explanation) {
  warning("Are you sure you want to calculate interaction strength of the GBM but skip other explanation metrics?")
}
# number of threads in xgb.cv (implicit parallelization)
setParam$fit$nThread <- 1

##### Random Forest #####
# tuning grid for random forests is a function
#   the tuning grid needs to change depending on the numbers of predictors 
#   the number of predictors depends on pTrash in the simulated condition
setParam$fit$setTuningGrid_RF <- function(nPred) {
  expand.grid(
    mtry = c(2, sqrt(nPred), nPred/3, nPred/2), # random predictors @ each node
    # splitrule = c("variance", "extratrees"), # splitting criterion
    splitrule = c("variance"), # splitting criterion
    min.node.size = c(5, 10, 20)) # min observations in end node
}
setParam$fit$numTreesRF <- 500


# organize data in outcome list objects
# number of list elements to save as outcome measures
#   ... for every model GBM, ENET, RF
setParam$fit$out <- 4 # performTrainStats, performTestStats, performPerSample, pvi
setParam$fit$outLabels <- c("performTrainStats", "performTestStats", "performPerSample", "pvi")
#   ... specifically for GBM 
setParam$fit$outGBM <- 3 # interStrength, selectionPerSample
setParam$fit$gbmLabels <- c("performCVtestStats", "interStrength", "selectionPerSample")
#   ... specifically for ENET 
setParam$fit$outENET <- 4 # estBeta, estBetaFull, varSelection, selectionPerSample 
setParam$fit$enetLabels <- c("estBeta", "estBetaFull", "varSelection", "selectionPerSample")
#   ... specifically for RF
setParam$fit$outRF <- 4 
setParam$fit$rfLabels <- c("performCVtestStats", "oobPredictions", "oobR2", "selectionPerSample")