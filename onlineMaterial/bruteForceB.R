# brute force algorithm to find beta coefficients with fixed Rsquared
# Gibbs Sampling to arrive at beta coefficients that lead to respective Rsquared 
#   for linear effects and interactions

# restrictions
#     - Var1 = Var2 = Var3 = Var4 = ?
#     - Var1:Var2 = Var1:Var4 = Var2:Var3 = Var3:Var4 = ?
#     - Var1:Var3 = Var2:Var4 = 0

# Rsquared either splits 20:80, 50:50 or 80:20 for interaction and linear effects respectively
# squared multiple correlation (r_y,yhat)^2 = multiple coefficient of determination
# multiple coefficient of determination = sum of semipartial determinantes (squared 
#     semipartial correlations) of increasingly higher order

# load packages
library(mvtnorm)
library(truncnorm)
library(parallel)

# load parameters & custom functions 
source("utils/setParameters.R") # parameter values
source("utils/simTools.R") # functions for data simulation

# generate folder for log files
logFolder = "log"
createFolder(logFolder)

pTrash <- setParam$bruteForceB$pTrash
N <- setParam$bruteForceB$N
reliability <- setParam$bruteForceB$reliability

P <- setParam$dgp$p + pTrash # total number of variables

X <- createPredictors(N = N, P = P, 
                      corMat = setParam$dgp$predictorCorMat[seq_len(P), seq_len(P)])

# # simulate X with uncorrelated predictors
# nullCorMat <- matrix(0, ncol = P, nrow = P)
# diag(nullCorMat) <- 1
# X <- createPredictors(N = N, P = P, 
#                       corMat = nullCorMat)

# add names to variables
colnames(X) <- paste0("Var", seq_len(P))

# create model formula (allows polynomial and interaction effects of any degree/depth)
popModel <- genModel(colnames(X), setParam$dgp$interDepth, setParam$dgp$poly)

# predictor matrix that allows for polynomials and interactions
X_int <- model.matrix(as.formula(popModel),data.frame(X))

# remove first degree polynomials from data (they are duplicates!)
if (setParam$dgp$poly > 0) {
  X_int <- rmDuplicatePoly(X_int)  
}

(empCor <- cor(X_int))
# round(empCor, 3)
cov(X_int)

# # for independent normally distributed random variables 
# mVar1 <- setParam$dgp$meanR
# sdVar1 <- sqrt(diag(setParam$dgp$predictorCorMat)[1])
# 
# mVar2 <- setParam$dgp$meanR
# sdVar2 <- sqrt(diag(setParam$dgp$predictorCorMat)[2])
# var12_ind <- (mVar2 * sdVar1)**2 + (mVar1 * sdVar2)**2 + sdVar1**2 * sdVar2**2

# the variance of the product of two correlated normally distributed random variables
rho <- setParam$dgp$predictorCorMat[1,2] # correlation between linear predictors
varProduct <- (1+rho**2)
sdProduct <- sqrt(varProduct)
# # without the assumption of standardized predictors (i.e., zero means and unit variance)
# var12 <- (mVar2 * sdVar1)**2 + (mVar1 * sdVar2)**2 + 
#   (sdVar1 * sdVar2)**2 * (1+rho**2) + 
#   2*mVar1*mVar2*sdVar1*sdVar2

# covariances and correlations between two products of random variables from a 
#   multivariate normal distribution
covCommon <- rho + (rho)**2
corCommon <- covCommon/(sdProduct * sdProduct)

covExclusive <- (rho)**2 + (rho)**2
corExclusive <- covExclusive/(sdProduct * sdProduct)
covExclusive/varProduct

# build correlation matrix as block structure with ...
#   ... block for linear predictors
#   ... block for inetraction terms (products)
corLin <- lavaan::lav_matrix_upper2full(rep(rho, setParam$dgp$p * (setParam$dgp$p-1) / 2), diagonal = F)
diag(corLin) <- 1

nInter <- setParam$dgp$p * (setParam$dgp$p-1) / 2
corInter <- lavaan::lav_matrix_upper2full(rep(corCommon, nInter * (nInter-1)/2), diagonal = F)
diag(corInter[,c(nInter:1)]) <- corExclusive
diag(corInter) <- 1

# full correlation matrix (+ labels)
theoCor <- as.matrix(Matrix::bdiag(corLin, corInter))
rownames(theoCor) = colnames(theoCor) = c(setParam$dgp$linEffects,
                                          "Var1:Var2", "Var1:Var3", "Var1:Var4",
                                          "Var2:Var3", "Var2:Var4", "Var3:Var4")

# # test formula in generall:
# #     for 50:50 and different fixed Rsquared
# setParam$dgp$Rsquared[3] * setParam$dgp$sigmaE / (1-setParam$dgp$Rsquared[3])
# t(c(rep(0.346, 8), 0, 0)) %*% empCor %*% c(rep(0.346, 8), 0, 0)
# 
# setParam$dgp$Rsquared[4] * setParam$dgp$sigmaE / (1-setParam$dgp$Rsquared[4])
# t(c(rep(0.693, 8), 0, 0)) %*% empCor %*% c(rep(0.693, 8), 0, 0)
# 
# # theta <- c(0.5009557, 0.4536119)
# beta <- vector(mode = "numeric", length = dim(empCor)[1])
# names(beta) <- colnames(empCor)
# beta[names(beta) %in% setParam$dgp$linEffects] <- 0.5009557
# beta[names(beta) %in% setParam$dgp$interEffects] <- 0.4536119
# t(beta) %*% empCor %*% beta

################################################################################
# functions
################################################################################
linOptim <- function(theta, R2, lin, inter, beta_inter) {
  # check if the argument values meet conditions  
  if (lin + inter != 1) {
    stop("lin and inter do not sum up to 1!")
  }
  
  ### linear effects and interaction effects
  # generate zero/empty vector of regression coefficients
  beta <- vector(mode = "numeric", length = dim(theoCor)[1])
  names(beta) <- colnames(theoCor) # name vector
  # only coefficients for simulated linear effects
  beta[names(beta) %in% setParam$dgp$linEffects] <- theta
  # add coefficients for interactions with simulated effects
  beta[names(beta) %in% setParam$dgp$interEffects] <- beta_inter
  
  ### only linear effects
  # generate zero/empty vector of regression coefficients
  betaVecLin <- vector(mode = "numeric", length = dim(theoCor)[1])
  names(betaVecLin) <- colnames(theoCor) # name vector
  # fill coefficients vector only with simulated linear effects
  betaVecLin[names(betaVecLin) %in% setParam$dgp$linEffects] <- theta 
  
  # calculate the explained variance for only linear effects
  R2_lin <- var(X_int %*% betaVecLin) / (var(X_int %*% beta) + setParam$dgp$sigmaE^2)
  
  # how much does the estimated R2_lin deviate from the target R2 for linear effects
  return(abs(R2_lin - (R2 * lin)))
}

interOptim <- function(theta, R2, lin, inter, beta_lin) {
  # check if the argument values meet conditions  
  if (lin + inter != 1) {
    stop("lin and inter do not sum up to 1!")
  }
  
  ### linear effects and interaction effects
  # generate zero/empty vector of regression coefficients
  beta <- vector(mode = "numeric", length = dim(theoCor)[1])
  names(beta) <- colnames(theoCor) # name vector
  # only coefficients for simulated linear effects
  beta[names(beta) %in% setParam$dgp$linEffects] <- beta_lin 
  beta[names(beta) %in% setParam$dgp$interEffects] <- theta
  
  ### only interaction effects
  # generate zero/empty vector of regression coefficients
  betaVecInter <- vector(mode = "numeric", length = dim(theoCor)[1])
  names(betaVecInter) <- colnames(theoCor) # name vector
  # fill coefficients vector only with simulated interaction effects
  betaVecInter[names(betaVecInter) %in% setParam$dgp$interEffects] <- theta
  
  # calculate the explained variance for only interaction effects
  R2_inter <- var(X_int %*% betaVecInter) / (var(X_int %*% beta) + setParam$dgp$sigmaE^2)
  
  # how much does the estimated R2_inter deviate from the target R2 for interaction effects
  return(abs(R2_inter - (R2 * inter)))
}

checkOptim <- function(betaLin, betaInter){
  # generate zero/empty vector of regression coefficients
  beta <- vector(mode = "numeric", length = dim(theoCor)[1])
  names(beta) <- colnames(theoCor) # name vector
  # add coefficients for simulated linear effects and for interactions with simulated effects
  beta[names(beta) %in% setParam$dgp$linEffects] <- betaLin
  beta[names(beta) %in% setParam$dgp$interEffects] <- betaInter
  
  betaVecLin <- vector(mode = "numeric", length = dim(theoCor)[1])
  names(betaVecLin) <- colnames(theoCor)
  betaVecLin[names(betaVecLin) %in% setParam$dgp$linEffects] <- betaLin
  
  betaVecInter <- vector(mode = "numeric", length = dim(theoCor)[1])
  names(betaVecInter) <- colnames(theoCor)
  betaVecInter[names(betaVecInter) %in% setParam$dgp$interEffects] <- betaInter
  
  R2_lin <- var(X_int %*% betaVecLin) / (var(X_int %*% beta) + setParam$dgp$sigmaE^2)
  R2_inter <- var(X_int %*% betaVecInter) / (var(X_int %*% beta) + setParam$dgp$sigmaE^2)
  
  R2_total <- getR2(X_int, beta, setParam$dgp$sigmaE)
  R2_add <- R2_lin + R2_inter
  
  betaOptim <- t(beta) %*% theoCor %*% c(beta)
  diff <- abs((setParam$dgp$Rsquared[4] * setParam$dgp$sigmaE / (1 - setParam$dgp$Rsquared[4])) -
                betaOptim)
  
  resVec <- c(R2_lin, R2_inter, R2_total, R2_add, diff)
  names(resVec) <- c("R2_lin", "R2_inter", "R2_total", "R2_add", "diff")
  return(resVec)
}

################################################################################

# two step optimization?
# 1. optim betaLin & betaInter via Gibbs Sampling

# parameters to estimate:
#     beta coefficient for all linear predictors
#     beta coefficient for all interactions

init <- c(1, 0.1)
names(init) <- c("betaLin", "betaInter")

condGrid <- expand.grid(R2 = setParam$dgp$Rsquared,
                        lin = setParam$dgp$percentLinear)
condGrid$inter <- 1 - condGrid$lin

optimBeta <- function(init, R2, lin, inter) {
  
  # estimate beta coefficient for linear effects conditioned on current value for
  #   interaction effects
  init["betaLin"] <- optim(par = init["betaLin"], # parameters and their initial value
                           fn = linOptim, # optimization criterion
                           R2 = R2, # fixed R2
                           lin = lin, inter = inter,
                           beta_inter = init["betaInter"],
                           method = "L-BFGS-B",
                           lower= setParam$fit$optimLowerLimit, # only positive beta coefficients
                           upper = setParam$fit$optimUpperLimit)$par
  
  # estimate beta coefficient for interaction effects conditioned on current value for
  #   linear effects
  init["betaInter"] <- optim(par = init["betaInter"], # parameters and their initial value
                             fn = interOptim, # optimization criterion
                             R2 = R2, # fixed R2
                             lin = lin, inter = inter,
                             beta_lin = init["betaLin"],
                             method = "L-BFGS-B",
                             lower= setParam$fit$optimLowerLimit, # only positive beta coefficients
                             upper = setParam$fit$optimUpperLimit)$par
  
  return(init)
}

# # test optim functions for betaLin and betaInter
# optimBeta(init, 
#           setParam$dgp$Rsquared[4], 
#           setParam$dgp$percentLinear[1],
#           setParam$dgp$percentInter[1])
# 
# # test check function
# checkOptim(init["betaLin"], init["betaInter"], setParam$dgp$Rsquared[4])

gibbsB <- function(init, R2, lin, inter) {
  
  repeat{
    
    # optimize beta for linear and interaction effects conditioned on the current other value
    tmp_init <- optimBeta(init, R2, lin, inter)
    print(tmp_init)
    
    # evaluate performance
    accEstBeta <- checkOptim(tmp_init["betaLin"], tmp_init["betaInter"])
    print(accEstBeta)
    
    # check that every R2 (lin, inter, total) is near enough around target R2
    acc <- c(abs((lin*R2) - accEstBeta["R2_lin"]), 
             abs((inter*R2) - accEstBeta["R2_inter"]), 
             # R2_total would be more accurate but does not necessarily converge 
             #    getting R2_total AND {R2_lin, R2_inter} right is more difficult
             abs(R2 - accEstBeta["R2_total"])) 
             # abs(R2 - accEstBeta["R2_add"])) 
    
    # repeat previous steps if criterions for convergence are not met
    #   a) R2_{lin, inter, add} near enough around target R2
    #   b) changes in beta{Lin, Inter} smaller than tolerance parameter
    if (all(acc < setParam$fit$optimTol) | all(abs(init - tmp_init) < setParam$fit$optimBetaTol)) {
      
      print(acc < setParam$fit$optimTol) 
      print(abs(init - tmp_init) < setParam$fit$optimBetaTol)
      
      return(list(init = tmp_init, 
                  checkAcc = accEstBeta))
      break
    } else {
      init <- tmp_init
    }
  }
}

# # test gibbs sampling including accuracy checks for target R2
# gibbsB(init, 
#        setParam$dgp$Rsquared[4], 
#        setParam$dgp$percentLinear[1],
#        setParam$dgp$percentInter[1])

# timeStampFolder <- format(Sys.time(), "%d%m%y_%H%M%S")
# nCoresSampling <- 6 # brute force beta estimation
# 
# # Initiate cluster; type = "FORK" only on Linux/MacOS: contains all environment variables automatically
# cl <- makeCluster(nCoresSampling, type = "FORK",
#                   outfile = paste0(logFolder, "/", "bruteForceBeta",
#                                    timeStampFolder, ".txt"))
# 
# # set seed that works for parallel processing
# set.seed(6723940)
# s <- .Random.seed
# clusterSetRNGStream(cl = cl, iseed = s)

# bruteForceB <- parLapply(cl, seq_len(dim(condGrid)[1]), function(iGrid) {
bruteForceB <- lapply(seq_len(dim(condGrid)[1]), function(iGrid) {
  
  # reset init
  init <- c(1, 0.1)
  names(init) <- c("betaLin", "betaInter")
  
  gibbsB(init, 
         condGrid[iGrid, "R2"], 
         condGrid[iGrid, "lin"], 
         condGrid[iGrid, "inter"])
})

# # close cluster to return resources (memory) back to OS
# stopCluster(cl)

names(bruteForceB) <- setParam$dgp$condLabels

betaCoef <- do.call(rbind, lapply(seq_along(bruteForceB), function(subList) {
  bruteForceB[[subList]][["init"]]
}))

checkAcc <- do.call(rbind, lapply(seq_along(bruteForceB), function(subList) {
  bruteForceB[[subList]][["checkAcc"]]
}))

# getOption("scipen") # 0
# options(scipen = 999)
# options(scipen = 0)

betaData <- cbind(condGrid, betaCoef, checkAcc)
round(betaData, 3)

# write.csv(betaData, "utils/bruteForceBcoeff_inter.csv", row.names=FALSE)
# bruteForceB <- read.table("utils/bruteForceBcoeff_inter.csv", header = T, sep = ",")
