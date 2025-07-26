# brute force algorithm to find beta coefficients with fixed Rsquared
# Gibbs Sampling to arrive at beta coefficients that lead to respective Rsquared 
#   for linear effects and nonlinear effects (stepwise)

# restrictions
#     - Var1 = Var2 = Var3 = Var4 = ?
#     - dumVar5.1 = dumVar6.1 = dumVar7.1 = ?

# Rsquared either splits 20:80, 50:50 or 80:20 for nonlinear and linear effects respectively
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

P <- setParam$dgp$p + setParam$dgp$pNL3 + pTrash # total number of variables
X <- createPredictors(N = N, P = P, 
                      corMat = setParam$dgp$predictorCorMat_pwl[seq_len(P), seq_len(P)])

# add names to variables
colnames(X) <- paste0("Var", seq_len(P))

# create model formula (allows polynomial and interaction effects of any degree/depth)
popModel <- genModel(colnames(X), setParam$dgp$interDepth, setParam$dgp$poly)

# predictor matrix that allows for polynomials and interactions
X_int <- model.matrix(as.formula(popModel), data.frame(X))

# remove first degree polynomials from data (they are duplicates!)
if (setParam$dgp$poly > 0) {
  X_int <- rmDuplicatePoly(X_int)  
}


X_int <- cbind(X_int, 
               dumVar5.1 = createDummy(X_int[, "Var5"], q = 0.5, effectCoding = T),
               dumVar6.1 = createDummy(X_int[, "Var6"], q = 0.5, effectCoding = T),
               dumVar7.1 = createDummy(X_int[, "Var7"], q = 0.5, effectCoding = T))
  
# remove every other variable from data
X_int <- X_int[,colnames(X_int) %in% c(setParam$dgp$linEffects, setParam$dgp$nonlinEffects3)]

(empCor <- cor(X_int))
# round(empCor, 3)
cov(X_int)

# the variance of the product of two correlated normally distributed random variables
rho <- setParam$dgp$predictorCorMat[1,2] # correlation between linear predictors

# build correlation matrix as block structure with ...
#   ... block for linear predictors
#   ... block for inetraction terms (products)
corLin <- lavaan::lav_matrix_upper2full(rep(rho, setParam$dgp$p * (setParam$dgp$p-1) / 2), diagonal = F)
diag(corLin) <- 1

nNL <- setParam$dgp$pNL3
corInter <- lavaan::lav_matrix_upper2full(rep(0, nNL * (nNL-1)/2), diagonal = F)
diag(corInter) <- 1

# full correlation matrix (+ labels)
theoCor <- as.matrix(Matrix::bdiag(corLin, corInter))
rownames(theoCor) = colnames(theoCor) = c(setParam$dgp$linEffects,
                                          setParam$dgp$nonlinEffects3)

################################################################################
# functions
################################################################################
linOptim <- function(theta, R2, lin, inter, beta_inter) {
  # check if the argument values meet conditions  
  if (lin + inter != 1) {
    stop("lin and inter do not sum up to 1!")
  }
  
  # generate zero/empty vector of regression coefficients
  beta <- betaVecLin <- vector(mode = "numeric", length = dim(theoCor)[1])
  names(beta) <- names(betaVecLin) <-  colnames(theoCor) # name vector
  linIdx <- names(beta) %in% setParam$dgp$linEffects
  beta[linIdx] <- betaVecLin[linIdx] <- theta 
  # beta[names(beta) %in% setParam$dgp$nonlinEffects] <- beta_inter
  beta[names(beta) %in% setParam$dgp$nonlinEffects3] <- beta_inter
  
  R2_lin <- var(X_int %*% betaVecLin) / (var(X_int %*% beta) + setParam$dgp$sigmaE^2)
  
  return(abs(R2_lin - (R2 * lin)))
}

interOptim <- function(theta, R2, lin, inter, beta_lin) {
  # check if the argument values meet conditions
  if (lin + inter != 1) {
    stop("lin and inter do not sum up to 1!")
  }
  
  beta <- betaVecInter <-  vector(mode = "numeric", length = dim(theoCor)[1])
  names(beta) <- names(betaVecInter) <- colnames(theoCor)
  # interIdx <- names(beta) %in% setParam$dgp$nonlinEffects
  interIdx <- names(beta) %in% setParam$dgp$nonlinEffects3
  beta[names(beta) %in% setParam$dgp$linEffects] <- beta_lin 
  beta[interIdx] <- betaVecInter[interIdx] <-  theta
  
  R2_inter <- var(X_int %*% betaVecInter) / (var(X_int %*% beta) + setParam$dgp$sigmaE^2)
  
  return(abs(R2_inter - (R2 * inter)))
}

checkOptim <- function(betaLin, betaInter){
  #
  beta <- betaVecLin <- betaVecInter <- vector(mode = "numeric", length = dim(theoCor)[1])
  names(beta) <- names(betaVecInter) <- names(betaVecLin) <- colnames(theoCor)
  linIdx <- names(beta) %in% setParam$dgp$linEffects
  # interIdx <- names(beta) %in% setParam$dgp$nonlinEffects
  interIdx <- names(beta) %in% setParam$dgp$nonlinEffects3
  beta[linIdx] <- betaLin
  beta[interIdx] <- betaInter
  
  betaVecLin[linIdx] <- betaLin
  betaVecInter[interIdx] <- betaInter
  
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

# write.csv(betaData, "utils/bruteForceBcoeff_nonlinear3.csv", row.names=FALSE)
# bruteForceB <- read.table("utils/bruteForceBcoeff_nonlinear3.csv", header = T, sep = ",")
