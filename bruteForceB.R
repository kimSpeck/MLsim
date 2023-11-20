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

# load parameters & custom functions 
source("setParameters.R") # parameter values
source("simTools.R") # functions for data simulation

setParam$dgp$interEffects <- c("Var1:Var2", "Var1:Var4", "Var2:Var3", "Var3:Var4")

pTrash <- 0
N <- 100000
reliability <- 1
setParam$dgp$poly <- 0

P <- setParam$dgp$p + pTrash # total number of variables

# to do: remove function and use rmvnorm directly now?!
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
  if (lin + inter != 1) {
    stop("lin and inter do not sum up to 1!")
  }
  
  beta <- vector(mode = "numeric", length = dim(empCor)[1])
  names(beta) <- colnames(empCor)
  beta[names(beta) %in% setParam$dgp$linEffects] <- theta 
  beta[names(beta) %in% setParam$dgp$interEffects] <- beta_inter
  
  # distribute effect size  
  Y <- calcDV(X = X_int, b = beta, sigmaE = setParam$dgp$sigmaE, N = N)
  
  
  resVar1 <- residuals(lm(Var1 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4, 
                          data = X_int_df))
  resVar2 <- residuals(lm(Var2 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4 + 
                            Var1, data = X_int_df))
  resVar3 <- residuals(lm(Var3 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4 + 
                            Var1 + Var2, data = X_int_df))
  resVar4 <- residuals(lm(Var4 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4 + 
                            Var1 + Var2 + Var3, data = X_int_df))
  (R2_lin <- cor(resVar1, Y)**2 + cor(resVar2, Y)**2 + cor(resVar3, Y)**2 + cor(resVar4, Y)**2)
  
  return(abs(R2_lin - (R2 * lin)))
}

interOptim <- function(theta, R2, lin, inter, beta_lin) {
  if (lin + inter != 1) {
    stop("lin and inter do not sum up to 1!")
  }
  
  beta <- vector(mode = "numeric", length = dim(empCor)[1])
  names(beta) <- colnames(empCor)
  beta[names(beta) %in% setParam$dgp$linEffects] <- beta_lin 
  beta[names(beta) %in% setParam$dgp$interEffects] <- theta
  
  # distribute effect size  
  Y <- calcDV(X = X_int, b = beta, sigmaE = setParam$dgp$sigmaE, N = N)
  
  #
  resVar12 <- residuals(lm(Var1.Var2 ~ Var1 + Var2 + Var3 + Var4,
                           data = X_int_df))
  resVar14 <- residuals(lm(Var1.Var4 ~ Var1 + Var2 + Var3 + Var4 +
                             Var1.Var2, data = X_int_df))
  resVar23 <- residuals(lm(Var2.Var3 ~ Var1 + Var2 + Var3 + Var4 +
                             Var1.Var2 + Var1.Var4, data = X_int_df))
  resVar34 <- residuals(lm(Var3.Var4 ~ Var1 + Var2 + Var3 + Var4 +
                             Var1.Var2 + Var1.Var4 + Var2.Var3, data = X_int_df))
  (R2_inter <- cor(resVar12, Y)**2 + cor(resVar14, Y)**2 + cor(resVar23, Y)**2 + cor(resVar34, Y)**2)
  
  return(abs(R2_inter - (R2 * inter)))
}

checkOptim <- function(betaLin, betaInter, R2){
  # 
  beta <- vector(mode = "numeric", length = dim(empCor)[1])
  names(beta) <- colnames(empCor)
  beta[names(beta) %in% setParam$dgp$linEffects] <- betaLin
  beta[names(beta) %in% setParam$dgp$interEffects] <- betaInter
  
  # distribute effect size
  Y <- calcDV(X = X_int, b = beta, sigmaE = setParam$dgp$sigmaE, N = N)
  
  resVar1 <- residuals(lm(Var1 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4,
                          data = X_int_df))
  resVar2 <- residuals(lm(Var2 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4 +
                            Var1, data = X_int_df))
  resVar3 <- residuals(lm(Var3 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4 +
                            Var1 + Var2, data = X_int_df))
  resVar4 <- residuals(lm(Var4 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4 +
                            Var1 + Var2 + Var3, data = X_int_df))
  R2_lin <- cor(resVar1, Y)**2 + cor(resVar2, Y)**2 + cor(resVar3, Y)**2 + cor(resVar4, Y)**2
  
  #
  resVar12 <- residuals(lm(Var1.Var2 ~ Var1 + Var2 + Var3 + Var4,
                           data = X_int_df))
  resVar14 <- residuals(lm(Var1.Var4 ~ Var1 + Var2 + Var3 + Var4 +
                             Var1.Var2, data = X_int_df))
  resVar23 <- residuals(lm(Var2.Var3 ~ Var1 + Var2 + Var3 + Var4 +
                             Var1.Var2 + Var1.Var4, data = X_int_df))
  resVar34 <- residuals(lm(Var3.Var4 ~ Var1 + Var2 + Var3 + Var4 +
                             Var1.Var2 + Var1.Var4 + Var2.Var3, data = X_int_df))
  R2_inter <- cor(resVar12, Y)**2 + cor(resVar14, Y)**2 + cor(resVar23, Y)**2 + cor(resVar34, Y)**2
  
  R2_total <- getR2(X_int, beta, setParam$dgp$sigmaE)
  R2_add <- R2_lin + R2_inter
  
  betaOptim <- t(beta) %*% empCor %*% c(beta)
  diff <- abs((R2 * setParam$dgp$sigmaE / (1 - R2)) -
                betaOptim)
  
  resVec <- c(R2_lin, R2_inter, R2_total, R2_add, diff)
  names(resVec) <- c("R2_lin", "R2_inter", "R2_total", "R2_add", "diff")
  return(resVec)
}

# #
# diffOptim <- function(theta, R2, lin, inter){
#   if (lin + inter != 1) {
#     stop("lin and inter do not sum up to 1!")
#   }
#   
#   # estimate beta for linear 
#   init <- 0.1
#   theta[1] <- optim(par = init, # parameters and their initial value
#                     fn = linOptim, # optimization criterion
#                     R2 = R2, # fixed R2
#                     lin = lin, inter = inter, 
#                     beta_inter = theta[2],
#                     method = "L-BFGS-B",
#                     lower= c(0), # only positive beta coefficients
#                     upper = c(2))$par
#   print(theta[1])
#   
#   init <- 0.1
#   theta[2] <- optim(par = init, # parameters and their initial value
#                     fn = interOptim, # optimization criterion
#                     R2 = R2, # fixed R2
#                     lin = lin, inter = inter, 
#                     beta_lin = theta[1],
#                     method = "L-BFGS-B",
#                     lower= c(0), # only positive beta coefficients
#                     upper = c(2))$par
#   print(theta[2])
#   
#   # 
#   beta <- vector(mode = "numeric", length = dim(empCor)[1])
#   names(beta) <- colnames(empCor)
#   beta[names(beta) %in% setParam$dgp$linEffects] <- theta[1] 
#   beta[names(beta) %in% setParam$dgp$interEffects] <- theta[2]
#   
#   # # distribute effect size  
#   # Y <- calcDV(X = X_int, b = beta, sigmaE = setParam$dgp$sigmaE, N = N)
#   # 
#   # resVar1 <- residuals(lm(Var1 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4, 
#   #                         data = X_int_df))
#   # resVar2 <- residuals(lm(Var2 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4 + 
#   #                           Var1, data = X_int_df))
#   # resVar3 <- residuals(lm(Var3 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4 + 
#   #                           Var1 + Var2, data = X_int_df))
#   # resVar4 <- residuals(lm(Var4 ~ Var1:Var2 + Var1:Var4 + Var2:Var3 + Var3:Var4 + 
#   #                           Var1 + Var2 + Var3, data = X_int_df))
#   # (R2_lin <- cor(resVar1, Y)**2 + cor(resVar2, Y)**2 + cor(resVar3, Y)**2 + cor(resVar4, Y)**2)  
#   # 
#   # #
#   # resVar12 <- residuals(lm(Var1.Var2 ~ Var1 + Var2 + Var3 + Var4,
#   #                          data = X_int_df))
#   # resVar14 <- residuals(lm(Var1.Var4 ~ Var1 + Var2 + Var3 + Var4 +
#   #                            Var1.Var2, data = X_int_df))
#   # resVar23 <- residuals(lm(Var2.Var3 ~ Var1 + Var2 + Var3 + Var4 +
#   #                            Var1.Var2 + Var1.Var4, data = X_int_df))
#   # resVar34 <- residuals(lm(Var3.Var4 ~ Var1 + Var2 + Var3 + Var4 +
#   #                            Var1.Var2 + Var1.Var4 + Var2.Var3, data = X_int_df))
#   # (R2_inter <- cor(resVar12, Y)**2 + cor(resVar14, Y)**2 + cor(resVar23, Y)**2 + cor(resVar34, Y)**2)
#   
#   # getR2(X_int, beta, setParam$dgp$sigmaE)
#   # R2_lin + R2_inter
#   
#   # # alternative: get R2 the "usual way" with only beta weights in linear effects
#   # #   calculating R2 this way highly overestimates R2
#   # beta_lin <- vector(mode = "numeric", length = dim(empCor)[1])
#   # names(beta_lin) <- colnames(empCor)
#   # beta_lin[names(beta_lin) %in% setParam$dgp$linEffects] <- init[1]
#   # getR2(X_int, beta_lin, setParam$dgp$sigmaE)
#   # 
#   # beta_inter <- vector(mode = "numeric", length = dim(empCor)[1])
#   # names(beta_inter) <- colnames(empCor)
#   # beta_inter[names(beta_inter) %in% setParam$dgp$interEffects] <- init[2]
#   # getR2(X_int, beta_inter, setParam$dgp$sigmaE)
#   
#   betaOptim <- t(beta) %*% empCor %*% c(beta)
#   diff <- abs((R2 * setParam$dgp$sigmaE / (1 - R2)) -
#                 betaOptim)
#   
#   print(diff)
#   return(diff)
# }
# 
# # example
# optim(par = init, # parameters and their initial value
#       fn = diffOptim, # optimization criterion
#       R2 = setParam$dgp$Rsquared[4], # fixed R2
#       lin = setParam$dgp$percentLinear[1], inter = setParam$dgp$percentInter[1], 
#       method = "L-BFGS-B",
#       lower= c(0, 0), # only positive beta coefficients
#       upper = c(2, 2))

################################################################################

# two step optimization?
# 1. optim betaLin & betaInter via Gibbs Sampling

# parameters to estimate:
#     beta coefficient for all linear predictors
#     beta coefficient for all interactions

X_int_df <- data.frame(X_int)


init <- c(1, 0.1)
names(init) <- c("betaLin", "betaInter")

condGrid <- expand.grid(R2 = setParam$dgp$Rsquared,
                        lin = setParam$dgp$percentLinear)
condGrid$inter <- 1 - condGrid$lin


setParam$fit$optimLowerLimit <- 0
setParam$fit$optimUpperLimit <- 2
setParam$fit$optimTol <- 5*1e-3
setParam$fit$optimBetaTol <- 1e-4

optimBeta <- function(init, R2, lin, inter) {
  
  init["betaLin"] <- optim(par = init["betaLin"], # parameters and their initial value
                           fn = linOptim, # optimization criterion
                           R2 = R2, # fixed R2
                           lin = lin, inter = inter,
                           beta_inter = init["betaInter"],
                           method = "L-BFGS-B",
                           lower= setParam$fit$optimLowerLimit, # only positive beta coefficients
                           upper = setParam$fit$optimUpperLimit)$par
  
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
    accEstBeta <- checkOptim(tmp_init["betaLin"], tmp_init["betaInter"], R2)
    print(accEstBeta)
    # check that every R2 (lin, inter, total) is near enough around target R2
    acc <- c(abs((lin*R2) - accEstBeta["R2_lin"]), 
             abs((inter*R2) - accEstBeta["R2_inter"]), 
             abs(R2 - accEstBeta["R2_add"])) # R2_total would be more accurate but does not converge as fast
    
    if (all(acc < setParam$fit$optimTol) | all(abs(init - tmp_init) < setParam$fit$optimBetaTol)) {
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

# lapply(seq_len(dim(condGrid)[1]), function(iGrid) {
lapply(seq_len(2), function(iGrid) {
  gibbsB(init, 
         condGrid[iGrid, "R2"], 
         condGrid[iGrid, "lin"], 
         condGrid[iGrid, "inter"])
})
