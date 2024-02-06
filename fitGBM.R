# gbm fitting

# fitted data in: (load this!)
#   testGBMresults.RData

# Vergleich von Enet-Results mit GBM-Results:
# Rsquared in test sample: the same for both models
# Rangreihe von extrahierten Prädiktoren ~ Variable importance measures
# welche Prädiktoren wurden ausgewählt?
#   gewählte Interaktionen in GBM?
# literature:
#   J.H. Friedman and B.E. Popescu (2005). “Predictive Learning via Rule Ensembles.” Section 8.1
#   H-statistics to quantify interactions
# package:
#   https://cran.r-project.org/web/packages/EIX/vignettes/EIX.html
# idea -> Interaktionen als Features im Sinne von feature engineering als Input zu GBM?
#   (wie geht GBM damit um? wählt GBM Interaktionen überhaupt?)
#   (pro: einfach zu quantifizieren, 
#   con: eigentlicher Vorteil des GBMs geht verloren + unklar, ob Interaktion überhaupt gewählt wird)

# for comparison with Enet without interaction we could constraint GBM to exclude 
#   all interactions 
# https://xgboost.readthedocs.io/en/stable/tutorials/feature_interaction_constraint.html

# params_constrained = params.copy()
# # Use nested list to define feature interaction constraints
# params_constrained['interaction_constraints'] = '[[0, 2], [1, 3, 4], [5, 6]]'
# # Features 0 and 2 are allowed to interact with each other but with no other feature
# # Features 1, 3, 4 are allowed to interact with one another but with no other feature
# # Features 5 and 6 are allowed to interact with each other but with no other feature

# to do: 
# was ist ein fairer Vergleich? Enet vs. GBM mit oder ohne Regularisierung irgendeiner Form?!
# add dropout to reduce overspecialization? (no: if it only increases Rsquared, yes: if it improves variable importance measures)
#   dropout: removing randomly selected trees from the decision tree ensemble;
#   i.e., instead of developing the next tree from the residual of all previous trees, 
#         develop the next tree from the residual of a sample of previous trees (thus, more information on residuals to fit next tree on)
#     -> Dropout Additive Regression Trees (DART) [argument: booster = "dart" in xgboost]
#     in xgboost the percentage of dropout as parameter that can be set in the tuning of the model.
#   (this is not equivalent to input-layer dropout: colsample_bytree and colsample_bylevel!)
#   with dropout regularization runtimes become longer!!!
#   Vinayak, R. K., & Gilad-Bachrach, R. (2015, February). Dart: Dropouts meet 
#   multiple additive regression trees. In Artificial Intelligence and Statistics (pp. 489-497). PMLR.

##### get one single data set to test code #####
# load parameters & custom functions 
source("setParameters.R") # parameter values
source("simTools.R")

# install.packages('xgboost')
# devtools :: install_github("ModelOriented/EIX")
library(xgboost)
library(iml)
# library(EIX)
# library(h2o)

dataFolder <- "data"
# condGrid <- expand.grid(N = setParam$dgp$N, 
#                         pTrash = setParam$dgp$pTrash)
# iSim = 1
# Probedatensatz mit 100 Menschen & 10 trash Variablen

# fileName <- paste0("simDataN", condGrid[iSim, "N"], "_pTrash", condGrid[iSim, "pTrash"], ".rda")
# load(paste0(dataFolder, "/", fileName))

load("data/simDataN1000_pTrash10_rel1_f0.rda")
iSample <- 1
colnames(data[[iSample]][["yMat"]])
iCond <- 3 # R² = 0.8; 50% lin vs. 50% inter

# data from single sample (with interaction as specified in elastic net)
XtrainWithInter <- as.matrix(data[[iSample]][["X_int"]])
# remove interaction to test if GBM is able to find interactions on its own
idxPoly <- stringr::str_detect(colnames(XtrainWithInter), "^poly")
idxInter <- stringr::str_detect(colnames(XtrainWithInter), ":")
Xtrain <- XtrainWithInter[, -c(which(idxPoly), which(idxInter))]

ytrain <- data[[iSample]][["yMat"]][,iCond]

#####

# original gbm package (https://github.com/cran/gbm)
#   is retired/no longer under active development
#   but functions from the gbm-package are used inside of caret
#   see line 68: https://github.com/topepo/caret/blob/master/models/files/gbm.R

# xgboost
#   approx. 10x faster than gbm (implemented in C++) + still under development
#   (+ options for regularised gbms to prevent overfitting)
#   https://github.com/dmlc/xgboost/tree/master/R-package

# parameter information
# https://xgboost.readthedocs.io/en/latest/parameter.html#learning-task-parameters

# only works for matrices that contain all numeric variables
#   dummy coding for categorical/binary variables beforehand!
# scaling and dummy coding is not necessary (so far) because we only have continuous
#     predictors on identical scales

# train model
# input data with xgboosts own class (recommended)
# https://xgboost.readthedocs.io/en/stable/R-package/xgboostPresentation.html
trainData <- xgboost::xgb.DMatrix(data = Xtrain, label = ytrain, missing = NA)
# testData <- xgboost::xgb.DMatrix(data = Xtest, label = ytest, missing = NA)
# xgboost::setinfo(x, "label", y)

# (tuning) parameter grid
# number of trees, depth of trees, learning rate, subsampling?
largeTuneGrid <- expand.grid(eta = seq(.001, .201, .02),         # shrinkage/learning rate
                        max_depth = c(1,2,3,4,5),                # tree depth; interaction depth hängt von simulation ab
                        min_child_weight = c(5,10,20,50),        # min node size
                        # subsample? # percent of training data to sample for each tree?/bag.fraction in gbm?
                        # colsample_bytree? # percent of columns to sample from for each tree
                        # nTrees = c(50,100,150,300,500,1000),     # max number of trees to run
                        # does max number of trees become irrelevant to tune as soon as we introduce early stopping?
                        optimalTrees = NA,                        # fill in cv with nTrees for each grid line    
                        minRMSE = NA)                             # fill in cv with RMSE for each grid line

# smaller grid that matches the grid for gbm (to compare run time)
tuneGrid <- expand.grid(max_depth = c(1,2,3), # tree depth (= interaction.depth in gbm)
                        min_child_weight = c(5,10), # end node size (= n.minobsinnode in gbm)
                        nTrees = c(50,100,150), # max number of trees (= n.trees in gbm)
                        eta = seq(.051, .201, .05), # shrinkage/learning rate (= shrinkage in gbm)
                        optimalTrees = NA,
                        minRMSE = NA)

# ! Ulrichs Grid hat mehr Parameter! müssten in cv angepasst werden! 
ulrichsGrid <- expand.grid(nrounds = 20, # max nTrees
                        max_depth = 2:3, # tree depth
                        eta = c(0.01, 0.001, 0.0001), # shrinkage
                        gamma = c(1, 2, 3), # ?
                        # percent of columns to sample from for each tree
                        colsample_bytree = c(0.4, 0.7, 1.0), 
                        min_child_weight = c(0.5, 1, 1.5), # min node size
                        # percent of training data to sample for each tree
                        subsample = 0.5)


tstart <- Sys.time()
# run tuning grid search 
for (iGline in seq_len(nrow(tuneGrid))) {

  # reproducibility
  set.seed(3829)
  
  # train model
  fit_cv <- xgb.cv(
    data = trainData,                 # X & y in xgboost data format
    nrounds = tuneGrid$nTrees[iGline],# max nTrees
    nfold = setParam$fit$nfolds,      # k-fold cross validation
    eta = tuneGrid$eta[iGline],
    max_depth = tuneGrid$max_depth[iGline],
    min_child_weight = tuneGrid$min_child_weight[iGline],
    # subsample = 1,                  # percent of training data to sample fpr each tree?/bag.fraction in gbm?
    # colsample_bytrees: percent of columns to sample from for each tree
    objective = "reg:squarederror",   # for regression models
    # eval_metric = "rmse",           # default for regression models 
    verbose = FALSE,                  # silent
    # booster = "dart",               # dropout regularisation (takes really long!)
    early_stopping_rounds = 10)       # stop running if the cross validated error 
                                      # does not improve for n continuous trees
  
  # identify the optimal number of trees and the corresponding minimum RMSE for the cv 
  # add number of Trees and corresponding training error to grid
  tuneGrid$optimalTrees[iGline] <- which.min(fit_cv$evaluation_log$test_rmse_mean)
  tuneGrid$minRMSE[iGline] <- min(fit_cv$evaluation_log$test_rmse_mean)
  
}
# tend <- Sys.time()
# difftime(tend, tstart)

# cross-validated tuning parameters
idxOptim <- which.min(tuneGrid$minRMSE)

# to do: save tuned parameters in output# Best parameters
tunedEta <- tuneGrid$eta[idxOptim]
tunedMax_depth <- tuneGrid$max_depth[idxOptim]
tunedMin_child_weight <- tuneGrid$min_child_weight[idxOptim]
tunedNrounds <- tuneGrid$optimalTrees[idxOptim]

# fit model
fit <- xgboost(
# fit <- h2o::h2o.xgboost(
  data = trainData,                 # X & y in xgboost data format
  eta = tuneGrid$eta[idxOptim],
  max_depth = tuneGrid$max_depth[idxOptim],
  min_child_weight = tuneGrid$min_child_weight[idxOptim],
  nrounds = tuneGrid$optimalTrees[idxOptim],
  objective = "reg:squarederror",   # for regression models
  # eval_metric = "rmse",           # default for regression models 
  verbose = FALSE)

tend <- Sys.time()
difftime(tend, tstart)

# performance test for train data
predTrain <- predict(fit, Xtrain)
performTrain <- evalPerformance(predTrain, ytrain)
performTrain <- matrix(performTrain, ncol = length(performTrain), nrow = 1)
colnames(performTrain) <- c("RMSE", "Rsquared", "MAE")

# test data
# performance test for test data
# get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model testing
performTest <- lapply(paste0("test", seq_len(setParam$dgp$nTest)), function(iTest) {
  XtestWithInter <- as.matrix(data[[iTest]][["X_int"]])
  idxPoly <- stringr::str_detect(colnames(XtestWithInter), "^poly")
  idxInter <- stringr::str_detect(colnames(XtestWithInter), ":")
  Xtest <- XtestWithInter[, -c(which(idxPoly), which(idxInter))]
  
  ytest <- data[[iTest]][["yMat"]][,iCond]
  pred <- predict(fit, Xtest)  
  
  evalPerformance(pred, ytest)
})
performTest <- do.call(rbind, performTest)

performTrain
performTest

# are true effects in importance? oder ranks somehow?
# create variable importance measure
# gain = average gain across all splits where feature was used
#   -> importance measures are model specific! 
#   -> use model agnostic permutation variable importance for Enet and GBM instead
# (importance <- xgb.importance(model = fit))
# importance$Feature
# importance$Gain

XtestWithInter <- as.matrix(data[["test1"]][["X_int"]])
idxPoly <- stringr::str_detect(colnames(XtestWithInter), "^poly")
idxInter <- stringr::str_detect(colnames(XtestWithInter), ":")
Xtest <- XtestWithInter[, -c(which(idxPoly), which(idxInter))]
ytest <- data[["test1"]][["yMat"]][,iCond]

# is there a function within xgboost or another package to evaluate the h statistic 
#   based on xgboost models?
# what exactly is the h statistic?
#   -> read paper! (Friedman & Popescu (2005))
#   -> https://christophm.github.io/interpretable-ml-book/interaction.html
#   -> https://christophm.github.io/interpretable-ml-book/pdp.html#pdp

# iml tutorial: https://uc-r.github.io/iml-pkg
# iml package: partial dependence plots and calculation of h-statistic based on 
#   partial dependence data; however: interaction functions are notably slow
# implemented in iml package: permutation-based approach for variable importance, 
#     which is model agnostic
# disadvanategs: 
#   The H-statistic interaction functions do not scale well to wide data (may predictor variables).
#   Only provides permutation-based variable importance scores (which become slow as number of features increase).

# iml does not support xgboost (only caret or randomForest); 
#   1. create data frame with the features
head(Xtrain)
featuresTrain <- data.frame(Xtrain)
#   2. create a vector with the actual responses
responseTrain <- ytrain
#   3. Create custom predict function that returns the predicted values as a vector 
predGBM <- function(model, newdata)  {
  # data format that works with xgboost models
  newData_x = xgb.DMatrix(data.matrix(newdata), missing = NA)
  results <- predict(model, newData_x)
  # results <- predict(model, newdata)
  return(results)
}

predict(fit, Xtest)
predGBM(fit, Xtest)

# calculation of partial dependence scores changes with every iteration
# why? where does randomness come in? 
# random sampling of predictors to sample dependent on density of values 
# increase grid.size to get more stable estimation of partial dependences in interactions
# prädiktorwert grid setzen?
# seed setzen? 
predictor.gbm <- Predictor$new(
  model = fit, 
  data = featuresTrain, 
  y = responseTrain, 
  predict.fun = predGBM,
  type = "prob"
)

# # check if predictor function works
# str(predictor.gbm)
# check <- predictor.gbm$predict(featuresTrain)

# feature importance measure
# increase of the model’s prediction error after permuting the feature
# Fisher et al., 2019; Pargent et al., 2023
#       increase n.repetitions? (How often should the shuffling of the feature 
#       be repeated? The higher the number of repetitions the more stable 
#       and accurate the results become.)
imp.gbm <- FeatureImp$new(predictor.gbm, loss = "mse")
pviRank <- imp.gbm$results$feature
pviValue <- imp.gbm$results$importance
(pIMP <- plot(imp.gbm) + ggtitle("GBM")) # plot this

# partial dependence (visualisation)
#   to compute two way interactions
# interactions with simulated effects
setParam$dgp$interEffects
pd_Var1Var2 <- Partial$new(predictor.gbm, c("Var1", "Var2"))
# partial dependences are identical for reruns (no randomness in calculation of pd values)
# pd_Var1Var2_v2 <- Partial$new(predictor.gbm, c("Var1", "Var2"))
# all(pd_Var1Var2$results$.value == pd_Var1Var2_v2$results$.value)

pd_Var1Var4 <- Partial$new(predictor.gbm, c("Var1", "Var4"))
pd_Var2Var3 <- Partial$new(predictor.gbm, c("Var2", "Var3"))
pd_Var3Var4 <- Partial$new(predictor.gbm, c("Var3", "Var4"))

# interactions without simulated effects
pd_Var11Var14 <- Partial$new(predictor.gbm, c("Var11", "Var14"))
pd_Var1Var14 <- Partial$new(predictor.gbm, c("Var1", "Var14"))

gridExtra::grid.arrange(pd_Var1Var2 %>% plot(), 
                        pd_Var1Var4 %>% plot(),
                        pd_Var2Var3 %>% plot(),
                        pd_Var3Var4 %>% plot(), 
                        pd_Var1Var14 %>% plot(),
                        pd_Var11Var14 %>% plot(), nrow = 3)

# measuring interactions with the h-statistic (Friedman & Popescu, 2005)
#   -> intereaction strength between 0 (no interaction) and 1 (all of variation of 
#   the predicted outcome depends on a given interaction) 
setParam$dgp$interEffects

# does feature interact with any other feature?
# f(x) = estimate predicted values with original model
# pd(x) = partial dependence of variable i
# pd(!x) = partial dependence of all features excluding i
# upper = sum(f(x) - pd(x) - pd(!x))
# lower = variance(f(x))
# rho = upper / lower

# caution!!!
#  -> these results change every time I run the code! why?
#   interactions are calculated based on subsample of the data
# https://github.com/christophM/iml/blob/a72863fe9fa42ac54c0dbfef16f7bcd512b17218/R/Interaction.R#L160C1-L160C77
#     To speed up the computation, we can sample from the n data points. This has 
#     the disadvantage of increasing the variance of the partial dependence estimates, 
#     which makes the H-statistic unstable. So if you are using sampling to reduce 
#     the computational burden, make sure to sample enough data points.
interact.main <- Interaction$new(predictor.gbm, grid.size = 100)
idxInter <- sort(interact.main$results$.interaction, decreasing = T)
featureImportance <- interact.main$results$.feature[match(idxInter, interact.main$results$.interaction)]
featureImportance
importanceMeasure <- interact.main$results$.interaction[match(idxInter, interact.main$results$.interaction)]
importanceMeasure

(p.interactMain <- Interaction$new(predictor.gbm) %>% plot() + ggtitle("GBM"))

#     - evaluate interaction strength with H-statistic
# evaluate two-way interactions:
# Molnar, C. (2020). Interpretable machine learning. Lulu. com.
# https://christophm.github.io/interpretable-ml-book/
# https://github.com/christophM/iml/blob/main/R/Interaction.R
#' The interaction strength between two features is the proportion of the
#' variance of the 2-dimensional partial dependence function that is not
#' explained by the sum of the two 1-dimensional partial dependence functions.
# pd(ij) = interaction partial dependence of variables i and j
# pd(i) = partial dependence of variable i
# pd(j) = partial dependence of variable j
# upper = sum(pd(ij) - pd(i) - pd(j))
# lower = variance(pd(ij))
# rho = upper / lower
# -> partial dependence of the interaction relative to partial dependence of the main effects
#    thus, for variables without simulated effects the overall interaction strength
#       might be as high as for variables with actually simulated interactions only because
#       these variables do have main effects as well 
#   see also Greenwell et al. (2018) and Henninger et al. (2023) 
# -> across variables comparison of overall interaction strength is meaningsless 

# caution!!!
#  -> these results change every time I run the code! why?
iVar1_v1 <- Interaction$new(predictor.gbm, feature = "Var1")
iVar1_v2 <- Interaction$new(predictor.gbm, feature = "Var1")
iVar1_v1_30$results$.interaction - iVar1_v2_30$results$.interaction

# grid size does not resolve the issue! 
iVar1_v1 <- Interaction$new(predictor.gbm, feature = "Var1", grid.size = 100)
iVar1_v2 <- Interaction$new(predictor.gbm, feature = "Var1", grid.size = 100)
iVar1_v1$results$.interaction - iVar1_v2$results$.interaction
all(iVar1_v1$results$.interaction == iVar1_v2$results$.interaction)

interStrength <- lapply(paste0("Var", seq_len(setParam$dgp$p)), function(iVar) {
  tmp <- Interaction$new(predictor.gbm, feature = iVar, grid.size = setParam$fit$nInterStrength)
  idxTmp <- sort(tmp$results$.interaction, decreasing = T)
  cbind(var = rep(iVar, length(idxTmp)),
        feature = tmp$results$.feature[match(idxTmp, tmp$results$.interaction)],
        interaction = tmp$results$.interaction[match(idxTmp, tmp$results$.interaction)])
  
})
interStrength <- do.call(rbind, interStrength)

(interact.Var1 <- Interaction$new(predictor.gbm, feature = "Var1") %>% plot())
(interact.Var2 <- Interaction$new(predictor.gbm, feature = "Var2") %>% plot())
(interact.Var3 <- Interaction$new(predictor.gbm, feature = "Var3") %>% plot())
(interact.Var4 <- Interaction$new(predictor.gbm, feature = "Var4") %>% plot())
(interact.Var9 <- Interaction$new(predictor.gbm, feature = "Var9") %>% plot())

# EIX package features for evaluation of interactions
# https://cran.r-project.org/web/packages/EIX/vignettes/EIX.html
# https://stats.stackexchange.com/questions/496876/interpreting-interaction-metrics-in-xgboost-model
#   sum gains of parent & child node
#   this only works if algorithm does not use a subset of variables at each splitting point! 
#   are these really interactions? 
# all pairs of variable, which occur in the model one above the other
# ! we cannot distinguish if pair of variables are real interaction or not
EIX::interactions(fit, trainData, option = "pairs")
# consider only these pairs of variables, where variable on the bottom (child) 
#   has higher gain than variable on the top (parent)
EIX::interactions(fit, trainData, option = "interactions")

# variable importance measures
EIX::importance(fit, trainData, option = "both")
EIX::importance(fit, trainData, option = "interactions")

# testData <- xgboost::xgb.DMatrix(data = Xtest, label = ytest, missing = NA)
# EIX::waterfall(fit, testData, trainData, option = "interactions")

# interaction evaluation based on gbm models
# maybe choose gbm instead of xgboost if there is h statistic only for gbm?
# https://stats.stackexchange.com/questions/141389/boosted-trees-and-variable-interactions

# # h2o package
# # h statistic for xgboost in h2o package?
# # https://github.com/h2oai/h2o-3/pull/5496
# # https://cran.r-project.org/web/packages/h2o/index.html
# # https://github.com/h2oai/h2o-3
# # this function does not work on xgboost models! 
# # run the entire simulation on h2o or switch to gbm due to h statistic?
# h2o::h2o.h(fit)

################################################################################
# test caret with gbm as estimation method
################################################################################
# evaluate time differences between caret/gbm and xgboost
# -> xgboost is faster than caret/gbm! thus, use xgboost
# try caret code
library(caret)
library(gbm)
# trainWithInter <- cbind(XtrainWithInter, ytrain)
# testWithInter <- cbind(XtestWithInter, ytest)
train <- cbind(Xtrain, ytrain)

# test data
XtestWithInter <- as.matrix(data[["test1"]][["X_int"]])
# remove interaction to test if GBM is able to find interactions on its own
idxPoly <- stringr::str_detect(colnames(XtestWithInter), "^poly")
idxInter <- stringr::str_detect(colnames(XtestWithInter), ":")
Xtest <- XtestWithInter[, -c(which(idxPoly), which(idxInter))]

ytest <- data[["test1"]][["yMat"]][,iCond]
test <- cbind(Xtest, ytest)

# # with default tuning grid
# model <-  train(ytrain ~ .,
#                 # train,
#                 trainWithInter,
#                 method = "gbm", #gradient boosting machines
#                 metric = "RMSE",
#                 preProc = "scale",
#                 trControl = trainControl(method="cv", number=10),
#                 verbose = FALSE)


# with our own tuning grid
# define tuning grid
smallerGrid <- expand.grid(interaction.depth = c(1,2,3),
                           n.minobsinnode = c(5, 10),        # end node size
                           n.trees = c(50,100,150),
                           shrinkage = seq(.051, .201, .05))

tstart <- Sys.time()
model <- train(ytrain ~ .,
               train,
               method = "gbm", #gradient boosting machines
               metric = "RMSE",
               trControl = trainControl(method="cv", number=10), #model training using ten-fold cross-validation
               tuneGrid = smallerGrid, #setting the tuning grid
               verbose = FALSE)
tend <- Sys.time()
difftime(tend, tstart)

#returns explained variance (R?), Root Mean Squared Error (RMSE), and Mean Absolute Error (MAE) for model training:
predictions <- predict(model, train) 
eval_train <- caret::postResample(pred = predictions, obs = train[,"ytrain"])  
# predictions <- predict(model, trainWithInter) 
# eval_train <- caret::postResample(pred = predictions, obs = trainWithInter[,"ytrain"])  

#returns explained variance (R?), Root Mean Squared Error (RMSE), and Mean Absolute Error (MAE) for model testing:
predictions <- predict(model, test) 
eval_test <- caret::postResample(pred = predictions, obs = test[,"ytest"])
# predictions <- predict(model, testWithInter)
# eval_test <- caret::postResample(pred = predictions, obs = testWithInter[,"ytest"])

eval_train
eval_test

# variable importance (only linear predictors)
vi <- t(varImp(model)$importance)  

# calculate interaction importance via h statistic
# https://rdrr.io/cran/gbm/man/interact.gbm.html
# look into source code and copy calculations for xgboost models?
# https://rdrr.io/cran/gbm/src/R/interact.gbm.R
interact.gbm(model, train, i.var = 1:14)
