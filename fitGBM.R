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
# library(EIX)

dataFolder <- "data"
# condGrid <- expand.grid(N = setParam$dgp$N, 
#                         pTrash = setParam$dgp$pTrash)
# iSim = 1
# Probedatensatz mit 100 Menschen & 10 trash Variablen

# fileName <- paste0("simDataN", condGrid[iSim, "N"], "_pTrash", condGrid[iSim, "pTrash"], ".rda")
# load(paste0(dataFolder, "/", fileName))

load("data/simDataN1000_pTrash10_rel1_f0.rda")

iCond <- 1
iSample <- 1

# data from single sample (with interaction as specified in elastic net)
XtrainWithInter <- as.matrix(data[[iSample]][["X_int"]])
# remove interaction to test if GBM is able to find interactions on its own
idxPoly <- stringr::str_detect(colnames(XtrainWithInter), "^poly")
idxInter <- stringr::str_detect(colnames(XtrainWithInter), ":")
Xtrain <- XtrainWithInter[, -c(which(idxPoly), which(idxInter))]

ytrain <- data[[iSample]][["yMat"]][,iCond]

# # test data 
# XtestWithInter <- as.matrix(data[["test1"]][["X_int"]])
# # remove interaction to test if GBM is able to find interactions on its own
# idxPoly <- stringr::str_detect(colnames(XtestWithInter), "^poly")
# idxInter <- stringr::str_detect(colnames(XtestWithInter), ":")
# Xtest <- XtestWithInter[, -c(which(idxPoly), which(idxInter))]
# 
# ytest <- data[["test1"]][["yMat"]][,iCond]
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

# are true effects in importance? oder ranks somehow?
# create importance matrix
(importance <- xgb.importance(model = fit))

# is there a function within xgboost or another package to evaluate the h statistic 
#   based on xgboost models?
# what exactly is the h statistic?
# read paper! (Friedman & Popescu (2005))
# https://christophm.github.io/interpretable-ml-book/interaction.html
# https://christophm.github.io/interpretable-ml-book/pdp.html#pdp

# h statistic for xgboost in h2o package?
# https://github.com/h2oai/h2o-3/pull/5496
# https://cran.r-project.org/web/packages/h2o/index.html
# https://github.com/h2oai/h2o-3

# interaction evaluation based on gbm models
# maybe choose gbm instead of xgboost if there is h statistic only for gbm?
# https://stats.stackexchange.com/questions/141389/boosted-trees-and-variable-interactions

# EIX package features for evaluation of interactions
# https://cran.r-project.org/web/packages/EIX/vignettes/EIX.html
# https://stats.stackexchange.com/questions/496876/interpreting-interaction-metrics-in-xgboost-model
#   sum gains of parent & child node
#   are those really interactions? 
# all pairs of variable, which occur in the model one above the other
# ! we cannot distinguish if pair of variables are real interaction or not
EIX::interactions(fit, trainData, option = "pairs")
# consider only these pairs of variables, where variable on the bottom (child) 
#   has higher gain than variable on the top (parent)
EIX::interactions(fit, trainData, option = "interactions")

# variable importance measures
EIX::importance(fit, trainData, option = "both")
EIX::importance(fit, trainData, option = "interactions")

XtestWithInter <- as.matrix(data[["test1"]][["X_int"]])
idxPoly <- stringr::str_detect(colnames(XtestWithInter), "^poly")
idxInter <- stringr::str_detect(colnames(XtestWithInter), ":")
Xtest <- XtestWithInter[, -c(which(idxPoly), which(idxInter))]
ytest <- data[["test1"]][["yMat"]][,iCond]

testData <- xgboost::xgb.DMatrix(data = Xtest, label = ytest, missing = NA)

EIX::waterfall(fit, testData, trainData, option = "interactions")
################################################################################
# evaluate time differences between caret/gbm and xgboost
# -> xgboost is faster than caret/gbm! thus, use xgboost
# try caret code
library(caret)
library(gbm)
# trainWithInter <- cbind(XtrainWithInter, ytrain)
# testWithInter <- cbind(XtestWithInter, ytest)
train <- cbind(Xtrain, ytrain)
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

# variable importance (only linear predictors)
vi <- t(varImp(model)$importance)  

# calculate interaction importance via h statistic
# https://rdrr.io/cran/gbm/man/interact.gbm.html
# look into source code and copy calculations for xgboost models?
# https://rdrr.io/cran/gbm/src/R/interact.gbm.R
interact.gbm(model, train)
