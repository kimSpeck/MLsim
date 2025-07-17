# full Data (list with nested conditions)
# full data list
#   now we have multiple DGPs, too {inter, nonlinear3, pwlinear} they live in different folders
#   model.{ENETw, ENETwo, GBM, RF} x N.{100, 300, 1000} x pTrash.{10, 50, 100} x rel{0.6, 0.8, 1}

# open all files and save data in one single rda-file but separate rda-files for each model type
# -> with all results included, these files would be too big (~ 5.8 GB for only one model) 
# -> reading in data and subsequently selecting variables would be really memory inefficient 

# parameter & general functions
source("utils/setParameters.R") # import parameter values

# save results to this folder
resFolder <- paste0("results")

# all results from the following combinations
dgpVec = c("inter", "pwlinear", "nonlinear3")
modelVec = c("ENETw", "ENETwo", "GBM", "RF")
condGrid <- expand.grid(data = dgpVec,
                        model = modelVec,
                        N = setParam$dgp$N, 
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability)

################################################################################
# extract and save performPerSample from all models for ANOVA
################################################################################

for (iDGP in dgpVec) {
  # generate folder to save new rda files
  depMeasures = paste0(resFolder, "/", iDGP, "/dependentMeasures")
  if (!file.exists(depMeasures)){
    dir.create(depMeasures)
  }
  
  for (iModel in modelVec) {
    subGrid <- condGrid[which(condGrid$data == iDGP & 
                                condGrid$model == iModel),]
    performPerSample <- vector(mode = "list", length = nrow(subGrid))
    for (iSim in seq_len(nrow(subGrid))) {
      resFileName <- paste0(resFolder, "/", iDGP, "/",  
                            "resultsModel", iModel, 
                            "_N", subGrid[iSim, "N"], 
                            "_pTrash", subGrid[iSim, "pTrash"], 
                            "_rel", subGrid[iSim, "reliability"], ".rds")
      
      print(resFileName)
      tmp <- readRDS(resFileName)[["performPerSample"]]
      
      tmp <- cbind(tmp, 
                   # add sample number as an additional variable
                   sample = rep(seq_len(setParam$dgp$nTrain),
                                dim(tmp)[1]/setParam$dgp$nTrain),
                   # add model variable (within factor in mixed ANOVA)
                   model = rep(iModel, dim(tmp)[1]),
                   # add dgp
                   dgp = rep(iDGP, dim(tmp)[1]),
                   # add N_pTrash_rel label
                   N_pTrash = rep(paste0("N", subGrid[iSim, "N"], 
                                         "_pTrash", subGrid[iSim, "pTrash"], 
                                         "_rel", subGrid[iSim, "reliability"]), 
                                  dim(tmp)[1]))
      
      performPerSample[[iSim]] <- tmp
    }
    
    names(performPerSample) <- paste0("N", subGrid$N, 
                                      "_pTrash", subGrid$pTrash,
                                      "_rel", subGrid$reliability)
    
    ppsFile <- paste0(depMeasures, "/performPerSample_", iDGP, "_", iModel, ".rda")
    save(performPerSample, file = ppsFile)
    print("done")
    gc()
  }
} 

################################################################################
# extract and save performTrainStats and performTestStats for basic plots
################################################################################

for (iDGP in dgpVec) {
  # generate folder to save new rda files
  depMeasures = paste0(resFolder, "/", iDGP, "/dependentMeasures")
  if (!file.exists(depMeasures)){
    dir.create(depMeasures)
  }
  
  for (iModel in modelVec) {
    subGrid <- condGrid[which(condGrid$data == iDGP & 
                                condGrid$model == iModel),]
    performTrain <- vector(mode = "list", length = nrow(subGrid))
    performTest <- vector(mode = "list", length = nrow(subGrid))
    for (iSim in seq_len(nrow(subGrid))) {
      resFileName <- paste0(resFolder, "/", iDGP, "/", "resultsModel", iModel, 
                            "_N", subGrid[iSim, "N"], 
                            "_pTrash", subGrid[iSim, "pTrash"], 
                            "_rel", subGrid[iSim, "reliability"], ".rds")
      # print(resFileName)
      tmpTrain <- readRDS(resFileName)[["performTrainStats"]]
      tmpTest <- readRDS(resFileName)[["performTestStats"]]
      
      tmpTrain <- cbind(tmpTrain, 
                        # add model variable (within factor in mixed ANOVA)
                        model = rep(iModel, dim(tmpTrain)[1]),
                        # add dgp
                        dgp = rep(iDGP, dim(tmpTrain)[1]),
                        # add N_pTrash_rel label
                        N_pTrash = rep(paste0("N", subGrid[iSim, "N"], 
                                              "_pTrash", subGrid[iSim, "pTrash"], 
                                              "_rel", subGrid[iSim, "reliability"]), 
                                       dim(tmpTrain)[1]))
      tmpTest <- cbind(tmpTest, 
                       # add model variable (within factor in mixed ANOVA)
                       model = rep(iModel, dim(tmpTest)[1]),
                       # add dgp
                       dgp = rep(iDGP, dim(tmpTest)[1]),
                       # add N_pTrash_rel label
                       N_pTrash = rep(paste0("N", subGrid[iSim, "N"], 
                                             "_pTrash", subGrid[iSim, "pTrash"], 
                                             "_rel", subGrid[iSim, "reliability"]), 
                                      dim(tmpTest)[1]))
      
      performTrain[[iSim]] <- tmpTrain
      performTest[[iSim]] <- tmpTest
    }
    
    names(performTrain) <- paste0("N", subGrid$N, 
                                  "_pTrash", subGrid$pTrash,
                                  "_rel", subGrid$reliability)
    names(performTest) <- paste0("N", subGrid$N, 
                                 "_pTrash", subGrid$pTrash,
                                 "_rel", subGrid$reliability)
    
    trainFile <- paste0(depMeasures, "/performTrainStats_", iDGP, "_", iModel, ".rda")
    testFile <- paste0(depMeasures, "/performTestStats_", iDGP, "_", iModel, ".rda")
    save(performTrain, file = trainFile)
    save(performTest, file = testFile)
    print("done")
    gc()
  } 
}

################################################################################
# extract hyperparameters (for all models)
################################################################################
# extract lambda and alpha to check where the weird results pattern in specificity 
#    values of the ENETwo comes from
# extract shrinkage, max_depth, min_child_weight and Nrounds to check where the 
#   weird results pattern in specificity values of the GBM comes from

for (iDGP in dgpVec) {
  # generate folder to save new rda files
  depMeasures = paste0(resFolder, "/", iDGP, "/dependentMeasures")
  if (!file.exists(depMeasures)){
    dir.create(depMeasures)
  }
  
  for (iModel in modelVec) {
    subGrid <- condGrid[which(condGrid$data == iDGP & 
                                condGrid$model == iModel),]
    hyper <- vector(mode = "list", length = nrow(subGrid))
    for (iSim in seq_len(nrow(subGrid))) {
      resFileName <- paste0(resFolder, "/", iDGP, "/", "resultsModel", iModel, 
                            "_N", subGrid[iSim, "N"], 
                            "_pTrash", subGrid[iSim, "pTrash"], 
                            "_rel", subGrid[iSim, "reliability"], ".rds")
      
      print(resFileName)
      tmp <- readRDS(resFileName)[["selectionPerSample"]]
      
      tmp <- cbind(tmp, 
                   sample = rep(seq_len(setParam$dgp$nTrain),
                                dim(tmp)[1]/setParam$dgp$nTrain),
                   # add model variable (within factor in mixed ANOVA)
                   model = rep(iModel, dim(tmp)[1]),
                   # add dgp
                   dgp = rep(iDGP, dim(tmp)[1]),
                   # add N_pTrash_rel label
                   N_pTrash = rep(paste0("N", subGrid[iSim, "N"], 
                                         "_pTrash", subGrid[iSim, "pTrash"], 
                                         "_rel", subGrid[iSim, "reliability"]), 
                                  dim(tmp)[1]))
      
      hyper[[iSim]] <- tmp
    }
    
    names(hyper) <- paste0("N", subGrid$N, 
                           "_pTrash", subGrid$pTrash,
                           "_rel", subGrid$reliability)
    
    hyperFile <- paste0(depMeasures, "/hyperParametersSample_", iDGP, "_", iModel, ".rda")
    save(hyper, file = hyperFile)
    print("done")
    gc()
  }
}
