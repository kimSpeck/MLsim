# full Data (list with nested conditions)
# full data list
#   model.{ENETw, ENETwo, GBM} x N.{100, 300, 1000} x pTrash.{10, 50, 100} x rel{0.6, 0.8, 1}

# open all files and save data in one single rda-file but separate rda-files for each model type
# -> with all results included, these files would be too big (~ 5.8 GB for only one model) 
# -> reading in data and subsequently selecting variables would be really memory inefficient 

# parameter & general functions
source("setParameters.R") # import parameter values

# save results to this folder
resFolder <- paste0("results/finalResults")

# generate folder to save new rda files
depMeasures = paste0(resFolder, "/dependentMeasures")
if (!file.exists(depMeasures)){
  dir.create(depMeasures)
}


modelVec = c("ENETw", "ENETwo", "GBM")
condGrid <- expand.grid(model = modelVec,
                        N = setParam$dgp$N, 
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability)

################################################################################
# extract and save performPerSample from all models for ANOVA
################################################################################
for (iModel in modelVec) {
  subGrid <- condGrid[which(condGrid$model == iModel),]
  performPerSample <- vector(mode = "list", length = nrow(subGrid))
  for (iSim in seq_len(nrow(subGrid))) {
    resFileName <- paste0(resFolder, "/", "resultsModel", iModel, 
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
  
  ppsFile <- paste0(depMeasures, "/performPerSample_", iModel, ".rda")
  save(performPerSample, file = ppsFile)
  print("done")
  gc()
} 

################################################################################
# extract and save performTrainStats and performTestStats for basic plots
################################################################################

for (iModel in modelVec) {
  subGrid <- condGrid[which(condGrid$model == iModel),]
  performTrain <- vector(mode = "list", length = nrow(subGrid))
  performTest <- vector(mode = "list", length = nrow(subGrid))
  for (iSim in seq_len(nrow(subGrid))) {
    resFileName <- paste0(resFolder, "/", "resultsModel", iModel, 
                          "_N", subGrid[iSim, "N"], 
                          "_pTrash", subGrid[iSim, "pTrash"], 
                          "_rel", subGrid[iSim, "reliability"], ".rds")
    # print(resFileName)
    tmpTrain <- readRDS(resFileName)[["performTrainStats"]]
    tmpTest <- readRDS(resFileName)[["performTestStats"]]
    
    tmpTrain <- cbind(tmpTrain, 
                 # add model variable (within factor in mixed ANOVA)
                 model = rep(iModel, dim(tmpTrain)[1]),
                 # add N_pTrash_rel label
                 N_pTrash = rep(paste0("N", subGrid[iSim, "N"], 
                                       "_pTrash", subGrid[iSim, "pTrash"], 
                                       "_rel", subGrid[iSim, "reliability"]), 
                                dim(tmpTrain)[1]))
    tmpTest <- cbind(tmpTest, 
                      # add model variable (within factor in mixed ANOVA)
                      model = rep(iModel, dim(tmpTest)[1]),
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
  
  trainFile <- paste0(depMeasures, "/performTrainStats_", iModel, ".rda")
  testFile <- paste0(depMeasures, "/performTestStats_", iModel, ".rda")
  save(performTrain, file = trainFile)
  save(performTest, file = testFile)
  print("done")
  gc()
} 

################################################################################
# extract and save pvi
################################################################################
for (iModel in modelVec) {
  subGrid <- condGrid[which(condGrid$model == iModel),]
  pvi <- vector(mode = "list", length = nrow(subGrid))
  for (iSim in seq_len(nrow(subGrid))) {
    resFileName <- paste0(resFolder, "/", "resultsModel", iModel, 
                          "_N", subGrid[iSim, "N"], 
                          "_pTrash", subGrid[iSim, "pTrash"], 
                          "_rel", subGrid[iSim, "reliability"], ".rds")
    print(resFileName)
    tmp <- data.frame(readRDS(resFileName)[["pvi"]])
    
    # remove entries with pvi value of 1 (pvi = 1: kein Einfluss)
    tmp <- tmp[tmp$pviValue != "1",]
    
    tmp <- cbind(tmp, 
                 # add model variable (within factor in mixed ANOVA)
                 model = rep(iModel, dim(tmp)[1]),
                 # add N_pTrash_rel label
                 N_pTrash = rep(paste0("N", subGrid[iSim, "N"], 
                                       "_pTrash", subGrid[iSim, "pTrash"], 
                                       "_rel", subGrid[iSim, "reliability"]), 
                                dim(tmp)[1]))
    
    pvi[[iSim]] <- tmp
  }
  
  names(pvi) <- paste0("N", subGrid$N, 
                       "_pTrash", subGrid$pTrash,
                       "_rel", subGrid$reliability)
  
  pviFile <- paste0(depMeasures, "/pvi_", iModel, ".rda")
  save(pvi, file = pviFile)
  print("done")
  gc()
} 

################################################################################
# extract interStrength (for GBM)
################################################################################
# this dependent measure only exists for GBM
subGrid <- condGrid[which(condGrid$model == "GBM"),]
interStrength <- vector(mode = "list", length = nrow(subGrid))
for (iSim in seq_len(nrow(subGrid))) {
  resFileName <- paste0(resFolder, "/", "resultsModel", "GBM", 
                        "_N", subGrid[iSim, "N"], 
                        "_pTrash", subGrid[iSim, "pTrash"], 
                        "_rel", subGrid[iSim, "reliability"], ".rds")
  print(resFileName)
  tmp <- data.frame(readRDS(resFileName)[["interStrength"]])
    
  # remove entries with pvi value of 1 (pvi = 1: kein Einfluss)
  tmp <- tmp[tmp$interaction != "0",]
    
  tmp <- cbind(tmp, 
               # add N_pTrash_rel label
                N_pTrash = rep(paste0("N", subGrid[iSim, "N"], 
                                      "_pTrash", subGrid[iSim, "pTrash"], 
                                      "_rel", subGrid[iSim, "reliability"]), 
                               dim(tmp)[1]))
    
  interStrength[[iSim]] <- tmp
}
  
names(interStrength) <- paste0("N", subGrid$N, 
                               "_pTrash", subGrid$pTrash,
                               "_rel", subGrid$reliability)

interStrengthFile <- paste0(depMeasures, "/interStrength_", "GBM", ".rda")
save(interStrength, file = interStrengthFile)

################################################################################
# extract beta coefficients (for ENET)
################################################################################
# extract beta coefficients for both ENET with and without interactions
modelVecSub <- c("ENETw", "ENETwo")
for (iModel in modelVecSub) {
  subGrid <- condGrid[which(condGrid$model == iModel),]
  estBeta <- vector(mode = "list", length = nrow(subGrid))
  for (iSim in seq_len(nrow(subGrid))) {
    resFileName <- paste0(resFolder, "/", "resultsModel", iModel, 
                          "_N", subGrid[iSim, "N"], 
                          "_pTrash", subGrid[iSim, "pTrash"], 
                          "_rel", subGrid[iSim, "reliability"], ".rds")
    
    print(resFileName)
    tmp <- readRDS(resFileName)[["estBetaFull"]]
    
    if (iModel == "ENETwo") {
      idxCondLabels <- rep(seq_along(setParam$dgp$condLabels), 
          each = length(setParam$dgp$linEffects) + subGrid[iSim, "pTrash"])  
    } else if (iModel == "ENETw") {
      tmpP = length(setParam$dgp$linEffects) + subGrid[iSim, "pTrash"]
      nrInter = choose(tmpP, setParam$dgp$interDepth)
      idxCondLabels <- rep(seq_along(setParam$dgp$condLabels), 
                           each = tmpP + nrInter)
    }
    
    
    tmp <- cbind(tmp, 
                 idxCondLabel = idxCondLabels,
                 # add model variable (within factor in mixed ANOVA)
                 model = rep(iModel, dim(tmp)[1]),
                 # add N_pTrash_rel label
                 N_pTrash = rep(paste0("N", subGrid[iSim, "N"], 
                                       "_pTrash", subGrid[iSim, "pTrash"], 
                                       "_rel", subGrid[iSim, "reliability"]), 
                                dim(tmp)[1]))
  
    estBeta[[iSim]] <- tmp
  }
  
  names(estBeta) <- paste0("N", subGrid$N, 
                           "_pTrash", subGrid$pTrash,
                           "_rel", subGrid$reliability)
  
  estBetaFile <- paste0(depMeasures, "/estBetaSample_", iModel, ".rda")
  save(estBeta, file = estBetaFile)
  print("done")
  gc()
}

################################################################################
# full data to one rda file
################################################################################
# # this creates huge rda files that would take a lot of time to read in
# # therefore this version of data management does not work for the final sample
# # read in data
# for (iModel in modelVec) {
#   subGrid <- condGrid[which(condGrid$model == iModel),]
#   fullData <- vector(mode = "list", length = nrow(subGrid))
#   for (iSim in seq_len(nrow(subGrid))) {
#     resFileName <- paste0(resFolder, "/", "resultsModel", iModel, 
#                           "_N", subGrid[iSim, "N"], 
#                           "_pTrash", subGrid[iSim, "pTrash"], 
#                           "_rel", subGrid[iSim, "reliability"], ".rds")
#     print(resFileName)
#     fullData[[iSim]] <- readRDS(resFileName)
#   }
#   
#   names(fullData) <- paste0("N", subGrid$N, 
#                             "_pTrash", subGrid$pTrash,
#                             "_rel", subGrid$reliability)
#   
#   fullDataFile <- paste0(resFolder, "/fullData_", iModel, ".rda")
#   save(fullData, file = fullDataFile)
#   print("done")
#   gc()
# } 

