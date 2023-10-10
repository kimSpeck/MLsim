# open all files and save data in one single rda-file
# full Data (list with nested conditions)
# full data list
#   N.{100, 300, 1000}xpTrash.{10, 50, 100}
#   {lambda, phi, eta, target, conv1step, conv2step, corLambda, cohensD, corDeriv, nFactorsPA}

# parameter & general functions
source("setParameters.R") # import parameter values

# save results to this folder
resFolder <- paste0("results/resultsServer")

condGrid <- expand.grid(N = setParam$dgp$N, 
                        pTrash = setParam$dgp$pTrash)

# read in data
fullData <- vector(mode = "list", length = nrow(condGrid))
for (iSim in seq_len(nrow(condGrid))) {
  resFileName <- paste0(resFolder, "/", "resultsN", condGrid[iSim, "N"], "_pTrash", condGrid[iSim, "pTrash"], ".rds")
  fullData[[iSim]] <- readRDS(resFileName)
}
names(fullData) <- paste0("N", condGrid$N, "_pTrash", condGrid$pTrash)

fullDataFile <- paste0(resFolder, "/fullData.rda")
save(fullData, file = fullDataFile)
