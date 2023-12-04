# open all files and save data in one single rda-file
# full Data (list with nested conditions)
# full data list
#   N.{100, 300, 1000}xpTrash.{10, 50, 100}xrel{0.6, 0.8, 1}
#   {lambda, phi, eta, target, conv1step, conv2step, corLambda, cohensD, corDeriv, nFactorsPA}

# parameter & general functions
source("setParameters.R") # import parameter values

# save results to this folder
# resFolder <- paste0("results/resultsFactorsWithoutInteraction")
resFolder <- paste0("results/resultsFactorsWithInteraction")

condGrid <- expand.grid(N = setParam$dgp$N, 
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability,
                        factors = c(TRUE, FALSE))

# read in data
fullData <- vector(mode = "list", length = nrow(condGrid))
for (iSim in seq_len(nrow(condGrid))) {
  resFileName <- paste0(resFolder, "/", "resultsN", condGrid[iSim, "N"], 
                        "_pTrash", condGrid[iSim, "pTrash"], 
                        "_rel", condGrid[iSim, "reliability"], 
                        "_f", ifelse(condGrid[iSim,"factors"], 1, 0), ".rds")
  fullData[[iSim]] <- readRDS(resFileName)
}
names(fullData) <- paste0("N", condGrid$N, 
                          "_pTrash", condGrid$pTrash,
                          "_rel", condGrid$reliability,
                          "_f", condGrid$factors)

fullDataFile <- paste0(resFolder, "/fullData.rda")
save(fullData, file = fullDataFile)
