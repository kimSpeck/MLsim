# load parameters & custom functions 
source("utils/setParameters.R") # parameter values
source("utils/simTools.R") # functions for data simulation

# generate folder for data files
dataFolder <- "data"
createFolder(dataFolder)

# grid to simulate data with mapply later
# simulation via mapply to easily simulate subsets of parameter combinations
# create factor structure for predictor variables to increase reliability
gridFull <- expand.grid(N = setParam$dgp$N,
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability)

for (iData in seq_len(dim(gridFull)[1])) {
  N <- gridFull[iData, "N"] 
  pTrash <- gridFull[iData, "pTrash"]
  reliability <- gridFull[iData, "reliability"]
  
  # create Folder for single sample files of respective simulated condition (N x pTrash x reliability)
  sampleFolder <- paste0("/simDataN", N, "_pTrash", pTrash, "_rel", reliability)
  createFolder(paste0(dataFolder, "/inter", sampleFolder))  
  
  # load data for respective simulated condition (N x pTrash x reliability)
  load(paste0(dataFolder, "/inter", sampleFolder, ".rda"))
  
  # split list of data in single sample files 
  for (iList in names(data)) {
    fileName <- paste0(dataFolder, "/inter", sampleFolder, paste0("/sample_", iList, ".rda"))  
    
    dataList <- data[[iList]]
    save(dataList, file = fileName)
  }
  
}

