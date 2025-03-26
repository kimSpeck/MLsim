# simulate data

# load packages
library(mvtnorm)
library(truncnorm)
library(parallel)

# load parameters & custom functions 
source("utils/setParameters.R") # parameter values
source("utils/simTools.R") # general functions for data simulation

# functions to generate data
source("utils/sampleInteractionData.R") # linear + interaction effects
source("utils/sampleNonlinearData.R") # linear + nonlinear effects

timeStampFolder <- format(Sys.time(), "%d%m%y_%H%M%S")

# generate folder for data files
dataFolder <- "data"
createFolder(dataFolder)

# generate folder for log files
logFolder = "log"
createFolder(logFolder)

# set number of cores to use in parallel computing
nCoresSampling <- 5 # simulate data

# grid to simulate data with mapply later
# simulation via mapply to easily simulate subsets of parameter combinations
# create factor structure for predictor variables to increase reliability
gridFull <- expand.grid(N = setParam$dgp$N,
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability)
# ! factors are interpreted as level numbers!; only character variables are interpreted by their name!
str(gridFull) 

# create seed number for parallel cluster (reproducibility of generated data)
set.seed(8967369)
seedNum <- sample(1:999999, dim(gridFull)[1], replace = FALSE) 
gridFull$sampleSeed <- seedNum[1:dim(gridFull)[1]]

# sample data in parallel
createData <- function(N, pTrash, reliability, sampleSeed){
  
  environment(sampleInteractionData) <- environment()
  
  # Initiate cluster; type = "FORK" only on Linux/MacOS: contains all environment variables automatically
  cl <- makeCluster(nCoresSampling, type = "FORK",
                    outfile = paste0(logFolder, "/", "simulateDataStatus",
                                     timeStampFolder, ".txt"))
  
  # set seed that works for parallel processing
  set.seed(sampleSeed)
  s <- .Random.seed
  clusterSetRNGStream(cl = cl, iseed = s)
  
  sampleInteractionData() # run function to actually create dataset
  
  # close cluster to return resources (memory) back to OS
  stopCluster(cl)
}

# # test it
# out <- do.call(mapply, c(FUN = createData, gridFull[1,]))

out <- do.call(mapply, c(FUN = createData, gridFull))

