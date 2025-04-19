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
source("utils/samplePiecewiseLinearData.R") # linear + nonlinear effects

timeStampFolder <- format(Sys.time(), "%d%m%y_%H%M%S")

# generate folder for data files
dataFolder <- "data"
createFolder(dataFolder)

dgpFolder <- "/inter"
createFolder(paste0(dataFolder, dgpFolder))
dgpFolder <- "/nonlinear"
createFolder(paste0(dataFolder, dgpFolder))
dgpFolder <- "/pwlinear"
createFolder(paste0(dataFolder, dgpFolder))

# generate folder for log files
logFolder = "log"
createFolder(logFolder)

# set number of cores to use in parallel computing
nCoresSampling <- 30 # simulate data

# grid to simulate data with mapply later
# simulation via mapply to easily simulate subsets of parameter combinations
# create factor structure for predictor variables to increase reliability
gridInter <- expand.grid(N = setParam$dgp$N,
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability)
# ! factors are interpreted as level numbers!; only character variables are interpreted by their name!
str(gridInter) 

# create seed number for parallel cluster (reproducibility of generated data)
#   keep seeds from simulating only linear effects and only add nonlinear conditions with seed
#   -> two step grid generating
set.seed(8967369)
seedNum <- sample(1:999999, dim(gridInter)[1], replace = FALSE) 
gridInter$sampleSeed <- seedNum[1:dim(gridInter)[1]]

# add dgp type column to the grid
gridNL <- cbind(data = "inter", gridInter)   

# add nonlinear grid
set.seed(4890920)
seedNum <- sample(1:999999, dim(gridInter)[1], replace = FALSE) 

gridNL <- rbind(gridNL, 
                cbind(data = "nonlinear", 
                      gridInter[,!colnames(gridInter) %in% "sampleSeed"], 
                      sampleSeed = seedNum))

# add piecewise linear grid
set.seed(7485936)
seedNum <- sample(1:999999, dim(gridInter)[1], replace = FALSE) 

gridFull <- rbind(gridNL, 
                  cbind(data = "pwlinear", 
                        gridInter[,!colnames(gridInter) %in% "sampleSeed"], 
                        sampleSeed = seedNum))


# sample data in parallel
createData <- function(data, N, pTrash, reliability, sampleSeed){
  
  if (data == "inter"){
    environment(sampleInteractionData) <- environment()  
  } else if (data == "nonlinear") {
    environment(sampleNonlinearData) <- environment()  
  } else if (data == "pwlinear") {
    environment(samplePiecewiseLinearData) <- environment()  
  } else {
    stop("We can only simulate inter, nonlinear or piecewise linear data!")
  }
  
  # Initiate cluster; type = "FORK" only on Linux/MacOS: contains all environment variables automatically
  cl <- makeCluster(nCoresSampling, type = "FORK",
                    outfile = paste0(logFolder, "/", "simulateDataStatus_", data, 
                                     timeStampFolder, ".txt"))
  
  # set seed that works for parallel processing
  set.seed(sampleSeed)
  s <- .Random.seed
  clusterSetRNGStream(cl = cl, iseed = s)
  
  if (data == "inter"){
    sampleInteractionData() # run function to actually create dataset
  } else if (data == "nonlinear") {
    sampleNonlinearData()
  } else if (data == "pwlinear") {
    samplePiecewiseLinearData()
  }
  
  # close cluster to return resources (memory) back to OS
  stopCluster(cl)
}

# # test simulating data for only one condition
# out <- do.call(mapply, c(FUN = createData, gridFull[1,])) # check inter
# out <- do.call(mapply, c(FUN = createData, gridFull[19,])) # check nonlinear

# # simulate only nonlinear data
# out <- do.call(mapply, c(FUN = createData, gridFull[gridFull$data == "nonlinear", ]))

# # simulate full data
# out <- do.call(mapply, c(FUN = createData, gridFull))

