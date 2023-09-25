setParam <- list()

# data generating process
setParam$dgp$nSamples <- 100

# setParam$dgp$N <- c(100, 300, 1000, 10000) # number of observations
setParam$dgp$N <- c(100, 300, 1000) # number of observations
setParam$dgp$p <- 4
setParam$dgp$pTrash <- c(10, 50, 100) # number of "trash" predictors

setParam$dgp$interDepth <- c(2) # depth of interactions (so far: only two-way interaction)
setParam$dgp$poly <- c(2) # degree of polynomials (so far: only quadratic effects)

# proportion of effect explained by linear effects vs. interaction
setParam$dgp$percentLinear <- c(0.5, 0.8, 0.2) 
setParam$dgp$percentInter <- c(0.5, 0.2, 0.8)
setParam$dgp$percentPoly <- c(0, 0, 0)

# check
if (!all(setParam$dgp$percentLinear+
         setParam$dgp$percentInter +
         setParam$dgp$percentPoly == 1)) stop("Proportion of different kinds of effects do not sum up to 1! Check values!")

# r squared
setParam$dgp$Rsquared <- c(.10, .30, .50, .80)

comboGrid <- expand.grid(setParam$dgp$Rsquared, 
                         paste(setParam$dgp$percentLinear, setParam$dgp$percentInter, sep = "_"))
setParam$dgp$condLabels <- sapply(seq_len(length(setParam$dgp$Rsquared) * length(setParam$dgp$percentLinear)), 
                                  function(x) paste0("R2", comboGrid$Var1[x], "lin_inter", comboGrid$Var2[x]))
rm(comboGrid)
# average correlations of predictors and their SD 
# ToDo: manipulate correlations between different variables (no separate simulated conditions)
setParam$dgp$meanR <- 0
setParam$dgp$sdR <- 0.1
# setParam$dgp$sdR <- 0.2 # original value

# error "variance" (= standard deviation)
setParam$dgp$sigmaE <- 1

# parameter for gbm or model fitting more general