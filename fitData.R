# Model fitting

library(glmnet)

# load parameters & custom functions 
source("setParameters.R") # parameter values
 
# N <- setParam$dgp$N[1]
# pTrash <- setParam$dgp$pTrash[1]
dataFolder <- "data"
resFolder <- "results"
if (!file.exists(resFolder)){
  dir.create(resFolder)
}


condGrid <- expand.grid(N = setParam$dgp$N, 
                        pTrash = setParam$dgp$pTrash)

res_estB <- lapply(seq_len(nrow(condGrid)), function(iSim) {
  
  fileName <- paste0("simDataN", condGrid[iSim, "N"], "_pTrash", condGrid[iSim, "pTrash"], ".rda")
  load(paste0(dataFolder, "/", fileName))
  
  # fitte eine regularisierte Regression
  tmp_estB <- lapply(seq_along(setParam$dgp$condLabels), function(iCond) {
    
    estB <- lapply(seq_len(setParam$dgp$nSamples), function(iSample) {
      fit_cv <- cv.glmnet(as.matrix(data[[iSample]][["X_int"]]), data[[iSample]][["yMat"]][,iCond]) 
      fit <- glmnet(as.matrix(data[[iSample]][["X_int"]]), data[[iSample]][["yMat"]][,iCond], 
                    lambda = fit_cv$lambda.1se)
      as.matrix(fit$beta)
    })
    
    estB <- do.call(cbind, estB)
    estB_M <- rowMeans(estB, na.rm = T)
    estB_SD <- apply(estB, MARGIN = 1, sd, na.rm = T)
    estB_SE <- estB_SD / sqrt(setParam$dgp$nSamples) # standard error of the mean 
    list(estB_M = estB_M, estB_SD = estB_SD, estB_SE = estB_SE)
  })
  
  names(tmp_estB) <- setParam$dgp$condLabels
  tmp_estB <- do.call(Map, c(f = cbind, tmp_estB))
  rm(data)
  gc()
  tmp_estB
})
# takes ~ 3.5 hours (without N = 10.000)

fileName = "estB_initialResults.rda"
save(res_estB, file = paste0(resFolder, "/", fileName))


