# parameter & general functions
source("utils/setParameters.R") # import parameter values
source("utils/simTools.R") 

# save results to this folder
dataFolder <- paste0("data/bigTestSamples")

listDir <- dir(dataFolder)

checkR2 <- vector(mode = "list", length = length(listDir))
names(checkR2) <- listDir

# iterate through both datasets
for (iData in seq_along(listDir)) {
  load(paste0(dataFolder, "/", listDir[iData]))
  print(dataList[["R2"]])
  
  Xtrain <- dataList[["X_int"]]
  if (stringr::str_detect(listDir[iData], "inter")) {
    Xtrain <- Xtrain[, colnames(Xtrain) %in% c(setParam$dgp$linEffects,
                                               setParam$dgp$interEffects)]  
    
    
    
  } else if (stringr::str_detect(listDir[iData], "pwlinear")) {
    Xtrain <- Xtrain[, colnames(Xtrain) %in% c(setParam$dgp$linEffects, 
                                               "Var5", "Var6", "Var7")] 
    
    Xtrain_lm <- Xtrain
    
    dumVar5 = createDummy(Xtrain[, "Var5"], q = 0.5, effectCoding = F) 
    dumVar6 = createDummy(Xtrain[, "Var6"], q = 0.5, effectCoding = F) 
    dumVar7 = createDummy(Xtrain[, "Var7"], q = 0.5, effectCoding = F) 
    Xtrain <- cbind(Xtrain, 
                    Var5.2nd = (Xtrain[, "Var5"] - quantile(Xtrain[, "Var5"], 0.5))*dumVar5,
                    Var6.2nd = (Xtrain[, "Var6"] - quantile(Xtrain[, "Var6"], 0.5))*dumVar6,
                    Var7.2nd = (Xtrain[, "Var7"] - quantile(Xtrain[, "Var7"], 0.5))*dumVar7)
    
    Xtrain_oracle <- Xtrain[, colnames(Xtrain) %in% c(setParam$dgp$linEffects, 
                                                      setParam$dgp$pwlinEffects)]
    
  } else if (stringr::str_detect(listDir[iData], "nonlinear3\\.")) {
    Xtrain <- Xtrain[, colnames(Xtrain) %in% c(setParam$dgp$linEffects, 
                                               "Var5", "Var6", "Var7")]  
    Xtrain_lm <- Xtrain
    
    Xtrain <- cbind(Xtrain, 
                    dumVar5.1 = createDummy(Xtrain[, "Var5"], q = 0.5, effectCoding = T),
                    dumVar6.1 = createDummy(Xtrain[, "Var6"], q = 0.5, effectCoding = T),
                    dumVar7.1 = createDummy(Xtrain[, "Var7"], q = 0.5, effectCoding = T))
    
    Xtrain <- Xtrain[, colnames(Xtrain) %in% c(setParam$dgp$linEffects, 
                                               setParam$dgp$nonlinEffects3)]  
    
  }
  
  
  # iterate through every condition in the ymat
  res <- lapply(seq_along(setParam$dgp$condLabels), function(iCond) {
    ytrain <- dataList[["yMat"]][,setParam$dgp$condLabels[iCond]]
    
    df <- data.frame(cbind(Xtrain, ytrain))
    
    if (stringr::str_detect(listDir[iData], "inter")) {
      colnames(df) <- c(setParam$dgp$linEffects, setParam$dgp$interEffects, "y")
      
    }  else if (stringr::str_detect(listDir[iData], "pwlinear")) {  
      # additionally the original variables as predictors 
      # does this affect the performance of the linear regression (due to multicollinearity)?
      colnames(df) <- c(setParam$dgp$linEffects, "Var5", "Var6", "Var7",
                        setParam$dgp$pwlinEffects, "y")
      
      # this is what the linear model gets
      df_lm <- data.frame(cbind(Xtrain_lm, ytrain))
      colnames(df_lm) <- c(setParam$dgp$linEffects, "Var5", "Var6", "Var7", "y")
      
      # this is what an oracle model would know
      df_oracle <- data.frame(cbind(Xtrain_oracle, ytrain))
      colnames(df_oracle) <- c(setParam$dgp$linEffects, setParam$dgp$pwlinEffects, "y")
      
    } else if (stringr::str_detect(listDir[iData], "nonlinear3\\.")) {
      colnames(df) <- c(setParam$dgp$linEffects, setParam$dgp$nonlinEffects3, "y")
      
      df_lm <- data.frame(cbind(Xtrain_lm, ytrain))
      colnames(df_lm) <- c(setParam$dgp$linEffects, "Var5", "Var6", "Var7", "y")
    }
    
    # fit a regression model on the true data
    #   - fit linear & other effect
    fit <- lm(y ~ ., data = df)
    R2_full <- summary(fit)$r.squared  
    
    if (!(stringr::str_detect(listDir[iData], "inter"))) {
      fit <- lm(y ~ ., data = df_lm)
      R2lm_full <- summary(fit)$r.squared  
    } 
    if (stringr::str_detect(listDir[iData], "pwlinear")) {
      fit <- lm(y ~ ., data = df_oracle)
      R2oracle_full <- summary(fit)$r.squared  
    } 
    
    
    
    #   - fit only linear effect
    df_lin <- df[,c(setParam$dgp$linEffects, "y")]
    fit <- lm(y ~ ., data = df_lin)
    R2_lin <- summary(fit)$r.squared  
    
    #   - fit only other effect
    if (stringr::str_detect(listDir[iData], "inter")) {
      
      df_other <- df[,c(setParam$dgp$interEffects, "y")] 
      
    } else if (stringr::str_detect(listDir[iData], "pwlinear")) {
      df_other <- df[,c(setParam$dgp$pwlinEffects, "y")] 
      
      # nonlinear effects variables as available in linear regression model
      fit <- lm(y ~ Var5 + Var6 + Var7, data = df_lm)
      R2lm_other <- summary(fit)$r.squared  
      
    } else if (stringr::str_detect(listDir[iData], "nonlinear3\\.")) {
      df_other <- df[,c(setParam$dgp$nonlinEffects3, "y")]
      
      # nonlinear effects variables as available in linear regression model
      fit <- lm(y ~ Var5 + Var6 + Var7, data = df_lm)
      R2lm_other <- summary(fit)$r.squared  
    }
    
    fit <- lm(y ~ ., data = df_other)
    R2_other <- summary(fit)$r.squared  
    
    trueR2 <- as.numeric(stringr::str_sub(setParam$dgp$condLabels[iCond], start = 3, end = 5))
    trueLinPC <- as.numeric(stringr::str_sub(setParam$dgp$condLabels[iCond], start = -7, end = -5)) 
    trueOtherPC <- as.numeric(stringr::str_sub(setParam$dgp$condLabels[iCond], start = -3, end = -1))
    
    trueLin <- trueLinPC * trueR2
    trueOther <- trueOtherPC * trueR2
    
    c(trueR2 = trueR2, R2_full = R2_full, 
      R2lm_full = if (stringr::str_detect(listDir[iData], "inter")) {NA} else {R2lm_full},
      # oracle model with only true predictors
      R2oracle_full = if (stringr::str_detect(listDir[iData], "pwlinear")) {R2oracle_full} else {NA},
      trueLin = trueLin, R2_lin = R2_lin, trueLinPC = trueLinPC,
      trueOther = trueOther, R2_other = R2_other, 
      R2lm_other = if (stringr::str_detect(listDir[iData], "inter")) {NA} else {R2lm_other},
      trueOtherPC = trueOtherPC)
  })
  
  checkR2[[iData]] <- do.call(rbind, res)
  
}

save(checkR2, file = paste0("results/dependentMeasures/checkR2manipulation.rda"))
