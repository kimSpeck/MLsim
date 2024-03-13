# calculate ANOVA to check for effects of manipulated variables in simulation 
#   interactions between manipulations 
library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²

source("setParameters.R")
source("analysisTools.R")

condGridGBM <- expand.grid(N = setParam$dgp$N,
                           pTrash = setParam$dgp$pTrash,
                           reliability = setParam$dgp$reliability)

condN_pTrashGBM <- paste0("N", condGridGBM$N, 
                          "_pTrash", condGridGBM$pTrash,
                          "_rel", condGridGBM$reliability)

condGridENET <- expand.grid(N = setParam$dgp$N,
                            pTrash = setParam$dgp$pTrash,
                            reliability = setParam$dgp$reliability,
                            factors = c(TRUE, FALSE))

# duplicate data with factors to run model on indicator data, too
condGridENET$indicators <- rep(FALSE, dim(condGridENET)[1])
addIndi <- condGridENET[which(condGridENET$factors == TRUE), ]
addIndi$indicators <- TRUE
condGridENET <- rbind(addIndi, condGridENET)

condN_pTrashENET <- paste0("N", condGridENET$N, 
                           "_pTrash", condGridENET$pTrash,
                           "_rel", condGridENET$reliability,
                           "_f", ifelse(condGridENET$factors, 1, 0),
                           "_ind", ifelse(condGridENET$indicators, 1, 0))

# GBM results
# resFolder <- "results/resultsGBM"
resFolder <- "results/resultsGBMwInter"

# to do: run fitData again with NA instead of 0
load(paste0(resFolder, "/fullData.rda"))
fullDataGBM <- fullData
rm(fullData)

# enet results
resFolderENET_woInt <- paste0("results/resultsWithoutInter_Indicators") # without interactions
resFolderENET_wInt <- paste0("results/resultsWithInter_Indicators") # with interactions

# to do: run fitData again with NA instead of 0
load(paste0(resFolderENET_woInt, "/fullData.rda"))
fullDataENET_wo <- fullData

load(paste0(resFolderENET_wInt, "/fullData.rda"))
fullDataENET_w <- fullData

rm(fullData)

################################################################################
# ANOVA - R² in test sample
################################################################################
# to do:
#   - für test R2 ist sample immer dasselber aber Hyperparameter kommen von anderem training

# pull data from nested list of all results (fullData)
rSquaredGBM <- rbindResults(fullDataGBM, "perfromPerSample")
rSquaredENET_w <- rbindResults(fullDataENET_w, "perfromPerSample")
rSquaredENET_wo <- rbindResults(fullDataENET_wo, "perfromPerSample")

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
rSquaredGBM <- idx2info(rSquaredGBM, condN_pTrashGBM, type = "gbm")
rSquaredENET_w <- idx2info(rSquaredENET_w, condN_pTrashENET, type = "enet")
rSquaredENET_wo<- idx2info(rSquaredENET_wo, condN_pTrashENET, type = "enet")

# remove data with simulated factor structure
rSquaredENET_w <- rSquaredENET_w[rSquaredENET_w$factor == 0, ]
rSquaredENET_wo <- rSquaredENET_wo[rSquaredENET_wo$factor == 0, ]

rSquaredENET_w[, c("factor", "indicators")] <- list(NULL)
rSquaredENET_wo[, c("factor", "indicators")] <- list(NULL)

# add sample number as an additional variable
# # check (should be true)
# dim(rSquaredGBM)[1]/setParam$dgp$nTrain == length(setParam$dgp$N) * length(setParam$dgp$pTrash) * 
#   length(setParam$dgp$Rsquared) * length(setParam$dgp$reliability) * length(setParam$dgp$percentInter) 
rSquaredGBM$sample <- rep(seq_len(setParam$dgp$nTrain), dim(rSquaredGBM)[1]/setParam$dgp$nTrain)
rSquaredENET_w$sample <- rep(seq_len(setParam$dgp$nTrain), dim(rSquaredENET_w)[1]/setParam$dgp$nTrain)
rSquaredENET_wo$sample <- rep(seq_len(setParam$dgp$nTrain), dim(rSquaredENET_wo)[1]/setParam$dgp$nTrain)

# add model variable (within factor in mixed ANOVA)
rSquaredGBM$model <- rep("GBM", dim(rSquaredGBM)[1])
rSquaredENET_w$model <- rep("ENET - mit", dim(rSquaredENET_w)[1])
rSquaredENET_wo$model <- rep("ENET - ohne", dim(rSquaredENET_wo)[1])

# in ANOVA sample will be the ID; but all samples run from 1:100 in all simulated 
#   conditions that represent independent samples (between factors); identical sample
#   names suggest within factor in ANOVA! 
#   -> introduce additional ID variable to indicate independent samples
# independent observations for samples {1:100} in different simulated conditions
#   that all have the same sample numbers
rSquaredGBM$ID <- seq_len(dim(rSquaredGBM)[1])
rSquaredENET_w$ID <- seq_len(dim(rSquaredENET_w)[1])
rSquaredENET_wo$ID <- seq_len(dim(rSquaredENET_wo)[1])

# check
colnames(rSquaredGBM) == colnames(rSquaredENET_w) 
colnames(rSquaredENET_w) == colnames(rSquaredENET_wo)

# # merge ENET_w, ENET_wo and GBM data; concatenate data for all models
rSquaredTest <- rbind(rSquaredGBM, rSquaredENET_w, rSquaredENET_wo)

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model")
rSquaredTest[col2fac] <- lapply(rSquaredTest[col2fac], factor)

# mixed ANOVA with ...
# ... id = sample but ID (independent samples between simulated conditions)
# ... dv = {R^2}
# ... between: 3 x 2 x 3 x 3 x 3
#   N (3) {100, 300, 1000}
#   pTrash (2) {10, 50}
#   R2 (3) {0.2, 0.5, 0.8}
#   rel (3) {0.6, 0.8, 1}
#   lin_inter (3) {0.2_0.8, 0.5_0.5, 0.8_0.2}
# ... within: 
#   model {Enet - mit; Enet - ohne, GBM}

# ANOVA: 
## mit ENET - ohne
anovaTestR2 <- aov_ez(id = "ID", 
                  dv = "Rsq_test",
                  data = rSquaredTest, 
                  between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"),
                  within = "model")
# summary(anovaTestR2)
# nice(anovaTestR2)

eta2TestR2 <- eta_squared(
  anovaTestR2, # fitted model
  partial = FALSE, # not partial!
  generalized = TRUE, # generalized eta squared
  ci = 0.95,
  verbose = TRUE)

# sort generalized eta-squared results 
# which higher order interactions do we need to illustrate to report simulation results?
(eta2TestR2.ordered <- eta2TestR2[order(eta2TestR2$Eta2_generalized, decreasing = T),])

pEta2TestR2.ordered <- print_html(
  eta2TestR2.ordered[eta2TestR2.ordered$Eta2_generalized >= 0.01,],
  digits = 2)

# pEtaFormat <- format(
#   eta2TestR2.ordered,
#   digits = 2,
#   output = c("text", "markdown", "html"))


## ohne ENET - ohne
anovaTestR2sub <- aov_ez(id = "ID", 
                      dv = "Rsq_test",
                      data = rSquaredTest[rSquaredTest$model != "ENET - ohne",], 
                      between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"),
                      within = "model")
# summary(anovaTestR2)
# nice(anovaTestR2)

eta2TestR2sub <- eta_squared(
  anovaTestR2sub, # fitted model
  partial = FALSE, # not partial!
  generalized = TRUE, # generalized eta squared
  ci = 0.95,
  verbose = TRUE)

# sort generalized eta-squared results 
# which higher order interactions do we need to illustrate to report simulation results?
(eta2TestR2sub.ordered <- eta2TestR2sub[order(eta2TestR2sub$Eta2_generalized, decreasing = T),])

##### Voraussetzungsprüfungen #####
# residual vs fitted plot (Punkte sollten sich vertikal unsystematisch und 
# über die gesamte x-Achse gleich streuend um die 0 verteilen)
plot(x = fitted(anovaTestR2), y = resid(anovaTestR2)) 

## Heteroskedastizität prüfen
# standardabweichungen pro Zelle deskriptiv prüfen
allSDs <- aggregate(Rsq_test ~ N * pTrash * R2 * rel * lin_inter * model, 
                    rSquaredTest,
                    function(x) sd(x)) 

plot(allSDs$Rsq_test)
allSDs[which(allSDs$Rsq_test > 0.07),]

################################################################################
# ANOVA - overfit (ENET vs. GBM)
################################################################################
rSquaredTest$Rsq_overfit <- rSquaredTest$Rsq_train - rSquaredTest$Rsq_test 

anovaTestoverfit <- aov_ez(id = "ID", 
                      dv = "Rsq_overfit",
                      data = rSquaredTest, 
                      between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"),
                      within = "model")
# summary(anovaTestR2)
# nice(anovaTestR2)

eta2Testoverfit <- eta_squared(
  anovaTestoverfit, # fitted model
  partial = FALSE, # not partial!
  generalized = TRUE, # generalized eta squared
  ci = 0.95,
  verbose = TRUE)

# sort generalized eta-squared results 
# which higher order interactions do we need to illustrate to report simulation results?
(eta2Testoverfit.ordered <- eta2Testoverfit[order(eta2Testoverfit$Eta2_generalized, 
                                                  decreasing = T),])
# -> depending on N and model the overfit varies

##### post-tests ##### 
# for model 
# averaged over the levels of: N, pTrash, R2, rel, lin_inter 
library(emmeans)
(postOverfitModel <- emmeans(anovaTestoverfit, specs = "model"))
# most overfit in GBM > ENET - mit > ENET - ohne

# # pairwise comparison for model + Bonferroni-Holm-correction for these 3 comparissons
# pairs(postOverfitModel, adjust="holm") 

################################################################################
# ANOVA - permutation variable importance (GBM)
################################################################################
# do not compare with ENET as different measures for variable selection are used 
#   -> comparison becomes interesting as soon as we apply permutation variable 
#       importance for the ENET as well 
# id = sample
# dv = {pvi}
# between: N {100, 300, 1000}, pTrash {10, 50}, R2 {0.2, 0.5, 0.8}, rel {0.6, 0.8, 1}, lin_inter {0.2_0.8, 0.5_0.5, 0.8_0.2}

