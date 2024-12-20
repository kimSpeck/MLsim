# result plots and ANOVA for H-statistic of the GBMs
# load parameters and helper functions 
source("utils/setParameters.R")
source("utils/analysisTools.R")

# plot H-statistic
library(ggplot2)
library(ggh4x)
library(patchwork)

# calculate ANOVA
library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²

# plot utils
colValuesR2 <- c('#db4a07', '#850c0c', '#3c1518')
colValuesInter <- c('#050440', '#181ff2', '#0eb2e8')
colValuesLin <- c('#0eb2e8', '#181ff2', '#050440')

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

# load results files
resFolder <- "results/finalResults/dependentMeasures" 

listDir <- dir(resFolder)
dataList <- listDir[stringr::str_detect(listDir, "^interStrength")]
assign("interStrength", loadRData(paste0(resFolder, "/", dataList)))

# 1. only look at interaction strength for variables with simulated linear effects 
#     measure is not reliable if variable does not have linear effects
#     (relation between interaction effect relative to linear effect)
#     var = {Var1, Var2, Var3, Var4} <- only H-Statistic for these variables extracted
# 2. remove interactions with variables without linear effects
#     only keep feature = {Var1:Var2, Var1:Var3, Var1:Var4,
#                           Var2:Var3, Var2:Var4, Var3:Var4}
#     ! any order of interacting variables is possible
# 3. sort these interactions for interactions with simulated effect vs. interactions
#     without simulated interaction effect
# 4. incorporate information of whether both linear effects were extracted based 
#     on PVI in the respective sample

# do 2. 
keepInters <- c("Var1:Var2", "Var1:Var3", "Var1:Var4", "Var2:Var3", "Var2:Var4", "Var3:Var4",
                "Var2:Var1", "Var3:Var1", "Var4:Var1", "Var3:Var2", "Var4:Var2", "Var4:Var3")
interSub <- lapply(interStrength, function(x) {
  x[x$feature %in% keepInters, ]
})

# join lists with different conditions to one data frame
interSub <- rbindSingleResults(interSub)    

# remove duplicate rows from interaction strength data
rmDuplicates <- interSub[duplicated(interSub),]
dim(interSub)
interSub <- interSub[!duplicated(interSub),]
dim(interSub)

# not all interactions do have estimated interaction strengths! fill up missing 
#   interStrength with 0!
# colnames(interSub)
countInters <- aggregate(feature ~ idxCondLabel + sample + var + N_pTrash, 
                         data = interSub, length)

# there are samples for which measures for interactions are missing
# write in data frame if variable was added by completion strategy
dim(countInters[countInters$feature != 3, ])[1]
nrMissing <- sum(3 - countInters[countInters$feature != 3, "feature"])

# not only single interactions are missing but whole conditions are missing
dim(interSub)[1] + nrMissing

# generate table with all possible interactions between variables with linear predictor effects
allCombos <- expand.grid(unique(interSub$N_pTrash), # 18 levels
                         1:4, 
                         unique(interSub$var), # 4 levels
                         unique(interSub$sample), # 1000 levels
                         unique(interSub$idxCondLabel)) # 9 levels
allCombos <- allCombos[,dim(allCombos)[2]:1] # reorder columns
colnames(allCombos) <- c("idxCondLabel", "sample", "var", "feature", "N_pTrash") # column names

# change label of possible interactions
allCombos$feature <- paste0("Var", allCombos$feature, ":", allCombos$var)
allCombos <- allCombos[allCombos$feature %in% keepInters,] # remove Var1:Var1, Var2:Var2, ...
allCombos$interaction <- 0 # give value to all rows and to keep for interactions that were not calculated

# manually added rows; manually added interaction value of 0
allCombos$added <- 1
interSub$added <- 0

col2fac <- c("idxCondLabel", "sample", "var", "feature", "N_pTrash")
interSub[col2fac] <- lapply(interSub[col2fac], factor) 
allCombos$feature <- factor(allCombos$feature)
interSub$interaction <- as.numeric(interSub$interaction)

str(allCombos)
str(interSub)

# unique(interSub$idxCondLabel) # 9 levels
# unique(interSub$sample) # 1000 levels
# unique(interSub$var) # 4 levels
# unique(interSub$feature) # 12 levels
# unique(interSub$N_pTrash) # 18 levels

# return rows in allCombos that are not present in interSub
notIn <- dplyr::anti_join(allCombos, interSub, by = col2fac)

dim(interSub)[1] + dim(notIn)[1] # there are 51 observations too much! 
# # check 
# # condLabal * sample * var * feature * N_pTrash
# 9 * 1000 * 4 * 3 * 18

# # return rows in interSub that are not present in allCombos 
# # there are no rows in interSub that are not in allCombos
# not <- dplyr::anti_join(interSub[,col2fac], allCombos[,col2fac], by = col2fac)

# # test where additional 51 rows came from (duplicates in data -> now already removed)
# interFull <- dplyr::full_join(interSub, allCombos, 
#                               by = c("idxCondLabel", "sample", "var", "feature", "N_pTrash")) 
# dim(unique(interFull[,col2fac]))
# interFull <- tidyr::unite(interFull[,col2fac], col = "ID", sep = "_", remove = F)
# tableTest <- table(interFull$ID)
# any(tableTest != 1)
# idx <- which(tableTest != 1)
# sum(tableTest[idx])
# length(tableTest[idx])
# 
# interSub[interSub$idxCondLabel == "7" &
#            interSub$sample == "156" &
#            interSub$var == "Var3" &
#            interSub$feature == "Var4:Var3" &
#            interSub$N_pTrash == "N100_pTrash50_rel1", ]
# 
# interSub$ID <- tidyr::unite(interSub[,col2fac], col = "ID", sep = "_", remove = T)
# interSub[duplicated(interSub),]
# dplyr::distinct(interSub, ID)
# interFull2 <- dplyr::right_join(interSub, allCombos, 
#                               by = c("idxCondLabel", "sample", "var", "feature", "N_pTrash")) 

interSub <- interSub[,c("idxCondLabel", "sample", "var", "feature", "N_pTrash", "interaction", "added")]
interSub <- rbind(interSub, notIn)

# which interaction has simulated effects which do not
simEffects <- c(setParam$dgp$interEffects, "Var2:Var1", "Var4:Var1", "Var3:Var2", "Var4:Var3")
interSub$simEffect <- ifelse(interSub$feature %in% simEffects, 1, 0)

# h-stats values with infinite value 
any(interSub$interaction == Inf)
idxInf <- which(interSub$interaction == Inf)
interSub[idxInf,]

# replace h-stats values with Inf by 0 (this affects 4 values)
interSub[idxInf,"interaction"] <- 0

# aggregate H-statistic effects
hStats <- aggregate(interaction ~ idxCondLabel + N_pTrash + simEffect, 
                    data = interSub, function(x) {
                      cbind(M = mean(x),
                            SD = sd(x),
                            nCond = length(x),
                            q025 = quantile(x, 0.025),
                            q975 = quantile(x, 0.975))})

hStats <- idx2infoNew(hStats)
hStats <- do.call(data.frame, hStats)
colnames(hStats) <- stringr::str_replace_all(colnames(hStats), "\\.1", "_M")
colnames(hStats) <- stringr::str_replace_all(colnames(hStats), "\\.2", "_SD")
colnames(hStats) <- stringr::str_replace_all(colnames(hStats), "\\.3", "_nCond")
colnames(hStats) <- stringr::str_replace_all(colnames(hStats), "\\.4", "_q025")
colnames(hStats) <- stringr::str_replace_all(colnames(hStats), "\\.5", "_q975")
colnames(hStats) <- stringr::str_replace_all(colnames(hStats), "^interaction_", "")

# # is count data correct?
# unique(hStats$nCond)  
# unique(interSub[interSub$simEffect == 0, "feature"]) # 4 different configurations
# unique(interSub[interSub$simEffect == 1, "feature"]) # 8 different configurations

# add standard error of the mean to the data
# ! for simulated effects the standard error of the mean might be smaller due to 
#   division by larger number of observations (8000 due to 8 configurations of the
#   interaction for 1000 samples vs. 4000 due to 4 configurations without effects)
hStats$SE <- hStats$SD / sqrt(hStats$nCond) 

# save(hStats, file = paste0(resFolder, "/hStatsPlottingData.rda"))

################################################################################
# between ANOVA for H-statistic
################################################################################
# calculate mean H-statistic within each sample for interactions with and without 
#     simulated effects

aovhStatsData <- aggregate(interaction ~ idxCondLabel + N_pTrash + sample + simEffect, 
                           data = interSub, mean)

aovhStatsData <- idx2infoNew(aovhStatsData)

# split, apply, combine to add IDs for between sample conditions
hStats_list <- split(aovhStatsData, f = aovhStatsData$simEffect)

hStats_out <- lapply(hStats_list, function(iList) {
  iList <- cbind(iList, ID = seq_len(dim(iList)[1]))   
})

aovhStatsData <- do.call(rbind, hStats_out)

colnames(aovhStatsData)
col2fac <- c("N", "pTrash", "rel", "simEffect", "R2", "lin_inter", "sample", "ID")
aovhStatsData[col2fac] <- lapply(aovhStatsData[col2fac], factor) 

# mixed ANOVA with ...
# ... id = sample (independent samples between simulated conditions)
# ... dv = {H-statistic}
# ... between: 3 x 2 x 3 x 3 x 3
#   N (3) {100, 300, 1000}
#   pTrash (2) {10, 50}
#   rel (3) {0.6, 0.8, 1}
#   R2 (3) {0.2, 0.5, 0.8}
#   lin_inter (3) {0.2_0.8, 0.5_0.5, 0.8_0.2}
# ... within: 
#   simEffect {interaction with effect vs. without effect}

aovHstats <- aov_ez(id = "ID",
                    dv = "interaction",
                    data = aovhStatsData,
                    between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"),
                    within = "simEffect")

# summary(aovHstats)

hStatsEta2 <- eta_squared(
  aovHstats, # fitted model
  partial = FALSE, # not partial!
  generalized = TRUE, # generalized eta squared
  ci = 0.95,
  verbose = TRUE)

hStatsEta2[order(hStatsEta2$Eta2_generalized, decreasing = T),]

hStatsEta2Table <- hStatsEta2[order(hStatsEta2$Eta2_generalized, decreasing = T),]

# print(xtable::xtable(hStatsEta2Table, type = "latex"), 
#       file = paste0(plotFolder, "/ANOVAresults/mixedANOVA_hStatistic.tex"))

(plotMixed_hStats <- plotEta2_4eachModel(hStatsEta2, 0.01, "#990000") +
  ggtitle("H-Statistic") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25)))

# ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_hStatistic.png"),
#                 plot = plotMixed_hStats,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

# between ANOVA for difference in H-Statistic between interactions with and without 
aovhStatsDiff <- tidyr::pivot_wider(aovhStatsData, names_from = simEffect, values_from = interaction)
aovhStatsDiff$diffH <- aovhStatsDiff$"1" - aovhStatsDiff$"0"

aovHdiffstats <- aov_ez(id = "ID",
                        dv = "diffH",
                        data = aovhStatsDiff,
                        between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"))

hDiffStatsEta2 <- eta_squared(
  aovHdiffstats, # fitted model
  partial = FALSE, # not partial!
  generalized = TRUE, # generalized eta squared
  ci = 0.95,
  verbose = TRUE)

hDiffStatsEta2[order(hDiffStatsEta2$Eta2_generalized, decreasing = T),]

hDiffStatsEta2Table <- hDiffStatsEta2[order(hDiffStatsEta2$Eta2_generalized, decreasing = T),]

# print(xtable::xtable(hDiffStatsEta2Table, type = "latex"), 
#       file = paste0(plotFolder, "/ANOVAresults/mixedANOVA_hDiffStatistic.tex"))

(plotMixed_hStatsDiff <- plotEta2_4eachModel(hDiffStatsEta2Table, 0.01, "#990000") +
    ggtitle("H-Statistic") + # for the main title
    theme(plot.title = element_text(hjust = 0.5, size = 25)))

# ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_hStatisticDifference.png"),
#                 plot = plotMixed_hStatsDiff,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")
