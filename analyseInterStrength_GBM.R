# load parameters and helper functions 
source("setParameters.R")
source("analysisTools.R")

library(ggplot2)

# plot utils
colValues <- c("green3", "darkblue", "darkmagenta")

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
# to do: is Var1:Var2 == Var2:Var1 for identical samples?
#     each interaction is estimated twice for each sample!
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
# to do: maybe if variables were not selected in GBM? 
#         if linear predictor was not extracted in GBM, interactions with this variable
#         cannot be calculated?
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

# replace h-stats values with Inf by 0
interSub[idxInf,"interaction"] <- 0

hStats <- aggregate(interaction ~ idxCondLabel + N_pTrash + simEffect, 
                    data = interSub, function(x) {
                      cbind(M = mean(x),
                            SE = sd(x),
                            q025 = quantile(x, 0.025),
                            q975 = quantile(x, 0.975))})

hStats <- idx2infoNew(hStats)
hStats <- do.call(data.frame, hStats)
colnames(hStats) <- stringr::str_replace_all(colnames(hStats), "\\.1", "_M")
colnames(hStats) <- stringr::str_replace_all(colnames(hStats), "\\.2", "_SE")
colnames(hStats) <- stringr::str_replace_all(colnames(hStats), "\\.3", "_q025")
colnames(hStats) <- stringr::str_replace_all(colnames(hStats), "\\.4", "_q975")
colnames(hStats) <- stringr::str_replace_all(colnames(hStats), "^interaction_", "")

str(hStats)
col2fac <- c("N", "pTrash", "rel", "simEffect", "R2", "lin_inter")
hStats[col2fac] <- lapply(hStats[col2fac], factor) 

hStats$N <- factor(hStats$N, levels = c(100, 300, 1000))

(hStats_Overview <- ggplot(hStats, aes(x = N, y = M, 
           group = interaction(R2, simEffect), colour = R2,
           linetype = simEffect, shape = simEffect)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_shape_manual(values = c(16, 1, 8)) +
  geom_errorbar(aes(ymin = q025, ymax = q975), 
                width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) +  
  scale_color_manual(values = colValues) +
  facet_grid(rel + pTrash ~ lin_inter, labeller = label_both) +
  ylab("sensitivity") +
  xlab("N") +
  ggtitle("sensitivity: TP / (TP + FN)") +
  theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
        panel.background = element_rect(color = "white", fill = "white"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15)))


# ggplot2::ggsave(filename = paste0(plotFolder, "/hStatsInter_Overview.png"),
#                 plot = hStats_Overview,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# pTrash does not make any difference
#     in general h-Stats numbers are smaller for large pTrash numbers but difference 
#     between interactions with and without effect do not really change
# große Stichprobengröße erforderlich, um Interaktionen mit Effekt überhaupt von 
#   Interaktionen ohne Effekt trennen zu können
# mit mehr R^2 für die Interaktionen ist Trennung auf Basis der H-Statistik ebenfalls besser
# mit höherer Reliabilität ist die Trennung von simulierten Interaktionen von Interaktionen ohne
#   echten Effekt besser

(hStats_pTrash10 <- ggplot(hStats[hStats$pTrash == 10, ], 
                           aes(x = N, y = M, 
                               group = interaction(R2, simEffect), colour = R2,
                               linetype = simEffect, shape = simEffect)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_shape_manual(values = c(16, 1, 8)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), 
                  width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) +  
    scale_color_manual(values = colValues) +
    facet_grid(rel ~ lin_inter, labeller = label_both) +
    ylab("sensitivity") +
    xlab("N") +
    ggtitle("sensitivity: TP / (TP + FN)") +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# ggplot2::ggsave(filename = paste0(plotFolder, "/hStatsInter_pTrash10.png"),
#                 plot = hStats_pTrash10,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")