# result plots and ANOVA for H-statistic of the GBMs
# load parameters and helper functions 
source("setParameters.R")
source("analysisTools.R")

# plot H-statistic
library(ggplot2)
library(ggh4x)
library(patchwork)

# calculate ANOVA
library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²

# plot utils
# colValues <- c("green3", "darkblue", "darkmagenta")
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

################################################################################
# plot h-Statistic effects
################################################################################
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
  scale_color_manual(values = colValuesR2) +
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

# save overview plot (with both trash variable conditions)
ggplot2::ggsave(filename = paste0(plotFolder, "/detectInteractions/hStatsInter_Overview.png"),
                plot = hStats_Overview,
                width = 13.08,
                height = 12.18,
                units = "in")

# result notes
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
    geom_point(position = position_dodge(width = 0.5), size = 2.2) +
    geom_line(position = position_dodge(width = 0.5), alpha = 0.2) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_shape_manual(values = c(16, 1, 8)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), 
                  width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) +  
    scale_color_manual(values = colValuesR2) +
    facet_grid(rel ~ lin_inter, labeller = label_both) +
    ylab("H-statistic") +
    xlab("N") +
    #ggtitle("sensitivity: TP / (TP + FN)") +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

ggplot2::ggsave(filename = paste0(plotFolder, "/detectInteractions/hStatsInter_pTrash10.png"),
                plot = hStats_pTrash10,
                width = 13.08,
                height = 12.18,
                units = "in")

######
hStats$lin_inter <- factor(hStats$lin_inter, 
                           levels = c("0.8_0.2", "0.5_0.5", "0.2_0.8"),
                           labels = c("80:20", "50:50", "20:80"))
# plot subset to display all relevant results 
(hStats_pTrash10rel1 <- ggplot(hStats[hStats$pTrash == 10 &
                                    hStats$rel == 1, ], 
                           aes(x = N, y = M, 
                               group = interaction(R2, simEffect), colour = R2,
                               linetype = simEffect, shape = simEffect)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_shape_manual(values = c(16, 17)) +
    # # plot 95% quantiles as they cannot reach impossible values
    # geom_errorbar(aes(ymin = q025, ymax = q975), 
    #               width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) + 
    # # plot two standard errors from the mean
    geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), 
                 width = 0.2, alpha = 0.4) +  
    scale_color_manual(name = "R²", values = colValuesR2) +
    # facet_grid(rel ~ lin_inter, labeller = label_both) +
    ylab("H-statistic") +
   facet_grid2(rel ~ lin_inter,
               strip = strip_themed(
                 background_x = list(element_rect(fill = alpha(colValuesLin[3], 0.4)),
                                     element_rect(fill = alpha(colValuesLin[2], 0.4)),
                                     element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
    xlab("N") +
    #ggtitle("sensitivity: TP / (TP + FN)") +
    theme(panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
          panel.background = element_rect(color = "white", fill = "white"),
          #plot.title = element_text(size = 30, face = "bold"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.position = c(.2, .85), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal"))

(hStats_pTrash10rel06 <- ggplot(hStats[hStats$pTrash == 10 &
                                        hStats$rel == 0.6 & 
                                        hStats$lin_inter == "20:80", ], 
                               aes(x = N, y = M, 
                                   group = interaction(R2, simEffect), colour = R2,
                                   linetype = simEffect, shape = simEffect)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_shape_manual(values = c(16, 17)) +
    guides(shape = "none", linetype = "none", colour = "none") +
    # # plot 95% quantiles as they cannot reach impossible values
    # geom_errorbar(aes(ymin = q025, ymax = q975), 
    #               width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) +  
    # # plot two standard errors from the mean
    geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), 
                  width = 0.2, alpha = 0.4) +  
    scale_color_manual(values = colValuesR2) +
    facet_grid2(rel ~ lin_inter,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
    ylim(0, 0.25) +
    ylab("H-statistic") +
    xlab("N") +
    #ggtitle("sensitivity: TP / (TP + FN)") +
    theme(panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
          panel.background = element_rect(color = "white", fill = "white"),
          #plot.title = element_text(size = 30, face = "bold"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.position = c(.2, .85), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal"))

# patch subplots 
(phStatsSub <- hStats_pTrash10rel06 + 
  hStats_pTrash10rel1 + 
  plot_layout(ncol = 2, widths = c(1, 3)))

ggplot2::ggsave(filename = paste0(plotFolder, "/detectInteractions/hStatsInterSubset_pTrash10_2SE.png"),
                plot = phStatsSub,
                width = 17.89,
                height = 5.88,
                units = "in")

(hStats_pTrash50rel1 <- ggplot(hStats[hStats$pTrash == 50 &
                                        hStats$rel == 1, ], 
                               aes(x = N, y = M, 
                                   group = interaction(R2, simEffect), colour = R2,
                                   linetype = simEffect, shape = simEffect)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_shape_manual(values = c(16, 17)) +
    # # plot 95% quantiles as they cannot reach impossible values
    # geom_errorbar(aes(ymin = q025, ymax = q975), 
    #               width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) + 
    # # plot two standard errors from the mean
    geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), 
                  width = 0.2, alpha = 0.4) +  
    scale_color_manual(name = "R²", values = colValuesR2) +
    # facet_grid(rel ~ lin_inter, labeller = label_both) +
    ylab("H-statistic") +
    facet_grid2(rel ~ lin_inter,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[3], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[2], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
    xlab("N") +
    #ggtitle("sensitivity: TP / (TP + FN)") +
    theme(panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
          panel.background = element_rect(color = "white", fill = "white"),
          #plot.title = element_text(size = 30, face = "bold"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.position = c(.2, .85), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal"))

(hStats_pTrash50rel06 <-  ggplot(hStats[hStats$pTrash == 50 &
                                          hStats$rel == 0.6 & 
                                          hStats$lin_inter == "20:80", ], 
                                 aes(x = N, y = M, 
                                     group = interaction(R2, simEffect), colour = R2,
                                     linetype = simEffect, shape = simEffect)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_shape_manual(values = c(16, 17)) +
    guides(shape = "none", linetype = "none", colour = "none") +
    # # plot 95% quantiles as they cannot reach impossible values
    # geom_errorbar(aes(ymin = q025, ymax = q975), 
    #               width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) +  
    # # plot two standard errors from the mean
    geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), 
                  width = 0.2, alpha = 0.4) +  
    scale_color_manual(values = colValuesR2) +
    facet_grid2(rel ~ lin_inter,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
    ylim(0, 0.25) +
    ylab("H-statistic") +
    xlab("N") +
    #ggtitle("sensitivity: TP / (TP + FN)") +
    theme(panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
          panel.background = element_rect(color = "white", fill = "white"),
          #plot.title = element_text(size = 30, face = "bold"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.position = c(.2, .85), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal"))

# patch subplots 
(phStatsSub_pTrash50 <- hStats_pTrash50rel06 + 
    hStats_pTrash50rel1 + 
    plot_layout(ncol = 2, widths = c(1, 3)))

ggplot2::ggsave(filename = paste0(plotFolder, "/detectInteractions/hStatsInterSubset_pTrash50_2SE.png"),
                plot = phStatsSub_pTrash50,
                width = 17.89,
                height = 5.88,
                units = "in")

################################################################################
# difference plots for H-statistic
################################################################################
hStatsDiff <- tidyr::pivot_wider(hStats[c("N", "pTrash", "rel", "R2", "lin_inter", "simEffect", "M")], 
                                 names_from = simEffect, 
                                 values_from =  M)

hStatsDiff$Hdiff <- hStatsDiff$'1' - hStatsDiff$'0'

(hStats_pTrash10 <- ggplot(hStatsDiff[hStatsDiff$pTrash == 10, ], 
                               aes(x = N, y = Hdiff, 
                                   group = interaction(R2, rel), colour = R2,
                                   linetype = rel, shape = rel)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(values = c(18, 16, 17)) +
    # # plot 95% quantiles as they cannot reach impossible values
    # geom_errorbar(aes(ymin = q025, ymax = q975), 
    #               width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) + 
    # # plot two standard errors from the mean
    # geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), 
    #               width = 0.2, alpha = 0.4) +  
    scale_color_manual(name = "R²", values = colValuesR2) +
    # facet_grid(rel ~ lin_inter, labeller = label_both) +
    ylab("H-statistic difference") +
    facet_grid2(~ lin_inter,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[3], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[2], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
    xlab("N") +
    #ggtitle("sensitivity: TP / (TP + FN)") +
    theme(panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
          panel.background = element_rect(color = "white", fill = "white"),
          #plot.title = element_text(size = 30, face = "bold"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.position = c(.2, .85), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal"))

(hStats_pTrash50 <- ggplot(hStatsDiff[hStatsDiff$pTrash == 50, ], 
                           aes(x = N, y = Hdiff, 
                               group = interaction(R2, rel), colour = R2,
                               linetype = rel, shape = rel)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(values = c(18, 16, 17)) +
    # # plot 95% quantiles as they cannot reach impossible values
    # geom_errorbar(aes(ymin = q025, ymax = q975), 
    #               width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) + 
    # # plot two standard errors from the mean
    # geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), 
    #               width = 0.2, alpha = 0.4) +  
    scale_color_manual(name = "R²", values = colValuesR2) +
    # facet_grid(rel ~ lin_inter, labeller = label_both) +
    ylab("H-statistic difference") +
    facet_grid2(~ lin_inter,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[3], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[2], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
    xlab("N") +
    #ggtitle("sensitivity: TP / (TP + FN)") +
    theme(panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
          panel.background = element_rect(color = "white", fill = "white"),
          #plot.title = element_text(size = 30, face = "bold"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.position = c(.2, .85), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal"))

ggplot2::ggsave(filename = paste0(plotFolder, "/detectInteractions/hStatsInterSubset_pTrash50_diff.png"),
                plot = hStats_pTrash50,
                width = 13.08,
                height = 5.88,
                units = "in")

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

print(xtable::xtable(hStatsEta2Table, type = "latex"), 
      file = paste0(plotFolder, "/ANOVAresults/mixedANOVA_hStatistic.tex"))

(plotMixed_hStats <- plotEta2_4eachModel(hStatsEta2, 0.01, "#990000") +
  ggtitle("H-Statistic") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25)))

ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_hStatistic.png"),
                plot = plotMixed_hStats,
                width = 17.52,
                height = 10.76,
                units = "in")

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

print(xtable::xtable(hDiffStatsEta2Table, type = "latex"), 
      file = paste0(plotFolder, "/ANOVAresults/mixedANOVA_hDiffStatistic.tex"))

(plotMixed_hStatsDiff <- plotEta2_4eachModel(hDiffStatsEta2Table, 0.01, "#990000") +
    ggtitle("H-Statistic") + # for the main title
    theme(plot.title = element_text(hjust = 0.5, size = 25)))

ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_hStatisticDifference.png"),
                plot = plotMixed_hStatsDiff,
                width = 17.52,
                height = 10.76,
                units = "in")
