# To Do: 
#     - plot positive predictive value (PPV) for variable selection
#       number of true positives (= selected variables that are actually linear/interaction effects)
#       divided by: number of true positives + number of false positives (= selected variables without effect)
#     - mean permutation variable importance for false positive predictors? 
library(ggplot2)

source("setParameters.R")
source("analysisTools.R")

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

condGrid <- expand.grid(N = setParam$dgp$N,
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability)

condN_pTrash <- paste0("N", condGrid$N, 
                       "_pTrash", condGrid$pTrash,
                       "_rel", condGrid$reliability)


resFolder <- "results/finalResults/dependentMeasures"

listDir <- dir(resFolder)
dataList <- listDir[stringr::str_detect(listDir, "^performT")]
models <- stringr::str_extract(dataList, "_[:alpha:]*.rda$")
models <- stringr::str_sub(models, start = 2L, end = -5)
tvst <- stringr::str_extract(dataList, "(Test|Train)") 
for (iData in seq_len(length(dataList))) {
  objectName <- paste0("p", tvst[iData], "_", models[iData])
  assign(objectName, loadRData(paste0(resFolder, "/", dataList[iData])))
}

################################################################################
# plot train and test performance
################################################################################
# change list of  matrices to data.frame
pTest_GBM <- rbindSingleResults(pTest_GBM)
pTrain_GBM <- rbindSingleResults(pTrain_GBM)

# variables from rownames to own column to work with variable information
pTrain_GBM <- rowNames2col(pTrain_GBM, "measure")
pTrain_GBM$measure <- paste0(pTrain_GBM$measure, "_train")

pTest_GBM <- rowNames2col(pTest_GBM, "measure")
pTest_GBM$measure[which(pTest_GBM$measure == "Rsq_test")] <- "Rsquared_test"

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
pTrain_GBM <- idx2infoNew(pTrain_GBM)
pTest_GBM <- idx2infoNew(pTest_GBM) 

# merge performance Train and performance Test
performanceStats <- rbind(pTrain_GBM, pTest_GBM)
performanceStats <- tidyr::separate(performanceStats, measure, c("measure", "trainTest"), sep = "_")

# change type of columns or specific entry details to prepare plotting  
performanceStats$N <- factor(performanceStats$N, levels = setParam$dgp$N)
performanceStats$pTrash <- factor(performanceStats$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
str(performanceStats)
chr2fac <- c("rel", "measure", "R2", "lin_inter")
performanceStats[chr2fac] <- lapply(performanceStats[chr2fac], factor)


# plot overfit instead of train as train - test
performanceStats <- tidyr::pivot_wider(performanceStats, 
                                       names_from = "trainTest", values_from = c(M, SE, SD))
chr2num <- c("M_train", "M_test", "SE_train", "SE_test", "SD_train", "SD_test")
performanceStats[chr2num] <- lapply(performanceStats[chr2num], as.numeric)

performanceStats$overfit <- performanceStats$M_train - performanceStats$M_test
performanceSub <- performanceStats[which(performanceStats$measure != "MAE" &
                                           performanceStats$measure != "RMSE"),]
performanceSub[,c("M_train", "SE_train", "SE_test", "SD_train", "SD_test")] <- list(NULL)
performanceSub <- tidyr::pivot_longer(performanceSub, c(M_test, overfit),
                                      names_to = "measures", values_to = "values")

# plot performance measures for train and test data
colValues <- c("green3", "darkblue", "darkmagenta")

(pPerformTrainVStestGBM <- ggplot(performanceSub,
                               aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                                   group = R2, colour = R2)) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
    # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(measures + rel ~ lin_inter, labeller = label_both) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "green3",
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "darkblue",
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "darkmagenta",
               alpha = 0.4) +
    ylab("") +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle("R^2: Training vs. Test performance") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTestGBM.eps"),
#                 plot = pPerformTrainVStestGBM,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTestGBM.png"),
#                 plot = pPerformTrainVStestGBM,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

################################################################################
# ENET results 
################################################################################
# change list of  matrices to data.frame
pTest_ENETw <- rbindSingleResults(pTest_ENETw)
pTest_ENETwo <- rbindSingleResults(pTest_ENETwo)
pTrain_ENETw <- rbindSingleResults(pTrain_ENETw)
pTrain_ENETwo <- rbindSingleResults(pTrain_ENETwo)

# variables from rownames to own column to work woth variable information
pTrain_ENETw <- rowNames2col(pTrain_ENETw, "measure")
pTrain_ENETw$measure <- paste0(pTrain_ENETw$measure, "_train")

pTrain_ENETwo <- rowNames2col(pTrain_ENETwo, "measure")
pTrain_ENETwo$measure <- paste0(pTrain_ENETwo$measure, "_train")

pTest_ENETw <- rowNames2col(pTest_ENETw, "measure")
pTest_ENETw$measure[which(pTest_ENETw$measure == "Rsq_test")] <- "Rsquared_test"

pTest_ENETwo <- rowNames2col(pTest_ENETwo, "measure")
pTest_ENETwo$measure[which(pTest_ENETwo$measure == "Rsq_test")] <- "Rsquared_test"

# join data
# pTrainENET <- rbind(pTrain_ENETw, pTrain_ENETwo)
# pTestENET <- rbind(pTest_ENETw, pTest_ENETwo)

# merge performance Train and performance Test
performanceStatsENET <- rbind(pTrain_ENETw, pTrain_ENETwo, pTest_ENETw, pTest_ENETwo)
performanceStatsENET <- tidyr::separate(performanceStatsENET, measure, c("measure", "trainTest"), sep = "_")
# unique(performanceStatsENET$trainTest)
# rm(pTrain_ENETw, pTrain_ENETwo, pTest_ENETw, pTest_ENETwo)

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
# pTrainENET <- idx2infoNew(pTrainENET)
# pTestENET <- idx2infoNew(pTestENET)
performanceStatsENET <- idx2infoNew(performanceStatsENET)

# code if interaction were fitted or not
unique(performanceStatsENET$model)
performanceStatsENET$fitInter <- ifelse(performanceStatsENET$model == "ENETw", 1, 0)

# change type of columns or specific entry details to prepare plotting  
str(performanceStatsENET)
performanceStatsENET$N <- factor(performanceStatsENET$N, levels = setParam$dgp$N)
performanceStatsENET$pTrash <- factor(performanceStatsENET$pTrash, 
                                      levels = sort(setParam$dgp$pTrash, decreasing = T))

# plot overfit instead of train as train - test
performanceStatsENET <- tidyr::pivot_wider(performanceStatsENET, 
                                           names_from = trainTest, values_from = c(M, SE, SD))

str(performanceStatsENET)
chr2fac <- c("model", "rel", "measure", "R2", "lin_inter", "fitInter")
performanceStatsENET[chr2fac] <- lapply(performanceStatsENET[chr2fac], factor)
chr2num <- c("M_train", "M_test", "SE_train", "SE_test", "SD_train", "SD_test")
performanceStatsENET[chr2num] <- lapply(performanceStatsENET[chr2num], as.numeric)

performanceStatsENET$overfit <- performanceStatsENET$M_train - performanceStatsENET$M_test
performanceSubENET <- performanceStatsENET[which(performanceStatsENET$measure != "MAE" &
                                               performanceStatsENET$measure != "RMSE"),]
performanceSubENET[,c("M_train", "SE_train", "SE_test", "SD_train", "SD_test")] <- list(NULL)
performanceSubENET <- tidyr::pivot_longer(performanceSubENET, c(M_test, overfit),
                                      names_to = "measures", values_to = "values")

head(performanceSubENET)
head(performanceSub)
performanceSub$fitInter <- rep(2, dim(performanceSub)[1])

colnames(performanceSub)
colnames(performanceSubENET)
performanceData <- rbind(performanceSub, performanceSubENET)
# performanceData$fitInter <- plyr::mapvalues(performanceData$fitInter, 
#                                             from=c(0,1,2), 
#                                             to=c("ENET - ohne","ENET - mit","GBM"))

unique(performanceData$fitInter)
performanceData$fit <- ifelse(performanceData$fitInter != 2, "enet", "gbm")
unique(performanceData$fit)
unique(performanceData$model)

# to use in analyseDataEnetGBM.R to calculate ANOVA
# save(performanceData, file = "results/finalResults/dependentMeasures/RsquaredData_stats.rda")

################################################################################
# hier einsteigen?
################################################################################
load("results/finalResults/dependentMeasures/RsquaredData_stats.rda")
load("results/finalResults/dependentMeasures/rSquaredData_eachSample.rda")


##### basic plot from exploratory tests #####
# plot performance measures for train and test data
colValues <- c("green3", "darkblue", "darkmagenta")

(pPerformTrainVStest <- ggplot(performanceData,
                               aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                               group = interaction(R2, model, fit), colour = R2,
                               linetype = model, shape = fit)) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_shape_manual(values = c(16, 8)) +
    # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(measures + rel ~ lin_inter, labeller = label_both) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "green3",
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "darkblue",
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "darkmagenta",
               alpha = 0.4) +
    ylab("") +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle("R^2: Training vs. Test performance") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTest_ENTvsGBM.eps"),
#                 plot = pPerformTrainVStest,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTest_ENTvsGBM.png"),
#                 plot = pPerformTrainVStest,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

##### basic plot from exploratory tests #####
# pTrash raus lassen (nicht mehr auf y-Achse; aber beide Varianten vergleichen) 

# x Achse: rel
# Farbe: R^2
# 3 Facetten für lin_inter

# Facetten für N zusätzlich ()
pR2_overview <- ggplot(performanceData[performanceData$pTrash == 50 &
                         performanceData$measures == "M_test", ],
       aes(x = rel, y = values, 
           group = interaction(R2, model, fit), colour = R2,
           linetype = model, shape = fit)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_shape_manual(values = c(16, 8)) +
  # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(N ~ lin_inter, labeller = label_both) +
  geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "green3",
             alpha = 0.4) +
  geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "darkblue",
             alpha = 0.4) +
  geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "darkmagenta",
             alpha = 0.4) +
  ylab("") +
  xlab("reliability of predictors") +
  ggtitle("R^2: Training vs. Test performance") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/R2_newOverview.png"),
#                 plot = pR2_overview,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

##### add Monte Carlo Error
# as standard deviation of the Monte Carlo estimate
# plot 2 Monte Carlo Errors as Error Bars  

getMCE <- aggregate(Rsq_test ~ model + N + pTrash + rel + R2 + lin_inter, 
          data = rSquaredTest, sd)
aggregate(Rsq_test ~ model + N + pTrash + rel + R2 + lin_inter, 
          data = rSquaredTest, mean)

pR2sub <- performanceData[performanceData$measure == "Rsquared" &
                            performanceData$measures == "M_test", ]

pR2sub <- merge(pR2sub, getMCE, by = c("model", "N", "pTrash", "rel", "R2", "lin_inter"))
colnames(pR2sub)[colnames(pR2sub) == 'Rsq_test'] <- 'MCE'

(pR2_N100 <- ggplot(pR2sub[pR2sub$pTrash == 50 &
                             pR2sub$N == 100, ],
                               aes(x = rel, y = values, 
                                   group = interaction(R2, model, fit), colour = R2,
                                   linetype = model, shape = fit)) +
   geom_point() +
   geom_line() +
   geom_errorbar(aes(ymin = values - 2* MCE, ymax = values + 2*MCE), 
                 width = 0.2, alpha = 0.4) +  
   scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
   scale_shape_manual(values = c(16, 8)) +
   scale_color_manual(values = colValues) +
   geom_hline(aes(yintercept = 0)) +
   facet_wrap(~ lin_inter) +
   geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "green3",
              alpha = 0.4) +
   geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "darkblue",
              alpha = 0.4) +
   geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "darkmagenta",
              alpha = 0.4) +
   ylab("test R²") +
   xlab("reliability of predictors") +
   ggtitle("R²: Training vs. Test performance {N = 100, pTrash = 50}") +
   theme(axis.text.y = element_text(size = 20),
         axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         strip.text.x = element_text(size = 15),
         strip.text.y = element_text(size = 15)))

(pR2_N100_0MCE <- ggplot(pR2sub[pR2sub$pTrash == 50 &
                             pR2sub$N == 100, ],
                    aes(x = rel, y = values, 
                        group = interaction(R2, model, fit), colour = R2,
                        linetype = model, shape = fit)) +
    geom_point() +
    geom_line() +
    #geom_errorbar(aes(ymin = values - 2* MCE, ymax = values + 2*MCE), width = 0.2) +  
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_shape_manual(values = c(16, 8)) +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_wrap(~ lin_inter) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "green3",
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "darkblue",
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "darkmagenta",
               alpha = 0.4) +
    ylab("test R²") +
    xlab("reliability of predictors") +
    ggtitle("R²: Training vs. Test performance {N = 100, pTrash = 50}") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

(pR2_N1000 <- ggplot(pR2sub[pR2sub$pTrash == 50 &
                             pR2sub$N == 1000, ],
                    aes(x = rel, y = values, 
                        group = interaction(R2, model, fit), colour = R2,
                        linetype = model, shape = fit)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = values - 2*MCE, ymax = values + 2*MCE), 
                  width = 0.2, alpha = 0.4) +  
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_shape_manual(values = c(16, 8)) +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_wrap(~ lin_inter) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "green3",
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "darkblue",
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "darkmagenta",
               alpha = 0.4) +
    ylab("test R²") +
    xlab("reliability of predictors") +
    ggtitle("R²: Training vs. Test performance {N = 1000, pTrash = 50}") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/R2_N100_2MCE.png"),
#                 plot = pR2_N100,
#                 width = 13.63,
#                 height = 8.59,
#                 units = "in")

# ggplot2::ggsave(filename = paste0(plotFolder, "/R2_N100_0MCE.png"),
#                 plot = pR2_N100_0MCE,
#                 width = 13.63,
#                 height = 8.59,
#                 units = "in")

# ggplot2::ggsave(filename = paste0(plotFolder, "/R2_N1000_2MCE.png"),
#                 plot = pR2_N1000,
#                 width = 13.63,
#                 height = 8.59,
#                 units = "in")



################################################################################
# permutation variable importance measure
# 
# stats calculation for ENET is reused in interaction strength plots
################################################################################
# here!
#### only GBM
# pull data from nested list of all results (fullData)
pvi <- rbindResults(fullDataGBM, "pvi")
head(pvi)

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
pvi <- idx2info(pvi, condN_pTrashGBM, type = "gbm")

# Wie häufig wurden lineare Effekte gefunden?
# hier hat GBM einen Vorteil, da in permutation variable importance der "linearen"
#   Effekte bereits die Effektstaerke, die auf Interaktionen simuliert wurde eingehen kann. 

idxPVI <- which(pvi$pviRank %in% c(setParam$dgp$linEffects))
pviEffects <- pvi[idxPVI,]

# pviValue > 1 == permutation der variable fuehrt zu mehr prediction error 
#   folglich hat Variable prädiktiven Wert
pviEffectsSub <- pviEffects[which(pviEffects$pviValue > 1),]

# count how often all linear effects were extracted per condition
pviEffectsSub <- tidyr::unite(pviEffectsSub, "N_pTrash_rel_R2_lin_inter", 
                              c(N, pTrash, rel, R2, lin_inter),
                              sep = "_", remove = FALSE) 
pviEffectsSub_list <- split(pviEffectsSub, 
                            f = pviEffectsSub$N_pTrash_rel_R2_lin_inter)

counter <- lapply(seq_along(pviEffectsSub_list), function(iSublist) {
  pviEffectsSub_list[[iSublist]]$pviValue <- as.numeric(pviEffectsSub_list[[iSublist]]$pviValue)
  # amount of extracted linear effects 
  nLin <- aggregate(pviValue ~ sample, 
            data = pviEffectsSub_list[[iSublist]], length)
  # mean permutation variable importance; not substantively interpretable but maybe 
  #     across different simulated condition?
  mLin <- mean(pviEffectsSub_list[[iSublist]]$pviValue)
  seLin <- sd(pviEffectsSub_list[[iSublist]]$pviValue)
  cbind(N_pTrash_rel_R2_lin_inter = names(pviEffectsSub_list)[iSublist],
        sample = nLin$sample,
        mLin = rep(mLin, length(nLin$sample)),
        seLin = rep(seLin, length(nLin$sample)),
        nLin = nLin$pviValue,
        allLin = ifelse(nLin$pviValue == length(setParam$dgp$linEffects), 1, 0))
})  

pviAllLin <- data.frame(do.call(rbind, counter))

pviAllLinCount <- aggregate(allLin ~ N_pTrash_rel_R2_lin_inter, 
                            data = pviAllLin[which(pviAllLin$allLin == 1),], length)

pviAllLinCount <- tidyr::separate(pviAllLinCount, N_pTrash_rel_R2_lin_inter, 
                             into = c("N", "pTrash", "rel", "R2", "lin", "inter"), sep = "_")
pviAllLinCount <- tidyr::unite(pviAllLinCount, "lin_inter", c(lin, inter), sep = "_")

pviAllLinCount$N <- factor(pviAllLinCount$N, levels = setParam$dgp$N)
pviAllLinCount$pTrash <- factor(pviAllLinCount$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
pviAllLinCount$rel <- factor(pviAllLinCount$rel, levels = setParam$dgp$rel)

colValues <- c("green3", "darkblue", "darkmagenta")

ggplot(pviAllLinCount,
       aes(x = interaction(pTrash, N, sep = " x "), y = allLin, 
           group = R2, colour = R2)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(rel ~ lin_inter, scales = "free_y", labeller = label_both, switch = "y") +
  ylab(paste0("n samples out of ", setParam$dgp$nTrain, " samples")) +
  xlab("pTrash (decreasing) x N (increasing)") +
  ggtitle("effect recovery across samples") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

# mean permutation variable importance and its standard error
str(pviAllLin)
pviAllLin$mLin <- as.numeric(pviAllLin$mLin)
pviAllLin$seLin <- as.numeric(pviAllLin$seLin)
pviAllLinStats <- aggregate(mLin ~ N_pTrash_rel_R2_lin_inter, 
                            data = pviAllLin, mean)
pviAllLinStats$seLin <- aggregate(seLin ~ N_pTrash_rel_R2_lin_inter, 
                            data = pviAllLin, mean)$seLin

pviAllLinStats <- tidyr::separate(pviAllLinStats, N_pTrash_rel_R2_lin_inter, 
                                  into = c("N", "pTrash", "rel", "R2", "lin", "inter"), sep = "_")
pviAllLinStats <- tidyr::unite(pviAllLinStats, "lin_inter", c(lin, inter), sep = "_")

pviAllLinStats$N <- factor(pviAllLinStats$N, levels = setParam$dgp$N)
pviAllLinStats$pTrash <- factor(pviAllLinStats$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
pviAllLinStats$rel <- factor(pviAllLinStats$rel, levels = setParam$dgp$rel)

pPVIstats <- ggplot(pviAllLinStats,
       aes(x = interaction(pTrash, N, sep = " x "), y = mLin, 
           group = R2, colour = R2)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 1)) +
  ylim(0, 5) +
  geom_errorbar(aes(ymin = mLin - seLin, ymax = mLin + seLin), 
                width = 0.2, alpha = 0.5) +
  # facet_grid(rel ~ lin_inter, scales = "free_y", labeller = label_both, switch = "y") +
  facet_grid(rel ~ lin_inter, labeller = label_both, switch = "y") +
  ylab(paste0("n samples out of ", setParam$dgp$nTrain, " samples")) +
  xlab("pTrash (decreasing) x N (increasing)") +
  ggtitle("effect recovery across samples") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/pPVIstatsGBM.eps"),
#                 plot = pPVIstats,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/pPVIstatsGBM.png"),
#                 plot = pPVIstats,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# standard error of the permutation variable importance measure increased for small sample
#   sizes, especially with large amount of trash variables 
# -> does this indicate that the permutation variable importance measure does not 
#     work?


#### for elastic net data
# pull data from nested list of all results (fullData)
selectedVarsENET_w <- rbindResults(fullDataENET_w, "selectSample")
selectedVarsENET_wo <- rbindResults(fullDataENET_wo, "selectSample")
selectedVarsENET_w$sample <- rep(seq_len(setParam$dgp$nTrain), 
                                 times = length(setParam$dgp$condLabels) * nrow(condGridENET))
selectedVarsENET_wo$sample <- rep(seq_len(setParam$dgp$nTrain), 
                                 times = length(setParam$dgp$condLabels) * nrow(condGridENET))

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
selectedVarsENET_wo <- idx2info(selectedVarsENET_wo, condN_pTrashENET, type = "enet")
selectedVarsENET_w <- idx2info(selectedVarsENET_w, condN_pTrashENET, type = "enet")

selectedVarsENET_wo <- selectedVarsENET_wo[which(selectedVarsENET_wo$factor == 0), ]
selectedVarsENET_w <- selectedVarsENET_w[which(selectedVarsENET_w$factor == 0), ]

selectedVarsENET_wo$fitInter <- rep(0, dim(selectedVarsENET_wo)[1])
selectedVarsENET_w$fitInter <- rep(1, dim(selectedVarsENET_w)[1])

selectedVars <- rbind(selectedVarsENET_w, selectedVarsENET_wo)
# unique(selectedVars$sample)
rm(selectedVarsENET_w, selectedVarsENET_wo)

# get number of unique predictors depending on pTrash
names(setParam$dgp$nModelPredictors) <- setParam$dgp$pTrash
selectedVars$nModelPredictors <- setParam$dgp$nModelPredictors[match(selectedVars$pTrash, 
                                                                     names(setParam$dgp$nModelPredictors))]

# samples that exactly recovered the simulated model (all simulated effects & no other predictors)
selectedVars$exactModel <- ifelse((selectedVars$all.T1F0 == 1) & 
                                    (selectedVars$nOthers == 0), 1, 0)

idxExactModel <- which(selectedVars$exactModel == 1)
selectedVars[idxExactModel,]

# split, apply, combine
# split data in different simulated conditions N x pTrash x R2 x lin_inter
selectedVars <- tidyr::unite(selectedVars, "N_pTrash_rel_R2_lin_inter_fitInter", 
                             c(N, pTrash, rel, R2, lin_inter, fitInter),
                             sep = "_", remove = FALSE) 
subLists <- split(selectedVars, 
                  f = selectedVars$N_pTrash_rel_R2_lin_inter_fitInter)

# identical computations for every subList (i.e., simulated condition)
# compute dependent measures to visualize in results plot
resSampleStats <- lapply(seq_along(subLists), function(iCond) {
  nLin_M = mean(subLists[[iCond]]$nLin, na.rm = T)
  nInter_M = mean(subLists[[iCond]]$nInter, na.rm = T)
  nAllLin = sum(subLists[[iCond]]$nLin == length(setParam$dgp$linEffects))
  nAllInter = sum(subLists[[iCond]]$nInter == length(setParam$dgp$interEffects))
  nAllEffects = sum(subLists[[iCond]]$all.T1F0)
  nOthers_M = mean(subLists[[iCond]]$nOthers, na.rm = T)
  nModelPredictors <- unique(subLists[[iCond]]$nModelPredictors)
  # to do: this ignores that the true predictors are not "others"; but true predictors are constant across conditions
  percentOthers = nOthers_M / nModelPredictors * 100 
  nExactModel_M = sum(subLists[[iCond]]$exactModel)
  cbind(N_pTrash_rel_R2_lin_inter_fitInter = names(subLists)[iCond], 
        nLin_M = nLin_M, nInter_M = nInter_M,
        nAllLin = nAllLin, nAllInter = nAllInter, nAllEffects = nAllEffects,
        nOthers_M = nOthers_M, percentOthers = percentOthers, nExactModel_M = nExactModel_M)
})
resSampleStats <- data.frame(do.call(rbind, resSampleStats))
resSampleStats <- tidyr::separate(resSampleStats, N_pTrash_rel_R2_lin_inter_fitInter, 
                                  into = c("N", "pTrash", "rel", "R2", "lin", "inter", "fitInter"), sep = "_")
resSampleStats <- tidyr::unite(resSampleStats, "lin_inter", c(lin, inter), sep = "_")
# resSampleStats <- tidyr::pivot_longer(resSampleStats, !c(N, pTrash, rel, R2, lin_inter), 
#                                       names_to = "measure", values_to = "values")
#resSampleStats$values <- as.numeric(resSampleStats$values)
resSampleStats$N <- factor(resSampleStats$N, levels = setParam$dgp$N)
resSampleStats$pTrash <- factor(resSampleStats$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
resSampleStats$rel <- factor(resSampleStats$rel, levels = setParam$dgp$reliability)
# resSampleStats$measure <- factor(resSampleStats$measure, 
#                                  levels = c("nAllLin", "nAllInter", "nAllEffects", 
#                                             "nLin_M", "nInter_M",
#                                             "nOthers_M", "percentOthers", "nExactModel_M"))

resSampleStatsLin <- resSampleStats[,c("N", "pTrash", "rel", "R2", "lin_inter", "fitInter", "nAllLin")]

# change type of columns or specific entry details to prepare plotting  
str(resSampleStatsLin)

names(resSampleStatsLin)[names(resSampleStatsLin) == 'nAllLin'] <- "allLin"

pviAllLinCount$fitInter <- rep(2, dim(pviAllLinCount)[1])

allLinData <- rbind(pviAllLinCount, resSampleStatsLin)
allLinData$fitInter <- plyr::mapvalues(allLinData$fitInter, 
                                       from=c(0,1,2), 
                                       to=c("ENET - ohne","ENET - mit","GBM"))
unique(allLinData$fitInter)
allLinData$fit <- ifelse(allLinData$fitInter != "GBM", "enet", "gbm")
allLinData$allLin <- as.numeric(allLinData$allLin)

plotNAllLin <- ggplot(allLinData,
       aes(x = interaction(pTrash, N, sep = " x "), y = allLin, 
           group = interaction(R2, fitInter, fit), colour = R2,
           linetype = fitInter, shape = fit)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_shape_manual(values = c(16, 8)) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(rel ~ lin_inter, scales = "free_y", labeller = label_both, switch = "y") +
  ylab(paste0("n samples out of ", setParam$dgp$nTrain, " samples")) +
  xlab("pTrash (decreasing) x N (increasing)") +
  ggtitle("effect recovery across samples") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/pCount_nAllLin_ENTvsGBM.eps"),
#                 plot = plotNAllLin,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/pCount_nAllLin_ENTvsGBM.png"),
#                 plot = plotNAllLin,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# measure: in wie vielen der 100 test samples wurden alle lineare Effekte gefunden?
# für starke Interaktionen (= Spalte auf der linken Seite) kann GBM bei kleiner 
#   sample size und hoher Anzahl von Trash-Variablen + einer Reliabilität von 0.8 
#   und 1 besser lineare Effekte finden als das ENET
#   -> dass GBM nicht zwischen Interaktion und linearem Effekt trennt wird zur 
#     Tugend?

################################################################################
# interaction strength measure 
################################################################################
#### only GBM
# pull data from nested list of all results (fullData)
interGBM <- rbindResults(fullDataGBM, "interStrength")
head(interGBM)

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
interGBM <- idx2info(interGBM, condN_pTrashGBM, type = "gbm")

# Wie häufig wurden simulierte/wahre Interaktionseffekte gefunden?
# wie gut kann GBM Interaktionen von kontinuierlichen Variablen abbilden?
idxInter <- which(interGBM$feature %in% c(setParam$dgp$interEffects))
interEffectsGBM <- interGBM[idxInter,]

# interaction > 0 == H-statistic (für genauer Interpretation der H-Statistic in 
#     Paper schauen)
interEffectsGBM <- interEffectsGBM[which(interEffectsGBM$interaction > 0),]

# count how often all linear effects were extracted per condition
interEffectsGBM <- tidyr::unite(interEffectsGBM, "N_pTrash_rel_R2_lin_inter", 
                              c(N, pTrash, rel, R2, lin_inter),
                              sep = "_", remove = FALSE) 
interEffectsGBM_list <- split(interEffectsGBM, 
                            f = interEffectsGBM$N_pTrash_rel_R2_lin_inter)

counter <- lapply(seq_along(interEffectsGBM_list), function(iSublist) {
  nInter <- aggregate(interaction ~ sample, 
                      data = interEffectsGBM_list[[iSublist]], length)
  cbind(N_pTrash_rel_R2_lin_inter = names(interEffectsGBM_list)[iSublist],
        sample = nInter$sample,
        nInter = nInter$interaction, 
        allInter = ifelse(nInter$interaction == length(setParam$dgp$interEffects), 1, 0))
})  

allInterGBM <- data.frame(do.call(rbind, counter))

allInterCount <- aggregate(allInter ~ N_pTrash_rel_R2_lin_inter, 
                            data = allInterGBM[which(allInterGBM$allInter == 1),], length)

allInterCount <- tidyr::separate(allInterCount, N_pTrash_rel_R2_lin_inter, 
                                  into = c("N", "pTrash", "rel", "R2", "lin", "inter"), sep = "_")
allInterCount <- tidyr::unite(allInterCount, "lin_inter", c(lin, inter), sep = "_")

allInterCount$N <- factor(allInterCount$N, levels = setParam$dgp$N)
allInterCount$pTrash <- factor(allInterCount$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
allInterCount$rel <- factor(allInterCount$rel, levels = setParam$dgp$rel)

colValues <- c("green3", "darkblue", "darkmagenta")

ggplot(allInterCount,
       aes(x = interaction(pTrash, N, sep = " x "), y = allInter, 
           group = R2, colour = R2)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(rel ~ lin_inter, scales = "free_y", labeller = label_both, switch = "y") +
  ylab(paste0("n samples out of ", setParam$dgp$nTrain, " samples")) +
  xlab("pTrash (decreasing) x N (increasing)") +
  ggtitle("effect recovery across samples") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

#### for elastic net data
# werden alle Interaktionseffekte extrahiert?
# use stats from plotting linear effects
#     -> already calculated number of samples with all interaction effects in each condition
resSampleStatsInter <- resSampleStats[which(resSampleStats$fitInter == 1), # only for ENET data with interactions!
                                      c("N", "pTrash", "rel", "R2", "lin_inter", "fitInter", "nAllInter")]

# change type of columns or specific entry details to prepare plotting  
str(resSampleStatsInter)

names(resSampleStatsInter)[names(resSampleStatsInter) == 'nAllInter'] <- "allInter"
resSampleStatsInter$fitInter <- NULL

allInterCount$fit <- rep("GBM", dim(allInterCount)[1])
resSampleStatsInter$fit <- rep("ENET", dim(resSampleStatsInter)[1])

allInterData <- rbind(allInterCount, resSampleStatsInter)

unique(allInterData$fit)
allInterData$allInter <- as.numeric(allInterData$allInter)

plotNAllInter <- ggplot(allInterData,
                      aes(x = interaction(pTrash, N, sep = " x "), y = allInter, 
                          group = interaction(R2, fit), colour = R2,
                          linetype = fit, shape = fit)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  scale_shape_manual(values = c(16, 8)) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(rel ~ lin_inter, scales = "free_y", labeller = label_both, switch = "y") +
  ylab(paste0("n samples out of ", setParam$dgp$nTrain, " samples")) +
  xlab("pTrash (decreasing) x N (increasing)") +
  ggtitle("effect recovery across samples") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/pCount_nAllInter_ENTvsGBM.eps"),
#                 plot = plotNAllInter,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/pCount_nAllInter_ENTvsGBM.png"),
#                 plot = plotNAllInter,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

################################################################################
# plot relative bias in estimated coefficients 
################################################################################
# pull data from nested list of all results (fullData)
estBeta <- rbindResults(fullData, "estBeta")
estBeta$Var <- rownames(estBeta)

# variables from rownames to own column to work woth variable information
estBeta$Var <- stringr::str_replace(estBeta$Var, "\\.[:digit:]{1,}$", "")
estBeta$Var <- stringr::str_replace(estBeta$Var, "\\.", ":")

# only keep variables with simulated effects 
idxEstBeta <- which(estBeta$Var %in% c(setParam$dgp$linEffects, setParam$dgp$interEffects))
# # check if idx really catches all effects
# length(idxEstBeta) == length(setParam$dgp$condLabels) *
#   length(condN_pTrash) * (length(setParam$dgp$linEffects) + length(setParam$dgp$interEffects))

estB_effects <- estBeta[idxEstBeta,]
estB_noEffects <- estBeta[-idxEstBeta,]

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter)
estB_effects <- idx2info(estB_effects)

# change type of columns or specific entry details to prepare plotting  
str(estB_effects)
estB_effects$N <- factor(estB_effects$N, levels = setParam$dgp$N)
estB_effects$pTrash <- factor(estB_effects$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
estB_effects$rel <- factor(estB_effects$rel, levels = setParam$dgp$rel)

estB_effects$effectType <- ifelse(estB_effects$Var %in% 
                                    unique(estB_effects$Var)[stringr::str_detect(unique(estB_effects$Var), ":")], "inter", "lin")

estB_effects$effectSize <- ifelse(estB_effects$effectType == "lin", 
                                  stringr::str_sub(estB_effects$lin_inter, start = 1L, end = 3L),
                                  stringr::str_sub(estB_effects$lin_inter, start = 5L, end = -1L))

# add true simulated effects to the data frame to be able to calculate biases in estimated coefficients
trueB_lin <- data.frame(effectType = rep("lin", times = length(setParam$dgp$Rsquared)),
                        R2 = setParam$dgp$Rsquared, 
                        setParam$dgp$trueEffects$lin)
trueB_inter <- data.frame(effectType = rep("inter", times = length(setParam$dgp$Rsquared)),
                          R2 = setParam$dgp$Rsquared, 
                          setParam$dgp$trueEffects$inter)
trueB <- rbind(trueB_lin, trueB_inter)
trueB <- tidyr::pivot_longer(trueB, cols = !c(R2, effectType), 
                             names_to = "effectSize", values_to = "trueB")

trueB$effectSize <- stringr::str_sub(trueB$effectSize, start = 2L)

trueB$R2 <- as.character(trueB$R2) # identical variable types to merge correctly
estB_effects <- merge(estB_effects, trueB, by = c("R2", "effectSize", "effectType"))

# calculate absolute & relative bias
# absolute bias biases interpretation of deviation in estimated coefficients due to different sizes of 
#     simulated coefficients depending on R2 budget size
estB_effects$absBiasB <- estB_effects$M - estB_effects$trueB
estB_effects$relBiasB <- (estB_effects$M - estB_effects$trueB) /  estB_effects$trueB

estB_effects$lin_inter <- factor(estB_effects$lin_inter, 
                                 levels = c("0.2_0.8", "0.5_0.5", "0.8_0.2"))
linInterTypes <- levels(estB_effects$lin_inter)
linIntLabels <- paste0("lin = ", 
                       stringr::str_sub(linInterTypes, start = 1L, end = 3L),
                       " + inter = ", stringr::str_sub(linInterTypes, start = -3L, end = -1L))
estB_effects$lin_inter <- factor(estB_effects$lin_inter, 
                                 levels = c("0.2_0.8", "0.5_0.5", "0.8_0.2"),
                                 labels = linIntLabels)

range(estB_effects$relBiasB)
# plot relative bias for every variable with simulated effects separately 
colValues <- c("green3", "darkcyan", "darkblue", "darkmagenta")

plotRelBias <- function(data) { # plot function
  ggplot(data,
         aes(x = interaction(pTrash, N, sep = " x "), y = relBiasB, group = R2, colour = R2)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = relBiasB - SE, ymax = relBiasB + SE), width=.2) +
    scale_y_continuous(limits = c(-1, 0.4), breaks = seq(-1, 0.4, 0.2)) +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(rel ~ lin_inter) +
    ylab("relative bias in est. \u03B2") +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle(paste0("relativ bias for ", iVar)) +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20))
}

# generate plot for every variable with simulated effects
for (iVar in c(setParam$dgp$linEffects, setParam$dgp$interEffects)) {
  subset <- estB_effects[which(estB_effects$Var == iVar), ]
  tmp_p <- plotRelBias(subset)
  
  plotName <- paste0("pRelBias_", sub(':','.', iVar))
  assign(plotName, tmp_p)
}


# # save all plots as files
# plotNames <- ls(pattern = "^pRelBias")
# for (iPlot in plotNames) {
#   pName <- paste0(iPlot, "_testResults")
# 
#   ggplot2::ggsave(filename = paste0(plotFolder, "/", pName, ".eps"),
#                     plot = eval(parse(text = iPlot)),
#                     device = cairo_ps,
#                     dpi = 300,
#                     width = 14.39,
#                     height = 10.83,
#                     units = "in")
# 
#   ggplot2::ggsave(filename = paste0(plotFolder, "/", pName, ".png"),
#                     plot = eval(parse(text = iPlot)),
#                     width = 14.39,
#                     height = 10.83,
#                     units = "in")
# }

################################################################################
# mean false positive coefficient
################################################################################
# for true negatives the coefficient is NA
# dependent measure for each simulated condition and separately for linear term, 
#     polynomial terms and interaction terms

# remove coefficients with simulated effects
estB_noEffects <- estBeta[-idxEstBeta,]

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter)
estB_noEffects <- idx2info(estB_noEffects)

# change type of columns or specific entry details to prepare plotting  
str(estB_effects)
estB_effects$N <- factor(estB_effects$N, levels = setParam$dgp$N)
estB_effects$pTrash <- factor(estB_effects$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
estB_effects$rel <- factor(estB_effects$rel, levels = setParam$dgp$rel)

# separate linear terms, polynomial terms and interaction terms
estB_noEffects$varType <- dplyr::case_when(
  stringr::str_detect(estB_noEffects$Var, "[:digit:]{1,}:Var") ~ "inter",
  stringr::str_detect(estB_noEffects$Var, "^poly:Var") ~ "poly",
  stringr::str_detect(estB_noEffects$Var, "^Var[:digit:]{1,}") ~ "lin")

# split - apply - combine
# split data in different simulated conditions N x pTrash x rel x R2 x lin_inter
estB_noEffects <- tidyr::unite(estB_noEffects, "N_pTrash_rel_R2_lin_inter", 
                               c(N, pTrash, rel, R2, lin_inter),
                               sep = "_", remove = FALSE) 
estB_noEffects_sub <- split(estB_noEffects, 
                            f = estB_noEffects$N_pTrash_rel_R2_lin_inter)

estB_falsePositive <- lapply(seq_along(estB_noEffects_sub), function(iSublist) {
  tmp <- t(aggregate(M ~ varType, 
                     data = estB_noEffects_sub[[iSublist]], range))
  # meanCoef <- tmp[2,]
  # names(meanCoef) <- tmp[1,]
  # as.numeric(meanCoef)
})

estB_falsePositive <- do.call(rbind, estB_falsePositive)



################################################################################
# plot selected variables
################################################################################
# pull data from nested list of all results (fullData)
selectedVars <- rbindResults(fullData, "selectSample")
selectedVars$sample <- rep(seq_len(setParam$dgp$nTrain), 
                           times = length(setParam$dgp$condLabels) * nrow(condGrid))

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
selectedVars <- idx2info(selectedVars)

# get number of unique predictors depending on pTrash
names(setParam$dgp$nModelPredictors) <- setParam$dgp$pTrash
selectedVars$nModelPredictors <- setParam$dgp$nModelPredictors[match(selectedVars$pTrash, 
                                                                     names(setParam$dgp$nModelPredictors))]

# samples that exactly recovered the simulated model (all simulated effects & no other predictors)
selectedVars$exactModel <- ifelse((selectedVars$all.T1F0 == 1) & 
                                    (selectedVars$nOthers == 0), 1, 0)

idxExactModel <- which(selectedVars$exactModel == 1)
selectedVars[idxExactModel,]

# split, apply, combine
# split data in different simulated conditions N x pTrash x R2 x lin_inter
selectedVars <- tidyr::unite(selectedVars, "N_pTrash_rel_R2_lin_inter", 
                             c(N, pTrash, rel, R2, lin_inter),
                             sep = "_", remove = FALSE) 
subLists <- split(selectedVars, 
                  f = selectedVars$N_pTrash_rel_R2_lin_inter)

# identical computations for every subList (i.e., simulated condition)
# compute dependent measures to visualize in results plot
resSampleStats <- lapply(seq_along(subLists), function(iCond) {
  nLin_M = mean(subLists[[iCond]]$nLin, na.rm = T)
  nInter_M = mean(subLists[[iCond]]$nInter, na.rm = T)
  nAllLin = sum(subLists[[iCond]]$nLin == length(setParam$dgp$linEffects))
  nAllInter = sum(subLists[[iCond]]$nInter == length(setParam$dgp$interEffects))
  nAllEffects = sum(subLists[[iCond]]$all.T1F0)
  nOthers_M = mean(subLists[[iCond]]$nOthers, na.rm = T)
  nModelPredictors <- unique(subLists[[iCond]]$nModelPredictors)
  # to do: this ignores that the true predictors are not "others"; but true predictors are constant across conditions
  percentOthers = nOthers_M / nModelPredictors * 100 
  nExactModel_M = sum(subLists[[iCond]]$exactModel)
  cbind(N_pTrash_rel_R2_lin_inter = names(subLists)[iCond], 
        nLin_M = nLin_M, nInter_M = nInter_M,
        nAllLin = nAllLin, nAllInter = nAllInter, nAllEffects = nAllEffects,
        nOthers_M = nOthers_M, percentOthers = percentOthers, nExactModel_M = nExactModel_M)
})
resSampleStats <- data.frame(do.call(rbind, resSampleStats))
resSampleStats <- tidyr::separate(resSampleStats, N_pTrash_rel_R2_lin_inter, 
                                  into = c("N", "pTrash", "rel", "R2", "lin", "inter"), sep = "_")
resSampleStats <- tidyr::unite(resSampleStats, "lin_inter", c(lin, inter), sep = "_")
resSampleStats <- tidyr::pivot_longer(resSampleStats, !c(N, pTrash, rel, R2, lin_inter), 
                                      names_to = "measure", values_to = "values")

# change type of columns or specific entry details to prepare plotting  
str(resSampleStats)
resSampleStats$values <- as.numeric(resSampleStats$values)
resSampleStats$N <- factor(resSampleStats$N, levels = setParam$dgp$N)
resSampleStats$pTrash <- factor(resSampleStats$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
resSampleStats$rel <- factor(resSampleStats$rel, levels = setParam$dgp$reliability)
resSampleStats$measure <- factor(resSampleStats$measure, 
                                 levels = c("nAllLin", "nAllInter", "nAllEffects", 
                                            "nLin_M", "nInter_M",
                                            "nOthers_M", "percentOthers", "nExactModel_M"))

# plot number of samples per simulated condition that recovers ...
#     ... all linear effects
#     ... all interaction effects
#     ... all linear and interaction effects

colValues <- c("green3", "darkcyan", "darkblue", "darkmagenta")

countPlot <- function(data) {
  ggplot(data,
         aes(x = interaction(pTrash, N, sep = " x "), y = values, 
             group = R2, colour = R2)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(rel ~ lin_inter, scales = "free_y", labeller = label_both, switch = "y") +
    ylab(paste0("n samples out of ", setParam$dgp$nTrain, " samples")) +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle("effect recovery across samples") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15))
}

countEffectType <- c("nAllLin", "nAllInter", "nAllEffects", "nExactModel_M")

for (iEffect in countEffectType) {
  subset <- resSampleStats[which(resSampleStats$measure %in% iEffect), ]
  tmp_p <- countPlot(subset)
  
  plotName <- paste0("pCount_", iEffect)
  assign(plotName, tmp_p)
}

# # save plots as files
# plotNames <- ls(pattern = "^pCount")
# for (iPlot in plotNames) {
#   pName <- paste0(iPlot)
# 
#   ggplot2::ggsave(filename = paste0(plotFolder, "/", pName, ".eps"),
#                     plot = eval(parse(text = iPlot)),
#                     device = cairo_ps,
#                     dpi = 300,
#                     width = 14.39,
#                     height = 10.83,
#                     units = "in")
# 
#   ggplot2::ggsave(filename = paste0(plotFolder, "/", pName, ".png"),
#                     plot = eval(parse(text = iPlot)),
#                     width = 14.39,
#                     height = 10.83,
#                     units = "in")
# }


# nAllEffects vs. nExactModel_M 
# barplot with number of models with all effects as transparent bars in the background
# with number of samples with extraction of exact models on top
pAllvsExactEffects <- ggplot(resSampleStats[which(resSampleStats$measure %in% c("nAllEffects")), ],
                             aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                                 group = R2, fill = R2)) +
  geom_col(position = "dodge", alpha = 0.4) +
  geom_bar(data = resSampleStats[which(resSampleStats$measure %in% c("nExactModel_M")), ],
           aes(x = interaction(pTrash, N, sep = " x "), y = values, 
               group = R2, fill = R2), stat = "identity", position = "dodge") +
  scale_fill_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(rel ~ lin_inter) +
  # facet_grid(measure ~ lin_inter, scales = "free_y", labeller = label_both, switch = "y") +
  ylab(paste0("n samples out of ", setParam$dgp$nTrain, " samples")) +
  xlab("pTrash (decreasing) x N (increasing)") +
  ggtitle("all effects recovered(transparent) vs. only simulated effects") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                        colour = "lightgrey"), 
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "lightgrey"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/nAllvsExactEffectsBarplot.eps"),
#                 plot = pAllvsExactEffects,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 15.33,
#                 height = 10.92,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/nAllvsExactEffectsBarplot.png"),
#                 plot = pAllvsExactEffects,
#                 width = 15.33,
#                 height = 10.92,
#                 units = "in")

# plot number of variables in model without simulated true effects
## in percent (nChosen/pTrash)*100
pPercentOthers <- ggplot(resSampleStats[which(resSampleStats$measure %in% c("percentOthers")), ],
                         aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                             group = R2, colour = R2)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(rel ~ lin_inter) +
  ylab(paste0("mean percentage of selected trash variables")) +
  xlab("pTrash (decreasing) x N (increasing)") +
  ggtitle("variables without simulated effect in model recovery") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/pPercentOthers.eps"),
#                 plot = pPercentOthers,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.95,
#                 height = 11.51,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/pPercentOthers.png"),
#                 plot = pPercentOthers,
#                 width = 13.95,
#                 height = 11.51,
#                 units = "in")

################################################################################
# plot hyperparameters (alpha & lambda)
################################################################################
hyperParamData <- tidyr::pivot_longer(selectedVars, cols = c(alpha, lambda), 
                                      names_to = "measures", values_to = "values")
hyperParamData$values <- factor(hyperParamData$values)

hyperParamData <- tidyr::unite(hyperParamData, "pTrashxN", c(pTrash, N), 
                               sep = " x ", remove = F)
pTrashxN_levels <- expand.grid(sort(setParam$dgp$pTrash, decreasing = T), setParam$dgp$N)
pTrashxN_levels <- as.vector(unlist(tidyr::unite(pTrashxN_levels, "", c(Var1, Var2), sep = " x ")))
hyperParamData$pTrashxN <- factor(hyperParamData$pTrashxN, levels = pTrashxN_levels)

# plot alpha & lambda
# turn factor back to numeric! 
hyperParamData$values <- as.numeric(as.character(hyperParamData$values))

# range of lambda (in case of warm start lambda grid is based on data)
range(hyperParamData[which(hyperParamData$measures == "lambda"), "values"])

tuneParamVec <- c("alpha", "lambda")
binwidthParam <- c(0.1, 0.5)

for (iParam in seq_along(tuneParamVec)) {
  for (iRel in setParam$dgp$reliability) {
    tmp_p <- ggplot(hyperParamData[which(hyperParamData$measures == tuneParamVec[iParam] &
                                           hyperParamData$rel == iRel), 
                                   c("sample", "pTrashxN", "R2", "lin_inter", "values")],
                    aes(x = values, group = R2, fill = R2)) +
      geom_histogram(stat = "bin", binwidth = binwidthParam[iParam], position = position_dodge2(width = 0.9, preserve = "single")) +
      scale_fill_manual(values = colValues) +
      scale_y_continuous(limits = c(0, setParam$dgp$nTrain), 
                         breaks = seq(0, setParam$dgp$nTrain, by = 30)) +
      facet_grid(pTrashxN ~ lin_inter) +
      ylab(paste0("n samples out of ", setParam$dgp$nTrain, " samples")) +
      xlab(paste0(tuneParamVec[iParam], " (tuning result)")) +
      ggtitle(tuneParamVec[iParam], " parameter value in (binned absolute) frequencies") +
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            linewidth = 0.5, linetype = "solid"),
            panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                            colour = "lightgrey"), 
            panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                            colour = "lightgrey"),
            axis.text.y = element_text(size = 20),
            axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            strip.text.x = element_text(size = 15))
    
    plotName <- paste0("p", tuneParamVec[iParam], "_", iRel)
    assign(plotName, tmp_p)
  }
}

# # save plots as files
# plotNames <- ls(pattern = "^palpha")
# plotNames <- c(plotNames, ls(pattern = "^plambda"))
# for (iPlot in plotNames) {
#   ggplot2::ggsave(filename = paste0(plotFolder, "/", iPlot, ".eps"),
#                     plot = eval(parse(text = iPlot)),
#                     device = cairo_ps,
#                     dpi = 300,
#                     width = 16.27,
#                     height = 11.38,
#                     units = "in")
# 
#   ggplot2::ggsave(filename = paste0(plotFolder, "/", iPlot, ".png"),
#                     plot = eval(parse(text = iPlot)),
#                     width = 16.27,
#                     height = 11.38,
#                     units = "in")
# }


