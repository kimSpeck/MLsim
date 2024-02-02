library(ggplot2)

source("setParameters.R")
source("analysisTools.R")

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

condGrid <- expand.grid(N = setParam$dgp$N,
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability,
                        factors = c(TRUE, FALSE))

# duplicate data with factors to run model on indicator data, too
condGrid$indicators <- rep(FALSE, dim(condGrid)[1])
addIndi <- condGrid[which(condGrid$factors == TRUE), ]
addIndi$indicators <- TRUE
condGrid <- rbind(addIndi, condGrid)

condN_pTrash <- paste0("N", condGrid$N, 
                       "_pTrash", condGrid$pTrash,
                       "_rel", condGrid$reliability,
                       "_f", ifelse(condGrid$factors, 1, 0),
                       "_ind", ifelse(condGrid$indicators, 1, 0))

# results data with no factor, 1 indicator; no factor, 5 indicators; factor, 5 indicators
# resFolder <- paste0("results/resultsWithoutInter_Indicators") # without interactions
resFolder <- paste0("results/resultsWithInter_Indicators") # with interactions
load(paste0(resFolder, "/fullData.rda"))

################################################################################
# plot train and test performance
################################################################################
# pull data from nested list of all results (fullData)
performanceTrain <- rbindResults(fullData, "performTrainStats")
performanceTest <- rbindResults(fullData, "performTestStats")

# variables from rownames to own column to work woth variable information
performanceTrain$measure <- rownames(performanceTrain)
performanceTrain$measure <- stringr::str_replace(performanceTrain$measure, "\\.[:digit:]{1,}$", "")
performanceTrain$measure <- paste0(performanceTrain$measure, "_train")

performanceTest$measure <- rownames(performanceTest)
performanceTest$measure <- stringr::str_replace(performanceTest$measure, "\\.[:digit:]{1,}$", "")
performanceTest$measure[which(performanceTest$measure == "Rsq_test")] <- "Rsquared_test"

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
performanceTrain <- idx2info(performanceTrain)
performanceTest <- idx2info(performanceTest)

# change type of columns or specific entry details to prepare plotting  
performanceTrain$N <- factor(performanceTrain$N, levels = setParam$dgp$N)
performanceTrain$pTrash <- factor(performanceTrain$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))

performanceTest$N <- factor(performanceTest$N, levels = setParam$dgp$N)
performanceTest$pTrash <- factor(performanceTest$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))

# merge performance Train and performance Test
performanceStats <- rbind(performanceTrain, performanceTest)
performanceStats <- tidyr::separate(performanceStats, measure, c("measure", "trainTest"), sep = "_")

performanceStats$factorIndicator <- paste0("f", performanceStats$factor, 
                                           "_ind", performanceStats$indicators) 

# plot overfit instead of train as train - test
performanceStats <- tidyr::pivot_wider(performanceStats, names_from = trainTest, values_from = c(M, SE, SD))
performanceStats$overfit <- performanceStats$M_train - performanceStats$M_test
performanceSub <- performanceStats[which(performanceStats$measure != "MAE" &
                                           performanceStats$measure != "RMSE"),]
performanceSub[,c("M_train", "SE_train", "SE_test", "SD_train", "SD_test")] <- list(NULL)
performanceSub <- tidyr::pivot_longer(performanceSub, c(M_test, overfit),
                                      names_to = "measures", values_to = "values")

# any(performanceStats[which(performanceStats$measure != "MAE" & 
#                              performanceStats$measure != "RMSE"), "overfit_test"] < 0)
# idx <- which(performanceStats[which(performanceStats$measure != "MAE" & 
#                              performanceStats$measure != "RMSE"), "overfit_test"] < 0)
# performanceStats[which(performanceStats$measure != "MAE" & 
#                          performanceStats$measure != "RMSE"), ][idx, ]

# plot performance measures for train and test data
# colValues <- c("green3", "darkcyan", "darkblue", "darkmagenta")
colValues <- c("green3", "darkblue", "darkmagenta")

(pPerformTrainVStest <- ggplot(performanceSub,
                               aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                                   group = interaction(R2, factorIndicator), colour = R2,
                                   linetype = factorIndicator, shape = factorIndicator)) +
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTest.eps"),
#                 plot = pPerformTrainVStest,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTest.png"),
#                 plot = pPerformTrainVStest,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

################################################################################
# plot selected variables
################################################################################
# werden alle linearen/Interaktionseffekte extrahiert?
# pull data from nested list of all results (fullData)
selectedVars <- rbindResults(fullData, "selectSample")
selectedVars$sample <- rep(seq_len(setParam$dgp$nTrain), 
                           times = length(setParam$dgp$condLabels) * nrow(condGrid))

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
selectedVars <- idx2info(selectedVars)

# get number of unique predictors depending on pTrash
names(setParam$dgp$nModelPredictors) <- setParam$dgp$pTrash
names(setParam$dgp$nModelPredictorsIndicators) <- setParam$dgp$pTrash

# split - apply - combine
#     separates data according to indicator values
#     {0 = factor/1 var, 1 = 5 indicators/factor} (in this exact order)
subLists <- split(selectedVars, 
                  f = selectedVars$indicators)

# check if data structure allows certain transformation
#   either data in nInd & nIndInter (if indicators = 1) OR in nLin and nInter (if indicators = 0) 
any(subLists[["0"]][,"nInd"] != 0) 
any(subLists[["1"]][,"nLin"] != 0) 
any(subLists[["1"]][,"nInter"] != 0) 
any(subLists[["1"]][,"all.T1F0"] != 0)

# identical computations for every subList (i.e., with and without indicators)
#     number of predictors without and with single indicators
nPredictors <- list(setParam$dgp$nModelPredictors, setParam$dgp$nModelPredictorsIndicators)
# compute dependent measures to visualize in results plot
resSampleStats <- lapply(seq_along(subLists), function(iSub) {
  # get number of predictors (different depending on indicators in model or not)
  subLists[[iSub]]$nModelPredictors <- nPredictors[[iSub]][match(
    subLists[[iSub]]$pTrash, names(nPredictors[[iSub]]))]
  
  # copy information from nInd & nIndInter to nLin & nInter (if indicators = 1)
  if (all(subLists[[iSub]]$indicators == 1)) {
    subLists[[iSub]]$nLin <- subLists[[iSub]]$nInd
    subLists[[iSub]]$nInter <- subLists[[iSub]]$nIndInter
    
    # !!! there is something wrong with nOthersInd - the count of nOthersInd is the 
    #     same as the total number of predictors! 
    #     nOthersInd = nOthers; thus, for the moment subtract nInd and nIndInter!
    #     change this later!!!!
    subLists[[iSub]]$nOthersInd <- subLists[[iSub]]$nOthersInd - subLists[[iSub]]$nInd - subLists[[iSub]]$nIndInter
    subLists[[iSub]]$nOthers <- subLists[[iSub]]$nOthersInd
    subLists[[iSub]]$all.T1F0 <- subLists[[iSub]]$all.Ind
  }
  
  subLists[[iSub]]
})

# here!
selectedVars <- data.frame(do.call(rbind, resSampleStats))
selectedVars[,c("nInd", "nIndInter", "nOthersInd", "all.Ind")]<- list(NULL)

# samples that exactly recovered the simulated model (all simulated effects & no other predictors)
selectedVars$exactModel <- ifelse((selectedVars$all.T1F0 == 1 & selectedVars$nOthers == 0) , 1, 0)

idxExactModel <- which(selectedVars$exactModel == 1)
selectedVars[idxExactModel,]

# split, apply, combine
# split data in different simulated conditions N x pTrash x R2 x lin_inter
selectedVars <- tidyr::unite(selectedVars, "N_pTrash_rel_f_ind_R2_lin_inter", 
                             c(N, pTrash, rel, factor, indicators, R2, lin_inter),
                             sep = "_", remove = FALSE) 
subLists <- split(selectedVars, 
                  f = selectedVars$N_pTrash_rel_f_ind_R2_lin_inter)

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
  cbind(N_pTrash_rel_f_ind_R2_lin_inter = names(subLists)[iCond], 
        nLin_M = nLin_M, nInter_M = nInter_M,
        nAllLin = nAllLin, nAllInter = nAllInter, nAllEffects = nAllEffects,
        nOthers_M = nOthers_M, percentOthers = percentOthers, nExactModel_M = nExactModel_M)
})
resSampleStats <- data.frame(do.call(rbind, resSampleStats))
resSampleStats <- tidyr::separate(resSampleStats, N_pTrash_rel_f_ind_R2_lin_inter, 
                                  into = c("N", "pTrash", "rel", "factor", "indicators", "R2", "lin", "inter"), sep = "_")
resSampleStats <- tidyr::unite(resSampleStats, "lin_inter", c(lin, inter), sep = "_")
resSampleStats <- tidyr::pivot_longer(resSampleStats, !c(N, pTrash, rel, factor, indicators, R2, lin_inter), 
                                      names_to = "measure", values_to = "values")

# change type of columns or specific entry details to prepare plotting  
str(resSampleStats)
resSampleStats$values <- as.numeric(resSampleStats$values)
resSampleStats$N <- factor(resSampleStats$N, levels = setParam$dgp$N)
resSampleStats$pTrash <- factor(resSampleStats$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
resSampleStats$rel <- factor(resSampleStats$rel, levels = setParam$dgp$reliability)
resSampleStats$factor <- factor(resSampleStats$factor)
resSampleStats$measure <- factor(resSampleStats$measure, 
                                 levels = c("nAllLin", "nAllInter", "nAllEffects", 
                                            "nLin_M", "nInter_M",
                                            "nOthers_M", "percentOthers", "nExactModel_M"))

resSampleStats$factorIndicator <- paste0("f", resSampleStats$factor, 
                                         "_ind", resSampleStats$indicators) 

# plot number of samples per simulated condition that recovers ...
#     ... all linear effects
#     ... all interaction effects
#     ... all linear and interaction effects

colValues <- c("green3", "darkblue", "darkmagenta")

countPlot <- function(data) {
  ggplot(data,
         aes(x = interaction(pTrash, N, sep = " x "), y = values, 
             group = interaction(R2, factorIndicator), colour = R2,
             linetype = factorIndicator, shape = factorIndicator)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = colValues) +
    scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
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

################################################################################
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
  facet_grid(rel ~ lin_inter + factorIndicator, labeller = label_both) +
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
#                 # with factors
#                 width = 19.96,
#                 height = 13.63,
#                 # # without factors
#                 # width = 15.33,
#                 # height = 10.92,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/nAllvsExactEffectsBarplot.png"),
#                 plot = pAllvsExactEffects,
#                 # with factors
#                 width = 19.96,
#                 height = 13.63,
#                 # # without factors
#                 # width = 15.33,
#                 # height = 10.92,
#                 units = "in")


# plot number of variables in model without simulated true effects
## in percent (nChosen/pTrash)*100
pPercentOthers <- ggplot(resSampleStats[which(resSampleStats$measure %in% c("percentOthers")), ],
                         aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                             group = interaction(R2, factorIndicator), colour = R2,
                             linetype = factorIndicator, shape = factorIndicator)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
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

# bisher: sind alle Effekte da oder nicht? sehr strenges Kriterium
#   - mindestens ein Indikator pro Faktor?
#   - wie viele Indikatoren pro Faktor?

estBetaFull <- rbindResults(fullData, "estBetaFull")
estBeta <- rbindResults(fullData, "estBeta")
estBetaFull$idxCondLabel <- estBeta$idxCondLabel
rm(estBeta)

head(estBetaFull)
estBetaFull <- idx2info(estBetaFull)

estBetaIndicator <- estBetaFull[which(estBetaFull$indicators == 1),]

# in wie vielen Samples ist mindestens 1 Indikator pro Faktor?
# in wie viele Indikatoren pro Faktor im Schnitt über alle Samples?
idxEstBetaCols <- stringr::str_which(colnames(estBetaIndicator), "^s0")
idxOtherCols <- colnames(estBetaIndicator)[-idxEstBetaCols]
estBetaIndicator <- cbind(estBetaIndicator[,idxEstBetaCols] != 0, estBetaIndicator[, idxOtherCols])

# variables from rownames to own column to work woth variable information
estBetaIndicator$Var <- rownames(estBetaIndicator)
estBetaIndicator$Var <- stringr::str_replace(estBetaIndicator$Var, "\\.[:digit:]{1,}$", "")
estBetaIndicator$Var <- stringr::str_replace(estBetaIndicator$Var, "\\.", ":")

# only keep variables with simulated effects 
idxEstBeta <- which(estBetaIndicator$Var %in% c(setParam$dgp$indEffects, setParam$dgp$indInterEffects))
estBetaInd_effects <- estBetaIndicator[idxEstBeta,]

estBetaInd_effects$effectType <- ifelse(estBetaInd_effects$Var %in% 
                                   unique(estBetaInd_effects$Var)[stringr::str_detect(unique(estBetaInd_effects$Var), ":")], "inter", "lin")

estBetaInd_effects$factorNum <- ifelse(estBetaInd_effects$effectType == "lin",
                                       stringr::str_sub(estBetaInd_effects$Var, start = 1L, end = 2L),
                                       paste0(stringr::str_sub(estBetaInd_effects$Var, start = 1L, end = 2L), 
                                              stringr::str_sub(estBetaInd_effects$Var, start = 5L, end = 7L)))

# split apply combine
# splitte data frames nach:
#   N, pTrash, rel, factor, R2, lin_inter, factorNum
#   indicator == 1 for all; effectType in factorNum

estBetaInd_effects <- tidyr::unite(estBetaInd_effects, "N_pTrash_rel_f_R2_lin_inter_fN", 
                             c(N, pTrash, rel, factor, R2, lin_inter, factorNum),
                             sep = "_", remove = FALSE) 
subLists <- split(estBetaInd_effects, 
                  f = estBetaInd_effects$N_pTrash_rel_f_R2_lin_inter_fN)

resSampleStats <- lapply(seq_along(subLists), function(iCond) {
  idxEstBetaCols <- stringr::str_which(colnames(subLists[[iCond]]), "^s0")
  indiCount_perSample <- colSums(subLists[[iCond]][,idxEstBetaCols], na.rm = TRUE)
  cbind(N_pTrash_rel_f_R2_lin_inter_fN = names(subLists)[iCond], 
        atLeast1 = sum(indiCount_perSample > 0)/setParam$dgp$nTrain,
        # counting does not work for interactions anymore (too many)
        # exactly1 = sum(indiCount_perSample == 1)/setParam$dgp$nTrain,
        # exactly2 = sum(indiCount_perSample == 2)/setParam$dgp$nTrain,
        # exactly3 = sum(indiCount_perSample == 3)/setParam$dgp$nTrain,
        # exactly4 = sum(indiCount_perSample == 4)/setParam$dgp$nTrain,
        # exactly5 = sum(indiCount_perSample == 5)/setParam$dgp$nTrain,
        meanIndiCount = mean(indiCount_perSample))
})

indEffectsCounter <- data.frame(do.call(rbind, resSampleStats))

indEffectsCounter <- tidyr::separate(indEffectsCounter, N_pTrash_rel_f_R2_lin_inter_fN, 
                                  into = c("N", "pTrash", "rel", "factor", "R2", "lin", "inter", "factorNum"), sep = "_")
indEffectsCounter <- tidyr::unite(indEffectsCounter, "lin_inter", c(lin, inter), sep = "_")

indEffectsCounter$effectType <- ifelse(indEffectsCounter$factorNum %in% 
                                          unique(indEffectsCounter$factorNum)[stringr::str_detect(unique(indEffectsCounter$factorNum), ":")], "inter", "lin")

library(tidyverse)

indEffectsCounter$atLeast1 <- as.numeric(indEffectsCounter$atLeast1)
indEffectsCounter$meanIndiCount <- as.numeric(indEffectsCounter$meanIndiCount)

indEffectsCounterM <- indEffectsCounter %>% 
  group_by(N, pTrash, rel, factor, R2, lin_inter, effectType) %>% 
  summarise(atLeast1M = mean(atLeast1, na.rm = T),
            meanIndiCountM = mean(meanIndiCount, na.rm = T))

indEffectsCounterM <- tidyr::pivot_longer(indEffectsCounterM, !c(N, pTrash, rel, factor, R2, lin_inter, effectType), 
                                         names_to = "measure", values_to = "values")

pAtLeast1 <- ggplot(indEffectsCounterM[which(indEffectsCounterM$measure == "atLeast1M"),],
       aes(x = interaction(pTrash, N, sep = " x "), y = values, 
           group = interaction(R2, effectType), colour = R2,
           linetype = effectType, shape = effectType)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  scale_linetype_manual(values = c("solid", "dashed")) +
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

pMeanIndicator_lin <- ggplot(indEffectsCounterM[which(indEffectsCounterM$measure == "meanIndiCountM" &
                                                    indEffectsCounterM$effectType == "lin"),],
       aes(x = interaction(pTrash, N, sep = " x "), y = values, 
           group = interaction(R2, effectType), colour = R2,
           linetype = effectType, shape = effectType)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  scale_linetype_manual(values = c("solid", "dashed")) +
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

pMeanIndicator_inter <- ggplot(indEffectsCounterM[which(indEffectsCounterM$measure == "meanIndiCountM" &
                                                        indEffectsCounterM$effectType == "inter"),],
                             aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                                 group = interaction(R2, effectType), colour = R2,
                                 linetype = effectType, shape = effectType)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  scale_linetype_manual(values = c("solid", "dashed")) +
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

# to do: plots
# proportion of selected Items auf der y-Achse
# Sortierung auf x-Achse fixen

# überprüfen, ob grouping effect der ridge regression zum Ergebnis führt oder nicht?
# Korrelation zwischen Items der Faktoren reduzieren? 
#  pro: alle Items der Faktoren werden in den meisten Samples aufgenommen
 
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
idxEstBeta <- which(estBeta$Var %in% c(setParam$dgp$linEffects, setParam$dgp$interEffects,
                                       setParam$dgp$indEffects, setParam$dgp$indInterEffects))
# # check if idx really catches all effects
# length(idxEstBeta) == length(setParam$dgp$condLabels) *
#   length(condN_pTrash) * (length(setParam$dgp$linEffects) + length(setParam$dgp$interEffects))


# 
estB_full <- idx2info(estBeta)

# # split - apply - combine
# estB_full <- tidyr::unite(estB_full, "N_pTrash_rel_f_R2_lin_inter", 
#                           c(N, pTrash, rel, factor, R2, lin_inter),
#                           sep = "_", remove = FALSE)
# estB_full_sub <- split(estB_full,
#                        f = estB_full$N_pTrash_rel_f_R2_lin_inter)
# 
# estB_full_arranged <- lapply(seq_along(estB_full_sub), function(iSublist) {
#   estB_full_sub[[iSublist]][order(-abs(estB_full_sub[[iSublist]]$M)),]
# })

estB_effects <- estBeta[idxEstBeta,]
estB_noEffects <- estBeta[-idxEstBeta,]

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter)
estB_effects <- idx2info(estB_effects)

# change type of columns or specific entry details to prepare plotting  
str(estB_effects)
estB_effects$N <- factor(estB_effects$N, levels = setParam$dgp$N)
estB_effects$pTrash <- factor(estB_effects$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
estB_effects$rel <- factor(estB_effects$rel, levels = setParam$dgp$rel)
estB_effects$factor <- factor(estB_effects$factor)
estB_effects$indicators <- factor(estB_effects$indicators)

estB_effects$factorIndicator <- paste0("f", estB_effects$factor, 
                                       "_ind", estB_effects$indicators)

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
# colValues <- c("green3", "darkcyan", "darkblue", "darkmagenta") # with 4 R^2 values
colValues <- c("green3", "darkblue", "darkmagenta") 

plotRelBias <- function(data) { # plot function
  ggplot(data,
         aes(x = interaction(pTrash, N, sep = " x "), y = relBiasB, 
             group = interaction(R2, factorIndicator), colour = R2,
             linetype = factorIndicator, shape = factorIndicator)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = relBiasB - SE, ymax = relBiasB + SE), width=.2) +
    scale_y_continuous(limits = c(-1, 0.65), breaks = seq(-1, 0.6, 0.2)) +
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


for (iVar in c(setParam$dgp$indEffects, setParam$dgp$indInterEffect)) {
  subset <- estB_effects[which(estB_effects$Var %in% c("F1_1", "F1_2", "F1_3", "F1_4", "F1_5")), ]
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
str(estB_noEffects)
estB_noEffects$N <- factor(estB_noEffects$N, levels = setParam$dgp$N)
estB_noEffects$pTrash <- factor(estB_noEffects$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
estB_noEffects$rel <- factor(estB_noEffects$rel, levels = setParam$dgp$rel)
estB_noEffects$factor <- factor(estB_noEffects$factor)

# separate linear terms, polynomial terms and interaction terms
estB_noEffects$varType <- dplyr::case_when(
  stringr::str_detect(estB_noEffects$Var, "[:digit:]{1,}:Var") ~ "inter",
  stringr::str_detect(estB_noEffects$Var, "^poly:Var") ~ "poly",
  stringr::str_detect(estB_noEffects$Var, "^Var[:digit:]{1,}") ~ "lin")

falsePosBoxplot <- function(data, varType) {
  ggplot(data[which(data$varType == varType),],
         aes(x = interaction(pTrash, N, sep = " x "), y = M, 
             group = interaction(pTrash, N, R2), colour = R2)) +
    geom_boxplot(outlier.colour="red", outlier.shape=1,
                 outlier.size=1, alpha = 0.5, linewidth = 0.25) +
    # geom_violin(trim=TRUE) +
    stat_summary(fun.y=mean, geom="point", shape=8, size=2) +
    scale_y_continuous(limits = c(-0.1, 0.1), breaks = seq(-0.1, 0.1, 0.05)) +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(rel ~ lin_inter + factor, labeller = label_both) +
    ylab(paste0("mean bias in est. \u03B2 across all 'false positive' ", varType, " predictors")) +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle(paste0("mean bias for 'false positive' ", varType, " predictors")) +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20))
  
}

varTypes <- unique(estB_noEffects$varType)
for (iVar in seq_along(varTypes)){
  tmp_p <- falsePosBoxplot(estB_noEffects, varTypes[iVar])
  
  plotName <- paste0("pMeanBias_falsePositive", varTypes[iVar])
  assign(plotName, tmp_p)
  
}

# # check
# # the number of missing values approximately matches the warning message regarding
# #   removed rows for each variable type
# sum(is.na(estB_noEffects[which(estB_noEffects$varType == "lin"), "M"]))
# sum(is.na(estB_noEffects[which(estB_noEffects$varType == "inter"), "M"]))
# sum(is.na(estB_noEffects[which(estB_noEffects$varType == "poly"), "M"]))

# # save all plots as files
# plotNames <- ls(pattern = "^pMeanBias_falsePositive")
# for (iPlot in plotNames) {
# 
#   ggplot2::ggsave(filename = paste0(plotFolder, "/", iPlot, ".eps"),
#                   plot = eval(parse(text = iPlot)),
#                   device = cairo_ps,
#                   dpi = 300,
#                   # with factors
#                   width = 19.96,
#                   height = 13.63,
#                   # # without factors
#                   # width = 15.33,
#                   # height = 10.92,
#                   units = "in")
# 
#   ggplot2::ggsave(filename = paste0(plotFolder, "/", iPlot, ".png"),
#                   plot = eval(parse(text = iPlot)),
#                   # with factors
#                   width = 19.96,
#                   height = 13.63,
#                   # # without factors
#                   # width = 15.33,
#                   # height = 10.92,
#                   units = "in")
# }

# # calculate values already illustrated by boxplots
# # split - apply - combine
# # split data in different simulated conditions N x pTrash x rel x R2 x lin_inter
# estB_noEffects <- tidyr::unite(estB_noEffects, "N_pTrash_rel_f_R2_lin_inter", 
#                                c(N, pTrash, rel, factor, R2, lin_inter),
#                                sep = "_", remove = FALSE) 
# estB_noEffects_sub <- split(estB_noEffects, 
#                             f = estB_noEffects$N_pTrash_rel_f_R2_lin_inter)
# 
# estB_falsePositive <- lapply(seq_along(estB_noEffects_sub), function(iSublist) {
#   # tmp_range <- t(aggregate(M ~ varType, 
#   #                    data = estB_noEffects_sub[[iSublist]], range))
#   
#   # use quantile instead of range (range included as 0 and 100)
#   tmp_quant <- t(aggregate(M ~ varType, 
#               data = estB_noEffects_sub[[iSublist]], quantile))
#   tmp_mean <- t(aggregate(M ~ varType, 
#               data = estB_noEffects_sub[[iSublist]], mean))
#   
#   # build zero matrix to fill up 
#   # in some conditions no additional linear predictor is chosen, therefore the number 
#   #   of columns does not match across simulated conditions
#   tmp <- matrix(0, ncol = 3, nrow = 6)
#   colnames(tmp) <- c("inter", "lin", "poly")
#   
#   # fill up zero matrix with actual values 
#   tmp[1, ] <- tmp_mean[2, match(colnames(tmp), tmp_mean[1,])]
#   tmp[2:nrow(tmp), ] <- tmp_quant[2:nrow(tmp_quant), match(colnames(tmp), tmp_quant[1,])]
# 
#   # replace NAs in tmp with zeros 
#   #   -> if no linear effect was wrongly chosen the bias is 0 not NA
#   tmp <- ifelse(is.na(tmp), 0, tmp)
# 
#   # # checks & debugging
#   # if (ncol(tmp_quant) != 3) {
#   #   print("! ncol differs from expected number of columns")
#   #   print(tmp)
#   # }
#   
#   rownames(tmp) <- c("mean", "M0", "M25", "M50", "M75", "M100")
#   cbind(N_pTrash_rel_f_R2_lin_inter = names(subLists)[iSublist], tmp)
# })
# 
# estB_falsePositive <- data.frame(do.call(rbind, estB_falsePositive))
# estB_falsePositive$measure <- rownames(estB_falsePositive)
# estB_falsePositive$measure <- stringr::str_replace(estB_falsePositive$measure, "\\.[:digit:]{1,}$", "")
# 
# estB_falsePositive <- tidyr::separate(estB_falsePositive, N_pTrash_rel_f_R2_lin_inter, 
#                                       into = c("N", "pTrash", "rel", "factor", "R2", "linPer", "interPer"), sep = "_")

