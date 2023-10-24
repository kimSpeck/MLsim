# to do: add measures in percent (now: duplicates due to nTrain = 100)

library(ggplot2)

source("setParameters.R")
source("analysisTools.R")

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

condGrid <- expand.grid(N = setParam$dgp$N,
                        pTrash = setParam$dgp$pTrash)
condN_pTrash <- paste0("N", condGrid$N, "_pTrash", condGrid$pTrash)

resFolder <- "results/resultsServer"
load(paste0(resFolder, "/fullData.rda"))

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

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter)
estB_effects <- idx2info(estB_effects)

# change type of columns or specific entry details to prepare plotting  
str(estB_effects)
estB_effects$N <- factor(estB_effects$N, levels = setParam$dgp$N)
estB_effects$pTrash <- factor(estB_effects$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))

estB_effects$effectType <- ifelse(estB_effects$Var %in% 
                                    unique(estB_effects$Var)[stringr::str_detect(unique(estB_effects$Var), ":")], "inter", "lin")

# add true simulated effects to the data frame to be able to calculate biases in estimated coefficients
trueB <- data.frame(R2 = setParam$dgp$Rsquared, setParam$dgp$trueEffects)
trueB <- tidyr::pivot_longer(trueB, cols = !R2, names_to = "lin", values_to = "trueB_lin")
trueB$lin <- stringr::str_sub(trueB$lin, start = 2L)
trueB$lin <- as.numeric(trueB$lin)
trueB$inter <- (1 - trueB$lin)

trueB2 <- data.frame(R2 = setParam$dgp$Rsquared, setParam$dgp$trueEffects)
trueB2 <- tidyr::pivot_longer(trueB2, cols = !R2, names_to = "inter", values_to = "trueB_inter")
trueB2$inter <- stringr::str_sub(trueB2$inter, start = 2L)

trueB <- merge(trueB, trueB2, by = c("R2", "inter"))
trueB <- tidyr::unite(trueB, "lin_inter", c(lin, inter), sep = "_")
colnames(trueB) <- c("R2", "lin_inter", "lin", "inter")
trueB <- tidyr::pivot_longer(trueB, cols = c(lin, inter), names_to = "effectType", values_to = "trueB")

trueB$R2 <- as.character(trueB$R2) # identical variable types to merge correctly
estB_effects <- merge(estB_effects, trueB, by = c("R2", "lin_inter", "effectType"))

# calculate absolute & relative bias
# absolute bias biases interpretation of deviation in estimated coefficients due to different sizes of 
#     simulated coefficients depending on R2 budget size
estB_effects$absBiasB <- estB_effects$M - estB_effects$trueB
estB_effects$relBiasB <- (estB_effects$M - estB_effects$trueB) /  estB_effects$trueB


linInterTypes <- unique(estB_effects$lin_inter)
linIntLabels <- paste0("lin = ", 
                       stringr::str_sub(linInterTypes, start = 1L, end = 3L),
                       " + inter = ", stringr::str_sub(linInterTypes, start = -3L, end = -1L))
estB_effects$lin_inter <- factor(estB_effects$lin_inter, 
                                 levels = c("0.2_0.8", "0.5_0.5", "0.8_0.2"),
                                 labels = linIntLabels)

# to do:
#   - pTrash/

# plot relative bias for every variable with simulated effects separately 
colValues <- c("green3", "darkcyan", "darkblue", "darkmagenta")

plotRelBias <- function(data) { # plot function
  ggplot(data,
         aes(x = interaction(pTrash, N, sep = " x "), y = relBiasB, group = R2, colour = R2)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = relBiasB - SE, ymax = relBiasB + SE), width=.2) +
    scale_y_continuous(limits = c(-1, 0), breaks = seq(-1, 0, 0.2)) +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_wrap(~ lin_inter) +
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
#                     width = 11.3,
#                     height = 8.22,
#                     units = "in")
#     
#   ggplot2::ggsave(filename = paste0(plotFolder, "/", pName, ".png"),
#                     plot = eval(parse(text = iPlot)),
#                     width = 11.3,
#                     height = 8.22,
#                     units = "in") 
# }

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

# plot performance measures for train and test data
colValues <- c("green3", "darkcyan", "darkblue", "darkmagenta")

pPerformTrainVStest <- ggplot(performanceStats[which(performanceStats$measure != "MAE"),],
       aes(x = interaction(pTrash, N, sep = " x "), y = M, 
           group = R2, colour = R2)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(measure + trainTest ~ lin_inter, scales = "free_y", labeller = label_both, switch = "y") +
  ylab("") +
  xlab("pTrash (decreasing) x N (increasing)") +
  ggtitle("Training vs. Test performance") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTest.eps"),
#                 plot = pPerformTrainVStest,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.95,
#                 height = 11.51,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTest.png"),
#                 plot = pPerformTrainVStest,
#                 width = 13.95,
#                 height = 11.51,
#                 units = "in") 

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
selectedVars <- tidyr::unite(selectedVars, "N_pTrash_R2_lin_inter", 
                             c(N, pTrash, R2, lin_inter),
                             sep = "_", remove = FALSE) 
subLists <- split(selectedVars, 
                  f = selectedVars$N_pTrash_R2_lin_inter)

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
  cbind(N_pTrash_R2_lin_inter = names(subLists)[iCond], 
        nLin_M = nLin_M, nInter_M = nInter_M,
        nAllLin = nAllLin, nAllInter = nAllInter, nAllEffects = nAllEffects,
        nOthers_M = nOthers_M, percentOthers = percentOthers, nExactModel_M = nExactModel_M)
})
resSampleStats <- data.frame(do.call(rbind, resSampleStats))
resSampleStats <- tidyr::separate(resSampleStats, N_pTrash_R2_lin_inter, 
                                  into = c("N", "pTrash", "R2", "lin", "inter"), sep = "_")
resSampleStats <- tidyr::unite(resSampleStats, "lin_inter", c(lin, inter), sep = "_")
resSampleStats <- tidyr::pivot_longer(resSampleStats, !c(N, pTrash, R2, lin_inter), 
                                      names_to = "measure", values_to = "values")

# change type of columns or specific entry details to prepare plotting  
str(resSampleStats)
resSampleStats$values <- as.numeric(resSampleStats$values)
resSampleStats$N <- factor(resSampleStats$N, levels = setParam$dgp$N)
resSampleStats$pTrash <- factor(resSampleStats$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
resSampleStats$measure <- factor(resSampleStats$measure, 
                                 levels = c("nAllLin", "nAllInter", "nAllEffects", 
                                            "nLin_M", "nInter_M",
                                            "nOthers_M", "percentOthers", "nExactModel_M"))

# plot number of samples per simulated condition that recovers ...
#     ... all linear effects
#     ... all interaction effects
#     ... all linear and interaction effects

colValues <- c("green3", "darkcyan", "darkblue", "darkmagenta")

pNsamplesEffect <- ggplot(resSampleStats[which(resSampleStats$measure %in% c("nAllLin", "nAllInter", "nAllEffects")), ],
       aes(x = interaction(pTrash, N, sep = " x "), y = values, 
           group = R2, colour = R2)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(measure ~ lin_inter, scales = "free_y", labeller = label_both, switch = "y") +
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/nSamplesWithEffects.eps"),
#                 plot = pNsamplesEffect,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.95,
#                 height = 11.51,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/nSamplesWithEffects.png"),
#                 plot = pNsamplesEffect,
#                 width = 13.95,
#                 height = 11.51,
#                 units = "in") 

# nAllEffects vs. nExactModel_M 
pAllvsExactEffects <- ggplot(resSampleStats[which(resSampleStats$measure %in% c("nAllEffects")), ],
       aes(x = interaction(pTrash, N, sep = " x "), y = values, 
           group = R2, fill = R2)) +
  geom_col(position = "dodge", alpha = 0.4) +
  geom_bar(data = resSampleStats[which(resSampleStats$measure %in% c("nExactModel_M")), ],
           aes(x = interaction(pTrash, N, sep = " x "), y = values, 
               group = R2, fill = R2), stat = "identity", position = "dodge") +
  scale_fill_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(~ lin_inter) +
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
#                 width = 24.24,
#                 height = 13.28,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/nAllvsExactEffectsBarplot.png"),
#                 plot = pAllvsExactEffects,
#                 width = 24.24,
#                 height = 13.28,
#                 units = "in") 
pAllvsExactEffectsLine <- ggplot(resSampleStats[which(resSampleStats$measure %in% c("nAllEffects", "nExactModel_M")), ],
                          aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                              group = R2, colour = R2)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(measure ~ lin_inter, scales = "free_y", labeller = label_both, switch = "y") +
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/nAllvsExactEffectsLineplot.eps"),
#                 plot = pAllvsExactEffectsLine,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.95,
#                 height = 11.51,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/nAllvsExactEffectsLineplot.png"),
#                 plot = pAllvsExactEffectsLine,
#                 width = 13.95,
#                 height = 11.51,
#                 units = "in") 

resSampleStats$measure <- factor(resSampleStats$measure, 
                                 levels = c("nExactModel_M", "percentOthers", "nAllEffects", 
                                            "nAllLin", "nAllInter", "nLin_M", "nInter_M",
                                            "nOthers_M"))

pExactModel <- ggplot(resSampleStats[which(resSampleStats$measure %in% c("nAllEffects", "nExactModel_M", "percentOthers")), ],
                                 aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                                     group = R2, colour = R2)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(measure ~ lin_inter) +
  ylab(paste0("n samples out of ", setParam$dgp$nTrain, " samples")) +
  xlab("pTrash (decreasing) x N (increasing)") +
  ggtitle("exact model recovery (other variables vs. all effect variables)") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/nExactModelAllvsOthers.eps"),
#                 plot = pExactModel,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.95,
#                 height = 11.51,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/nExactModelAllvsOthers.png"),
#                 plot = pExactModel,
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

pAlpha <- ggplot(hyperParamData[which(hyperParamData$measures == "alpha"), 
                      c("sample", "pTrashxN", "R2", "lin_inter", "values")],
       aes(x = values, group = R2, fill = R2)) +
  geom_histogram(stat = "count", position = position_dodge2(width = 0.9, preserve = "single")) +
  scale_fill_manual(values = colValues) +
  scale_y_continuous(limits = c(0, setParam$dgp$nTrain), 
                     breaks = seq(0, setParam$dgp$nTrain, by = 50)) +
  facet_grid(pTrashxN ~ lin_inter) +
  ylab(paste0("n samples out of ", setParam$dgp$nTrain, " samples")) +
  xlab(paste0("alpha (tuning result)")) +
  ggtitle("alpha parameter value in (absolute) frequencies") +
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/plotTunedAlpha.eps"),
#                 plot = pAlpha,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 16.27,
#                 height = 11.38,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/plotTunedAlpha.png"),
#                 plot = pAlpha,
#                 width = 16.27,
#                 height = 11.38,
#                 units = "in")

