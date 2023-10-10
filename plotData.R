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

# save all plots as files
plotNames <- ls(pattern = "^pRelBias")
for (iPlot in plotNames) {
  pName <- paste0(iPlot, "_testResults")
    
  ggplot2::ggsave(filename = paste0(plotFolder, "/", pName, ".eps"),
                    plot = eval(parse(text = iPlot)),
                    device = cairo_ps,
                    dpi = 300,
                    width = 11.3,
                    height = 8.22,
                    units = "in")
    
  ggplot2::ggsave(filename = paste0(plotFolder, "/", pName, ".png"),
                    plot = eval(parse(text = iPlot)),
                    width = 11.3,
                    height = 8.22,
                    units = "in") 
}

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
        axis.title.y = element_text(size = 20))

# save plots as files
ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTest.eps"),
                plot = pPerformTrainVStest,
                device = cairo_ps,
                dpi = 300,
                width = 13.95,
                height = 11.51,
                units = "in")

ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTest.png"),
                plot = pPerformTrainVStest,
                width = 13.95,
                height = 11.51,
                units = "in") 
