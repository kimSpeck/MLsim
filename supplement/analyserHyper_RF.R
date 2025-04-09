# load parameters and helper functions 
source("utils/setParameters.R")
source("utils/analysisTools.R")

# # ANOVA
# library(afex) # für aov_ez()
# library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²

# plot results
library(ggplot2)
library(ggh4x)
library(patchwork)
# library(see)

# colValues <- c("green3", "darkblue", "darkmagenta")
colValuesR2 <- c('#db4a07', '#850c0c', '#3c1518')
colValuesInter <- c('#050440', '#181ff2', '#0eb2e8')
colValuesLin <- c('#0eb2e8', '#181ff2', '#050440')

colValuesRel <- c("#81C784", "#388E3C", "#1B5E20")

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

data <- "nonlinear"

# load results files
# resFolder <- "results/finalResults/dependentMeasures" 
# resFolder <- "results/dependentMeasures" 
# resFolder <- "results/inter/oldDataRF/dependentMeasures" 
resFolder <- paste0("results/", data, "/dependentMeasures")

listDir <- dir(resFolder)
hyperRF <- loadRData(paste0(resFolder, "/hyperParametersSample_RF.rda"))

# pull data from nested list of all results (fullData)
hyperRF <- rbindSingleResults(hyperRF)

colnames(hyperRF)

hyperRF <- idx2infoNew(hyperRF)

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model", "sample")
hyperRF[col2fac] <- lapply(hyperRF[col2fac], factor)
str(hyperRF)

################################################################################
# plot count data for hyper parameter values (how often were they chosen?)
################################################################################
library(tidyverse)

plotHyperBars <- function(data, x, measureLabel) {
  # ggplot(data, aes(x = factor(x), y = n, fill = pTrash)) +
  ggplot(data, aes(x = factor(x), y = n, fill = rel)) +
    geom_col(position = position_dodge(width = 0.7)) +
    facet_grid2(N + lin_inter ~ R2,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesR2[1], 0.4)),
                                      element_rect(fill = alpha(colValuesR2[2], 0.4)),
                                      element_rect(fill = alpha(colValuesR2[3], 0.4))))) +
    #guides(fill = "none") + 
    scale_fill_manual(values = colValuesRel) +
    xlab(measureLabel) +
    ylab("Count per level") +
    theme(panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
          panel.background = element_rect(color = "white", fill = "white"),
          plot.title = element_text(size = 30, face = "bold"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.position = c(.2, .85), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal")
}


mTryCount <- hyperRF %>% 
  group_by(N, pTrash, rel, R2, lin_inter) %>% 
  count(mtry) %>% 
  ungroup()
mTryCount$N <- factor(mTryCount$N, levels = c(100, 300, 1000))
mTryCount$mtry <- factor(mTryCount$mtry, levels = unique(mTryCount$mtry))

mTryCount_pTrash10 <- mTryCount[mTryCount$pTrash == 10,]
(pMtryBar_pTrash10 <- plotHyperBars(mTryCount_pTrash10, mTryCount_pTrash10$mtry, "mtry"))
mTryCount_pTrash50 <- mTryCount[mTryCount$pTrash == 50,]
(pMtryBar_pTrash50 <- plotHyperBars(mTryCount_pTrash50, mTryCount_pTrash50$mtry, "mtry"))

# ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_mTryCount_pTrash10_", data, ".png"),
#                 plot = pMtryBar_pTrash10,
#                 width = 17.78,
#                 height = 9.5,
#                 units = "in")
# ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_mTryCount_pTrash50_ ", data, ".png"),
#                 plot = pMtryBar_pTrash50,
#                 width = 17.78,
#                 height = 9.5,
#                 units = "in")

minNodeCount <- hyperRF %>% 
  group_by(N, pTrash, rel, R2, lin_inter) %>% 
  count(minNode) %>% 
  ungroup()
minNodeCount$N <- factor(minNodeCount$N, levels = c(100, 300, 1000))
minNodeCount$minNode <- factor(minNodeCount$minNode, levels = c(5, 10, 20))

minNodeCount_pTrash10 <- minNodeCount[minNodeCount$pTrash == 10,]
(pMinNodeBar_pTrash10 <- plotHyperBars(minNodeCount_pTrash10, minNodeCount_pTrash10$minNode, "min obs in endnode"))
minNodeCount_pTrash50 <- minNodeCount[minNodeCount$pTrash == 50,]
(pMinNodeBar_pTrash50 <- plotHyperBars(minNodeCount_pTrash50, minNodeCount_pTrash50$minNode, "min obs in endnode"))

# ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_minNodeCount_pTrash10.png"),
#                 plot = pMinNodeBar_pTrash10,
#                 width = 17.78,
#                 height = 9.5,
#                 units = "in")
# ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_minNodeCount_pTrash50.png"),
#                 plot = pMinNodeBar_pTrash50,
#                 width = 17.78,
#                 height = 9.5,
#                 units = "in")

splitRuleCount <- hyperRF %>% 
  group_by(N, pTrash, rel, R2, lin_inter) %>% 
  count(splitRule) %>% 
  ungroup()
splitRuleCount$N <- factor(splitRuleCount$N, levels = c(100, 300, 1000))
splitRuleCount$splitRule <- factor(splitRuleCount$splitRule, levels = c("variance", "extratrees"))

splitRuleCount_pTrash10 <- splitRuleCount[splitRuleCount$pTrash == 10,]
(pSplitruleBar_pTrash10 <- plotHyperBars(splitRuleCount_pTrash10, splitRuleCount_pTrash10$splitRule, "splitrule"))
splitRuleCount_pTrash50 <- splitRuleCount[splitRuleCount$pTrash == 50,]
(pSplitruleBar_pTrash50 <- plotHyperBars(splitRuleCount_pTrash50, splitRuleCount_pTrash50$splitRule, "splitrule"))

# ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_splitRuleCount_pTrash10.png"),
#                 plot = pSplitruleBar_pTrash10,
#                 width = 17.78,
#                 height = 9.5,
#                 units = "in")
# ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_splitRuleCount_pTrash50.png"),
#                 plot = pSplitruleBar_pTrash50,
#                 width = 17.78,
#                 height = 9.5,
#                 units = "in")
