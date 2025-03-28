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

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

# load results files
# resFolder <- "results/finalResults/dependentMeasures" 
resFolder <- "results/dependentMeasures" 

listDir <- dir(resFolder)
# grid <- "mitNpred_500trees"
grid <- "ohneNpred_1000trees"
hyperRF <- loadRData(paste0(resFolder, "/hyperParametersSample_RF_", grid, ".rda"))

# pull data from nested list of all results (fullData)
hyperRF <- rbindSingleResults(hyperRF)

colnames(hyperRF)

hyperRF <- idx2infoNew(hyperRF)
hyperRF$sample <- rep(1:100, times = 6*9) 

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model", "sample")
hyperRF[col2fac] <- lapply(hyperRF[col2fac], factor)
str(hyperRF)

################################################################################
# plot count data for hyper parameter values (how often were they chosen?)
################################################################################
library(tidyverse)

plotHyperBars <- function(data, x, measureLabel) {
  ggplot(data, aes(x = factor(x), y = n, fill = pTrash)) +
    geom_col(position = position_dodge(width = 0.7)) +
    facet_grid2(N + lin_inter ~ R2,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesR2[1], 0.4)),
                                      element_rect(fill = alpha(colValuesR2[2], 0.4)),
                                      element_rect(fill = alpha(colValuesR2[3], 0.4))))) +
    guides(fill = "none") + 
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

(pMtryBar <- plotHyperBars(mTryCount, mTryCount$mtry, "mtry"))

ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_mTryCount_rel06_", grid,".png"),
                plot = pMtryBar,
                width = 17.78,
                height = 9.5,
                units = "in")

minNodeCount <- hyperRF %>% 
  group_by(N, pTrash, rel, R2, lin_inter) %>% 
  count(minNode) %>% 
  ungroup()
minNodeCount$N <- factor(minNodeCount$N, levels = c(100, 300, 1000))
minNodeCount$minNode <- factor(minNodeCount$minNode, levels = c(5, 10, 20))

(pMinNodeBar <- plotHyperBars(minNodeCount, minNodeCount$minNode, "min obs in endnode"))

ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_minNodeCount_rel06_", grid,".png"),
                plot = pMinNodeBar,
                width = 17.78,
                height = 9.5,
                units = "in")

splitRuleCount <- hyperRF %>% 
  group_by(N, pTrash, rel, R2, lin_inter) %>% 
  count(splitRule) %>% 
  ungroup()
splitRuleCount$N <- factor(splitRuleCount$N, levels = c(100, 300, 1000))

(pSplitruleBar <- plotHyperBars(splitRuleCount, splitRuleCount$splitRule, "splitrule"))

ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_splitRukeCount_rel06_", grid,".png"),
                plot = pSplitruleBar,
                width = 17.78,
                height = 9.5,
                units = "in")
