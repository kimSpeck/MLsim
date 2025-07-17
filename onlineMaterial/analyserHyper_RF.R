# load parameters and helper functions 
source("utils/setParameters.R")
source("utils/analysisTools.R")

# plot results
library(tidyverse)
library(ggplot2)
library(ggh4x)
library(patchwork)

colValuesR2 <- c('#db4a07', '#850c0c', '#3c1518')
colValuesInter <- c('#050440', '#181ff2', '#0eb2e8')
colValuesLin <- c('#0eb2e8', '#181ff2', '#050440')

colValuesRel <- c("#999999", "#444444", "#222222")
colValuesN <-  c("#81C784", "#388E3C", "#1B5E20")

plotFolder <- "plots/hyperParameter"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

################################################################################
# plot utility functions
################################################################################
plotHyperBars <- function(data, x, measureLabel) {
  ggplot(data, aes(x = factor(x), y = n, fill = rel)) +
    geom_col(position = position_dodge(width = 0.7)) +
    scale_y_continuous(limits = c(-100, 1000), breaks = seq(0, 1000, 500)) +
    facet_grid2(facetY ~ N,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesN[1], 0.7)),
                                      element_rect(fill = alpha(colValuesN[2], 0.7)),
                                      element_rect(fill = alpha(colValuesN[3], 0.7))),
                  background_y = list(element_rect(fill = alpha(colValuesR2[1], 0.5)),
                                      element_rect(fill = alpha(colValuesR2[1], 0.5)),
                                      element_rect(fill = alpha(colValuesR2[2], 0.5)),
                                      element_rect(fill = alpha(colValuesR2[2], 0.5)),
                                      element_rect(fill = alpha(colValuesR2[3], 0.5)),
                                      element_rect(fill = alpha(colValuesR2[3], 0.5))),
                  text_y = elem_list_text(angle = 0))) +
    #guides(fill = "none") + 
    scale_fill_manual(name = "Reliability", values = colValuesRel) +
    geom_line(aes(y = 0, color = R2), linewidth = 8, alpha = 0.5) +
    scale_color_manual(name = "R²", values = c("0.2" = colValuesR2[1], 
                                               "0.5" = colValuesR2[2], 
                                               "0.8" = colValuesR2[3])) +
    xlab(measureLabel) +
    ylab("Count per level") +
    theme(panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
          panel.background = element_rect(color = "white", fill = "white"),
          plot.title = element_text(size = 30, face = "bold"),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          # legend.position = c(.2, .85), 
          legend.position = "bottom", 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal")
}

################################################################################
# load data
################################################################################
# load results files
dgpVec <- c("inter", "nonlinear3", "pwlinear")
resFolderVec <- paste0("results/", dgpVec, "/dependentMeasures")


listDir <- dir(resFolderVec)
listDir <- stringr::str_subset(listDir, "hyperParametersSample_[:alnum:]+_RF.rda")

################################################################################
# plot hyper parameters
################################################################################
for (iDGP in seq_along(listDir)){
  
  hyperRF <- loadRData(paste0("results/", dgpVec[iDGP], "/dependentMeasures/", listDir[iDGP]))
  
  # pull data from nested list of all results (fullData)
  hyperRF <- rbindSingleResults(hyperRF)
  
  colnames(hyperRF)
  
  hyperRF <- idx2infoNew(hyperRF)
  
  # change variables to factors
  col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model", "sample")
  hyperRF[col2fac] <- lapply(hyperRF[col2fac], factor)
  # str(hyperRF)
  
  ################################################################################
  # plot count data for hyper parameter values (how often were they chosen?)
  ################################################################################
  hyperRF <- hyperRF[which(hyperRF$lin_inter != "0.5_0.5"),]
  
  mTryCount <- hyperRF %>% 
    group_by(N, pTrash, rel, R2, lin_inter) %>% 
    count(mtry) %>% 
    ungroup()
  mTryCount$N <- factor(mTryCount$N, levels = c(100, 300, 1000),
                        labels = c("N = 100", "N = 300", "N = 1000"))
  mTryCount$mtry <- factor(mTryCount$mtry, levels = unique(mTryCount$mtry))
  mTryCount$facetY <- interaction(mTryCount$lin_inter, mTryCount$R2, sep = "x")
  mTryCount$facetY <- factor(mTryCount$facetY, labels = 
                                    c("R² = 0.2\nlinear 20%", "R² = 0.2\nlinear 80%", 
                                      "R² = 0.5\nlinear 20%", "R² = 0.5\nlinear 80%",
                                      "R² = 0.8\nlinear 20%", "R² = 0.8\nlinear 80%"))
  
  
  mTryCount_pTrash10 <- mTryCount[mTryCount$pTrash == 10,]
  pMtryBar_pTrash10 <- plotHyperBars(mTryCount_pTrash10, mTryCount_pTrash10$mtry, "# Predictors for Split (mtry)")
  mTryCount_pTrash50 <- mTryCount[mTryCount$pTrash == 50,]
  pMtryBar_pTrash50 <- plotHyperBars(mTryCount_pTrash50, mTryCount_pTrash50$mtry, "# Predictors for Split (mtry)")
  
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_mTryCount_pTrash10_", dgpVec[iDGP], ".png"),
                  plot = pMtryBar_pTrash10,
                  width = 17.78,
                  height = 9.5,
                  units = "in")
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_mTryCount_pTrash50_", dgpVec[iDGP], ".png"),
                  plot = pMtryBar_pTrash50,
                  width = 17.78,
                  height = 9.5,
                  units = "in")
  
  minNodeCount <- hyperRF %>% 
    group_by(N, pTrash, rel, R2, lin_inter) %>% 
    count(minNode) %>% 
    ungroup()
  minNodeCount$N <- factor(minNodeCount$N, levels = c(100, 300, 1000),
                           labels = c("N = 100", "N = 300", "N = 1000"))
  minNodeCount$minNode <- factor(minNodeCount$minNode, levels = c(5, 10, 20))
  minNodeCount$facetY <- interaction(minNodeCount$lin_inter, minNodeCount$R2, sep = "x")
  minNodeCount$facetY <- factor(minNodeCount$facetY, labels = 
                               c("R² = 0.2\nlinear 20%", "R² = 0.2\nlinear 80%", 
                                 "R² = 0.5\nlinear 20%", "R² = 0.5\nlinear 80%",
                                 "R² = 0.8\nlinear 20%", "R² = 0.8\nlinear 80%"))
  
  minNodeCount_pTrash10 <- minNodeCount[minNodeCount$pTrash == 10,]
  pMinNodeBar_pTrash10 <- plotHyperBars(minNodeCount_pTrash10, minNodeCount_pTrash10$minNode, "min. End Node Size")
  minNodeCount_pTrash50 <- minNodeCount[minNodeCount$pTrash == 50,]
  pMinNodeBar_pTrash50 <- plotHyperBars(minNodeCount_pTrash50, minNodeCount_pTrash50$minNode, "min. End Node Size")
  
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_minNodeCount_pTrash10_", dgpVec[iDGP], ".png"),
                  plot = pMinNodeBar_pTrash10,
                  width = 17.78,
                  height = 9.5,
                  units = "in")
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_minNodeCount_pTrash50_", dgpVec[iDGP], ".png"),
                  plot = pMinNodeBar_pTrash50,
                  width = 17.78,
                  height = 9.5,
                  units = "in")
  
  splitRuleCount <- hyperRF %>% 
    group_by(N, pTrash, rel, R2, lin_inter) %>% 
    count(splitRule) %>% 
    ungroup()
  splitRuleCount$N <- factor(splitRuleCount$N, levels = c(100, 300, 1000),
                             labels = c("N = 100", "N = 300", "N = 1000"))
  splitRuleCount$splitRule <- factor(splitRuleCount$splitRule, levels = c("variance", "extratrees"))
  splitRuleCount$facetY <- interaction(splitRuleCount$lin_inter, splitRuleCount$R2, sep = "x")
  splitRuleCount$facetY <- factor(splitRuleCount$facetY, labels = 
                                  c("R² = 0.2\nlinear 20%", "R² = 0.2\nlinear 80%", 
                                    "R² = 0.5\nlinear 20%", "R² = 0.5\nlinear 80%",
                                    "R² = 0.8\nlinear 20%", "R² = 0.8\nlinear 80%"))
  
  splitRuleCount_pTrash10 <- splitRuleCount[splitRuleCount$pTrash == 10,]
  pSplitruleBar_pTrash10 <- plotHyperBars(splitRuleCount_pTrash10, splitRuleCount_pTrash10$splitRule, "Splitrule")
  splitRuleCount_pTrash50 <- splitRuleCount[splitRuleCount$pTrash == 50,]
  pSplitruleBar_pTrash50 <- plotHyperBars(splitRuleCount_pTrash50, splitRuleCount_pTrash50$splitRule, "Splitrule")
  
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_splitRuleCount_pTrash10_", dgpVec[iDGP], ".png"),
                  plot = pSplitruleBar_pTrash10,
                  width = 17.78,
                  height = 9.5,
                  units = "in")
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperRF_splitRuleCount_pTrash50_", dgpVec[iDGP], ".png"),
                  plot = pSplitruleBar_pTrash50,
                  width = 17.78,
                  height = 9.5,
                  units = "in")
}
