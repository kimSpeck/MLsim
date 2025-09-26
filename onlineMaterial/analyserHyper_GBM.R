# load parameters and helper functions 
source("utils/setParameters.R")
source("utils/analysisTools.R")

# ANOVA
library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²

# plot results
library(tidyverse)
library(ggplot2)
library(ggh4x)
library(patchwork)
library(see)

# colValues <- c("green3", "darkblue", "darkmagenta")
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
    geom_line(aes(color = R2), linewidth = 8, alpha = 0.5) +
    scale_color_manual(name = "R²", values = c("0.2" = colValuesR2[1], 
                                              "0.5" = colValuesR2[2], 
                                              "0.8" = colValuesR2[3])) +
    xlab(measureLabel) +
    ylab("Count per level")
}

themeFunction <- function(plotObject){
  plotObject + theme(
    panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
    panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.background = element_rect(color = "white", fill = "white"),
    plot.title = element_text(size = 30, face = "bold"),
    axis.text.y = element_text(size = 15),
    # axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
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
listDir <- stringr::str_subset(listDir, "hyperParametersSample_[[:alnum:]]+_GBM.rda")

################################################################################
# anovas and plot hyper parameters
################################################################################
for (iDGP in seq_along(listDir)){
  
  hyperGBM <- loadRData(paste0("results/", dgpVec[iDGP], "/dependentMeasures/", listDir[iDGP]))
  
  # pull data from nested list of all results (fullData)
  hyperGBM <- rbindSingleResults(hyperGBM)
  
  colnames(hyperGBM)
  
  hyperGBM <- idx2infoNew(hyperGBM)
  
  # change variables to factors
  col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model", "sample")
  hyperGBM[col2fac] <- lapply(hyperGBM[col2fac], factor)
  
  col2num <- c("shrinkage", "max_depth", "min_child_weight", "Nrounds")
  hyperGBM[col2num] <- lapply(hyperGBM[col2num], as.numeric)
  # str(hyperGBM)
  ################################################################################
  hyperGBM$ID <- seq_len(dim(hyperGBM)[1])
  
  dvVec <- c("shrinkage", "max_depth", "min_child_weight", "Nrounds")
  
  for (iDV in dvVec){
    # only between-sample ANOVA
    anovaObjectName <- paste0("anova", "GBM", "_", iDV)
    
    tmp <- aov_ez(id = "ID",
                  dv = iDV,
                  data = hyperGBM,
                  between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"))
    assign(anovaObjectName, tmp)
    
    # calculate generalized eta squared as effect size measure
    tmpEta2 <- eta_squared(
      tmp, # fitted model
      partial = FALSE, # not partial!
      generalized = TRUE, # generalized eta squared
      ci = 0.95,
      verbose = TRUE)
    
    # sort gen. eta squared in decreasing order
    tmpEta2 <- tmpEta2[order(tmpEta2$Eta2_generalized, decreasing = T),]
    tmpEta2 <- cbind(tmpEta2, dv = iDV)
    
    # which higher order interactions do we need to illustrate to report simulation results?
    etaObjectName <- paste0("eta2", "GBM", "_", iDV)
    assign(etaObjectName, tmpEta2) 
  }
  
  eta2Thresh <- 0.01
  
  
  
  (pMaxDepth_GBM <- plotEta2_4eachModel(eta2GBM_max_depth, eta2Thresh, "#990000") + 
      ggtitle("max tree depth") + # for the main title
      theme(plot.title = element_text(hjust = 0.5, size = 25)))
  
  (pMinChildWeight_GBM <- plotEta2_4eachModel(eta2GBM_min_child_weight, eta2Thresh, "#990000") + 
      ggtitle("min child weight") + # for the main title
      theme(plot.title = element_text(hjust = 0.5, size = 25)))
  
  
  (pNrounds_GBM <- plotEta2_4eachModel(eta2GBM_Nrounds, eta2Thresh, "#990000") + 
      ggtitle("N rounds") + # for the main title
      theme(plot.title = element_text(hjust = 0.5, size = 25)))
  
  (pShrinkage_GBM <- plotEta2_4eachModel(eta2GBM_shrinkage, eta2Thresh, "#990000") + 
      ggtitle("shrinkage") + # for the main title
      theme(plot.title = element_text(hjust = 0.5, size = 25)))
  
  (pHyper <- pMaxDepth_GBM + pShrinkage_GBM + 
      pNrounds_GBM + pMinChildWeight_GBM +
      plot_layout(ncol = 2))
  
  
  # save plot as files
  ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_hyperGBM_", dgpVec[iDGP], ".png"),
                  plot = pHyper,
                  width = 13.68,
                  height = 15,
                  units = "in")
  
  rm(list=ls(pattern="eta2GBM_"))
  rm(list=ls(pattern="anovaGBM_"))
  
  
  names(hyperGBM)[names(hyperGBM) == 'max_depth'] <- 'maxDepth'
  names(hyperGBM)[names(hyperGBM) == 'min_child_weight'] <- 'minChildWeight'
  
  ################################################################################
  # plot hyper parameter (barplots)
  ################################################################################
  hyperGBM <- hyperGBM[which(hyperGBM$lin_inter != "0.5_0.5"),]
  
  ##### shrinkage #####
  shrinkageCount <- hyperGBM %>% 
    group_by(N, pTrash, rel, R2, lin_inter) %>% 
    count(shrinkage) %>% 
    ungroup()
  shrinkageCount$N <- factor(shrinkageCount$N, levels = c(100, 300, 1000),
                             labels = c("N = 100", "N = 300", "N = 1000"))
  shrinkageCount$shrinkage <- factor(shrinkageCount$shrinkage, levels = unique(shrinkageCount$shrinkage), 
                                     labels = sub("^0\\.", ".", as.character(unique(shrinkageCount$shrinkage))))
  shrinkageCount$facetY <- interaction(shrinkageCount$lin_inter, shrinkageCount$R2, sep = "x")
  shrinkageCount$facetY <- factor(shrinkageCount$facetY, labels = 
                                    c("R² = 0.2\nlinear 20%", "R² = 0.2\nlinear 80%", 
                                      "R² = 0.5\nlinear 20%", "R² = 0.5\nlinear 80%",
                                      "R² = 0.8\nlinear 20%", "R² = 0.8\nlinear 80%"))
  
  shrinkageCount_pTrash10 <- shrinkageCount[shrinkageCount$pTrash == 10,]
  pSkrinkBar_pTrash10 <- plotHyperBars(shrinkageCount_pTrash10, shrinkageCount_pTrash10$shrinkage, "Shrinkage")
  (pSkrinkBar_pTrash10 <- themeFunction(pSkrinkBar_pTrash10))
  shrinkageCount_pTrash50 <- shrinkageCount[shrinkageCount$pTrash == 50,]
  (pSkrinkBar_pTrash50 <- plotHyperBars(shrinkageCount_pTrash50, shrinkageCount_pTrash50$shrinkage, "Shrinkage"))
  pSkrinkBar_pTrash50 <- themeFunction(pSkrinkBar_pTrash50)
  
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperGBM_shrinkCount_pTrash10_", dgpVec[iDGP], ".png"),
                  plot = pSkrinkBar_pTrash10,
                  width = 17.78,
                  height = 9.5,
                  units = "in")
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperGBM_shrinkCount_pTrash50_ ", dgpVec[iDGP], ".png"),
                  plot = pSkrinkBar_pTrash50,
                  width = 17.78,
                  height = 9.5,
                  units = "in")
  
  ##### maxDepth #####
  treeDepthCount <- hyperGBM %>% 
    group_by(N, pTrash, rel, R2, lin_inter) %>% 
    count(maxDepth) %>% 
    ungroup()
  treeDepthCount$N <- factor(treeDepthCount$N, levels = c(100, 300, 1000),
                             labels = c("N = 100", "N = 300", "N = 1000"))
  treeDepthCount$maxDepth <- factor(treeDepthCount$maxDepth, levels = unique(treeDepthCount$maxDepth))
  treeDepthCount$facetY <- interaction(treeDepthCount$lin_inter, treeDepthCount$R2, sep = "x")
  treeDepthCount$facetY <- factor(treeDepthCount$facetY, labels = 
                                    c("R² = 0.2\nlinear 20%", "R² = 0.2\nlinear 80%", 
                                      "R² = 0.5\nlinear 20%", "R² = 0.5\nlinear 80%",
                                      "R² = 0.8\nlinear 20%", "R² = 0.8\nlinear 80%"))
  
  
  treeDepthCount_pTrash10 <- treeDepthCount[treeDepthCount$pTrash == 10,]
  (pDepthBar_pTrash10 <- plotHyperBars(treeDepthCount_pTrash10, treeDepthCount_pTrash10$maxDepth, "max tree depth"))
  pDepthBar_pTrash10 <- themeFunction(pDepthBar_pTrash10)
  treeDepthCount_pTrash50 <- treeDepthCount[treeDepthCount$pTrash == 50,]
  (pDepthBar_pTrash50 <- plotHyperBars(treeDepthCount_pTrash50, treeDepthCount_pTrash50$maxDepth, "max tree depth"))
  pDepthBar_pTrash50 <- themeFunction(pDepthBar_pTrash50)
  
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperGBM_treeDepthCount_pTrash10_", dgpVec[iDGP], ".png"),
                  plot = pDepthBar_pTrash10,
                  width = 17.78,
                  height = 9.5,
                  units = "in")
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperGBM_treeDepthCount_pTrash50_ ", dgpVec[iDGP], ".png"),
                  plot = pDepthBar_pTrash50,
                  width = 17.78,
                  height = 9.5,
                  units = "in")
  
  ##### Ntrees #####
  nTreesCount <- hyperGBM %>% 
    group_by(N, pTrash, rel, R2, lin_inter) %>% 
    count(Nrounds) %>% 
    ungroup()
  nTreesCount$N <- factor(nTreesCount$N, levels = c(100, 300, 1000),
                          labels = c("N = 100", "N = 300", "N = 1000"))
  nTreesCount$Nrounds <- factor(nTreesCount$Nrounds, levels = sort(unique(nTreesCount$Nrounds)))
  nTreesCount$facetY <- interaction(nTreesCount$lin_inter, nTreesCount$R2, sep = "x")
  nTreesCount$facetY <- factor(nTreesCount$facetY, labels = 
                                    c("R² = 0.2\nlinear 20%", "R² = 0.2\nlinear 80%", 
                                      "R² = 0.5\nlinear 20%", "R² = 0.5\nlinear 80%",
                                      "R² = 0.8\nlinear 20%", "R² = 0.8\nlinear 80%"))
  
  
  nTreesCount_pTrash10 <- nTreesCount[nTreesCount$pTrash == 10,]
  pNtreesBar_pTrash10 <- plotHyperBars(nTreesCount_pTrash10, nTreesCount_pTrash10$Nrounds, "N trees")
  (pNtreesBar_pTrash10 <- themeFunction(pNtreesBar_pTrash10) +
      theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1)))
  
  nTreesCount_pTrash50 <- nTreesCount[nTreesCount$pTrash == 50,]
  pNtreesBar_pTrash50 <- plotHyperBars(nTreesCount_pTrash50, nTreesCount_pTrash50$Nrounds, "N trees")
  (pNtreesBar_pTrash50 <- themeFunction(pNtreesBar_pTrash50) + 
      theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1)))
  
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperGBM_nTreesCount_pTrash10_", dgpVec[iDGP], ".png"),
                  plot = pNtreesBar_pTrash10,
                  width = 17.78,
                  height = 9.5,
                  units = "in")
  ggplot2::ggsave(filename = paste0(plotFolder, "/hyperGBM_nTreesCount_pTrash50_ ", dgpVec[iDGP], ".png"),
                  plot = pNtreesBar_pTrash50,
                  width = 17.78,
                  height = 9.5,
                  units = "in")

}




