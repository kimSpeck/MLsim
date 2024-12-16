# load parameters and helper functions 
source("setParameters.R")
source("analysisTools.R")

# ANOVA
library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²

# plot results
library(ggplot2)
library(ggh4x)
library(patchwork)
library(see)

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
hyperGBM <- loadRData(paste0(resFolder, "/hyperParametersSample_GBM.rda"))

# pull data from nested list of all results (fullData)
hyperGBM <- rbindSingleResults(hyperGBM)

colnames(hyperGBM)

hyperGBM <- idx2infoNew(hyperGBM)

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model", "sample")
hyperGBM[col2fac] <- lapply(hyperGBM[col2fac], factor)

col2num <- c("shrinkage", "max_depth", "min_child_weight", "Nrounds")
hyperGBM[col2num] <- lapply(hyperGBM[col2num], as.numeric)
str(hyperGBM)
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

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_hyperGBM.png"),
#                 plot = pHyper,
#                 width = 13.68,
#                 height = 15,
#                 units = "in")

################################################################################
# plot hyper parameter values (alpha and lambda for the ENET)
################################################################################
names(hyperGBM)[names(hyperGBM) == 'max_depth'] <- 'maxDepth'
names(hyperGBM)[names(hyperGBM) == 'min_child_weight'] <- 'minChildWeight'

plotHyper <- aggregate(cbind(shrinkage, maxDepth, minChildWeight, Nrounds) ~ 
                         N + pTrash + rel + R2 + lin_inter, 
                       data = hyperGBM,  
                       function(x) {cbind(mean(x), 
                                          sd(x),
                                          quantile(x, 0.025),
                                          quantile(x, 0.975))})

# turn matrices inside data frame into columns
plotHyper <- do.call(data.frame, plotHyper)
colnames(plotHyper) <- stringr::str_replace_all(colnames(plotHyper), "\\.1", "_M")
colnames(plotHyper) <- stringr::str_replace_all(colnames(plotHyper), "\\.2", "_SE")
colnames(plotHyper) <- stringr::str_replace_all(colnames(plotHyper), "\\.3", "_q025")
colnames(plotHyper) <- stringr::str_replace_all(colnames(plotHyper), "\\.4", "_q975")

plotHyper$N <- factor(plotHyper$N, levels = c(100, 300, 1000))

plotHyper <- tidyr::pivot_longer(plotHyper, 
                                 cols = !c(N, pTrash, rel, R2, lin_inter),
                                 names_to = c("DV", "measure"), 
                                 names_sep = "_",
                                 values_to = "values")

plotHyper <- tidyr::pivot_wider(plotHyper,
                                names_from = measure, 
                                values_from = values)

plotHyper <- tidyr::separate(plotHyper, "lin_inter", 
                             into = c("lin", "inter"), sep = "_", remove = F) 

# N and R2 (+ their interaction) influence the max tree depth and the shrinkage
(pTreeDepth_rel08 <- ggplot(plotHyper[plotHyper$DV == "maxDepth" &
                                                       plotHyper$rel == 0.8,],
                                           aes(x = N, y = M, 
                                               group = interaction(lin_inter, pTrash), colour = lin_inter,
                                               linetype = lin_inter, shape = pTrash)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(name = "Lin:Inter", values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(name = "Noise", values = c(16, 17)) +
    scale_color_manual(name = "Lin:Inter", values = colValuesLin) +
    facet_grid2(~ R2,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesR2[1], 0.4)),
                                      element_rect(fill = alpha(colValuesR2[2], 0.4)),
                                      element_rect(fill = alpha(colValuesR2[3], 0.4))))) + 
    ylab("max tree depth") +
    xlab("N") +
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
          legend.box = "horizontal"))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/hyperGBM_treeDepth_R2facette_rel08.png"),
#                 plot = pTreeDepth_rel08,
#                 width = 13.08,
#                 height = 6.68,
#                 units = "in")

(pShrinkage_rel08 <- ggplot(plotHyper[plotHyper$DV == "shrinkage" &
                                        plotHyper$rel == 0.8,],
                            aes(x = N, y = M, 
                                group = interaction(lin_inter, pTrash), colour = lin_inter,
                                linetype = lin_inter, shape = pTrash)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(name = "Lin:Inter", values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(name = "Noise", values = c(16, 17)) +
    scale_color_manual(name = "Lin:Inter", values = colValuesLin) +
    facet_grid2(~ R2,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesR2[1], 0.4)),
                                      element_rect(fill = alpha(colValuesR2[2], 0.4)),
                                      element_rect(fill = alpha(colValuesR2[3], 0.4))))) + 
    ylab("shrinkage") +
    xlab("N") +
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
          legend.box = "horizontal"))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/hyperGBM_shrinkage_R2facette_rel08.png"),
#                 plot = pShrinkage_rel08,
#                 width = 13.08,
#                 height = 6.68,
#                 units = "in")