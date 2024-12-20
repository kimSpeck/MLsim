# load parameters and helper functions 
source("utils/setParameters.R")
source("utils/analysisTools.R")

# plot
library(ggplot2)
library(ggh4x)
library(patchwork)
library(see)

# plot utils
colValuesR2 <- c('#db4a07', '#850c0c', '#3c1518')
colValuesInter <- c('#050440', '#181ff2', '#0eb2e8')
colValuesLin <- c('#0eb2e8', '#181ff2', '#050440')

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

# load results files
resFolder <- "results/finalResults/dependentMeasures" 
################################################################################
# plot specificity, sensitivity and (balanced) accuracy 
################################################################################
# calculate M, SE, 2.5% Quantile, 97.5% Quantile for dependent measures:
#   ... linear TP effects
#   ... linear FP effects
#   ... positive predictive value
#   ... accuracy
#   ... specificity
#   ... sensitivity
#   ... balanced accuracy

load(paste0(resFolder, "/relFrequencyMeasures.rda"))

# average sensitivity across all simulated conditions (for linear effects in the ENETinter)
mean(linENETw$sensitivity)

plotDiags <- aggregate(cbind(linTP, linFP, PPV, ACC, specificity, sensitivity, 
                             balACC, linTN, linFN) ~ 
                         model + N + pTrash + rel + R2 + lin_inter, 
                       data = rbind(linENETwo, linENETw, linGBM), 
                       function(x) {cbind(mean(x), 
                                          sd(x),
                                          quantile(x, 0.025),
                                          quantile(x, 0.975))})

# turn matrices inside data frame into columns
plotDiags <- do.call(data.frame, plotDiags)
colnames(plotDiags) <- stringr::str_replace_all(colnames(plotDiags), "\\.1", "_M")
colnames(plotDiags) <- stringr::str_replace_all(colnames(plotDiags), "\\.2", "_SD")
colnames(plotDiags) <- stringr::str_replace_all(colnames(plotDiags), "\\.3", "_q025")
colnames(plotDiags) <- stringr::str_replace_all(colnames(plotDiags), "\\.4", "_q975")

plotDiags$N <- factor(plotDiags$N, levels = c(100, 300, 1000))
plotDiags$model <- factor(plotDiags$model, 
                          levels = c("ENETw", "ENETwo", "GBM"))

plotDiags <- tidyr::pivot_longer(plotDiags, 
                                 cols = !c(model, N, pTrash, rel, R2, lin_inter),
                                 names_to = c("DV", "measure"), 
                                 names_sep = "_",
                                 values_to = "values")

plotDiags <- tidyr::pivot_wider(plotDiags,
                                names_from = measure, 
                                values_from = values)

# generic plot function does not work
# as indicated by ANOVA results, different manipulations affected different kind of measures

################################################################################
# plot sensitivity and specificity
################################################################################
plotDiags$lin_inter <- factor(plotDiags$lin_inter, 
                              levels = c("0.2_0.8", "0.5_0.5", "0.8_0.2"),
                              labels = c("20:80", "50:50", "80:20"))

plotSensSpec <- function(data, DV, model, rel, 
                     quantiles = F, guides = T, 
                     title = "", yTitle = "") {
  pTmp <- ggplot(data[data$DV == DV &
                     data$model == model &
                     data$rel == rel,],
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
    ylim(c(0, 1)) +
    ylab(yTitle) +
    xlab("N") +
    ggtitle(title) +
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
          legend.position = c(.85, .15), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal")
  if (quantiles == T) {
    pTmp <- pTmp + geom_errorbar(aes(ymin = q025, ymax = q975), 
                                 width = 0.2, alpha = 0.4, 
                                 position = position_dodge(width = 0.5)) 
  }
  if (guides == F) {
    pTmp <- pTmp + guides(color = "none", shape = "none", linetype = "none")
  }
  pTmp
}

(pSensitivityENETw_R2facette_rel08 <- plotSensSpec(plotDiags, DV = "sensitivity", 
                                                   model = "ENETw", rel = 0.8,
                                                   quantiles = F, yTitle = "Sensitivity", title = "A"))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/detectMains/mainSensitivity_pviENETwR2facette_rel08.png"),
#                 plot = pSensitivityENETw_R2facette_rel08,
#                 width = 13.08,
#                 height = 6.68,
#                 units = "in")

(pSensitivityGBM_R2facette_rel08 <- plotSensSpec(plotDiags, DV = "sensitivity", 
             model = "GBM", rel = 0.8,
             quantiles = F, yTitle = "Sensitivity", title = "A"))
 
# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/detectMains/mainSensitivity_pviGBMR2facette_rel08.png"),
#                 plot = pSensitivityGBM_R2facette_rel08,
#                 width = 13.08,
#                 height = 6.68,
#                 units = "in")

(pSpecificityENETw_R2facette_rel08 <- plotSensSpec(plotDiags, DV = "specificity", 
                                                 model = "ENETw", rel = 0.8,
                                                 guides = F, quantiles = F, yTitle = "Specificity", title = "B"))
# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/detectMains/mainSpecificity_pviENETwR2facette_rel08.png"),
#                 plot = pSpecificityENETw_R2facette_rel08,
#                 width = 13.08,
#                 height = 6.68,
#                 units = "in")

(pSpecificityGBM_R2facette_rel08 <- plotSensSpec(plotDiags, DV = "specificity", 
                                                   model = "GBM", rel = 0.8,
                                                   guides = F, quantiles = F, yTitle = "Specificity", title = "B"))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/detectMains/mainSpecificity_pviGBMR2facette_rel08.png"),
#                 plot = pSpecificityGBM_R2facette_rel08,
#                 width = 13.08,
#                 height = 6.68,
#                 units = "in")