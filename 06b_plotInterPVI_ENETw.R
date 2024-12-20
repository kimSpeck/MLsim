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

load(paste0(resFolder, "/interENETinterMeasures.rda"))
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

plotInterPVI <- aggregate(cbind(interTP, interFP, interPPV, interACC,
                                interSpecificity, interSensitivity, interBalACC) ~ 
                            model + N + pTrash + rel + R2 + lin_inter, 
                          data = interENETw, 
                          function(x) {cbind(mean(x), 
                                             sd(x),
                                             quantile(x, 0.025),
                                             quantile(x, 0.975))})

plotInterPVI <- do.call(data.frame, plotInterPVI)
colnames(plotInterPVI) <- stringr::str_replace_all(colnames(plotInterPVI), "\\.1", "_M")
colnames(plotInterPVI) <- stringr::str_replace_all(colnames(plotInterPVI), "\\.2", "_SE")
colnames(plotInterPVI) <- stringr::str_replace_all(colnames(plotInterPVI), "\\.3", "_q025")
colnames(plotInterPVI) <- stringr::str_replace_all(colnames(plotInterPVI), "\\.4", "_q975")

plotInterPVI$N <- factor(plotInterPVI$N, levels = c(100, 300, 1000))

plotInterPVI <- tidyr::pivot_longer(plotInterPVI, 
                                    cols = !c(model, N, pTrash, rel, R2, lin_inter),
                                    names_to = c("DV", "measure"), 
                                    names_sep = "_",
                                    values_to = "values")

plotInterPVI <- tidyr::pivot_wider(plotInterPVI,
                                   names_from = measure, 
                                   values_from = values)

################################################################################
# plots 
################################################################################

plotInterPVI$lin_inter <- factor(plotInterPVI$lin_inter, 
                                 levels = c("0.2_0.8", "0.5_0.5", "0.8_0.2"),
                                 labels = c("20:80", "50:50", "80:20"))

plotInterENETw <- function(data, DV, rel, guides = T, yTitle = "", title = "") {
  pTmp <- ggplot(data[data$DV == DV & data$rel == rel,],
         aes(x = N, y = M, 
             group = interaction(lin_inter, pTrash), colour = lin_inter,
             linetype = lin_inter, shape = pTrash)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(name = "Lin:Inter", values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(name = "Noise", values = c(16, 17)) +
    scale_color_manual(name = "Lin:Inter", values = colValuesInter) +
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
  if (guides == F) {
    pTmp <- pTmp + guides(color = "none", shape = "none", linetype = "none")
  }
  pTmp
}

(pInterSensitivity_pviR2facette_rel08 <- plotInterENETw(plotInterPVI, DV = "interSensitivity", rel = 0.8, 
               yTitle = "Sensitivity", title = "A", guides = T))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/detectInteractions/interSensitivity_pviR2facette_rel08.png"),
#                 plot = pInterSensitivity_pviR2facette_rel08,
#                 width = 13.08,
#                 height = 6.68,
#                 units = "in")

(pInterSpecificity_pviR2facette_rel08 <- plotInterENETw(plotInterPVI, DV = "interSpecificity", rel = 0.8, 
               yTitle = "Specificity", title = "B", guides = F))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/detectInteractions/interSpecificity_pviR2facette_rel08.png"),
#                 plot = pInterSpecificity_pviR2facette_rel08,
#                 width = 13.08,
#                 height = 6.68,
#                 units = "in")