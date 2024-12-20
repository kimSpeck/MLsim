# result plots for H-statistic of the GBMs
# load parameters and helper functions 
source("utils/setParameters.R")
source("utils/analysisTools.R")

# plot H-statistic
library(ggplot2)
library(ggh4x)
library(patchwork)

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
load(paste0(resFolder, "/hStatsPlottingData.rda"))

################################################################################
# plot h-Statistic effects
################################################################################
str(hStats)
col2fac <- c("N", "pTrash", "rel", "simEffect", "R2", "lin_inter")
hStats[col2fac] <- lapply(hStats[col2fac], factor) 

hStats$N <- factor(hStats$N, levels = c(100, 300, 1000))

######
hStats$lin_inter <- factor(hStats$lin_inter, 
                           levels = c("0.8_0.2", "0.5_0.5", "0.2_0.8"),
                           labels = c("80:20", "50:50", "20:80"))

(hStats_pTrash50rel1 <- ggplot(hStats[hStats$pTrash == 50 &
                                        hStats$rel == 1, ], 
                               aes(x = N, y = M, 
                                   group = interaction(R2, simEffect), colour = R2,
                                   linetype = simEffect, shape = simEffect)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_shape_manual(values = c(16, 17)) +
    # plot two standard errors from the mean
    geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), 
                  width = 0.2, alpha = 0.4) +  
    scale_color_manual(name = "R²", values = colValuesR2) +
    ylab("H-statistic") +
    facet_grid2(rel ~ lin_inter,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[3], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[2], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
    xlab("N") +
    theme(panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
          panel.background = element_rect(color = "white", fill = "white"),
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

(hStats_pTrash50rel06 <-  ggplot(hStats[hStats$pTrash == 50 &
                                          hStats$rel == 0.6 & 
                                          hStats$lin_inter == "20:80", ], 
                                 aes(x = N, y = M, 
                                     group = interaction(R2, simEffect), colour = R2,
                                     linetype = simEffect, shape = simEffect)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_shape_manual(values = c(16, 17)) +
    guides(shape = "none", linetype = "none", colour = "none") +
    # # plot two standard errors from the mean
    geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), 
                  width = 0.2, alpha = 0.4) +  
    scale_color_manual(values = colValuesR2) +
    facet_grid2(rel ~ lin_inter,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
    ylim(0, 0.25) +
    ylab("H-statistic") +
    xlab("N") +
    theme(panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
          panel.background = element_rect(color = "white", fill = "white"),
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

# patch subplots 
(phStatsSub_pTrash50 <- hStats_pTrash50rel06 + 
    hStats_pTrash50rel1 + 
    plot_layout(ncol = 2, widths = c(1, 3)))

ggplot2::ggsave(filename = paste0(plotFolder, "/detectInteractions/hStatsInterSubset_pTrash50_2SE.png"),
                plot = phStatsSub_pTrash50,
                width = 17.89,
                height = 5.88,
                units = "in")

################################################################################
# difference plots for H-statistic
################################################################################
hStatsDiff <- tidyr::pivot_wider(hStats[c("N", "pTrash", "rel", "R2", "lin_inter", "simEffect", "M")], 
                                 names_from = simEffect, 
                                 values_from =  M)

hStatsDiff$Hdiff <- hStatsDiff$'1' - hStatsDiff$'0'

(hStats_pTrash50 <- ggplot(hStatsDiff[hStatsDiff$pTrash == 50, ], 
                           aes(x = N, y = Hdiff, 
                               group = interaction(R2, rel), colour = R2,
                               linetype = rel, shape = rel)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(name = "Reliability", values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(name = "Reliability", values = c(18, 16, 17)) +
    scale_color_manual(name = "R²", values = colValuesR2) +
    ylab("H-statistic difference") +
    facet_grid2(~ lin_inter,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[3], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[2], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
    xlab("N") +
    theme(panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
          panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
          panel.background = element_rect(color = "white", fill = "white"),
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

ggplot2::ggsave(filename = paste0(plotFolder, "/detectInteractions/hStatsInterSubset_pTrash50_diff.png"),
                plot = hStats_pTrash50,
                width = 13.08,
                height = 5.88,
                units = "in")
