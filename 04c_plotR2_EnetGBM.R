library(ggplot2)
library(ggh4x)

source("utils/setParameters.R")
source("utils/analysisTools.R")

colValuesR2 <- c('#db4a07', '#850c0c', '#3c1518')
colValuesInter <- c('#050440', '#181ff2', '#0eb2e8')
colValuesLin <- c('#0eb2e8', '#181ff2', '#050440')

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

resFolder <- "results/finalResults/dependentMeasures"

################################################################################
# hier einsteigen?
################################################################################
load("results/finalResults/dependentMeasures/rSquaredData_stats.rda")
load("results/finalResults/dependentMeasures/rSquaredData_eachSample.rda")


##### basic plot from exploratory tests #####
# plot performance measures for train and test data

(pPerformTrainVStest <- ggplot(performanceData,
                               aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                                   group = interaction(R2, model, fit), colour = R2,
                                   linetype = model, shape = fit)) +
   geom_point() +
   geom_line() +
   scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
   scale_shape_manual(values = c(16, 8)) +
   # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
   scale_color_manual(values = colValuesR2) +
   geom_hline(aes(yintercept = 0)) +
   facet_grid(measures + rel ~ lin_inter, labeller = label_both) +
   geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "green3",
              alpha = 0.4) +
   geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "darkblue",
              alpha = 0.4) +
   geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "darkmagenta",
              alpha = 0.4) +
   ylab("") +
   xlab("pTrash (decreasing) x N (increasing)") +
   ggtitle("R^2: Training vs. Test performance") +
   theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "lightgrey"), 
         panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "lightgrey"),
         panel.background = element_rect(color = "white", fill = "white"),
         axis.text.y = element_text(size = 20),
         # axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
         axis.text.x = element_text(size = 15),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         strip.text.x = element_text(size = 15),
         strip.text.y = element_text(size = 15)))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTest_ENTvsGBM.png"),
#                 plot = pPerformTrainVStest,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

##### basic plot from exploratory tests #####
# pTrash raus lassen (nicht mehr auf y-Achse; aber beide Varianten vergleichen) 

# x Achse: rel
# Farbe: R^2
# 3 Facetten für lin_inter

# Facetten für N zusätzlich ()
(pR2_overview <- ggplot(performanceData[performanceData$pTrash == 50 &
                                          performanceData$measures == "M_test", ],
                        aes(x = rel, y = values, 
                            group = interaction(R2, model, fit), colour = R2,
                            linetype = model, shape = fit)) +
   geom_point() +
   geom_line() +
   scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
   scale_shape_manual(values = c(16, 8)) +
   # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
   scale_color_manual(values = colValuesR2) +
   geom_hline(aes(yintercept = 0)) +
   facet_grid(N ~ lin_inter, labeller = label_both) +
   geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "green3",
              alpha = 0.4) +
   geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "darkblue",
              alpha = 0.4) +
   geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "darkmagenta",
              alpha = 0.4) +
   ylab("") +
   xlab("reliability of predictors") +
   ggtitle("R^2: Training vs. Test performance") +
   theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "lightgrey"), 
         panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "lightgrey"),
         panel.background = element_rect(color = "white", fill = "white"),
         axis.text.y = element_text(size = 20),
         # axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
         axis.text.x = element_text(size = 15),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         strip.text.x = element_text(size = 15),
         strip.text.y = element_text(size = 15)))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/R2_newOverview.png"),
#                 plot = pR2_overview,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

##### add Monte Carlo Error
# as standard deviation of the Monte Carlo estimate
# plot 2 Monte Carlo Errors as Error Bars  
str(rSquaredTest)
# rSquaredTest$testR2_relational <- rSquaredTest$Rsq_test / as.numeric(as.character(rSquaredTest$R2))


getSD <- aggregate(Rsq_test ~ model + N + pTrash + rel + R2 + lin_inter, 
                   data = rSquaredTest, sd)
getN <- aggregate(Rsq_test ~ model + N + pTrash + rel + R2 + lin_inter, 
                  data = rSquaredTest, length)
# unique(getN$Rsq_test) # always based on 1000 sample
aggregate(Rsq_test ~ model + N + pTrash + rel + R2 + lin_inter, 
          data = rSquaredTest, mean)

pR2sub <- performanceData[performanceData$measure == "Rsquared" &
                            performanceData$measures == "M_test", ]

pR2sub <- merge(pR2sub, getSD, by = c("model", "N", "pTrash", "rel", "R2", "lin_inter"))
colnames(pR2sub)[colnames(pR2sub) == 'Rsq_test'] <- 'SD'

pR2sub$SE <- pR2sub$SD / sqrt(setParam$dgp$nTrain)
pR2sub$lin_inter <- factor(pR2sub$lin_inter,
                           levels = c("0.8_0.2", "0.5_0.5", "0.2_0.8"),
                           labels = c("80:20", "50:50", "20:80"))

pR2sub$model <- factor(pR2sub$model,
                       levels = c("ENETw", "ENETwo", "GBM"),
                       labels = c("ENETinter", "ENETlin", "GBM"))

(pR2_N100 <- ggplot(pR2sub[pR2sub$pTrash == 50 &
                             pR2sub$N == 100, ],
                    aes(x = rel, y = values, 
                        group = interaction(R2, model), colour = R2,
                        linetype = model, shape = model)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    # geom_errorbar(aes(ymin = values - 2*SD, ymax = values + 2*SD),
    #                width = 0.2, alpha = 0.4) +
    geom_errorbar(aes(ymin = values - 2*SE, ymax = values + 2*SE),
                  width = 0.2, alpha = 0.4) +
    scale_linetype_manual(name = "ML model", values = c("solid", "dashed", "dotted")) +
    scale_shape_manual(name = "ML model", values = c(16, 16, 17)) +
    scale_color_manual(name = "R²", values = colValuesR2) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid2(~ lin_inter,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[3], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[2], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = colValuesR2[1],
               alpha = 0.6, linewidth = 1) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = colValuesR2[2],
               alpha = 0.6, linewidth = 1) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = colValuesR2[3],
               alpha = 0.6, linewidth = 1) +
    ylab("Test R²") +
    xlab("Reliability of predictors") +
    ggtitle("B") +
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
          legend.position = c(.85, .85), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal"))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/R2_N100_2SE.png"),
#                 plot = pR2_N100,
#                 width = 13.63,
#                 height = 8.59,
#                 units = "in")

# (pR2_N100_0MCE <- ggplot(pR2sub[pR2sub$pTrash == 50 &
#                              pR2sub$N == 100, ],
#                     aes(x = rel, y = values, 
#                         group = interaction(R2, model, fit), colour = R2,
#                         linetype = model, shape = fit)) +
#     geom_point() +
#     geom_line() +
#     #geom_errorbar(aes(ymin = values - 2* MCE, ymax = values + 2*MCE), width = 0.2) +  
#     scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
#     scale_shape_manual(values = c(16, 8)) +
#     scale_color_manual(values = colValues) +
#     geom_hline(aes(yintercept = 0)) +
#     facet_wrap(~ lin_inter) +
#     geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "green3",
#                alpha = 0.4) +
#     geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "darkblue",
#                alpha = 0.4) +
#     geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "darkmagenta",
#                alpha = 0.4) +
#     ylab("test R²") +
#     xlab("reliability of predictors") +
#     ggtitle("R²: Training vs. Test performance {N = 100, pTrash = 50}") +
#     theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "lightgrey"), 
#           panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "lightgrey"),
#           panel.background = element_rect(color = "white", fill = "white"),
#           axis.text.y = element_text(size = 20),
#           # axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
#           axis.text.x = element_text(size = 20),
#           axis.title.x = element_text(size = 20),
#           axis.title.y = element_text(size = 20),
#           strip.text.x = element_text(size = 15),
#           strip.text.y = element_text(size = 15)))

(pR2_N1000 <- ggplot(pR2sub[pR2sub$pTrash == 50 &
                              pR2sub$N == 1000, ],
                     aes(x = rel, y = values, 
                         group = interaction(R2, model, model), colour = R2,
                         linetype = model, shape = model)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    # geom_errorbar(aes(ymin = values - 2*SD, ymax = values + 2*SD),
    #               width = 0.2, alpha = 0.4) +
    geom_errorbar(aes(ymin = values - 2*SE, ymax = values + 2*SE),
                  width = 0.2, alpha = 0.4) +
    scale_linetype_manual(name = "ML model", values = c("solid", "dashed", "dotted")) +
    scale_shape_manual(name = "ML model", values = c(16, 16, 17)) +
    scale_color_manual(name = "R²", values = colValuesR2) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid2(~ lin_inter,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[3], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[2], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = colValuesR2[1],
               alpha = 0.6, linewidth = 1) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = colValuesR2[2],
               alpha = 0.6, linewidth = 1) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = colValuesR2[3],
               alpha = 0.6, linewidth = 1) +
    ylab("Test R²") +
    xlab("Reliability of predictors") +
    ggtitle("A") +
    guides(color = "none", 
           shape = "none",
           linetype = "none") +
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
          legend.position = c(.85, .85), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal"))

# ggplot2::ggsave(filename = paste0(plotFolder, "/R2_N1000_2SE.png"),
#                 plot = pR2_N1000,
#                 width = 13.63,
#                 height = 8.59,
#                 units = "in")

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/R2_N100_2SE.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/R2_N100_2SD.png"),
#                 plot = pR2_N100,
#                 width = 13.63,
#                 height = 8.59,
#                 units = "in")
# ggplot2::ggsave(filename = paste0(plotFolder, "/R2_N1000_2SE.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/R2_N1000_2SD.png"),
#                 plot = pR2_N1000,
#                 width = 13.63,
#                 height = 8.59,
#                 units = "in")

# ggplot2::ggsave(filename = paste0(plotFolder, "/ENETvsGBM/R2_N100_0MCE.png"),
#                 plot = pR2_N100_0MCE,
#                 width = 13.63,
#                 height = 8.59,
#                 units = "in")