# load packages for plotting
library(ggplot2)
library(ggh4x)

# get parameter values and utility functions 
source("utils/setParameters.R")
source("utils/analysisTools.R")

# plotting colors across final paper plots
colValuesR2 <- c('#db4a07', '#850c0c', '#3c1518')
colValuesInter <- c('#050440', '#181ff2', '#0eb2e8')
colValuesLin <- c('#0eb2e8', '#181ff2', '#050440')

# path to plot folder
plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

################################################################################
# load and merge new and old GBM results for dgp == "inter"
################################################################################

# load old data
load("results/finalResults/dependentMeasures/rSquaredData_stats.rda")
load("results/finalResults/dependentMeasures/rSquaredData_eachSample.rda")

performanceGBM <- performanceData[performanceData$model == "GBM",]
rSquaredTestGBM <- rSquaredTest[rSquaredTest$model == "GBM",]

# load new GBM data
load("results/inter/dependentMeasures/rSquaredData_stats.rda")
load("results/dependentMeasures/rSquaredData_eachSample.rda")

newPerformanceGBM <- performanceData[performanceData$model == "GBM" &
                                       performanceData$dgp == "inter",]

newPerformanceGBM$dgp <- NULL
colnames(newPerformanceGBM) == colnames(performanceGBM)
performanceGBM$grid <- "oldGrid"
newPerformanceGBM$grid <- "newGrid"

performanceData <- rbind(performanceGBM, newPerformanceGBM)

newRSquaredTest <- rSquaredTest[rSquaredTest$model == "GBM" &
                                  rSquaredTest$dgp == "inter", ]
newRSquaredTest$dgp <- NULL
rSquaredTestGBM$ID <- NULL

colnames(newRSquaredTest) == colnames(rSquaredTestGBM)
rSquaredTestGBM$grid <- "oldGrid"
newRSquaredTest$grid <- "newGrid"

rSquaredTest <- rbind(rSquaredTestGBM, newRSquaredTest)

rm(newPerformanceGBM, newRSquaredTest, performanceGBM, rSquaredTestGBM)
################################################################################
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
pR2sub$grid <- factor(pR2sub$grid)

plotR2sub <- function(data, pTrash, N, title = "", guides = T) {
  pTmp <- ggplot(data[data$pTrash == pTrash &
                        data$N == N, ],
                 aes(x = rel, y = values, 
                     group = interaction(R2, grid), colour = R2,
                     linetype = grid, shape = grid)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    geom_errorbar(aes(ymin = values - 2*SE, ymax = values + 2*SE),
                  width = 0.2, alpha = 0.4) +
    scale_linetype_manual(name = "GBM grid", values = c("solid", "dashed", "dotted", "dotdash")) +
    scale_shape_manual(name = "GBM grid", values = c(16, 17)) +
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
          legend.position = "bottom", 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal")
  if (guides == F) {
    pTmp <- pTmp + guides(color = "none", shape = "none", linetype = "none")
  }
  return(pTmp)
}

pR2_N100_pTrash50 <- plotR2sub(pR2sub, pTrash = 50, N = 100, title = "B", guides = T)
pR2_N1000_pTrash50 <-plotR2sub(pR2sub, pTrash = 50, N = 1000, title = "A", guides = T)
pR2_N300_pTrash50 <-plotR2sub(pR2sub, pTrash = 50, N = 300, title = "", guides = F)

pR2_N100_pTrash10 <- plotR2sub(pR2sub, pTrash = 10, N = 100, title = "B", guides = T)
pR2_N1000_pTrash10 <-plotR2sub(pR2sub, pTrash = 10, N = 1000, title = "A", guides = T)
pR2_N300_pTrash10 <-plotR2sub(pR2sub, pTrash = 10, N = 300, title = "", guides = F)

# save plot as files
ggplot2::ggsave(filename = paste0(plotFolder, "/gridGBM_R2_N100_pTrash10.png"),
                plot = pR2_N100_pTrash10,
                width = 13.63,
                height = 8.59,
                units = "in")
ggplot2::ggsave(filename = paste0(plotFolder, "/gridGBM_R2_N300_pTrash10.png"),
                plot = pR2_N300_pTrash10,
                width = 13.63,
                height = 8.59,
                units = "in")
ggplot2::ggsave(filename = paste0(plotFolder, "/gridGBM_R2_N1000_pTrash10.png"),
                plot = pR2_N1000_pTrash10,
                width = 13.63,
                height = 8.59,
                units = "in")

# save plot as files
ggplot2::ggsave(filename = paste0(plotFolder, "/gridGBM_R2_N100_pTrash50.png"),
                plot = pR2_N100_pTrash50,
                width = 13.63,
                height = 8.59,
                units = "in")
ggplot2::ggsave(filename = paste0(plotFolder, "/gridGBM_R2_N300_pTrash50.png"),
                plot = pR2_N300_pTrash50,
                width = 13.63,
                height = 8.59,
                units = "in")
ggplot2::ggsave(filename = paste0(plotFolder, "/gridGBM_R2_N1000_pTrash50.png"),
                plot = pR2_N1000_pTrash50,
                width = 13.63,
                height = 8.59,
                units = "in")
