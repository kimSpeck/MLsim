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

condGrid <- expand.grid(N = setParam$dgp$N,
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability)

condN_pTrash <- paste0("N", condGrid$N, 
                       "_pTrash", condGrid$pTrash,
                       "_rel", condGrid$reliability)

################################################################################
# plot utility functions
################################################################################
themeFunction <- function(plotObject) {
  plotObject + theme(
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "lightgrey"), 
    panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "lightgrey"),
    panel.background = element_rect(color = "white", fill = "white"),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15))
}

themeFunctionPaper <- function(plotObject, guides = T){
  pTmp <- plotObject + theme(
    panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
    panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.background = element_rect(color = "white", fill = "white"),
    plot.title = element_text(size = 30, face = "bold"),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    #legend.position = c(.85, .39), 
    legend.position = "bottom", 
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.key.width = unit(2, "cm"),
    legend.box = "horizontal")
  if (guides == F) {
    pTmp <- pTmp + guides(color = "none", shape = "none", linetype = "none")
  }
  return(pTmp)
}
################################################################################
# load data
################################################################################
dgpVec <- c("inter", "pwlinear", "nonlinear3")
resFolderVec <- paste0("results/", dgpVec, "/dependentMeasures")
resFolder <- "results/"

# read in data
listDir <- dir(resFolderVec)
dataList <- listDir[stringr::str_detect(listDir, "^performT")]

dgps <- stringr::str_extract(dataList, "_[:alpha:]*3{0,1}_")
dgps <- stringr::str_sub(dgps, start = 2L, end = -2)

models <- stringr::str_extract(dataList, "_[:alpha:]*.rda$")
models <- stringr::str_sub(models, start = 2L, end = -5)

tvst <- stringr::str_extract(dataList, "(Test|Train)") 

for (iData in seq_len(length(dataList))) {
  objectName <- paste0("p", tvst[iData], "_", models[iData], "_", dgps[iData])
  assign(objectName, loadRData(paste0(resFolder, "/", dgps[iData], 
                                      "/dependentMeasures/", dataList[iData])))
}

################################################################################
# prepare data
################################################################################
trainList <- sapply(grep("^pTrain", ls(), value = TRUE), get, simplify = FALSE)
testList <- sapply(grep("^pTest", ls(), value = TRUE), get, simplify = FALSE)

trainList <- lapply(trainList, function(iTrain) {
  # pull data from nested list of all results (fullData)
  tmp <- rbindSingleResults(iTrain)
  
  # variables from rownames to own column to work with variable information
  tmp <- rowNames2col(tmp, "measure")
  tmp$measure <- paste0(tmp$measure, "_train")
  
  # get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
  idx2infoNew(tmp)
})

testList <- lapply(testList, function(iTest) {
  # pull data from nested list of all results (fullData)
  tmp <- rbindSingleResults(iTest)
  
  # variables from rownames to own column to work with variable information
  tmp <- rowNames2col(tmp, "measure")
  tmp$measure[which(tmp$measure == "Rsq_test")] <- "Rsquared_test"
  
  # get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
  idx2infoNew(tmp) 
})

# merge performance Train and performance Test
trainList <- do.call(rbind, trainList)
testList <- do.call(rbind, testList)

# # check
# colnames(trainList)
# colnames(testList)

performanceStats <- rbind(trainList, testList)
performanceStats <- tidyr::separate(performanceStats, measure, c("measure", "trainTest"), sep = "_")

# change type of columns or specific entry details to prepare plotting  
performanceStats$N <- factor(performanceStats$N, levels = setParam$dgp$N)
performanceStats$pTrash <- factor(performanceStats$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
str(performanceStats)
chr2fac <- c("rel", "measure", "R2", "lin_inter")
performanceStats[chr2fac] <- lapply(performanceStats[chr2fac], factor)

performanceStats$lin_inter <- factor(performanceStats$lin_inter,
                           levels = c("0.8_0.2", "0.5_0.5", "0.2_0.8"),
                           labels = c("80:20", "50:50", "20:80"))

performanceStats$model <- factor(performanceStats$model,
                                     levels = c("ENETw", "ENETwo", "GBM", "RF"),
                                     labels = c("◆ ENETint", "○ ENETlin", "+ GBM", "x RF"))

performanceStats$dgp <- factor(performanceStats$dgp, 
                               levels = c("inter", "pwlinear", "nonlinear3"), 
                               labels = c("Interaction", "Piecewise", "Stepwise"))

# plot overfit instead of train as train - test
performanceStats <- tidyr::pivot_wider(performanceStats, 
                                       names_from = "trainTest", values_from = c(M, SE, SD))
chr2num <- c("M_train", "M_test", "SE_train", "SE_test", "SD_train", "SD_test")
performanceStats[chr2num] <- lapply(performanceStats[chr2num], as.numeric)

performanceStats$overfit <- performanceStats$M_train - performanceStats$M_test

rm(list = ls(pattern = "^pT"))

################################################################################
# plot results for paper (subset of the data)
################################################################################
# ENETinter, ENETlin, GBM, RF
# modelLineColors <- c("#0a2463", "#2a9d8f", "#8338ec", "#a44200")
modelLineColors <- c("#A92525", "#F57F17", "#2a9d8f", "#1A237E")
modelPointShapes <- c(18, 1, 3, 4)

# interaction, piecewise, stepwise
dgpPanelColors <- c("#F44336", "#9C27B0", "#1565C0")
# relColValues <- c("#7BD580", "#1B5E20") # rel = 0.8, rel = 1; greenish colors
relColValues <- c("#999999", "#444444") # rel = 0.8, rel = 1
linetypeVec <- c("dotted", "dashed", "solid")


plotPaperR2 <- function(data, plotMeasure, title = "", yLabel = "",
                        yMin = NULL, yMax = NULL){
  ggplot(data,
         aes(x = N, y = plotMeasure, 
             group = interaction(R2, model), linetype = R2, color = model)) +
    geom_line(linewidth = 1, alpha = 0.4) +
    geom_point(aes(shape = model), size = 3) +
    scale_linetype_manual(name = expression(R[sim]^2), values = linetypeVec,
                          guide = guide_legend(override.aes = list(size = 2, alpha = 0.8))) +
    scale_color_manual(values = modelLineColors) +
    scale_shape_manual(values = modelPointShapes) +
    scale_y_continuous(limits = c(yMin, yMax), breaks=seq(
      ifelse(yMin %% 0.2 == 0, yMin, yMin + 0.1), yMax, 0.2)) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid2(rel ~ dgp,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(dgpPanelColors[1], 0.4)),
                                      element_rect(fill = alpha(dgpPanelColors[2], 0.4)),
                                      element_rect(fill = alpha(dgpPanelColors[3], 0.4))),
                  background_y = list(element_rect(fill = alpha(relColValues[2], 0.4)),
                                      element_rect(fill = alpha(relColValues[1], 0.4))))) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "black",
               alpha = 0.4, linetype = linetypeVec[1]) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "black",
               alpha = 0.4, linetype = linetypeVec[2]) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "black",
               alpha = 0.4, linetype = linetypeVec[3]) +
    ylab(yLabel) +
    xlab("Sample Size (N)") +
    ggtitle(title) +
    guides(colour = "none", shape = "none")
}

subFig1 <- performanceStats[which(performanceStats$measure == "Rsquared" & 
                                    performanceStats$lin_inter == "20:80" & 
                                    performanceStats$rel != 0.6 &
                                    performanceStats$pTrash == 50),]

subFig2 <- performanceStats[which(performanceStats$measure == "Rsquared" & 
                                   performanceStats$lin_inter == "80:20" & 
                                   performanceStats$rel != 0.6 &
                                   performanceStats$pTrash == 50),]
subFig1$rel <- factor(subFig1$rel, levels = c(1, 0.8))
subFig2$rel <- factor(subFig2$rel, levels = c(1, 0.8))

fig1 <- plotPaperR2(subFig1, plotMeasure = subFig1$M_test, title = "", yLabel = expression(R[test]^2),
                    yMin = 0, yMax = 0.9)
(fig1 <- themeFunctionPaper(fig1))
(fig1 <- fig1 + geom_text(data = subFig1[which(subFig1$R2 == 0.8 &
                                                 subFig1$N == 300),],
                         aes(x = N, y = 0.85, group = interaction(R2, model),
                             color = model, label = model), position = position_dodge(3)))

fig2 <- plotPaperR2(subFig2, plotMeasure = subFig2$M_test, title = "", yLabel = expression(R[test]^2),
                    yMin = 0, yMax = 0.9)
fig2 <- themeFunctionPaper(fig2, guides = T)
(fig2 <- fig2 + geom_text(data = subFig2[which(subFig2$R2 == 0.8 &
                                                 subFig2$N == 300),],
                          aes(x = N, y = 0.85, group = interaction(R2, model),
                              color = model, label = model), position = position_dodge(3)))

ggplot2::ggsave(filename = paste0(plotFolder, "/R2_20:80_pTrash50_paper.png"),
                plot = fig1,
                width = 13.63,
                height = 12.07,
                units = "in")
ggplot2::ggsave(filename = paste0(plotFolder, "/R2_80:20_pTrash50_paper.png"),
                plot = fig2,
                width = 13.63,
                height = 12.07,
                units = "in")

################################################################################
# plot overfitting
################################################################################
subOverfitFig1 <- performanceStats[which(performanceStats$measure == "Rsquared" & 
                                    performanceStats$lin_inter == "20:80" & 
                                    performanceStats$rel != 0.6 &
                                    performanceStats$pTrash == 50),]

subOverfitFig1$rel <- factor(subOverfitFig1$rel, levels = c(1, 0.8))
overfitFig1 <- plotPaperR2(subOverfitFig1, plotMeasure = subOverfitFig1$overfit, title = "", 
                           yLabel = expression(R[train]^2 - R[test]^2),
                           yMin = -0.1, yMax = 1.1)
(overfitFig1 <- themeFunctionPaper(overfitFig1))
(overfitFig1 <- overfitFig1 + geom_text(data = subOverfitFig1[which(subOverfitFig1$R2 == 0.8 &
                                                                   subOverfitFig1$N == 300),],
                          aes(x = N, y = 1.05, group = interaction(R2, model),
                              color = model, label = model), position = position_dodge(3)))

subOverfitFig2 <- performanceStats[which(performanceStats$measure == "Rsquared" & 
                                    performanceStats$lin_inter == "80:20" & 
                                    performanceStats$rel != 0.6 &
                                    performanceStats$pTrash == 50),]

subOverfitFig2$rel <- factor(subOverfitFig2$rel, levels = c(1, 0.8))
overfitFig2 <- plotPaperR2(subOverfitFig2, plotMeasure = subOverfitFig2$overfit, title = "", 
                           yLabel = expression(R[train]^2 - R[test]^2),
                           yMin = -0.1, yMax = 1.1)
(overfitFig2 <- themeFunctionPaper(overfitFig2, guides = T))
(overfitFig2 <- overfitFig2 + geom_text(data = subOverfitFig2[which(subOverfitFig2$R2 == 0.8 &
                                                                      subOverfitFig2$N == 300),],
                                        aes(x = N, y = 1.05, group = interaction(R2, model),
                                            color = model, label = model), position = position_dodge(3)))

# ggplot2::ggsave(filename = paste0(plotFolder, "/overfit_20:80_pTrash50_paper.png"),
#                 plot = overfitFig1,
#                 width = 13.63,
#                 height = 12.07,
#                 units = "in")
# ggplot2::ggsave(filename = paste0(plotFolder, "/overfit_80:20_pTrash50_paper.png"),
#                 plot = overfitFig2,
#                 width = 13.63,
#                 height = 12.07,
#                 units = "in")

################################################################################
# plot results for additional online material (full data)
################################################################################
relColValues <- c("white","#999999", "#444444") # rel = 0.8, rel = 1

plotOnlineR2 <- function(data, plotMeasure, title = "", yLabel = "",
                         yMin = NULL, yMax = NULL){
  ggplot(data,
         aes(x = N, y = plotMeasure, 
             group = interaction(R2, model), linetype = R2, color = model)) +
    geom_line(linewidth = 1, alpha = 0.4) +
    geom_point(aes(shape = model), size = 3) +
    scale_linetype_manual(name = expression(R[sim]^2), values = linetypeVec,
                          guide = guide_legend(override.aes = list(size = 2, alpha = 0.8))) +
    scale_color_manual(values = modelLineColors) +
    scale_shape_manual(values = modelPointShapes) +
    scale_y_continuous(limits = c(yMin, yMax), breaks=seq(
      ifelse(yMin %% 0.2 == 0, yMin, yMin + 0.1), yMax, 0.2)) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid2(rel ~ dgp,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(dgpPanelColors[1], 0.4)),
                                      element_rect(fill = alpha(dgpPanelColors[2], 0.4)),
                                      element_rect(fill = alpha(dgpPanelColors[3], 0.4))),
                  background_y = list(element_rect(fill = alpha(relColValues[3], 0.4)),
                                      element_rect(fill = alpha(relColValues[2], 0.4)),
                                      element_rect(fill = alpha(relColValues[1], 0.4))))) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "black",
               alpha = 0.4, linetype = linetypeVec[1]) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "black",
               alpha = 0.4, linetype = linetypeVec[2]) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "black",
               alpha = 0.4, linetype = linetypeVec[3]) +
    ylab(yLabel) +
    xlab("Sample Size (N)") +
    ggtitle(title) +
    guides(colour = "none", shape = "none")
}


linInterVec <- unique(performanceStats$lin_inter)
pTrashVec <- unique(performanceStats$pTrash)

for (iLinInter in seq_along(linInterVec)) {
  for (iTrash in seq_along(pTrashVec)) {
    
    # plot R²
    subFig1_full <- performanceStats[which(performanceStats$measure == "Rsquared" & 
                                             performanceStats$lin_inter == linInterVec[iLinInter] & 
                                             performanceStats$pTrash == pTrashVec[iTrash]),]
    
    subFig1_full$rel <- factor(subFig1_full$rel, levels = c(1, 0.8, 0.6))
    
    fig1 <- plotOnlineR2(subFig1_full, plotMeasure = subFig1_full$M_test, title = "", yLabel = expression(R[test]^2),
                         yMin = 0, yMax = 0.9)
    (fig1 <- themeFunctionPaper(fig1))
    (fig1 <- fig1 + geom_text(data = subFig1_full[which(subFig1_full$R2 == 0.8 &
                                                          subFig1_full$N == 300),],
                              aes(x = N, y = 0.85, group = interaction(R2, model),
                                  color = model, label = model), position = position_dodge(3)))
    
    ggplot2::ggsave(filename = paste0(plotFolder, "/R2_", 
                                      linInterVec[iLinInter], "_pTrash", 
                                      pTrashVec[iTrash],"_full.png"),
                    plot = fig1,
                    width = 13.63,
                    height = 12.07,
                    units = "in")
    
    # plot Overfit
    overfitFig1 <- plotOnlineR2(subFig1_full, plotMeasure = subFig1_full$overfit, title = "", 
                               yLabel = expression(R[train]^2 - R[test]^2),
                               yMin = -0.1, yMax = 1.1)
    (overfitFig1 <- themeFunctionPaper(overfitFig1))
    (overfitFig1 <- overfitFig1 + geom_text(data = subOverfitFig1[which(subOverfitFig1$R2 == 0.8 &
                                                                          subOverfitFig1$N == 300),],
                                            aes(x = N, y = 1.05, group = interaction(R2, model),
                                                color = model, label = model), position = position_dodge(3)))
    
    ggplot2::ggsave(filename = paste0(plotFolder, "/overfit_",
                                      linInterVec[iLinInter], "_pTrash", 
                                      pTrashVec[iTrash],"_full.png"),
                    plot = overfitFig1,
                    width = 13.63,
                    height = 12.07,
                    units = "in")
    
  }
}
