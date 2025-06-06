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
    # axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15))
}

################################################################################
# load data
################################################################################
dgpVec <- c("inter", "nonlinear", "pwlinear", "nonlinear3")
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
                                     labels = c("ENET_inter", "ENET_lin", "GBM", "RF"))

# plot overfit instead of train as train - test
performanceStats <- tidyr::pivot_wider(performanceStats, 
                                       names_from = "trainTest", values_from = c(M, SE, SD))
chr2num <- c("M_train", "M_test", "SE_train", "SE_test", "SD_train", "SD_test")
performanceStats[chr2num] <- lapply(performanceStats[chr2num], as.numeric)

performanceStats$overfit <- performanceStats$M_train - performanceStats$M_test

rm(list = ls(pattern = "^pT"))

################################################################################
# plot results for each model (compare DGPs)
################################################################################

modelVec <- unique(performanceStats$model)
pTrashVec <- unique(performanceStats$pTrash)

dgpColors <- c("#254441", "#FF6F59", "darkred", "#43AA8B")
linetypeVec <- c("dotted", "dashed", "solid")

plotDGPcomparison <- function(data, model, pTrash){
  ggplot(data,
         aes(x = N, y = M_test, 
             group = interaction(R2, dgp), linetype = R2, color = dgp)) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = linetypeVec) +
    scale_color_manual(values = dgpColors) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid2(rel ~ lin_inter + pTrash,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[3], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[2], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[1], 0.4))))) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "black",
               alpha = 0.4, linetype = linetypeVec[1]) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "black",
               alpha = 0.4, linetype = linetypeVec[2]) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "black",
               alpha = 0.4, linetype = linetypeVec[3]) +
    ylab("test R²") +
    xlab("N (increasing)") +
    ggtitle(paste0("Model: ", model, ", # Noise Variables: ", pTrash))
}

for (iModel in seq_along(modelVec)) {
  for (iNoise in seq_along(pTrashVec)) {
    
    performanceSub <- performanceStats[which(performanceStats$measure == "Rsquared" & 
                                               performanceStats$model == modelVec[iModel] & 
                                               performanceStats$pTrash == pTrashVec[iNoise]),]
    
    pTMP <- plotDGPcomparison(performanceSub, modelVec[iModel], pTrashVec[iNoise])
    pTMP <- themeFunction(pTMP)
    
    # save plots as files
    ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTest_", 
                                      modelVec[iModel], "_pTrash", 
                                      pTrashVec[iNoise], ".png"),
                    plot = pTMP,
                    width = 13.08,
                    height = 12.18,
                    units = "in")
        
  }
} 

################################################################################
# plot results for each model (compare models)
################################################################################

dgpVec <- unique(performanceStats$dgp)
pTrashVec <- unique(performanceStats$pTrash)

modelColors <- c("#0a2463", "#2a9d8f", "#8338ec", "#a44200")
linetypeVec <- c("dotted", "dashed", "solid")

plotDGPcomparison <- function(data, dgp, pTrash){
  ggplot(data,
         aes(x = N, y = M_test, 
             group = interaction(R2, model), linetype = R2, color = model)) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = linetypeVec) +
    scale_color_manual(values = modelColors) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid2(rel ~ lin_inter + pTrash,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesLin[3], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[2], 0.4)),
                                      element_rect(fill = alpha(colValuesLin[1], 0.4))))) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = "black",
               alpha = 0.4, linetype = linetypeVec[1]) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = "black",
               alpha = 0.4, linetype = linetypeVec[2]) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = "black",
               alpha = 0.4, linetype = linetypeVec[3]) +
    ylab("test R²") +
    xlab("N (increasing)") +
    ggtitle(paste0("DGP: ", dgp, ", # Noise Variables: ", pTrash))
}

for (iDGP in seq_along(dgpVec)) {
  for (iNoise in seq_along(pTrashVec)) {
    
    performanceSub <- performanceStats[which(performanceStats$measure == "Rsquared" & 
                                               performanceStats$dgp == dgpVec[iDGP] & 
                                               performanceStats$pTrash == pTrashVec[iNoise]),]
    
    pTMP <- plotDGPcomparison(performanceSub, dgpVec[iDGP], pTrashVec[iNoise])
    pTMP <- themeFunction(pTMP)
    
    # save plots as files
    ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTest_", 
                                      dgpVec[iDGP], "_pTrash", 
                                      pTrashVec[iNoise], ".png"),
                    plot = pTMP,
                    width = 13.08,
                    height = 12.18,
                    units = "in")
    
  }
} 


################################################################################
# (overview) plot train and test performance (RF)
################################################################################



performanceSubRF <- performanceStats[which(performanceStats$measure != "MAE" &
                                           performanceStats$measure != "RMSE"),]
performanceSubRF[,c("M_train", "SE_train", "SE_test", "SD_train", "SD_test")] <- list(NULL)
performanceSubRF <- tidyr::pivot_longer(performanceSubRF, c(M_test, overfit),
                                      names_to = "measures", values_to = "values")

# performanceSubRF <- rbind(performanceSubRF_sc, performanceSubRF)
# plot performance measures for train and test data
(pPerformTrainVStestRF <- ggplot(performanceSubRF,
# (pPerformTrainVStestRF <- ggplot(performanceSubRF[performanceSubRF$measures == "M_test", ],
                                  aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                                      group = R2, colour = R2)) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
    # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
    scale_color_manual(values = colValuesR2) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(measures + rel ~ lin_inter, labeller = label_both) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = colValuesR2[1],
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = colValuesR2[2],
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = colValuesR2[3],
               alpha = 0.4) +
    ylab("") +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle("R^2: Training vs. Test performance (RF)"))
 
# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTestRF_", data, "_sanityCheck.png"),
#                 plot = pPerformTrainVStestRF,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTestRF_", data, "_sanityCheck.png"),
#                 plot = pPerformTrainVStestRF,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")


################################################################################
# (overview) plot train and test performance (GBM)
################################################################################


performanceSubGBM <- performanceStats[which(performanceStats$measure != "MAE" &
                                           performanceStats$measure != "RMSE"),]
performanceSubGBM[,c("M_train", "SE_train", "SE_test", "SD_train", "SD_test")] <- list(NULL)
performanceSubGBM <- tidyr::pivot_longer(performanceSubGBM, c(M_test, overfit),
                                      names_to = "measures", values_to = "values")

# plot performance measures for train and test data
(pPerformTrainVStestGBM <- ggplot(performanceSubGBM,
                               aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                                   group = R2, colour = R2)) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
    # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
    scale_color_manual(values = colValuesR2) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(measures + rel ~ lin_inter, labeller = label_both) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = colValuesR2[1],
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = colValuesR2[2],
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = colValuesR2[3],
               alpha = 0.4) +
    ylab("") +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle("R^2: Training vs. Test performance (GBM)"))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTestGBM.png"),
#                 plot = pPerformTrainVStestGBM,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

################################################################################
# ENET results 
################################################################################



chr2fac <- c("model", "rel", "measure", "R2", "lin_inter", "fitInter")
performanceStatsENET[chr2fac] <- lapply(performanceStatsENET[chr2fac], factor)
chr2num <- c("M_train", "M_test", "SE_train", "SE_test", "SD_train", "SD_test")
performanceStatsENET[chr2num] <- lapply(performanceStatsENET[chr2num], as.numeric)

performanceSubENET <- performanceStatsENET[which(performanceStatsENET$measure != "MAE" &
                                               performanceStatsENET$measure != "RMSE"),]
performanceSubENET[,c("M_train", "SE_train", "SE_test", "SD_train", "SD_test")] <- list(NULL)
performanceSubENET <- tidyr::pivot_longer(performanceSubENET, c(M_test, overfit),
                                      names_to = "measures", values_to = "values")

# (pPerformTrainVStestENETwo_sc <- ggplot(performanceSubENET[performanceSubENET$measures == "M_test",],
#                                   aes(x = interaction(pTrash, N, sep = " x "), y = values, 
#                                       group = R2, colour = R2)) +
#     geom_point() +
#     geom_line() +
#     scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
#     # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
#     scale_color_manual(values = colValuesR2) +
#     geom_hline(aes(yintercept = 0)) +
#     facet_grid(measures + rel ~ lin_inter, labeller = label_both) +
#     geom_hline(yintercept = setParam$dgp$Rsquared[1], col = colValuesR2[1],
#                alpha = 0.4) +
#     geom_hline(yintercept = setParam$dgp$Rsquared[2], col = colValuesR2[2],
#                alpha = 0.4) +
#     geom_hline(yintercept = setParam$dgp$Rsquared[3], col = colValuesR2[3],
#                alpha = 0.4) +
#     ylab("") +
#     xlab("pTrash (decreasing) x N (increasing)") +
#     ggtitle("R^2: Training vs. Test performance (GBM)") +
#     theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "lightgrey"), 
#           panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "lightgrey"),
#           panel.background = element_rect(color = "white", fill = "white"),
#           axis.text.y = element_text(size = 20),
#           # axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
#           axis.text.x = element_text(size = 15),
#           axis.title.x = element_text(size = 20),
#           axis.title.y = element_text(size = 20),
#           strip.text.x = element_text(size = 15),
#           strip.text.y = element_text(size = 15)))

head(performanceSubENET)
head(performanceSubGBM)
head(performanceSubRF)
# code if interaction were fitted or not
unique(performanceStatsENET$model)
performanceStatsENET$fitInter <- ifelse(performanceStatsENET$model == "ENETw", 1, 0)
performanceSubGBM$fitInter <- rep(2, dim(performanceSubGBM)[1])
performanceSubRF$fitInter <- rep(3, dim(performanceSubRF)[1])

colnames(performanceSubRF)
colnames(performanceSubGBM)
colnames(performanceSubENET)

performanceData <- rbind(performanceSubRF, performanceSubGBM, performanceSubENET)
# performanceData$fitInter <- plyr::mapvalues(performanceData$fitInter, 
#                                             from=c(0,1,2), 
#                                             to=c("ENET - lin","ENET - inter","GBM"))

unique(performanceData$fitInter)
performanceData$fit <- ifelse(performanceData$fitInter >= 2, "treeBased", "reg")
unique(performanceData$fit)
unique(performanceData$model)

# save(performanceData, file = paste0("results/", data, "/dependentMeasures/rSquaredData_stats.rda"))

###############################################################################
# (overview) plot train and test performance (GBM + ENET)
################################################################################
# plot performance measures for train and test data

(pPerformTrainVStest <- ggplot(performanceData,
                               aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                                   group = interaction(R2, model, fit), colour = R2,
                                   linetype = model, shape = fit)) +
   geom_point() +
   geom_line() +
   scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
   scale_shape_manual(values = c(16, 8)) +
   # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
   scale_color_manual(values = colValuesR2) +
   geom_hline(aes(yintercept = 0)) +
   facet_grid(measures + rel ~ lin_inter, labeller = label_both) +
   geom_hline(yintercept = setParam$dgp$Rsquared[1], col = colValuesR2[1],
              alpha = 0.4) +
   geom_hline(yintercept = setParam$dgp$Rsquared[2], col = colValuesR2[2],
              alpha = 0.4) +
   geom_hline(yintercept = setParam$dgp$Rsquared[3], col = colValuesR2[3],
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTest_ENTvsGBMvsRF.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTest_ENTvsGBMvsRF_nonlinear.png"),
#                 plot = pPerformTrainVStest,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

(pPerformTest <- ggplot(performanceData[performanceData$measures == "M_test",],
                               # aes(x = interaction(pTrash, N, sep = " x "), y = values, 
                                   aes(x = N, y = values, 
                                   group = interaction(R2, model, fit), colour = R2,
                                   linetype = model, shape = fit)) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
    scale_shape_manual(values = c(16, 8)) +
    # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
    scale_color_manual(values = colValuesR2) +
    geom_hline(aes(yintercept = 0)) +
    # facet_grid(rel ~ lin_inter, labeller = label_both) +
    facet_grid(rel ~ pTrash + lin_inter, labeller = label_both) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = colValuesR2[1],
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = colValuesR2[2],
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = colValuesR2[3],
               alpha = 0.4) +
    ylab("") +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle("R^2 test performance") +
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

# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTest_ENTvsGBMvsRF.png"),
ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTest_ENTvsGBMvsRF_pwlinear.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTest_ENTvsGBMvsRF_nonlinear.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTest_ENTvsGBMvsRF_inter.png"),
                plot = pPerformTest,
                width = 13.08,
                height = 12.18,
                units = "in")

###############################################################################
# (overview) plot train and test performance
################################################################################
interData <- loadRData("results/inter/dependentMeasures/rSquaredData_stats.rda")
nlData <- loadRData("results/nonlinear/dependentMeasures/rSquaredData_stats.rda")

fullData <- rbind(interData, nlData)
(pPerformTest <- ggplot(fullData[fullData$measures == "M_test" & 
                                   fullData$fit == "treeBased",],
                        aes(x = N, y = values, 
                            group = interaction(R2, dgp, model), colour = R2,
                            linetype = dgp, shape = model)) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
    scale_shape_manual(values = c(16, 8)) +
    # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
    scale_color_manual(values = colValuesR2) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(rel + pTrash ~ lin_inter, labeller = label_both) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = colValuesR2[1],
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = colValuesR2[2],
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = colValuesR2[3],
               alpha = 0.4) +
    ylab("") +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle("R^2 test performance") +
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

(pPerformTest <- ggplot(fullData[fullData$measures == "M_test" & 
                                   fullData$fit == "reg",],
                        aes(x = N, y = values, 
                            group = interaction(R2, dgp, model), colour = R2,
                            linetype = dgp, shape = model)) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
    scale_shape_manual(values = c(16, 8)) +
    # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
    scale_color_manual(values = colValuesR2) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(rel + pTrash ~ lin_inter, labeller = label_both) +
    geom_hline(yintercept = setParam$dgp$Rsquared[1], col = colValuesR2[1],
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[2], col = colValuesR2[2],
               alpha = 0.4) +
    geom_hline(yintercept = setParam$dgp$Rsquared[3], col = colValuesR2[3],
               alpha = 0.4) +
    ylab("") +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle("R^2 test performance") +
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
