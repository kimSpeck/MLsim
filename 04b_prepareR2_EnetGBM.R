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


resFolder <- "results/finalResults/dependentMeasures"

# read in data
listDir <- dir(resFolder)
dataList <- listDir[stringr::str_detect(listDir, "^performT")]
models <- stringr::str_extract(dataList, "_[:alpha:]*.rda$")
models <- stringr::str_sub(models, start = 2L, end = -5)
tvst <- stringr::str_extract(dataList, "(Test|Train)") 
for (iData in seq_len(length(dataList))) {
  objectName <- paste0("p", tvst[iData], "_", models[iData])
  assign(objectName, loadRData(paste0(resFolder, "/", dataList[iData])))
}

################################################################################
# (overview) plot train and test performance
################################################################################
# change list of  matrices to data.frame
pTest_GBM <- rbindSingleResults(pTest_GBM)
pTrain_GBM <- rbindSingleResults(pTrain_GBM)

# variables from rownames to own column to work with variable information
pTrain_GBM <- rowNames2col(pTrain_GBM, "measure")
pTrain_GBM$measure <- paste0(pTrain_GBM$measure, "_train")

pTest_GBM <- rowNames2col(pTest_GBM, "measure")
pTest_GBM$measure[which(pTest_GBM$measure == "Rsq_test")] <- "Rsquared_test"

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
pTrain_GBM <- idx2infoNew(pTrain_GBM)
pTest_GBM <- idx2infoNew(pTest_GBM) 

# merge performance Train and performance Test
performanceStats <- rbind(pTrain_GBM, pTest_GBM)
performanceStats <- tidyr::separate(performanceStats, measure, c("measure", "trainTest"), sep = "_")

# change type of columns or specific entry details to prepare plotting  
performanceStats$N <- factor(performanceStats$N, levels = setParam$dgp$N)
performanceStats$pTrash <- factor(performanceStats$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))
str(performanceStats)
chr2fac <- c("rel", "measure", "R2", "lin_inter")
performanceStats[chr2fac] <- lapply(performanceStats[chr2fac], factor)


# plot overfit instead of train as train - test
performanceStats <- tidyr::pivot_wider(performanceStats, 
                                       names_from = "trainTest", values_from = c(M, SE, SD))
chr2num <- c("M_train", "M_test", "SE_train", "SE_test", "SD_train", "SD_test")
performanceStats[chr2num] <- lapply(performanceStats[chr2num], as.numeric)

performanceStats$overfit <- performanceStats$M_train - performanceStats$M_test
performanceSub <- performanceStats[which(performanceStats$measure != "MAE" &
                                           performanceStats$measure != "RMSE"),]
performanceSub[,c("M_train", "SE_train", "SE_test", "SD_train", "SD_test")] <- list(NULL)
performanceSub <- tidyr::pivot_longer(performanceSub, c(M_test, overfit),
                                      names_to = "measures", values_to = "values")

# plot performance measures for train and test data
(pPerformTrainVStestGBM <- ggplot(performanceSub,
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/performanceTrainTestGBM.png"),
#                 plot = pPerformTrainVStestGBM,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

################################################################################
# ENET results 
################################################################################
# change list of  matrices to data.frame
pTest_ENETw <- rbindSingleResults(pTest_ENETw)
pTest_ENETwo <- rbindSingleResults(pTest_ENETwo)
pTrain_ENETw <- rbindSingleResults(pTrain_ENETw)
pTrain_ENETwo <- rbindSingleResults(pTrain_ENETwo)

# variables from rownames to own column to work woth variable information
pTrain_ENETw <- rowNames2col(pTrain_ENETw, "measure")
pTrain_ENETw$measure <- paste0(pTrain_ENETw$measure, "_train")

pTrain_ENETwo <- rowNames2col(pTrain_ENETwo, "measure")
pTrain_ENETwo$measure <- paste0(pTrain_ENETwo$measure, "_train")

pTest_ENETw <- rowNames2col(pTest_ENETw, "measure")
pTest_ENETw$measure[which(pTest_ENETw$measure == "Rsq_test")] <- "Rsquared_test"

pTest_ENETwo <- rowNames2col(pTest_ENETwo, "measure")
pTest_ENETwo$measure[which(pTest_ENETwo$measure == "Rsq_test")] <- "Rsquared_test"

# merge performance Train and performance Test
performanceStatsENET <- rbind(pTrain_ENETw, pTrain_ENETwo, pTest_ENETw, pTest_ENETwo)
performanceStatsENET <- tidyr::separate(performanceStatsENET, measure, c("measure", "trainTest"), sep = "_")
# unique(performanceStatsENET$trainTest)
# rm(pTrain_ENETw, pTrain_ENETwo, pTest_ENETw, pTest_ENETwo)

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
performanceStatsENET <- idx2infoNew(performanceStatsENET)

# code if interaction were fitted or not
unique(performanceStatsENET$model)
performanceStatsENET$fitInter <- ifelse(performanceStatsENET$model == "ENETw", 1, 0)

# change type of columns or specific entry details to prepare plotting  
str(performanceStatsENET)
performanceStatsENET$N <- factor(performanceStatsENET$N, levels = setParam$dgp$N)
performanceStatsENET$pTrash <- factor(performanceStatsENET$pTrash, 
                                      levels = sort(setParam$dgp$pTrash, decreasing = T))

# plot overfit instead of train as train - test
performanceStatsENET <- tidyr::pivot_wider(performanceStatsENET, 
                                           names_from = trainTest, values_from = c(M, SE, SD))

str(performanceStatsENET)
chr2fac <- c("model", "rel", "measure", "R2", "lin_inter", "fitInter")
performanceStatsENET[chr2fac] <- lapply(performanceStatsENET[chr2fac], factor)
chr2num <- c("M_train", "M_test", "SE_train", "SE_test", "SD_train", "SD_test")
performanceStatsENET[chr2num] <- lapply(performanceStatsENET[chr2num], as.numeric)

performanceStatsENET$overfit <- performanceStatsENET$M_train - performanceStatsENET$M_test
performanceSubENET <- performanceStatsENET[which(performanceStatsENET$measure != "MAE" &
                                               performanceStatsENET$measure != "RMSE"),]
performanceSubENET[,c("M_train", "SE_train", "SE_test", "SD_train", "SD_test")] <- list(NULL)
performanceSubENET <- tidyr::pivot_longer(performanceSubENET, c(M_test, overfit),
                                      names_to = "measures", values_to = "values")

head(performanceSubENET)
head(performanceSub)
performanceSub$fitInter <- rep(2, dim(performanceSub)[1])

colnames(performanceSub)
colnames(performanceSubENET)
performanceData <- rbind(performanceSub, performanceSubENET)
# performanceData$fitInter <- plyr::mapvalues(performanceData$fitInter, 
#                                             from=c(0,1,2), 
#                                             to=c("ENET - lin","ENET - inter","GBM"))

unique(performanceData$fitInter)
performanceData$fit <- ifelse(performanceData$fitInter != 2, "enet", "gbm")
unique(performanceData$fit)
unique(performanceData$model)

# save(performanceData, file = "results/finalResults/dependentMeasures/rSquaredData_stats.rda")

