# compare selected variables between GBM and ENET by using PVI (model agnostic)
# for linear effects: 
#     comparison between ENETw and ENETwo is fair
#     GBM might benefit from interaction effects (confounded in "linear effect")

# load parameters and helper functions 
source("utils/setParameters.R")
source("utils/analysisTools.R")

# ANOVA
library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²

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

# interactions written as Var1.Var2 instead of Var1:Var2
interEffectsDot <- gsub("\\:", "\\.", setParam$dgp$interEffects)

# read in data
listDir <- dir(resFolder)
dataList <- listDir[stringr::str_detect(listDir, "^pvi")]
models <- stringr::str_extract(dataList, "_[:alpha:]*.rda$")
models <- stringr::str_sub(models, start = 2L, end = -5)
for (iData in seq_len(length(dataList))){
  objectName <- paste0("pvi", models[iData])
  assign(objectName, loadRData(paste0(resFolder, "/", dataList[iData])))
}

# pull data from nested list of all results (fullData)
pviENETw <- rbindSingleResults(pviENETw)
pviENETwo <- rbindSingleResults(pviENETwo)
pviGBM <- rbindSingleResults(pviGBM)
################################################################################
# evaluate frequency stats
################################################################################
# how often are predictors "selected" in model? 
#   -> i.e., how often is a certain variable relevant for prediction (pvi > 1)
#       all variables with pvi <= 1 are already removed
# ! pvi for ENETw includes interactions
#   otherwise {ENETwo & GBM} pvi does not inlcude interactions

# 1. calculate number of true positives for each sample
# how many (frequently do) variables **with** simulated linear effects (setParam$dgp$linEffects) 
#     affect (pvi > 1) the prediction? 
#     -> = true positive (linear)
pviENETwo$linTP <- ifelse(pviENETwo$pviRank %in% setParam$dgp$linEffects, 1, 0)
pviENETw$linTP <- ifelse(pviENETw$pviRank %in% setParam$dgp$linEffects, 1, 0)
pviGBM$linTP <- ifelse(pviGBM$pviRank %in% setParam$dgp$linEffects, 1, 0)

# 2. calculate number of false positives for each sample
## extract false positive main effects
# how many (frequently do) variables **without** simulated main effects (setParam$dgp$linEffects) 
#     affect (pvi > 1) the prediction? 
#     -> = false positive (linear)
# ENETwo and GBM without interactions in PVI, therefore only one version of FP necessary
pviENETwo$linFP <- ifelse(pviENETwo$pviRank %in% setParam$dgp$linEffects, 0, 1)
pviGBM$linFP <- ifelse(pviGBM$pviRank %in% setParam$dgp$linEffects, 0, 1)

# ENETw with interactions in PVI, therefore test if exclusion of interactions does 
#     make a difference
# get all possible interactions to exclude them from FP calculation
# check if all possible interactions are in data to avoid generating them 
#   -> should be TRUE
length(grep("\\.", unique(pviENETw$pviRank))) == 
  choose(max(setParam$dgp$pTrash) + length(setParam$dgp$linEffects), setParam$dgp$interDepth)

allInterVars <- unique(pviENETw$pviRank)[grep("\\.", unique(pviENETw$pviRank))]

# every PVI > 1 that is ...
#   ... a variable with a simulated main effect or an interaction **with- or without** simulated effects receives a 0
#   ... a variable without a simulated main effect and is not an interaction receives 1
# -> false positives (linear) with any interactions excluded
pviENETw$linFP <- ifelse(pviENETw$pviRank %in% setParam$dgp$linEffects |
                             pviENETw$pviRank %in% allInterVars, 0, 1)

# aggregate number of true positives and false negatives within each sample
linENETwo <- aggregate(cbind(linTP, linFP) ~ sample + idxCondLabel + model + N_pTrash, 
          data = pviENETwo, sum)
linENETw <- aggregate(cbind(linTP, linFP) ~ sample + idxCondLabel + model + N_pTrash, 
          data = pviENETw, sum)
linGBM <- aggregate(cbind(linTP, linFP) ~ sample + idxCondLabel + model + N_pTrash, 
          data = pviGBM, sum)

# id 
# independent observations for samples {1:100} in different simulated conditions
#   that all have the same sample numbers
# add ID variable for models separately since models are within factor in within ANOVA
linENETwo$ID <- seq_len(dim(linENETwo)[1])
linENETw$ID <- seq_len(dim(linENETw)[1])
linGBM$ID <- seq_len(dim(linGBM)[1])

# 3. calculate positive predictive value (TP / (TP + FP)) for every sample
# calculate positive predictive value from TP and FP
#     TP / (TP + FP)
linENETwo$PPV <- linENETwo$linTP / (linENETwo$linTP + linENETwo$linFP) 
linENETw$PPV <- linENETw$linTP / (linENETw$linTP + linENETw$linFP) 
linGBM$PPV <- linGBM$linTP / (linGBM$linTP + linGBM$linFP) 

# separate N_pTrash
linENETwo <- idx2infoNew(linENETwo)
linENETw <- idx2infoNew(linENETw)
linGBM <- idx2infoNew(linGBM)

# to do: 
# 1. draw matrix with TP, FP, TN, FN for all models in terms of PVI measure
# ENETw & GBM:
#   TP = extracted and simulated effect in DGP
#   FP = extracted and no simulated effect in DGP 
#   TN = not extracted and no simulated effect in DGP; Anzahl an pTrash - FP
#   FN = not extracted but simulated in DGP; Anzahl linearer Effekte - TP
# ENETwo
#   TP = extracted and simulated (linear) effect in DGP
#   FP = extracted and no simulated effect in DGP
#   FP_woInter = extracted and no simulated effect in DGP (true interaction are excluded) 
#   TN = not extracted and no simulated effect in DGP; 
#         (Anzahl an pTrash + alle Interaktionen) - FP_woInter (true interactions are excluded)
#   FN = not extracted but simulated linear effects in DGP; Anzahl linearer Effekte - TP

# 2a. add measures...
#   ... linTN = trash variable (or interaction effect) that is correctly not extracted

# ENETwo and GBM without interactions in PVI, therefore only one version of FP necessary
linENETwo$linTN <- as.numeric(linENETwo$pTrash) - linENETwo$linFP
linGBM$linTN <- as.numeric(linGBM$pTrash) - linGBM$linFP

# all possible variables without simulated main effects - false positive main effects 
#    i.e., ignore interactions all together 
#   as this only counts main effects as true negative options
linENETw$linTN <- as.numeric(linENETw$pTrash) - # all possible interactions
  linENETw$linFP

#   ... linFN = linear effect that is not extracted
linENETwo$linFN <- length(setParam$dgp$linEffects) - linENETwo$linTP
linGBM$linFN <- length(setParam$dgp$linEffects) - linGBM$linTP
linENETw$linFN <- length(setParam$dgp$linEffects) - linENETw$linTP

# 2b. add measures...
#   ... accuracy = (TP + TN) / (TP + TN + FP + FN)
# getAccuracy <- function(data) {
#   data$ACC <- (data$linTP + data$linTN) / (data$linTP + data$linTN + data$linFP + data$linFN)
#   return(data)
# }
# linENETwo <- getAccuracy(linENETwo)
# linGBM <- getAccuracy(linGBM)
linENETwo$ACC <- (linENETwo$linTP + linENETwo$linTN) / (linENETwo$linTP + linENETwo$linTN + 
                                                          linENETwo$linFP + linENETwo$linFN)
linGBM$ACC <- (linGBM$linTP + linGBM$linTN) / (linGBM$linTP + linGBM$linTN + 
                                                 linGBM$linFP + linGBM$linFN)

linENETw$ACC <- (linENETw$linTP + linENETw$linTN) / (linENETw$linTP + linENETw$linTN + 
                                                          linENETw$linFP + linENETw$linFN)

#   ... sensitivity = TP / (TP + FN)
linENETwo$sensitivity <- linENETwo$linTP / (linENETwo$linTP + linENETwo$linFN)
linGBM$sensitivity <- linGBM$linTP / (linGBM$linTP + linGBM$linFN)

linENETw$sensitivity <- linENETw$linTP / (linENETw$linTP + linENETw$linFN)

#   ... specificity = TN / (TN + FP)
linENETwo$specificity <- linENETwo$linTN / (linENETwo$linTN + linENETwo$linFP)
linGBM$specificity <- linGBM$linTN / (linGBM$linTN + linGBM$linFP)
linENETw$specificity <- linENETw$linTN / (linENETw$linTN + linENETw$linFP)

#   ... balanced accuracy = (sensitivity + specificity) / 2 
linENETwo$balACC <- (linENETwo$specificity + linENETwo$sensitivity) / 2
linGBM$balACC <- (linGBM$specificity + linGBM$sensitivity) / 2
linENETw$balACC <- (linENETw$specificity + linENETw$sensitivity) / 2

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model")
linENETwo[col2fac] <- lapply(linENETwo[col2fac], factor)
linENETw[col2fac] <- lapply(linENETw[col2fac], factor)
linGBM[col2fac] <- lapply(linGBM[col2fac], factor)

# # check
# range(linENETw$specificity)

# # save data for sensitivity and specificity result plots
# save(linENETw, linENETwo, linGBM, file = paste0(resFolder, "/relFrequencyMeasures.rda"))
################################################################################
# run ANOVAs for ...
#     ... different models {ENETwo, ENETw, GBM}
#     ... different dependent variables {PPV, linTP, linFP_woTInter}
################################################################################
# run mixed anova 
#     ... id = sample but ID (independent samples between simulated conditions)
#     ... dv = {ppv}
#     ... between: 3 x 2 x 3 x 3 x 3
#         N (3) {100, 300, 1000}
#         pTrash (2) {10, 50}
#         rel (3) {0.6, 0.8, 1}
#         R2 (3) {0.2, 0.5, 0.8}
#         lin_inter (3) {0.2_0.8, 0.5_0.5, 0.8_0.2}
#     ... within: 
#         model {Enetw; Enetwo, GBM}
# merge data
linPPV <- rbind(linENETw, linENETwo, linGBM)
str(linPPV)

dvVec <- c("linTP", "linFP", "PPV", "ACC", 
           "sensitivity", "specificity", "balACC")
modVec <- c("ENETwo", "ENETw", "GBM")
# modVec <- "GBM"
for (iModel in modVec) {
  for (iDV in dvVec) {
    # only between factor ANOVA (only within factor was model)
    # between-sample ANOVA for every model separately
    anovaObjectName <- paste0("anova", iModel, "_", iDV)
    tmp <- aov_ez(id = "ID",
                  dv = iDV,
                  data = linPPV[linPPV$model == iModel,],
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
    tmpEta2 <- cbind(tmpEta2,
                     model = iModel,
                     dv = iDV)
    
    # which higher order interactions do we need to illustrate to report simulation results?
    etaObjectName <- paste0("eta2", iModel, "_", iDV)
    #assign(etaObjectName, tmpEta2[order(tmpEta2$Eta2_generalized, decreasing = T),]) 
    assign(etaObjectName, tmpEta2) 
  }
}

eta2Thresh <- 0.01

##### plot all models separately and join into one plot #####
plotEta2_4eachModel <- function(data, eta2Thresh, fillVal = "grey") {
  
  # sort parameters by importance
  sortIdx <- order(data$Eta2_generalized, decreasing = T)
  data$Parameter <- factor(data$Parameter, 
                                levels = data$Parameter[sortIdx])
  
  # plot generalized eta2 for every model
  tmp <- ggplot(data[data$Eta2_generalized > eta2Thresh, ], 
         aes(x = Parameter, y = Eta2_generalized)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single"),
             fill = fillVal) +
    geom_text(aes(label=round(Eta2_generalized, 2)), 
              angle = 90, hjust = 1.5, vjust=0.5, 
              #angle = 0, vjust=1.5, 
              position = position_dodge(width = .9), 
              color="black", size=3.5)+
    ylab("generalisiertes eta^2") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15))
  
  # plotName <- paste("pEta2_", model)
  # assign(plotName, tmp)
  return(tmp)
}

plotEta2 <- function(dataENETw, dataENETwo, dataGBM, eta2Thresh) {
  pEta2_ENETw <- plotEta2_4eachModel(dataENETw, eta2Thresh, "#006600") + 
    ggtitle("ENETw") + # for the main title
    theme(plot.title = element_text(hjust = 0.5, size = 25))
  pEta2_ENETwo <- plotEta2_4eachModel(dataENETwo, eta2Thresh, "#009999") +
    ggtitle("ENETwo") + # for the main title
    theme(plot.title = element_text(hjust = 0.5, size = 25)) 
  pEta2_GBM <- plotEta2_4eachModel(dataGBM, eta2Thresh, "#990000") +
    ggtitle("GBM") + # for the main title
    theme(plot.title = element_text(hjust = 0.5, size = 25))
  
  (tmpP <- pEta2_GBM | pEta2_ENETw | pEta2_ENETwo)
  # plotName <- paste0("pEta2_", DV)
  # assign(plotName, tmpP) 
  return(tmpP)
}

eta2Thresh <- 0.01

# sensitivity does not change! depending on calculation
pEta2_sensitivity <- plotEta2_4eachModel(eta2ENETw_sensitivity, eta2Thresh, "#006600") + 
  ggtitle("sensitivity") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pEta2_ACC <- plotEta2_4eachModel(eta2ENETw_ACC, eta2Thresh, "#006600") + 
  ggtitle("ACC") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pEta2_balACC <- plotEta2_4eachModel(eta2ENETw_balACC, eta2Thresh, "#006600") + 
  ggtitle("balACC") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pEta2_specificity <- plotEta2_4eachModel(eta2ENETw_specificity, eta2Thresh, "#006600") + 
  ggtitle("specificity") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

(pCheck_part1 <- pEta2_specificity + pEta2_balACC + pEta2_ACC +
  plot_layout(ncol = 2))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/detectMains/checkpviMainMeasures_part1.png"),
#                 plot = pCheck_part1,
#                 width = 13.68,
#                 height = 15,
#                 units = "in")

pEta2_TP <- plotEta2_4eachModel(eta2ENETw_linTP, eta2Thresh, "#006600") + 
  ggtitle("TP") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pEta2_FP <- plotEta2_4eachModel(eta2ENETw_linFP, eta2Thresh, "#006600") + 
  ggtitle("FP") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pEta2_PPV <- plotEta2_4eachModel(eta2ENETw_PPV, eta2Thresh, "#006600") + 
  ggtitle("PPV") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

(pCheck_part2 <- pEta2_TP + pEta2_sensitivity + pEta2_FP + pEta2_PPV + 
  plot_layout(ncol = 2))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/detectMains/checkpviMainMeasures_part2.png"),
#                 plot = pCheck_part2,
#                 width = 13.68,
#                 height = 15,
#                 units = "in")

# sensitivity
pEta2_sensitivity <- plotEta2(eta2ENETw_sensitivity, eta2ENETwo_sensitivity, eta2GBM_sensitivity, eta2Thresh)

# ggplot2::ggsave(filename = paste0(plotFolder, "/detectMains/betweenANOVA_pviMainSensitivity.png"),
#                 plot = pEta2_sensitivity,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

# specificity
pEta2_specificity <- plotEta2(eta2ENETw_specificity, eta2ENETwo_specificity, eta2GBM_specificity, eta2Thresh)

# ggplot2::ggsave(filename = paste0(plotFolder, "/detectMains/betweenANOVA_pviMainSpecificity.png"),
#                 plot = pEta2_specificity,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

# remove saved plots
rm(list = c(ls(pat = "^pEta2_")))
# remove ANOVA result data
rm(list = c(ls(pat="^anova")))
# remove eta2 data
rm(list = c(ls(pat="^eta2")))
gc()
