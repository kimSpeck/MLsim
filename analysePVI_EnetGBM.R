# load parameters and helper functions 
source("setParameters.R")
source("analysisTools.R")

library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²
library(ggplot2)
library(patchwork)
library(see)

colValues <- c("green3", "darkblue", "darkmagenta")

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

# load results files
resFolder <- "results/finalResults/dependentMeasures" 

# interactions written as Var1.Var2 instead of Var1:Var2
interEffectsDot <- gsub("\\:", "\\.", setParam$dgp$interEffects)

################################################################################
# ANOVA - permutation variable importance
################################################################################
# compare selected variables between GBM and ENET by using PVI (model agnostic)
# for linear effects: 
#     comparison between ENETw and ENETwo is fair
#     GBM might benefit from interaction effects (confounded in "linear effect")
#     how to count interactions in ENETw? 
#       -> two versions: ...
#                       ... exclude true interactions, but include false interactions
#                       ... remove both, true and false positive, interactions 
#     thus, maybe switch to between ANOVA for each model
# for interaction effects:
#     exclude ENETwo; ENETw has interactions in pvi; GBM has interactions in interStrength 
#     run ANOVA for each model?

# 1. calculate number of true positives for each sample
# 2. calculate number of false positives for each sample
# 3. calculate positive predictive value (TP / (TP + FP)) for every sample
# 4. run mixed anova 
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
# how often are predictors "selected" in model? 
#   -> i.e., how often is a certain variable relevant for prediction (pvi > 1)
#       all variables with pvi <= 1 are already removed
# ! pvi for ENETw includes interactions
#   otherwise {ENETwo & GBM} pvi does not inlcude interactions

# how many (frequently do) variables **with** simulated linear effects (setParam$dgp$linEffects) 
#     affect (pvi > 1) the prediction? 
#     -> = true positive (linear)
pviENETwo$linTP <- ifelse(pviENETwo$pviRank %in% setParam$dgp$linEffects, 1, 0)
pviENETw$linTP <- ifelse(pviENETw$pviRank %in% setParam$dgp$linEffects, 1, 0)
pviGBM$linTP <- ifelse(pviGBM$pviRank %in% setParam$dgp$linEffects, 1, 0)

# how many (frequently do) variables **without** simulated linear effects (setParam$dgp$linEffects) 
#     affect (pvi > 1) the prediction? 
#     -> = false positive (linear)
pviENETwo$linFP <- ifelse(pviENETwo$pviRank %in% setParam$dgp$linEffects, 0, 1)
pviENETw$linFP <- ifelse(pviENETw$pviRank %in% setParam$dgp$linEffects, 0, 1)
pviGBM$linFP <- ifelse(pviGBM$pviRank %in% setParam$dgp$linEffects, 0, 1)

# false positives (linear) with true interactions excluded
# with ENETwo and GBM not having PVI measures for interactions:
#   linFP_woTInter is exactly the same as linFP for ENETwo and GBM
# pviENETwo$linFP_woTInter <- ifelse(pviENETwo$pviRank %in% setParam$dgp$linEffects & 
#                                      !(pviENETwo$pviRank %in% interEffectsDot), 0, 1)
pviENETwo$linFP_woTInter <- pviENETwo$linFP 
pviGBM$linFP_woTInter <- pviGBM$linFP
pviENETw$linFP_woTInter <- ifelse(pviENETw$pviRank %in% setParam$dgp$linEffects |
                                    pviENETw$pviRank %in% interEffectsDot, 0, 1)


# aggregate number of true positives and false negatives
linENETwo <- aggregate(cbind(linTP, linFP, linFP_woTInter) ~ sample + idxCondLabel + model + N_pTrash, 
          data = pviENETwo, sum)
linENETw <- aggregate(cbind(linTP, linFP, linFP_woTInter) ~ sample + idxCondLabel + model + N_pTrash, 
          data = pviENETw, sum)
linGBM <- aggregate(cbind(linTP, linFP, linFP_woTInter) ~ sample + idxCondLabel + model + N_pTrash, 
          data = pviGBM, sum)

# predictors, trash variables and all possible interactions
setParam$dgp$nModelPredictors

# id 
# independent observations for samples {1:100} in different simulated conditions
#   that all have the same sample numbers
# add ID variable for models separately since models are within factor in within ANOVA
linENETwo$ID <- seq_len(dim(linENETwo)[1])
linENETw$ID <- seq_len(dim(linENETw)[1])
linGBM$ID <- seq_len(dim(linGBM)[1])

# merge data
linPPV <- rbind(linENETw, linENETwo, linGBM)

# calculate positive predictive value from TP and FP
#     TP / (TP + FP)
linPPV$PPV <- linPPV$linTP / (linPPV$linTP + linPPV$linFP) 

linENETwo$PPV <- linENETwo$linTP / (linENETwo$linTP + linENETwo$linFP) 
linENETw$PPV <- linENETw$linTP / (linENETw$linTP + linENETw$linFP_woTInter) 
linGBM$PPV <- linGBM$linTP / (linGBM$linTP + linGBM$linFP) 

# separate N_pTrash
linPPV <- idx2infoNew(linPPV)

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
linENETwo$linTN <- as.numeric(linENETwo$pTrash) - linENETwo$linFP
linGBM$linTN <- as.numeric(linGBM$pTrash) - linGBM$linFP

linENETw$linTN <- (as.numeric(linENETw$pTrash) + # trash variables
                     choose(as.numeric(linENETw$pTrash) + length(setParam$dgp$linEffects), 
                            setParam$dgp$interDepth)) - # all possible interactions
  length(setParam$dgp$linEffects) -  # remove true interactions
  linENETw$linFP_woTInter

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

################################################################################
# run ANOVAs for ...
#     ... different models {ENETwo, ENETw, GBM}
#     ... different dependent variables {PPV, linTP, linFP_woTInter}
################################################################################
str(linPPV)
col2fac <- c("N", "pTrash", "R2", "rel", "lin_inter")
linPPV[col2fac] <- lapply(linPPV[col2fac], factor)

dvVec <- c("PPV", "linTP", "linFP_woTInter")
modVec <- c("ENETwo", "ENETw", "GBM")
for (iModel in modVec) {
  for (iDV in dvVec) {
    # only between factor ANOVA (only within factor was model)
    anovaObjectName <- paste0("anovaPPV", iModel, "_", iDV)
    tmp <- aov_ez(id = "ID",
                  dv = iDV,
                  # dv = "PPV",
                  # dv = "linTP",
                  # dv = "linFP_woTInter",
                  data = linPPV[linPPV$model == iModel,],
                  between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"))
    assign(anovaObjectName, tmp)
    
    tmpEta2 <- eta_squared(
      tmp, # fitted model
      partial = FALSE, # not partial!
      generalized = TRUE, # generalized eta squared
      ci = 0.95,
      verbose = TRUE)
    
    tmpEta2 <- tmpEta2[order(tmpEta2$Eta2_generalized, decreasing = T),]
    tmpEta2 <- cbind(tmpEta2,
                     model = iModel,
                     dv = iDV)
    
    # sort generalized eta-squared results
    # which higher order interactions do we need to illustrate to report simulation results?
    etaObjectName <- paste0("eta2", iModel, "_", iDV)
    #assign(etaObjectName, tmpEta2[order(tmpEta2$Eta2_generalized, decreasing = T),]) 
    assign(etaObjectName, tmpEta2) 
  }
}

# Datensätze zusammenführen
eta2Thresh <- 0.01
# dim(eta2ENETw[eta2ENETw$Eta2_generalized >= eta2Thresh,])
# dim(eta2ENETwo[eta2ENETwo$Eta2_generalized >= eta2Thresh,])
# dim(eta2GBM[eta2GBM$Eta2_generalized >= eta2Thresh,])

eta2_PPV <- rbind(eta2ENETw_PPV[eta2ENETw_PPV$Eta2_generalized >= eta2Thresh,], 
              eta2ENETwo_PPV[eta2ENETwo_PPV$Eta2_generalized >= eta2Thresh,],
              eta2GBM_PPV[eta2GBM_PPV$Eta2_generalized >= eta2Thresh,])

# sort variables according to total eta2 across all models
eta2Sums <- aggregate(Eta2_generalized ~ Parameter, data = eta2_PPV, sum)
sortIdx <- order(eta2Sums$Eta2_generalized, decreasing = T)
eta2_PPV$Parameter <- factor(eta2_PPV$Parameter, 
                         levels = eta2Sums$Parameter[sortIdx])

eta2_PPV$model <- factor(eta2_PPV$model, 
                         levels = c("GBM", "ENETw", "ENETwo")) # obvious order of models

# als Balkendiagramme mit Cut-Off bei gen. eta² von .1 plotten
#   y-Achse unterschiedliche Variablen; Farben unterschiedliche Modelle
(pEta2Model_PPV <- ggplot(eta2_PPV, aes(x = Parameter, y = Eta2_generalized,
                                group = model, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    geom_text(aes(label=round(Eta2_generalized, 2), group = model), 
              angle = 90, hjust = 1.5, vjust=0.5, 
              #angle = 0, vjust=1.5, 
              position = position_dodge(width = .9), 
              color="black", size=3.5)+
    ylab("generalisiertes eta^2") +
    scale_fill_manual(values = c("#990000", "#006600", "#009999")) +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_PVI.png"),
#                 plot = pEta2Model_PPV,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

# these two are not replicable right now
# however, facet versions of generalised eta² are better anyway!
# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_linTP.png"),
#                 plot = pEta2Model,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_linFP.png"),
#                 plot = pEta2Model,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")


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

pEta2_ENETw <- plotEta2_4eachModel(eta2ENETw_PPV, eta2Thresh, "#006600") + 
  ggtitle("ENETw") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))
pEta2_ENETwo <- plotEta2_4eachModel(eta2ENETwo_PPV, eta2Thresh, "#009999") +
  ggtitle("ENETwo") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))
pEta2_GBM <- plotEta2_4eachModel(eta2GBM_PPV, eta2Thresh, "#990000") +
  ggtitle("GBM") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

(pEta2_PVI <- pEta2_GBM | pEta2_ENETw | pEta2_ENETwo)

# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_PVI_facetten.png"),
#                 plot = pEta2_PVI,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

pEta2_ENETw <- plotEta2_4eachModel(eta2ENETw_linTP, eta2Thresh, "#006600") + 
  ggtitle("ENETw") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))
pEta2_ENETwo <- plotEta2_4eachModel(eta2ENETwo_linTP, eta2Thresh, "#009999") +
  ggtitle("ENETwo") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))
pEta2_GBM <- plotEta2_4eachModel(eta2GBM_linTP, eta2Thresh, "#990000") +
  ggtitle("GBM") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

(pEta2_linTP <- pEta2_GBM | pEta2_ENETw | pEta2_ENETwo)

# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_linTP_facetten.png"),
#                 plot = pEta2_linTP,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

pEta2_ENETw <- plotEta2_4eachModel(eta2ENETw_linFP_woTInter, eta2Thresh, "#006600") + 
  ggtitle("ENETw") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))
pEta2_ENETwo <- plotEta2_4eachModel(eta2ENETwo_linFP_woTInter, eta2Thresh, "#009999") +
  ggtitle("ENETwo") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))
pEta2_GBM <- plotEta2_4eachModel(eta2GBM_linFP_woTInter, eta2Thresh, "#990000") +
  ggtitle("GBM") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

(pEta2_linFP_woTInter <- pEta2_GBM | pEta2_ENETw | pEta2_ENETwo)

# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_linFP_facetten.png"),
#                 plot = pEta2_PVI,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

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

plotDiags <- aggregate(cbind(linTP, linFP_woTInter, PPV, ACC,
                             specificity, sensitivity, balACC) ~ 
                          model + N + pTrash + rel + R2 + lin_inter, 
                        data = rbind(linENETwo, linENETw, linGBM), 
                        function(x) {cbind(mean(x), 
                                           sd(x),
                                           quantile(x, 0.025),
                                           quantile(x, 0.975))})

get(plotDiags)
for (iCol in paste0("df",1:3)){
  d=get(i)
  colnames(d)[10]=i
  assign(i,d)
}

# names(plotDiags) <- "M_balACC"
# 
# plotBalACC$M <- plotBalACC[,"balACC"][,1]
# plotBalACC$SE <- plotBalACC[,"balACC"][,2]
# plotBalACC$balACC <- NULL

plotDiags$N <- factor(plotDiags$N, levels = c(100, 300, 1000))
plotDiags$model <- factor(plotDiags$model, 
                        levels = c("ENETw", "ENETwo", "GBM"))


ggplot(plotDiags[plotDiags$pTrash == 50 &
                   , ],
       aes(x = rel, y = M, 
           group = interaction(R2, model), colour = R2,
           linetype = model)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  # scale_shape_manual(values = c(16, 8)) +
  geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), 
                width = 0.2, alpha = 0.4) +  
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 1)) +
  facet_grid(N~ lin_inter, labeller = label_both) +
  ylab("PPV") +
  xlab("reliability of predictors") +
  ggtitle("PPV: TP / (TP + FP)") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

################################################################################
# plot sensitivity and specificity
################################################################################
plotSensitivity <- aggregate(sensitivity ~ model + N + pTrash + rel + R2 + lin_inter, 
                        data = rbind(linENETwo, linENETw, linGBM), 
                        function(x) {cbind(mean(x), sd(x))})

plotSensitivity$M <- plotSensitivity[,"sensitivity"][,1]
plotSensitivity$SE <- plotSensitivity[,"sensitivity"][,2]
plotSensitivity$sensitivity <- NULL

str(plotSensitivity)
plotSensitivity$N <- factor(plotSensitivity$N, levels = c(100, 300, 1000))
plotSensitivity$model <- factor(plotSensitivity$model, 
                           levels = c("ENETw", "ENETwo", "GBM"))


ggplot(plotSensitivity[plotSensitivity$pTrash == 10, ],
       aes(x = rel, y = M, 
           group = interaction(R2, model), colour = R2,
           linetype = model)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  # scale_shape_manual(values = c(16, 8)) +
  geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), 
                width = 0.2, alpha = 0.4) +  
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 1)) +
  facet_grid(N~ lin_inter, labeller = label_both) +
  ylab("PPV") +
  xlab("reliability of predictors") +
  ggtitle("PPV: TP / (TP + FP)") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

plotSpecificity <- aggregate(specificity ~ model + N + pTrash + rel + R2 + lin_inter, 
                             data = rbind(linENETwo, linENETw, linGBM), 
                             function(x) {cbind(mean(x), sd(x))})

plotSpecificity$M <- plotSpecificity[,"sensitivity"][,1]
plotSpecificity$SE <- plotSpecificity[,"speci"][,2]
plotSpecificity$sensitivity <- NULL

str(plotSensitivity)
plotSensitivity$N <- factor(plotSensitivity$N, levels = c(100, 300, 1000))
plotSensitivity$model <- factor(plotSensitivity$model, 
                                levels = c("ENETw", "ENETwo", "GBM"))


ggplot(plotSpeci[plotSensitivity$pTrash == 10, ],
       aes(x = rel, y = M, 
           group = interaction(R2, model), colour = R2,
           linetype = model)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  # scale_shape_manual(values = c(16, 8)) +
  geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), 
                width = 0.2, alpha = 0.4) +  
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 1)) +
  facet_grid(N~ lin_inter, labeller = label_both) +
  ylab("PPV") +
  xlab("reliability of predictors") +
  ggtitle("PPV: TP / (TP + FP)") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))


################################################################################
# plot PVI
################################################################################
plotPPV <- aggregate(PPV ~ model + N + pTrash + rel + R2 + lin_inter, 
                     data = rbind(linENETwo, linENETw, linGBM), 
                     function(x) {cbind(mean(x), sd(x))})

plotPPV$M_PPV <- plotPPV[,"PPV"][,1]
plotPPV$MCE_PPV <- plotPPV[,"PPV"][,2]
plotPPV$PPV <- NULL

str(plotPPV)
plotPPV$N <- factor(plotPPV$N, levels = c(100, 300, 1000))
plotPPV$model <- factor(plotPPV$model, 
                        levels = c("ENETw", "ENETwo", "GBM"))

# show ppv separately for ENETs and GBM 
# -> dadurch, dass PPV sowieso durch andere Faktoren beeinflusst wird in GBM und ENETs,
#     und daher die Abbildungen getrennt werden, vielleicht doch alle Interaktionen aus 
#     ENETw rausnehmen?
# ... as can be seen in ANOVA: the relevant simulated conditions for performance 
#     in terms of PPV are really different for models 



##### ENETs #####
# R2, lin_inter, N
# p Trash mactht hier quasi keinen Unterschied
(plotPPV_overview <- ggplot(plotPPV[plotPPV$pTrash == 50 &
                                      plotPPV$model != "GBM", ],
                            aes(x = rel, y = M_PPV, 
                                group = interaction(R2, model), colour = R2,
                                linetype = model)) +
   geom_point() +
   geom_line() +
   scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
   # scale_shape_manual(values = c(16, 8)) +
   geom_errorbar(aes(ymin = M_PPV - 2*MCE_PPV, ymax = M_PPV + 2*MCE_PPV), 
                 width = 0.2, alpha = 0.4) +  
   scale_color_manual(values = colValues) +
   geom_hline(aes(yintercept = 1)) +
   facet_grid(N~ lin_inter, labeller = label_both) +
   ylab("PPV") +
   xlab("reliability of predictors") +
   ggtitle("PPV: TP / (TP + FP)") +
   theme(axis.text.y = element_text(size = 20),
         axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         strip.text.x = element_text(size = 15),
         strip.text.y = element_text(size = 15)))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/PPV_ENETs.png"),
#                 plot = plotPPV_overview,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

##### GBM #####
## only pTrash and N:R2 really relevant (lin_inter... maybe)
(plotPPV_GBM <- ggplot(plotPPV[plotPPV$model == "GBM" &
                                 plotPPV$rel == 1, ],
                            aes(x = N, y = M_PPV, 
                                group = interaction(R2, pTrash), colour = R2,
                                linetype = pTrash)) +
   geom_point() +
   geom_line() +
   scale_linetype_manual(values = c("solid", "dotted")) +
   # scale_shape_manual(values = c(16, 8)) +
   geom_errorbar(aes(ymin = M_PPV - 2*MCE_PPV, ymax = M_PPV + 2*MCE_PPV), 
                 width = 0.2, alpha = 0.4) +  
   scale_color_manual(values = colValues) +
   geom_hline(aes(yintercept = 1)) +
   facet_wrap(~ lin_inter) +
   ylab("") +
   xlab("reliability of predictors") +
   ggtitle("PPV: TP / (TP + FP)") +
   theme(axis.text.y = element_text(size = 20),
         axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         strip.text.x = element_text(size = 15),
         strip.text.y = element_text(size = 15)))

# ggplot2::ggsave(filename = paste0(plotFolder, "/PPV_GBM.png"),
#                 plot = plotPPV_GBM,
#                 width = 13.63,
#                 height = 8.59,
#                 units = "in")

# Fazit: GBMs extrahieren wahrscheinlich viel Trash
# -> wie sieht TP Score aus? 
# -> wie sieht Rangreihe der PVIs aus? profitiert GBM hier von Konfundierung zwischen
#     linearen und Interaktionseffekten im GBM?

################################################################################
# plot lin TP & FP
################################################################################
plotLinTP <- aggregate(linTP ~ model + N + pTrash + rel + R2 + lin_inter, 
                       data = rbind(linENETwo, linENETw, linGBM), 
                       function(x) {cbind(M_TP = mean(x), 
                                          SE_TP = sd(x),
                                          q25 = quantile(x, 0.25),
                                          q75 = quantile(x, 0.75))})

plotLinTP$M_TP <- plotLinTP[,"linTP"][,1]
plotLinTP$SE_TP <- plotLinTP[,"linTP"][,2]
plotLinTP$q25 <- plotLinTP[,"linTP"][,3]
plotLinTP$q75 <- plotLinTP[,"linTP"][,4]
plotLinTP$linTP <- NULL

str(plotLinTP)
plotLinTP$N <- factor(plotLinTP$N, levels = c(100, 300, 1000))
plotLinTP$model <- factor(plotLinTP$model, 
                          levels = c("ENETw", "ENETwo", "GBM"))

plotLinFP <- aggregate(linFP_woTInter ~ model + N + pTrash + rel + R2 + lin_inter, 
                       data = rbind(linENETwo, linENETw, linGBM), 
                       function(x) {cbind(mean(x), 
                                          sd(x),
                                          q25 = quantile(x, 0.25),
                                          q75 = quantile(x, 0.75))})

plotLinFP$M_FP <- plotLinFP[,"linFP_woTInter"][,1]
plotLinFP$SE_FP <- plotLinFP[,"linFP_woTInter"][,2]
plotLinFP$q25 <- plotLinFP[,"linFP_woTInter"][,3]
plotLinFP$q75 <- plotLinFP[,"linFP_woTInter"][,4]
plotLinFP$linFP_woTInter <- NULL

str(plotLinFP)
plotLinFP$N <- factor(plotLinFP$N, levels = c(100, 300, 1000))
plotLinFP$model <- factor(plotLinFP$model, 
                          levels = c("ENETw", "ENETwo", "GBM"))

colValues <- c("green3", "darkblue", "darkmagenta")

##### GBM #####
# je größer N/je größer R2, desto eher werden TP gefunden
# je höher R2, desto weniger entscheidend ist N für TP = bei höherem R2 genügt 
#     kleinere Stichprobe, um alle TP zu finden
(pTP_GBM <- ggplot(plotLinTP[plotLinTP$model == "GBM", ],
        aes(x = N, y = M_TP, 
            group = interaction(R2, pTrash), colour = R2,
            linetype = pTrash)) +
   geom_point() +
   geom_line() +
   scale_linetype_manual(values = c("dotted", "solid")) +
   # scale_shape_manual(values = c(16, 8)) +
   # geom_errorbar(aes(ymin = M_TP - 2*SE_TP, ymax = M_TP + 2*SE_TP), 
   geom_errorbar(aes(ymin = q25, ymax = q75), 
                 width = 0.2, alpha = 0.4) +  
   scale_color_manual(values = colValues) +
   geom_hline(aes(yintercept = 1)) +
   facet_grid(rel ~ lin_inter, labeller = label_both) +
   ylab("number of true positive predictors") +
   xlab("reliability of predictors") +
   ggtitle("TP (GBM)") +
   theme(axis.text.y = element_text(size = 20),
         axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         strip.text.x = element_text(size = 15),
         strip.text.y = element_text(size = 15)))

# basically consistently all 10 trash variables extracted; 
# je größer R2 und je größer N, desto mehr Trash wird mit extrahiert
(pFP_GBM <- ggplot(plotLinFP[plotLinFP$model == "GBM", ],
        aes(x = N, y = M_FP, 
            group = interaction(R2, pTrash), colour = R2,
            linetype = pTrash)) +
    geom_point() +
    geom_line() +
    scale_linetype_manual(values = c("dotted", "solid")) +
    # scale_shape_manual(values = c(16, 8)) +
    #geom_errorbar(aes(ymin = M_FP - 2*SE_FP, ymax = M_FP + 2*SE_FP), 
    geom_errorbar(aes(ymin = q25, ymax = q75), 
                  width = 0.2, alpha = 0.4) +  
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 1)) +
    facet_grid(rel ~ lin_inter, labeller = label_both) +
    ylab("PPV") +
    xlab("reliability of predictors") +
    ggtitle("PPV: TP / (TP + FP)") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

##### ENETs #####
(pTP_ENETs <- ggplot(plotLinTP[plotLinTP$model != "GBM", ],
        aes(x = N, y = M_TP, 
            group = interaction(R2, model, pTrash), colour = R2,
            linetype = model, alpha = pTrash)) +
   geom_point() +
   geom_line() +
   scale_linetype_manual(values = c("solid", "dashed")) +
   scale_alpha_manual(values = c(1, 0.3)) +
   # scale_shape_manual(values = c(16, 8)) +
   #geom_errorbar(aes(ymin = M_TP - 2*SE_TP, ymax = M_TP + 2*SE_TP), 
   geom_errorbar(aes(ymin = q25, ymax = q75), 
                 width = 0.2, alpha = 0.4) +  
   scale_color_manual(values = colValues) +
   geom_hline(aes(yintercept = 4)) +
   facet_grid(rel ~ lin_inter, labeller = label_both) +
   ylab("PPV") +
   xlab("reliability of predictors") +
   ggtitle("PPV: TP / (TP + FP)") +
   theme(axis.text.y = element_text(size = 20),
         axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         strip.text.x = element_text(size = 15),
         strip.text.y = element_text(size = 15)))

(pFP_ENETs <- ggplot(plotLinFP[plotLinFP$model != "GBM", ],
                     aes(x = N, y = M_FP, 
                         group = interaction(R2, model, pTrash), colour = R2,
                         linetype = model, alpha = pTrash)) +
                geom_point() +
                geom_line() +
                scale_linetype_manual(values = c("solid", "dashed")) +
                scale_alpha_manual(values = c(1, 0.3)) +
                # scale_shape_manual(values = c(16, 8)) +
                # geom_errorbar(aes(ymin = M_FP - 2*SE_FP, ymax = M_FP + 2*SE_FP), 
                geom_errorbar(aes(ymin = q25, ymax = q75), 
                              width = 0.2, alpha = 0.4) +  
                scale_color_manual(values = colValues) +
                geom_hline(aes(yintercept = 1)) +
                facet_grid(rel ~ lin_inter, labeller = label_both) +
                ylab("PPV") +
                xlab("reliability of predictors") +
                ggtitle("PPV: TP / (TP + FP)") +
                theme(axis.text.y = element_text(size = 20),
                      axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
                      axis.title.x = element_text(size = 20),
                      axis.title.y = element_text(size = 20),
                      strip.text.x = element_text(size = 15),
                      strip.text.y = element_text(size = 15)))


################################################################################
# calculate rank correlations
################################################################################
# to do: issues with calculating rank correlations
#     - how do we get variance in the rank of variables of the dgm?
#       (a) random ranks for variables with effects (near zero correlations)
#       (b) cheat it: result ranks for linear effect variables in dgm
#       issue: multiple true linear effects might have the same rank therefore reducing variance again
#       issue: multiple samples; ranks are therefore used nTrain times 
#     - what if only a subset of variables is in the model? 
#       issue: if only one of 4 TP variables are extracted but nothing else the correlation will be 1
#             -> more conservative models benefit in terms of rank correlation

# alternative Maße: 
#     - Wie oft ist Variable mit simuliertem linearen Effekt in Top 4?
#     - Wie viele der Prädiktoren sind in Top 4?
#     - Sind die ersten 4 extrahierten Variablen die wahren Prädiktoren?

# to get rank values by group
library(tidyverse)

# pvi values to numeric
pviGBM$pviValue <- as.numeric(pviGBM$pviValue) 
pviENETw$pviValue <- as.numeric(pviENETw$pviValue) 
pviENETwo$pviValue <- as.numeric(pviENETwo$pviValue)

# ranks to variables
# ENETwo and GBM do only have linear predictors 
# ENETw has interactions as well; interactions on top ranks are correctly extracted, too
#   -> interactions that are in setParam$dgp$interEffects need to be ignored! 
pviENETwo <- pviENETwo %>% 
  group_by(idxCondLabel, sample, model, N_pTrash) %>% 
  # negative pvi values as ranks are given in ascending instead of descending order
  # 1 = 1, 1.001 = 2, ... instead of 1.9 = 1, 1.8 = 2, ...
  mutate(rankModel = rank(pviValue * (-1), ties.method = 'average'))

pviGBM <- pviGBM %>% 
  group_by(idxCondLabel, sample, model, N_pTrash) %>% 
  mutate(rankModel = rank(pviValue * (-1), ties.method = 'average'))

pviENETw <- pviENETw %>% 
  group_by(idxCondLabel, sample, model, N_pTrash) %>% 
  # remove true interactions from pvi list
  # thereby avoid ENETw's disadvantage for selecting simulated interactions 
  filter(!(pviRank %in% interEffectsDot)) %>% 
  mutate(rankModel = rank(pviValue * (-1), ties.method = 'average'))

# remove every variable that does not have simulated effects in dgm
pviENETwo_rank <- pviENETwo[pviENETwo$pviRank %in% setParam$dgp$linEffects, ]
pviENETw_rank <- pviENETw[pviENETw$pviRank %in% setParam$dgp$linEffects, ]
pviGBM_rank <- pviGBM[pviGBM$pviRank %in% setParam$dgp$linEffects, ]

rankDGM <- expand.grid(idxCondLabel = unique(pviENETwo$idxCondLabel), 
                       sample = 1:setParam$dgp$nTrain,
                       N_pTrash = unique(pviENETwo$N_pTrash),
                       pviRank = setParam$dgp$linEffects)

rankDGM <- rankDGM %>%  
  group_by(idxCondLabel, sample, N_pTrash) %>% 
  mutate(rankDGM = sample(1:4, 4, replace = FALSE))

pviENETwo_rank <- merge(rankDGM, pviENETwo_rank, all = T, 
                        by = c("idxCondLabel", "sample", "N_pTrash", "pviRank"))
pviENETw_rank <- merge(rankDGM, pviENETw, all = T, 
                        by = c("idxCondLabel", "sample", "N_pTrash", "pviRank"))
pviGBM_rank <- merge(rankDGM, pviGBM_rank, all = T, 
                        by = c("idxCondLabel", "sample", "N_pTrash", "pviRank"))
#####
# # correlation with random rankDGM
# # > range(pviGBM_rankCor$rankCor, na.rm = T)
# # [1] -0.0468  0.0620
# # > range(pviENETwo_rankCor$rankCor, na.rm = T)
# # [1] -0.0526  0.0490
# # > range(pviENETw_rankCor$rankCor, na.rm = T)
# # [1] -0.0530018  0.0438000
# pviENETwo_rankCor <- pviENETwo_rank %>% 
#   group_by(idxCondLabel, model, N_pTrash) %>% 
#   # negative pvi values as ranks are given in ascending instead of descending order
#   # 1 = 1, 1.001 = 2, ... instead of 1.9 = 1, 1.8 = 2, ...
#   mutate(rankCor = cor(x = rankDGM, y = rankModel, method = "spearman", 
#                        use = "pairwise.complete.obs"))
# 
# pviENETw_rankCor <- pviENETw_rank %>% 
#   group_by(idxCondLabel, model, N_pTrash) %>% 
#   # negative pvi values as ranks are given in ascending instead of descending order
#   # 1 = 1, 1.001 = 2, ... instead of 1.9 = 1, 1.8 = 2, ...
#   mutate(rankCor = cor(x = rankDGM, y = rankModel, method = "spearman", 
#                        use = "pairwise.complete.obs"))
# 
# pviGBM_rankCor <- pviGBM_rank %>% 
#   group_by(idxCondLabel, model, N_pTrash) %>% 
#   # negative pvi values as ranks are given in ascending instead of descending order
#   # 1 = 1, 1.001 = 2, ... instead of 1.9 = 1, 1.8 = 2, ...
#   mutate(rankCor = cor(x = rankDGM, y = rankModel, method = "spearman", 
#                        use = "pairwise.complete.obs"))
# 
# range(pviENETwo_rankCor$rankCor, na.rm = T)
# range(pviENETw_rankCor$rankCor, na.rm = T)
# range(pviGBM_rankCor$rankCor, na.rm = T)
#####

# push correlations by using the rank of the model als rankDGM if rankModel in 1:4
pviENETwo_rankCor$rankDGMpush <- ifelse(pviENETwo_rankCor$rankModel %in% 1:4, 
                                        pviENETwo_rankCor$rankModel, pviENETwo_rankCor$rankDGM)
pviENETw_rankCor$rankDGMpush <- ifelse(pviENETw_rankCor$rankModel %in% 1:4, 
                                        pviENETw_rankCor$rankModel, pviENETw_rankCor$rankDGM)
pviGBM_rankCor$rankDGMpush <- ifelse(pviGBM_rankCor$rankModel %in% 1:4, 
                                     pviGBM_rankCor$rankModel, pviGBM_rankCor$rankDGM)

pviENETwo_rankCor <- pviENETwo_rankCor %>% 
  group_by(idxCondLabel, model, N_pTrash) %>% 
  # negative pvi values as ranks are given in ascending instead of descending order
  # 1 = 1, 1.001 = 2, ... instead of 1.9 = 1, 1.8 = 2, ...
  mutate(rankCorPush = cor(x = rankDGMpush, y = rankModel, method = "spearman", 
                       use = "pairwise.complete.obs"))

pviENETw_rankCor <- pviENETw_rankCor %>% 
  group_by(idxCondLabel, model, N_pTrash) %>% 
  # negative pvi values as ranks are given in ascending instead of descending order
  # 1 = 1, 1.001 = 2, ... instead of 1.9 = 1, 1.8 = 2, ...
  mutate(rankCorPush = cor(x = rankDGMpush, y = rankModel, method = "spearman", 
                           use = "pairwise.complete.obs"))

pviGBM_rankCor <- pviGBM_rankCor %>% 
  group_by(idxCondLabel, model, N_pTrash) %>% 
  # negative pvi values as ranks are given in ascending instead of descending order
  # 1 = 1, 1.001 = 2, ... instead of 1.9 = 1, 1.8 = 2, ...
  mutate(rankCorPush = cor(x = rankDGMpush, y = rankModel, method = "spearman", 
                           use = "pairwise.complete.obs"))

# range: 0.4432967 1.0000000
range(pviENETwo_rankCor$rankCorPush, na.rm = T)
range(pviENETw_rankCor$rankCorPush, na.rm = T)
range(pviGBM_rankCor$rankCorPush, na.rm = T)

str(pviENETwo_rankCor)
str(pviENETw_rankCor)
str(pviGBM_rankCor)
pviENETw_rankCor$sample <- as.integer(pviENETw_rankCor$sample)

pviRankCor <- rbind(pviENETwo_rankCor, pviENETw_rankCor, pviGBM_rankCor)
plotRankCor <- aggregate(rankCorPush ~ idxCondLabel + N_pTrash + model, 
          data = pviRankCor, mean) 

plotRankCor_ENETwo <- aggregate(rankCorPush ~ idxCondLabel + N_pTrash + model, 
                         data = pviENETwo_rankCor, mean) 
plotRankCor_ENETw <- aggregate(rankCorPush ~ idxCondLabel + N_pTrash + model, 
                                data = pviENETw_rankCor, mean) 
plotRankCor_GBM <- aggregate(rankCorPush ~ idxCondLabel + N_pTrash + model, 
                                data = pviGBM_rankCor, mean) 

plotRankCor <- idx2infoNew(plotRankCor)
plotRankCor_ENETwo <- idx2infoNew(plotRankCor_ENETwo)
plotRankCor_ENETw <- idx2infoNew(plotRankCor_ENETw)
plotRankCor_GBM <- idx2infoNew(plotRankCor_GBM)

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model")
plotRankCor[col2fac] <- lapply(plotRankCor[col2fac], factor)
plotRankCor$N <- factor(plotRankCor$N, 
                        levels = c(100, 300, 1000))

plotRankCor_ENETwo[col2fac] <- lapply(plotRankCor_ENETwo[col2fac], factor)
plotRankCor_ENETwo$N <- factor(plotRankCor_ENETwo$N, 
                        levels = c(100, 300, 1000))
plotRankCor_ENETw[col2fac] <- lapply(plotRankCor_ENETw[col2fac], factor)
plotRankCor_ENETw$N <- factor(plotRankCor_ENETw$N, 
                               levels = c(100, 300, 1000))
plotRankCor_GBM[col2fac] <- lapply(plotRankCor_GBM[col2fac], factor)
plotRankCor_GBM$N <- factor(plotRankCor_GBM$N, 
                               levels = c(100, 300, 1000))

pRankCor_ENETwo <- ggplot(plotRankCor_ENETwo,
       aes(x = rel, y = rankCorPush, 
           group = interaction(R2, pTrash), colour = R2,
           linetype = pTrash)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  #scale_shape_manual(values = c(16, 8)) +
  # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(N ~ lin_inter, labeller = label_both) +
  ylab("rank correlation (Spearman)") +
  xlab("reliability of predictors") +
  ggtitle("rank correlation: ENETwo") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

pRankCor_ENETw <- ggplot(plotRankCor_ENETw,
                          aes(x = rel, y = rankCorPush, 
                              group = interaction(R2, pTrash), colour = R2,
                              linetype = pTrash)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  #scale_shape_manual(values = c(16, 8)) +
  # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(N ~ lin_inter, labeller = label_both) +
  ylab("rank correlation (Spearman)") +
  xlab("reliability of predictors") +
  ggtitle("rank correlation: ENETwo") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

pRankCor_GBM <- ggplot(plotRankCor_GBM,
                          aes(x = rel, y = rankCorPush, 
                              group = interaction(R2, pTrash), colour = R2,
                              linetype = pTrash)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  #scale_shape_manual(values = c(16, 8)) +
  # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(N ~ lin_inter, labeller = label_both) +
  ylab("rank correlation (Spearman)") +
  xlab("reliability of predictors") +
  ggtitle("rank correlation: ENETwo") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

pRankCor_ENETwo
pRankCor_ENETw
pRankCor_GBM

# ggplot2::ggsave(filename = paste0(plotFolder, "/rankCorrelation_ENETwo.png"),
#                 plot = pRankCor_ENETwo,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/rankCorrelation_ENETw.png"),
#                 plot = pRankCor_ENETw,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/rankCorrelation_GBM.png"),
#                 plot = pRankCor_GBM,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

pRankCor <- ggplot(plotRankCor,
       aes(x = rel, y = rankCorPush, 
           group = interaction(R2, model, pTrash), colour = R2,
           linetype = model, alpha = pTrash, shape = pTrash)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_shape_manual(values = c(16, 8)) +
  scale_alpha_manual(values = c(1, 0.2)) +
  # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
  scale_color_manual(values = colValues) +
  geom_hline(aes(yintercept = 0)) +
  facet_grid(N ~ lin_inter, labeller = label_both) +
  ylab("rank correlation (Spearman)") +
  xlab("reliability of predictors") +
  ggtitle("rank correlation: ENETwo") +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

# ggplot2::ggsave(filename = paste0(plotFolder, "/rankCorrelation.png"),
#                 plot = pRankCor,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")


################################################################################


# how many (frequently do) variables with simulated interaction effects (setParam$dgp$interEffects) 
#     affect (pvi > 1 for ENETw & interStrength for GBM) the prediction?
#     -> = true positive (interaction)
pviENETw$interTP <- ifelse(pviENETw$pviRank %in% interEffectsDot, 1, 0)

# how many (frequently do) variables with simulated interaction effects (setParam$dgp$interEffects) 
#     affect (pvi > 1 for ENETw & interStrength for GBM) the prediction?
#     -> = true positive (interaction)
# calculate positive predictive value: TP / (TP + FP)


# are the linear effects in the top 4 of simulated effects (for ENETwo and GBM) or 
#     in the top 8 of simulated effects (for ENETw)?
# -> for GBM interaction and linear effects are conflated in variable selection and
#     in pvi. How does this effect the identification of linear (and interaction) effects?
# -> rank correlation with Kendalls Tau? 
# are only true (linear) predictors in model 
#     for ENETw: 8 predictors with simulated effects and every other variable == 0
#     for GBM & ENETwo: 4 predictors
# how many of the linear effects are recovered?
# how many of the interaction effects are recovered? 
#     ENETw: pvi
#     GBM: interStrength
#     ENETwo: no interactions extracted
# all simulated effects selected in model?
# only simulated effects selected (i.e., every other predictor is not selected!)
