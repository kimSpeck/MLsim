# analyse variable selection for interactions
# compare results from coefficients to PVI results for ENET with interactions (ENETw)

# confusion matrix: extracted in ENET {yes, no} vs. PVI > 1 {yes, no}

library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²
library(patchwork)

# plot results
library(ggplot2)

source("setParameters.R")
source("analysisTools.R")

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

# load results files
resFolder <- "results/finalResults/dependentMeasures" 

listDir <- dir(resFolder)

################################################################################
# ... for beta coefficients
################################################################################
estBetaSampleENETw <- loadRData(paste0(resFolder, "/estBetaSample_ENETw.rda"))

# write list of matrices in one single data.frame
estBetaSampleENETw <- rbindSingleResults(estBetaSampleENETw) 

# write variables in columns instead of rownames
estBetaSampleENETw <- rowNames2col(estBetaSampleENETw)

colnames(estBetaSampleENETw)[colnames(estBetaSampleENETw) == 'rowNames'] <- 'var'

# # check if transformations were successfull
# estBetaSampleENETw[1:6, 998:1004]
# unique(estBetaSampleENETw$var)

estBetaSampleENETw <- tidyr::pivot_longer(estBetaSampleENETw, 
                                          cols = !c(idxCondLabel, model, N_pTrash, var), 
                                          names_to = "sample",
                                          values_to = "betaEst")

# # remove all rows with variables that were not selected by ENET
# #     if variable is not selected by ENET, the value in betaEst is NA
# #     therefore remove the NAs from the ENET data
# rmNAidx <- is.na(estBetaSampleENETw$betaEst)
# sum(rmNAidx) # 121258924
# dim(estBetaSampleENETw)[1] # 128790000
# dim(estBetaSampleENETw)[1] - sum(rmNAidx) # 7531076

estBetaSampleENETw <- na.omit(estBetaSampleENETw)
gc()

# 1. calculate number of true positives for each sample (interactions)
## extract true interaction effects
# # check interaction labels and number of unique labels
# length(unique(estBetaSampleENETw$var))
interLabels <- stringr::str_replace_all(setParam$dgp$interEffects, "\\:", "\\.") 
estBetaSampleENETw$interTP <- ifelse(estBetaSampleENETw$var %in% interLabels, 1, 0)

# 2. calculate number of false positives for each sample ()
## extract false psotive interaction effects
## it does not matter if main trash predictors are subtracted or not 
##    -> number either 10 or 50 but with more main trash predictors even more interactions 
##        do arrive and therefore this does not make any difference for the results
##    -> maybe easier to report would be false positive out of every interaction?
##        instead of false positive out of every interaction and main effects without simulated effects
allLinVars <- paste0("Var", seq_len(max(setParam$dgp$pTrash)))
# every beta coefficient != 0 that is ...
#   ... an interaction with simulated effects or is a main effect receives 0
#   ... an interaction without simulated effects and is not a main effect receives 1
estBetaSampleENETw$interFP.a <- ifelse(estBetaSampleENETw$var %in% interLabels |
                                       estBetaSampleENETw$var %in% allLinVars, 0, 1)

# # every beta coefficient != 0 that is ...
# #   ... an interaction with simulated effects or is a main effect with simulated effect receives 0
# #   ... an interaction without simulated effects and is not a main effect with simualted effect receives 1
# estBetaSampleENETw$interFP.b <- ifelse(estBetaSampleENETw$var %in% interLabels |
#                                        estBetaSampleENETw$var %in% setParam$dgp$linEffects, 0, 1)

# sum up number of true positives and false negatives within each sample
estBetaENETw <- aggregate(cbind(interTP, interFP.a, interFP.b) ~ sample + idxCondLabel + model + N_pTrash, 
                          data = estBetaSampleENETw, sum)
rm(estBetaSampleENETw)
gc()

# 3. calculate positive predictive value (TP / (TP + FP)) for every sample
estBetaENETw$interPPV.a <- estBetaENETw$interTP / (estBetaENETw$interTP + estBetaENETw$interFP.a) 
# estBetaENETw$interPPV.b <- estBetaENETw$interTP / (estBetaENETw$interTP + estBetaENETw$interFP.b) 

estBetaENETw <- idx2infoNew(estBetaENETw)

# two variants ... 
#   ... .a only out of interactions vs. 
#   ... .b out of interactions and main predictors without effect)

# 2a. add measures... 
#   ... interTN =  interactions without simulated effect and which are not extracted 

# all possible interactions - true interactions (with simulated effects) - false positive interactions
# .a as this only counts interactions as true negative options
estBetaENETw$interTN.a <- choose(as.numeric(estBetaENETw$pTrash) + length(setParam$dgp$linEffects), 
                               setParam$dgp$interDepth) - 
  length(setParam$dgp$interEffects) - estBetaENETw$interFP.a

# # all possible interactions + linear predictors without effect - true interactions (with simulated effects) - false positive interactions
# # .b as this counts interactions as true negative options and pTrash variables without simulated effects
# # (linear predictors with simulated effects are ignored alltogether here)
# estBetaENETw$interTN.b <- choose(as.numeric(estBetaENETw$pTrash) + length(setParam$dgp$linEffects), 
#                                  setParam$dgp$interDepth) + 
#   as.numeric(estBetaENETw$pTrash) - 
#   length(setParam$dgp$interEffects) - estBetaENETw$interFP.b 

#   ... interFN = simulated interaction effect that is not extracted
estBetaENETw$interFN <- length(setParam$dgp$interEffects) - estBetaENETw$interTP

# 2b. add measures...
#   ... accuracy = (TP + TN) / (TP + TN + FP + FN)
estBetaENETw$interACC.a <- (estBetaENETw$interTP + estBetaENETw$interTN.a) / 
  (estBetaENETw$interTP + estBetaENETw$interTN.a + estBetaENETw$interFP.a + estBetaENETw$interFN)

# estBetaENETw$interACC.b <- (estBetaENETw$interTP + estBetaENETw$interTN.b) / 
#   (estBetaENETw$interTP + estBetaENETw$interTN.b + estBetaENETw$interFP.b + estBetaENETw$interFN)

#   ... sensitivity = TP / (TP + FN)
estBetaENETw$interSensitivity <- estBetaENETw$interTP / (estBetaENETw$interTP + estBetaENETw$interFN)

#   ... specificity = TN / (TN + FP)
estBetaENETw$interSpecificity.a <- estBetaENETw$interTN.a / (estBetaENETw$interTN.a + estBetaENETw$interFP.a)
# estBetaENETw$interSpecificity.b <- estBetaENETw$interTN.b / (estBetaENETw$interTN.b + estBetaENETw$interFP.b)

#   ... balanced accuracy = (sensitivity + specificity) / 2 
estBetaENETw$interBalACC.a <- (estBetaENETw$interSpecificity.a + estBetaENETw$interSensitivity) / 2
# estBetaENETw$interBalACC.b <- (estBetaENETw$interSpecificity.b + estBetaENETw$interSensitivity) / 2

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter")
estBetaENETw[col2fac] <- lapply(estBetaENETw[col2fac], factor)
str(estBetaENETw)
# how often are predictors "selected" in model? beta-coefficient > 0 (i.e., NA)
# ! pvi for ENETw includes interactions
#   otherwise {ENETwo & GBM} pvi does not inlcude interactions

################################################################################
# run ANOVAs for ...
#     ... different dependent variables {interTP, interFP.a, ...}
################################################################################
# there are only between-sample factors for the ANOVA
# generate ID along data set
estBetaENETw$ID <- seq_len(dim(estBetaENETw)[1])

dvVec <- c("interTP", "interFP.a", "interFP.b", "interPPV.a", "interPPV.b", 
           "interACC.a", "interACC.b", "interSensitivity", 
           "interSpecificity.a", "interSpecificity.b", "interBalACC.a", "interBalACC.b")
for (iDV in dvVec) {
  # only between-sample ANOVA
  anovaObjectName <- paste0("anova", iDV)
  tmp <- aov_ez(id = "ID",
                dv = iDV,
                data = estBetaENETw,
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
  tmpEta2 <- cbind(tmpEta2, dv = iDV)
    
  # which higher order interactions do we need to illustrate to report simulation results?
  etaObjectName <- paste0("eta2", iDV)
  assign(etaObjectName, tmpEta2) 
}

eta2Thresh <- 0.01

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

# there seems to be no difference between ...
# ... specificity.a and specificity.b 
# ... balACC.a and balACC.b 
# ... ACC.a and ACC.b 
# therefore remove all b-versions 
pSensitivity <- plotEta2_4eachModel(eta2interSensitivity, eta2Thresh, "#006600") + 
  ggtitle("Sensitivity") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pSpecificity.a <- plotEta2_4eachModel(eta2interSpecificity.a, eta2Thresh, "#006600") + 
  ggtitle("Specificity") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

# pSpecificity.b <- plotEta2_4eachModel(eta2interSpecificity.b, eta2Thresh, "#006600") + 
#   ggtitle("Specificity.b") + # for the main title
#   theme(plot.title = element_text(hjust = 0.5, size = 25))

pBalACC.a <- plotEta2_4eachModel(eta2interBalACC.a, eta2Thresh, "#006600") + 
  ggtitle("BalACC") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

# pBalACC.b <- plotEta2_4eachModel(eta2interBalACC.b, eta2Thresh, "#006600") + 
#   ggtitle("BalACC.b") + # for the main title
#   theme(plot.title = element_text(hjust = 0.5, size = 25))

pACC.a <- plotEta2_4eachModel(eta2interACC.a, eta2Thresh, "#006600") + 
  ggtitle("ACC") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

# pACC.b <- plotEta2_4eachModel(eta2interACC.b, eta2Thresh, "#006600") + 
#   ggtitle("ACC.b") + # for the main title
#   theme(plot.title = element_text(hjust = 0.5, size = 25))

# pSensitivity + pSpecificity.a + pSpecificity.b + pBalACC.a + pBalACC.b + pACC.a + pACC.b

pTP <- plotEta2_4eachModel(eta2interTP, eta2Thresh, "#006600") + 
  ggtitle("TP") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pFP.a <- plotEta2_4eachModel(eta2interFP.a, eta2Thresh, "#006600") + 
  ggtitle("FP") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

# pFP.b <- plotEta2_4eachModel(eta2interFP.b, eta2Thresh, "#006600") + 
#   ggtitle("FP.b") + # for the main title
#   theme(plot.title = element_text(hjust = 0.5, size = 25))

pPPV.a <- plotEta2_4eachModel(eta2interPPV.a, eta2Thresh, "#006600") + 
  ggtitle("PPV") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

# pPPV.b <- plotEta2_4eachModel(eta2interPPV.b, eta2Thresh, "#006600") + 
#   ggtitle("PPV.b") + # for the main title
#   theme(plot.title = element_text(hjust = 0.5, size = 25))

# pTP + pFP.a + pFP.b + pPPV.a + pPPV.b

pbetaInter <- pSensitivity + pSpecificity.a + pBalACC.a + 
  pTP + pFP.a + pPPV.a + plot_layout(ncol = 3)

# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_betaInterMeasures.png"),
#                 plot = pbetaInter,
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

colnames(estBetaENETw) <- stringr::str_replace_all(colnames(estBetaENETw), "\\.a", "")

plotInterBeta <- aggregate(cbind(interTP, interFP, interPPV, interACC,
                                 interSpecificity, interSensitivity, interBalACC) ~ 
                             model + N + pTrash + rel + R2 + lin_inter, 
                           data = estBetaENETw, 
                           function(x) {cbind(mean(x), 
                                              sd(x),
                                              quantile(x, 0.025),
                                              quantile(x, 0.975))})

plotInterBeta <- do.call(data.frame, plotInterBeta)
colnames(plotInterBeta) <- stringr::str_replace_all(colnames(plotInterBeta), "\\.1", "_M")
colnames(plotInterBeta) <- stringr::str_replace_all(colnames(plotInterBeta), "\\.2", "_SE")
colnames(plotInterBeta) <- stringr::str_replace_all(colnames(plotInterBeta), "\\.3", "_q025")
colnames(plotInterBeta) <- stringr::str_replace_all(colnames(plotInterBeta), "\\.4", "_q975")

plotInterBeta$N <- factor(plotInterBeta$N, levels = c(100, 300, 1000))

plotInterBeta <- tidyr::pivot_longer(plotInterBeta, 
                                     cols = !c(model, N, pTrash, rel, R2, lin_inter),
                                     names_to = c("DV", "measure"), 
                                     names_sep = "_",
                                     values_to = "values")

plotInterBeta <- tidyr::pivot_wider(plotInterBeta,
                                    names_from = measure, 
                                    values_from = values)

colValues <- c("green3", "darkblue", "darkmagenta")


(pInterSensitivity <- ggplot(plotInterBeta[plotInterBeta$DV == "interSensitivity",],
                             aes(x = N, y = M, 
                                 group = interaction(R2, rel), colour = R2,
                                 linetype = rel, shape = rel)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) +
    scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(values = c(16, 1, 8)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), 
                  width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) +  
    scale_color_manual(values = colValues) +
    facet_grid(pTrash ~ lin_inter, labeller = label_both) +
    ylab("sensitivity") +
    xlab("N") +
    ggtitle("sensitivity: TP / (TP + FN)") +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/interSensitivity_betaCoef.png"),
#                 plot = pInterSensitivity,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")


(pInterSpecificity <- ggplot(plotInterBeta[plotInterBeta$DV == "interSpecificity",],
                             aes(x = N, y = M, 
                                 group = interaction(R2, rel), colour = R2,
                                 linetype = rel, shape = rel)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) +
    scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(values = c(16, 1, 8)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), 
                  width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) +  
    scale_color_manual(values = colValues) +
    facet_grid(pTrash ~ lin_inter, labeller = label_both) +
    ylab("specificity") +
    xlab("N") +
    ggtitle("specificity: TN / (TN + FP)") +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/interSpecificity_betaCoef.png"),
#                 plot = pInterSpecificity,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

(pInterBalACC <- ggplot(plotInterBeta[plotInterBeta$DV == "interBalACC",],
                             aes(x = N, y = M, 
                                 group = interaction(R2, rel), colour = R2,
                                 linetype = rel, shape = rel)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) +
    scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(values = c(16, 1, 8)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), 
                  width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) +  
    scale_color_manual(values = colValues) +
    facet_grid(pTrash ~ lin_inter, labeller = label_both) +
    ylab("balanced accuracy") +
    xlab("N") +
    ggtitle("balanced accuracy: (sensitivity + specificity) / 2") +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/interBalAcc_betaCoef.png"),
#                 plot = pInterBalACC,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
