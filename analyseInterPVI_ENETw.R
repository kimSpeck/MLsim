# load parameters and helper functions 
source("setParameters.R")
source("analysisTools.R")

# ANOVA
library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²

# plot results
library(ggplot2)
library(patchwork)
library(see)

colValues <- c("green3", "darkblue", "darkmagenta")

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

# load results files
resFolder <- "results/finalResults/dependentMeasures" 

# interactions written as Var1.Var2 instead of Var1:Var2
interEffectsDot <- gsub("\\:", "\\.", setParam$dgp$interEffects)

listDir <- dir(resFolder)
pviENETw <- loadRData(paste0(resFolder, "/pvi_ENETw.rda"))

# pull data from nested list of all results (fullData)
pviENETw <- rbindSingleResults(pviENETw)

################################################################################
# evaluate frequency stats
################################################################################
# how often are predictors "selected" in model? 
#   -> i.e., how often is a certain variable relevant for prediction (pvi > 1)
#       all variables with pvi <= 1 are already removed
# ! pvi for ENETw includes interactions

# 1. calculate number of true positives for each sample (interactions)
## extract true interaction effects
interLabels <- stringr::str_replace_all(setParam$dgp$interEffects, "\\:", "\\.") 
pviENETw$interTP <- ifelse(pviENETw$pviRank %in% interLabels, 1, 0)

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
pviENETw$interFP.a <- ifelse(pviENETw$pviRank %in% interLabels |
                               pviENETw$pviRank %in% allLinVars, 0, 1)

# every beta coefficient != 0 that is ...
#   ... an interaction with simulated effects or is a main effect with simulated effect receives 0
#   ... an interaction without simulated effects and is not a main effect with simualted effect receives 1
pviENETw$interFP.b <- ifelse(pviENETw$pviRank %in% interLabels |
                               pviENETw$pviRank %in% setParam$dgp$linEffects, 0, 1)

# sum up number of true positives and false negatives within each sample
interENETw <- aggregate(cbind(interTP, interFP.a, interFP.b) ~ sample + idxCondLabel + model + N_pTrash, 
                          data = pviENETw, sum)

# 3. calculate positive predictive value (TP / (TP + FP)) for every sample
interENETw$interPPV.a <- interENETw$interTP / (interENETw$interTP + interENETw$interFP.a) 
interENETw$interPPV.b <- interENETw$interTP / (interENETw$interTP + interENETw$interFP.b)

interENETw <- idx2infoNew(interENETw)

# two variants ... 
#   ... .a only out of interactions vs. 
#   ... .b out of interactions and main predictors without effect)

# 2a. add measures... 
#   ... interTN =  interactions without simulated effect and which are not extracted 

# all possible interactions - true interactions (with simulated effects) - false positive interactions
# .a as this only counts interactions as true negative options
interENETw$interTN.a <- choose(as.numeric(interENETw$pTrash) + length(setParam$dgp$linEffects), 
                                 setParam$dgp$interDepth) - 
  length(setParam$dgp$interEffects) - interENETw$interFP.a

# all possible interactions + linear predictors without effect - true interactions (with simulated effects) - false positive interactions
# .b as this counts interactions as true negative options and pTrash variables without simulated effects
# (linear predictors with simulated effects are ignored alltogether here)
interENETw$interTN.b <- choose(as.numeric(interENETw$pTrash) + length(setParam$dgp$linEffects),
                                 setParam$dgp$interDepth) +
  as.numeric(interENETw$pTrash) -
  length(setParam$dgp$interEffects) - interENETw$interFP.b

#   ... interFN = simulated interaction effect that is not extracted
interENETw$interFN <- length(setParam$dgp$interEffects) - interENETw$interTP

# 2b. add measures...
#   ... accuracy = (TP + TN) / (TP + TN + FP + FN)
interENETw$interACC.a <- (interENETw$interTP + interENETw$interTN.a) / 
  (interENETw$interTP + interENETw$interTN.a + interENETw$interFP.a + interENETw$interFN)

interENETw$interACC.b <- (interENETw$interTP + interENETw$interTN.b) /
  (interENETw$interTP + interENETw$interTN.b + interENETw$interFP.b + interENETw$interFN)

#   ... sensitivity = TP / (TP + FN)
interENETw$interSensitivity <- interENETw$interTP / (interENETw$interTP + interENETw$interFN)

#   ... specificity = TN / (TN + FP)
interENETw$interSpecificity.a <- interENETw$interTN.a / (interENETw$interTN.a + interENETw$interFP.a)
interENETw$interSpecificity.b <- interENETw$interTN.b / (interENETw$interTN.b + interENETw$interFP.b)

#   ... balanced accuracy = (sensitivity + specificity) / 2 
interENETw$interBalACC.a <- (interENETw$interSpecificity.a + interENETw$interSensitivity) / 2
interENETw$interBalACC.b <- (interENETw$interSpecificity.b + interENETw$interSensitivity) / 2

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter")
interENETw[col2fac] <- lapply(interENETw[col2fac], factor)
str(interENETw)

################################################################################
# run ANOVAs for ...
#     ... different dependent variables {interTP, interFP.a, ...}
################################################################################
# there are only between-sample factors for the ANOVA
# generate ID along data set
interENETw$ID <- seq_len(dim(interENETw)[1])

dvVec <- c("interTP", "interFP.a", "interFP.b", "interPPV.a", "interPPV.b", 
           "interACC.a", "interACC.b", "interSensitivity", 
           "interSpecificity.a", "interSpecificity.b", "interBalACC.a", "interBalACC.b")
for (iDV in dvVec) {
  # only between-sample ANOVA
  anovaObjectName <- paste0("anova", iDV)
  tmp <- aov_ez(id = "ID",
                dv = iDV,
                data = interENETw,
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

pPviInter <- pSensitivity + pSpecificity.a + pBalACC.a + 
  pTP + pFP.a + pPPV.a + plot_layout(ncol = 3)

# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_pviInterMeasures.png"),
#                 plot = pPviInter,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

# results for pvi and beta are exactly the same

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

colnames(interENETw) <- stringr::str_replace_all(colnames(interENETw), "\\.a", "")

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

colValues <- c("green3", "darkblue", "darkmagenta")


(pInterSensitivity <- ggplot(plotInterPVI[plotInterPVI$DV == "interSensitivity",],
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/interSensitivity_pvi.png"),
#                 plot = pInterSensitivity,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")


(pInterSpecificity <- ggplot(plotInterPVI[plotInterPVI$DV == "interSpecificity",],
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/interSpecificity_pvi.png"),
#                 plot = pInterSpecificity,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

(pInterBalACC <- ggplot(plotInterPVI[plotInterPVI$DV == "interBalACC",],
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/interBalAcc_pvi.png"),
#                 plot = pInterBalACC,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
