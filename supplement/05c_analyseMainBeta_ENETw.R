# analyse variable selection for interactions
# compare results from coefficients to PVI results for ENET with interactions (ENETw)
library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²
library(patchwork)

# plot results
library(ggplot2)
library(ggh4x)

# plot utils
colValuesR2 <- c('#db4a07', '#850c0c', '#3c1518')
colValuesInter <- c('#050440', '#181ff2', '#0eb2e8')
colValuesLin <- c('#0eb2e8', '#181ff2', '#050440')


source("utils/setParameters.R")
source("utils/analysisTools.R")

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

# 1. calculate number of true positives for each sample (main effects)
## extract true main effects
# length(unique(estBetaSampleENETw$var))
estBetaSampleENETw$mainTP <- ifelse(estBetaSampleENETw$var %in% setParam$dgp$linEffects, 1, 0)

# 2. calculate number of false positives for each sample ()
## extract false positive linear/main effects
allInterVars <- stringr::str_subset(estBetaSampleENETw$var, "\\.")
interLabels <- stringr::str_replace_all(setParam$dgp$interEffects, "\\:", "\\.") 

# every beta coefficient != 0 that is ...
#   ... a main effect with simulated effects or is an interaction of any kind receives 0
#   ... a main effect without simulated effects and is not an interaction of any kind receives 1
estBetaSampleENETw$mainFP <- ifelse(estBetaSampleENETw$var %in% setParam$dgp$linEffects |
                                         estBetaSampleENETw$var %in% allInterVars, 0, 1)

# sum up number of true positives and false negatives within each sample
estBetaENETw <- aggregate(cbind(mainTP, mainFP) ~ sample + idxCondLabel + model + N_pTrash, 
                          data = estBetaSampleENETw, sum)
rm(estBetaSampleENETw)
gc()

# 3. calculate positive predictive value (TP / (TP + FP)) for every sample
estBetaENETw$mainPPV <- estBetaENETw$mainTP / (estBetaENETw$mainTP + estBetaENETw$mainFP) 

estBetaENETw <- idx2infoNew(estBetaENETw)

# 2a. add measures... 
#   ... mainTN =  main effects without simulated effect and which are not extracted 

# all possible main effects (without effect) - false positive main effects
# as this only counts main effects as true negative options
estBetaENETw$mainTN <- as.numeric(as.character(estBetaENETw$pTrash)) - estBetaENETw$mainFP

#   ... mainFN = simulated main effect that is not extracted
estBetaENETw$mainFN <- length(setParam$dgp$linEffects) - estBetaENETw$mainTP

# 2b. add measures...
#   ... accuracy = (TP + TN) / (TP + TN + FP + FN)
estBetaENETw$mainACC <- (estBetaENETw$mainTP + estBetaENETw$mainTN) / 
  (estBetaENETw$mainTP + estBetaENETw$mainTN + estBetaENETw$mainFP + estBetaENETw$mainFN)

#   ... sensitivity = TP / (TP + FN)
estBetaENETw$mainSensitivity <- estBetaENETw$mainTP / (estBetaENETw$mainTP + estBetaENETw$mainFN)

#   ... specificity = TN / (TN + FP)
estBetaENETw$mainSpecificity <- estBetaENETw$mainTN / (estBetaENETw$mainTN + estBetaENETw$mainFP)

#   ... balanced accuracy = (sensitivity + specificity) / 2 
estBetaENETw$mainBalACC <- (estBetaENETw$mainSpecificity + estBetaENETw$mainSensitivity) / 2

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter")
estBetaENETw[col2fac] <- lapply(estBetaENETw[col2fac], factor)
str(estBetaENETw)
# how often are predictors "selected" in model? beta-coefficient > 0 (i.e., NA)
# ! pvi for ENETw includes interactions
#   otherwise {ENETwo & GBM} pvi does not inlcude interactions

# save(estBetaENETw, file = "relFrequencyMeasuresBeta.rda")
################################################################################
# run ANOVAs for ...
#     ... different dependent variables {interTP, interFP, ...}
################################################################################
# there are only between-sample factors for the ANOVA
# generate ID along data set
estBetaENETw$ID <- seq_len(dim(estBetaENETw)[1])

dvVec <- c("mainTP", "mainFP", "mainPPV", 
           "mainACC", "mainSensitivity",
           "mainSpecificity", "mainBalACC")
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

pSensitivity <- plotEta2_4eachModel(eta2mainSensitivity, eta2Thresh, "#006600") + 
  ggtitle("Sensitivity") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pSpecificity <- plotEta2_4eachModel(eta2mainSpecificity, eta2Thresh, "#006600") + 
  ggtitle("Specificity") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))


pBalACC <- plotEta2_4eachModel(eta2mainBalACC, eta2Thresh, "#006600") + 
  ggtitle("BalACC") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pACC <- plotEta2_4eachModel(eta2mainACC, eta2Thresh, "#006600") + 
  ggtitle("ACC") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pSensitivity + pSpecificity + pBalACC + pACC

pTP <- plotEta2_4eachModel(eta2mainTP, eta2Thresh, "#006600") + 
  ggtitle("TP") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pFP <- plotEta2_4eachModel(eta2mainFP, eta2Thresh, "#006600") + 
  ggtitle("FP") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pPPV <- plotEta2_4eachModel(eta2mainPPV, eta2Thresh, "#006600") + 
  ggtitle("PPV") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

# pTP + pFP + pPPV 

pbetaInter <- pSensitivity + pSpecificity + pBalACC + 
  pTP + pFP + pPPV + plot_layout(ncol = 3)

# ggplot2::ggsave(filename = paste0(plotFolder, "/detectMains/betweenANOVA_betaMainMeasures.png"),
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

load(paste0(resFolder, "/relFrequencyMeasuresBeta.rda"))

plotMainBeta <- aggregate(cbind(mainTP, mainFP, mainPPV, mainACC,
                                mainSpecificity, mainSensitivity, mainBalACC) ~ 
                             model + N + pTrash + rel + R2 + lin_inter, 
                           data = estBetaENETw, 
                           function(x) {cbind(mean(x), 
                                              sd(x),
                                              quantile(x, 0.025),
                                              quantile(x, 0.975))})

plotMainBeta <- do.call(data.frame, plotMainBeta)
colnames(plotMainBeta) <- stringr::str_replace_all(colnames(plotMainBeta), "\\.1", "_M")
colnames(plotMainBeta) <- stringr::str_replace_all(colnames(plotMainBeta), "\\.2", "_SE")
colnames(plotMainBeta) <- stringr::str_replace_all(colnames(plotMainBeta), "\\.3", "_q025")
colnames(plotMainBeta) <- stringr::str_replace_all(colnames(plotMainBeta), "\\.4", "_q975")

plotMainBeta$N <- factor(plotMainBeta$N, levels = c(100, 300, 1000))

plotMainBeta <- tidyr::pivot_longer(plotMainBeta, 
                                     cols = !c(model, N, pTrash, rel, R2, lin_inter),
                                     names_to = c("DV", "measure"), 
                                     names_sep = "_",
                                     values_to = "values")

plotMainBeta <- tidyr::pivot_wider(plotMainBeta,
                                    names_from = measure, 
                                    values_from = values)

plotMainBeta$lin_inter <- factor(plotMainBeta$lin_inter, 
                              levels = c("0.2_0.8", "0.5_0.5", "0.8_0.2"),
                              labels = c("20:80", "50:50", "80:20"))

plotSensSpec <- function(data, DV, model, rel, 
                         quantiles = F, guides = T, 
                         title = "", yTitle = "") {
  pTmp <- ggplot(data[data$DV == DV &
                        data$model == model &
                        data$rel == rel,],
                 aes(x = N, y = M, 
                     group = interaction(lin_inter, pTrash), colour = lin_inter,
                     linetype = lin_inter, shape = pTrash)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.75) +
    scale_linetype_manual(name = "Lin:Inter", values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(name = "Noise", values = c(16, 17)) +
    scale_color_manual(name = "Lin:Inter", values = colValuesLin) +
    facet_grid2(~ R2,
                strip = strip_themed(
                  background_x = list(element_rect(fill = alpha(colValuesR2[1], 0.4)),
                                      element_rect(fill = alpha(colValuesR2[2], 0.4)),
                                      element_rect(fill = alpha(colValuesR2[3], 0.4))))) + 
    ylim(c(0, 1)) +
    ylab(yTitle) +
    xlab("N") +
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
          legend.position = c(.85, .15), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal")
  if (quantiles == T) {
    pTmp <- pTmp + geom_errorbar(aes(ymin = q025, ymax = q975), 
                                 width = 0.2, alpha = 0.4, 
                                 position = position_dodge(width = 0.5)) 
  }
  if (guides == F) {
    pTmp <- pTmp + guides(color = "none", shape = "none", linetype = "none")
  }
  pTmp
}

(pSensitivityENETw_beta <- plotSensSpec(plotMainBeta, DV = "mainSensitivity", 
                                        model = "ENETw", rel = 0.8,
                                        quantiles = F, yTitle = "Sensitivity", title = "A"))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/detectMains/mainSensitivity_betaENETwR2facette_rel08.png"),
#                 plot = pSensitivityENETw_beta,
#                 width = 13.08,
#                 height = 6.68,
#                 units = "in")

(pSpecificityENETw_beta <- plotSensSpec(plotMainBeta, DV = "mainSpecificity", 
                                        model = "ENETw", rel = 0.8,
                                        guides = F, quantiles = F, yTitle = "Specificity", title = "B"))
# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/detectMains/mainSpecificity_betaENETwR2facette_rel08.png"),
#                 plot = pSpecificityENETw_beta,
#                 width = 13.08,
#                 height = 6.68,
#                 units = "in")