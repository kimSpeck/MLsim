# compare PVI and ENET predictors
# compare results from coefficients to PVI results for ENETs
# is PVI > 1 equivalent to predictor selection based on regularization?

# can we somehow check this for GBMs as well?

# 1. read in PVI data and ENET model data
# 2. compare if variables are extracted or not?
# confusion matrix: extracted in ENET {yes, no} vs. PVI > 1 {yes, no}

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
dataList <- listDir[stringr::str_detect(listDir, "^estBetaSample")]
models <- stringr::str_extract(dataList, "_[:alpha:]*.rda$")
models <- stringr::str_sub(models, start = 2L, end = -5)
for (iData in seq_len(length(dataList))){
  objectName <- paste0("estBetaSample", models[iData])
  assign(objectName, loadRData(paste0(resFolder, "/", dataList[iData])))
}

# write list of matrices in one single data.frame
estBetaSampleENETw <- rbindSingleResults(estBetaSampleENETw) 
estBetaSampleENETwo <- rbindSingleResults(estBetaSampleENETwo) 

# write variables in columns instead of rownames
estBetaSampleENETwo <- rowNames2col(estBetaSampleENETwo)
estBetaSampleENETw <- rowNames2col(estBetaSampleENETw)

colnames(estBetaSampleENETwo)[colnames(estBetaSampleENETwo) == 'rowNames'] <- 'var'
colnames(estBetaSampleENETw)[colnames(estBetaSampleENETw) == 'rowNames'] <- 'var'

# # check if transformations were successfull
# estBetaSampleENETwo[1:6, 998:1004]
# unique(estBetaSampleENETwo$rowNames)

estBetaSampleENETwo <- tidyr::pivot_longer(estBetaSampleENETwo, 
                                           cols = !c(idxCondLabel, model, N_pTrash, var), 
                                           names_to = "sample",
                                           values_to = "betaEst")

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
estBetaSampleENETwo <- na.omit(estBetaSampleENETwo)
gc()

# 1. calculate number of true positives for each sample (linear and interactions)
## extract true linear effects
estBetaSampleENETwo$linTP <- ifelse(estBetaSampleENETwo$var %in% setParam$dgp$linEffects, 1, 0)
estBetaSampleENETw$linTP <- ifelse(estBetaSampleENETw$var %in% setParam$dgp$linEffects, 1, 0)

## extract true inter effects
# # check interaction labels and number of unique labels
# length(unique(estBetaSampleENETw$var))
interLabels <- stringr::str_replace_all(setParam$dgp$interEffects, "\\:", "\\.") 
estBetaSampleENETw$interTP <- ifelse(estBetaSampleENETw$var %in% interLabels, 1, 0)

# 2. calculate number of false positives for each sample ()
## extract false positive linear effects
estBetaSampleENETw$linFP <- ifelse(estBetaSampleENETw$var %in% setParam$dgp$linEffects |
                                     estBetaSampleENETw$var %in% interLabels, 0, 1)

## extract false psotive inter effects
estBetaSampleENETw$interFP <- ifelse(estBetaSampleENETw$var %in% interLabels |
                                       estBetaSampleENETw$var %in% setParam$dgp$linEffects, 0, 1)

# aggregate number of true positives and false negatives
estBetaENETw <- aggregate(cbind(linTP, linFP, interTP, interFP) ~ sample + idxCondLabel + model + N_pTrash, 
                      data = estBetaSampleENETw, sum)
rm(estBetaSampleENETw)
gc()

# 3. calculate positive predictive value (TP / (TP + FP)) for every sample
estBetaENETw$linPPV <- estBetaENETw$linTP / (estBetaENETw$linTP + estBetaENETw$linFP) 
estBetaENETw$interPPV <- estBetaENETw$interTP / (estBetaENETw$interTP + estBetaENETw$interFP) 

estBetaENETw <- idx2infoNew(estBetaENETw)

# 2a. add measures...
#   ... linTN = trash variable (or interaction effect) that is correctly not extracted
estBetaENETw$linTN <- (as.numeric(estBetaENETw$pTrash) + # trash variables
    choose(as.numeric(estBetaENETw$pTrash) + length(setParam$dgp$linEffects), 
           setParam$dgp$interDepth)) - # all possible interactions
  length(setParam$dgp$linEffects) -  # remove true interactions
  estBetaENETw$linFP

# estBetaENETw$linTN.b <- as.numeric(estBetaENETw$pTrash) +
#   length(setParam$dgp$linEffects) - estBetaENETw$linFP.b

# all possible interactions
estBetaENETw$interTN <- choose(as.numeric(estBetaENETw$pTrash) + length(setParam$dgp$linEffects), 
                                setParam$dgp$interDepth) - estBetaENETw$interFP

#   ... linFN = linear effect that is not extracted
estBetaENETw$linFN <- length(setParam$dgp$linEffects) - estBetaENETw$linTP
estBetaENETw$interFN <- length(setParam$dgp$interEffects) - estBetaENETw$interTP

# 2b. add measures...
#   ... accuracy = (TP + TN) / (TP + TN + FP + FN)
estBetaENETw$linACC <- (estBetaENETw$linTP + estBetaENETw$linTN) / (estBetaENETw$linTP + estBetaENETw$linTN + 
                                                                   estBetaENETw$linFP + estBetaENETw$linFN)

estBetaENETw$interACC <- (estBetaENETw$interTP + estBetaENETw$interTN) / (estBetaENETw$interTP + estBetaENETw$interTN + 
                                                                      estBetaENETw$interFP + estBetaENETw$interFN)

#   ... sensitivity = TP / (TP + FN)
estBetaENETw$linSensitivity <- estBetaENETw$linTP / (estBetaENETw$linTP + estBetaENETw$linFN)

estBetaENETw$interSensitivity <- estBetaENETw$interTP / (estBetaENETw$interTP + estBetaENETw$interFN)

#   ... specificity = TN / (TN + FP)
estBetaENETw$linSpecificity <- estBetaENETw$linTN / (estBetaENETw$linTN + estBetaENETw$linFP)

estBetaENETw$interSpecificity <- estBetaENETw$interTN / (estBetaENETw$interTN + estBetaENETw$interFP)

#   ... balanced accuracy = (sensitivity + specificity) / 2 
estBetaENETw$linBalACC <- (estBetaENETw$linSpecificity + estBetaENETw$linSensitivity) / 2

estBetaENETw$interBalACC <- (estBetaENETw$interSpecificity + estBetaENETw$interSensitivity) / 2

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model")
estBetaENETw[col2fac] <- lapply(estBetaENETw[col2fac], factor)

# to do: think about what counts as false positive linear or interaction effect

# how often are predictors "selected" in model? beta-coefficient > 0 (i.e., NA)
# ! pvi for ENETw includes interactions
#   otherwise {ENETwo & GBM} pvi does not inlcude interactions

################################################################################
# run ANOVAs for ...
#     ... different models {ENETwo, ENETw, GBM}
#     ... different dependent variables {PPV, linTP, linFP_woTInter}
################################################################################

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

(pInterSensitivity <- ggplot(plotInterBeta[plotInterBeta$DV == "interSensitivity" &
                                             plotInterBeta$rel == 1,],
                        aes(x = N, y = M, 
                            group = interaction(R2, model), colour = R2,
                            linetype = model, shape = model)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
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


(pInterSensitivity <- ggplot(plotInterBeta[plotInterBeta$DV == "interSensitivity",],
                        aes(x = rel, y = M, 
                            group = interaction(R2, N), colour = R2,
                            linetype = N, shape = N)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
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


plotLinBeta <- aggregate(cbind(linTP, linFP, linPPV, linACC,
                             linSpecificity, linSensitivity, linBalACC) ~ 
                         model + N + pTrash + rel + R2 + lin_inter, 
                       data = estBetaENETw, 
                       function(x) {cbind(mean(x), 
                                          sd(x),
                                          quantile(x, 0.025),
                                          quantile(x, 0.975))})

plotLinBeta <- do.call(data.frame, plotLinBeta)
colnames(plotLinBeta) <- stringr::str_replace_all(colnames(plotLinBeta), "\\.1", "_M")
colnames(plotLinBeta) <- stringr::str_replace_all(colnames(plotLinBeta), "\\.2", "_SE")
colnames(plotLinBeta) <- stringr::str_replace_all(colnames(plotLinBeta), "\\.3", "_q025")
colnames(plotLinBeta) <- stringr::str_replace_all(colnames(plotLinBeta), "\\.4", "_q975")

plotLinBeta$N <- factor(plotLinBeta$N, levels = c(100, 300, 1000))


plotLinBeta <- tidyr::pivot_longer(plotLinBeta, 
                                     cols = !c(model, N, pTrash, rel, R2, lin_inter),
                                     names_to = c("DV", "measure"), 
                                     names_sep = "_",
                                     values_to = "values")

plotLinBeta <- tidyr::pivot_wider(plotLinBeta,
                                    names_from = measure, 
                                    values_from = values)

(pLinSensitivity <- ggplot(plotLinBeta[plotLinBeta$DV == "linSensitivity" &
                                             plotLinBeta$rel == 0.6,],
                             aes(x = N, y = M, 
                                 group = interaction(R2, model), colour = R2,
                                 linetype = model, shape = model)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
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
