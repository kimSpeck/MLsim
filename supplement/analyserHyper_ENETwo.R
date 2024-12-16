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

listDir <- dir(resFolder)
hyperENETwo <- loadRData(paste0(resFolder, "/hyperParametersSample_ENETwo.rda"))
hyperENETw <- loadRData(paste0(resFolder, "/hyperParametersSample_ENETw.rda"))

# pull data from nested list of all results (fullData)
hyperENETwo <- rbindSingleResults(hyperENETwo)
hyperENETw <- rbindSingleResults(hyperENETw)

colnames(hyperENETwo)
hyperENETwo[,c("nLin", "nInter", "all.T1F0", "nOthers")] <- list(NULL)
hyperENETw[,c("nLin", "nInter", "all.T1F0", "nOthers")] <- list(NULL)

hyperENETwo <- idx2infoNew(hyperENETwo)
hyperENETw <- idx2infoNew(hyperENETw)

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model", "sample")
hyperENETwo[col2fac] <- lapply(hyperENETwo[col2fac], factor)
hyperENETw[col2fac] <- lapply(hyperENETw[col2fac], factor)

col2num <- c("alpha", "lambda")
hyperENETwo[col2num] <- lapply(hyperENETwo[col2num], as.numeric)
hyperENETw[col2num] <- lapply(hyperENETw[col2num], as.numeric)
str(hyperENETwo)
str(hyperENETw)
################################################################################
hyperENETwo$ID <- seq_len(dim(hyperENETwo)[1])
hyperENETw$ID <- seq_len(dim(hyperENETw)[1])

modVec <- c("ENETw", "ENETwo")
dvVec <- c("alpha", "lambda") 

for (iModel in modVec) {
  for (iDV in dvVec){
    # only between-sample ANOVA
    anovaObjectName <- paste0("anova", iModel, "_", iDV)
    
    tmp <- aov_ez(id = "ID",
                  dv = iDV,
                  data = hyperENETwo,
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
    etaObjectName <- paste0("eta2", iModel, "_", iDV)
    assign(etaObjectName, tmpEta2) 
  }
}

eta2Thresh <- 0.01

# alpha changes with R2:lin_inter which is the weird disordinal interaction effect
#   for specificity in the ENETwo
# but alpha shows the same pattern of results for the ENETw
# thus, the hyperparameter selection is affected by identical simulated conditions
# -> therefore, check the range of hyperparameters in both models
#     only differing ranges might be able to explain the differences in result patterns
(pAlpha_ENETwo <- plotEta2_4eachModel(eta2ENETwo_alpha, eta2Thresh, "#009999") + 
  ggtitle("alpha") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25)))

(pAlpha_ENETw <- plotEta2_4eachModel(eta2ENETw_alpha, eta2Thresh, "#006600") + 
    ggtitle("alpha") + # for the main title
    theme(plot.title = element_text(hjust = 0.5, size = 25)))


(pLambda_ENETwo <- plotEta2_4eachModel(eta2ENETwo_lambda, eta2Thresh, "#009999") + 
    ggtitle("lambda") + # for the main title
    theme(plot.title = element_text(hjust = 0.5, size = 25)))

(pLambda_ENETw <- plotEta2_4eachModel(eta2ENETw_lambda, eta2Thresh, "#006600") + 
    ggtitle("lambda") + # for the main title
    theme(plot.title = element_text(hjust = 0.5, size = 25)))

(pHyper <- pAlpha_ENETw + pAlpha_ENETwo + 
    pLambda_ENETw + pLambda_ENETwo +
    plot_layout(ncol = 2))

# # save plot as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_hyperENETs.png"),
#                 plot = pHyper,
#                 width = 13.68,
#                 height = 15,
#                 units = "in")

################################################################################
# plot hyper parameter values (alpha and lambda for the ENET)
################################################################################

plotHyper <- aggregate(cbind(alpha, lambda) ~ 
                         model + N + pTrash + rel + R2 + lin_inter, 
                       data = rbind(hyperENETwo, hyperENETw), 
                       function(x) {cbind(mean(x), 
                                          sd(x),
                                          quantile(x, 0.025),
                                          quantile(x, 0.975))})

# turn matrices inside data frame into columns
plotHyper <- do.call(data.frame, plotHyper)
colnames(plotHyper) <- stringr::str_replace_all(colnames(plotHyper), "\\.1", "_M")
colnames(plotHyper) <- stringr::str_replace_all(colnames(plotHyper), "\\.2", "_SE")
colnames(plotHyper) <- stringr::str_replace_all(colnames(plotHyper), "\\.3", "_q025")
colnames(plotHyper) <- stringr::str_replace_all(colnames(plotHyper), "\\.4", "_q975")

plotHyper$N <- factor(plotHyper$N, levels = c(100, 300, 1000))
plotHyper$model <- factor(plotHyper$model, 
                          levels = c("ENETw", "ENETwo"))

plotHyper <- tidyr::pivot_longer(plotHyper, 
                                 cols = !c(model, N, pTrash, rel, R2, lin_inter),
                                 names_to = c("DV", "measure"), 
                                 names_sep = "_",
                                 values_to = "values")

plotHyper <- tidyr::pivot_wider(plotHyper,
                                names_from = measure, 
                                values_from = values)

plotHyper <- tidyr::separate(plotHyper, "lin_inter", 
                into = c("lin", "inter"), sep = "_", remove = F) 
plotHyper$modelR2 <- as.numeric(as.character(plotHyper$R2)) * as.numeric(plotHyper$lin)
plotHyper$modelR2 <- ifelse(plotHyper$model == "ENETw", 
                            as.numeric(as.character(plotHyper$R2)), 
                            as.numeric(as.character(plotHyper$R2)) * as.numeric(plotHyper$lin))
  
# pTrash does not make any meaningful difference
# reliability does not even appear in ANOVA analysis of alpha or lambda

# alpha = 0 -> ridge regression (no variable selection)
# alpha = 1 -> lasso regression (variable selection)
(pAlpha_rel1pTrash50 <- ggplot(plotHyper[plotHyper$DV == "alpha" &
                         plotHyper$rel == 1 & 
                         plotHyper$pTrash == 50,],
                            aes(x = N, y = M, 
                                group = interaction(R2, model), colour = R2,
                                linetype = model, shape = model)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_line(position = position_dodge(width = 0.5), alpha = 0.6) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_shape_manual(values = c(16, 1, 8)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), 
                  width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) +  
    scale_color_manual(values = colValues) +
    facet_grid(pTrash ~ lin_inter, labeller = label_both) +
    ylab("alpha") +
    xlab("N") +
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/hyperENET_alphaRel1pTrash50.png"),
#                 plot = pAlpha_rel1pTrash50,
#                 width = 13.08,
#                 height = 7.18,
#                 units = "in")

(pAlpha_effectiveR2 <- ggplot(plotHyper[plotHyper$DV == "alpha" &
                                           plotHyper$rel == 1 & 
                                           plotHyper$pTrash == 50 & 
                                           plotHyper$model == "ENETwo",],
                               aes(x = modelR2, y = M, 
                                   group = interaction(N, R2), colour = R2,
                                   linetype = N, shape = N)) +
    geom_point(size = 2) +
    geom_line(alpha = 0.6) +
    scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(values = c(16, 1, 8)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), 
                  width = 0.03, alpha = 0.4) +  
    scale_color_manual(values = colValues) +
    ylab("alpha") +
    xlab("effective R²") +
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/hyperENET_alphaEffectiveR2.png"),
#                 plot = pAlpha_effectiveR2,
#                 width = 13.08,
#                 height = 7.18,
#                 units = "in")

(pAlpha_effectiveR2 <- ggplot(plotHyper[plotHyper$DV == "alpha" &
                                          plotHyper$rel == 1 & 
                                          plotHyper$pTrash == 50 & 
                                          plotHyper$model == "ENETw",],
                              aes(x = modelR2, y = M, 
                                  group = interaction(N, R2), colour = R2,
                                  linetype = N, shape = N)) +
    geom_point(size = 2) +
    geom_line(alpha = 0.6) +
    scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(values = c(16, 1, 8)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), 
                  width = 0.03, alpha = 0.4) +  
    scale_color_manual(values = colValues) +
    ylab("alpha") +
    xlab("effective R²") +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))


(pLambda <- ggplot(plotHyper[plotHyper$DV == "lambda" &
                               plotHyper$rel == 1 & 
                               plotHyper$pTrash == 50,],
                               aes(x = N, y = M, 
                                   group = interaction(R2, model), colour = R2,
                                   linetype = model, shape = model)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_line(position = position_dodge(width = 0.5), alpha = 0.6) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_shape_manual(values = c(16, 1, 8)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), 
                  width = 0.2, alpha = 0.4, position = position_dodge(width = 0.5)) +  
    scale_color_manual(values = colValues) +
    facet_grid(pTrash ~ lin_inter, labeller = label_both) +
    ylab("lambda") +
    xlab("N") +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

(pLambda_effectiveR2 <- ggplot(plotHyper[plotHyper$DV == "lambda" &
                                          plotHyper$rel == 1 & 
                                          plotHyper$pTrash == 50 & 
                                          plotHyper$model == "ENETwo",],
                              aes(x = modelR2, y = M, 
                                  group = interaction(N, R2), colour = R2,
                                  linetype = N, shape = N)) +
    geom_point(size = 2) +
    geom_line(alpha = 0.6) +
    scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
    scale_shape_manual(values = c(16, 1, 8)) +
    geom_errorbar(aes(ymin = q025, ymax = q975), 
                  width = 0.03, alpha = 0.4) +  
    scale_color_manual(values = colValues) +
    ylab("lambda") +
    xlab("effective R²") +
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
# ggplot2::ggsave(filename = paste0(plotFolder, "/hyperENET_lambdaEffectiveR2.png"),
#                 plot = pLambda_effectiveR2,
#                 width = 13.08,
#                 height = 7.18,
#                 units = "in")