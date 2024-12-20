# calculate ANOVA for R² to check for effects of manipulated variables in simulation 
#   interactions between manipulations 
library(afex) # for aov_ez()
library(effectsize) # for efect size calculation; generalized eta²
library(emmeans)
library(ggplot2)

# load parameters and helper functions 
source("utils/setParameters.R")
source("utils/analysisTools.R")

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

# load results files
resFolder <- "results/finalResults/dependentMeasures" 
################################################################################
# ANOVA - R² in test sample
################################################################################
# load perform per sample data
listDir <- dir(resFolder)
dataList <- listDir[stringr::str_detect(listDir, "^performPerSample")]
models <- stringr::str_extract(dataList, "_[:alpha:]*.rda$")
models <- stringr::str_sub(models, start = 2L, end = -5)
for (iData in seq_len(length(dataList))){
  objectName <- paste0("pps", models[iData])
  assign(objectName, loadRData(paste0(resFolder, "/", dataList[iData])))
}

# pull data from nested list of all results (fullData)
ppsENETw <- rbindSingleResults(ppsENETw)
ppsENETwo <- rbindSingleResults(ppsENETwo)
ppsGBM <- rbindSingleResults(ppsGBM)

# get informative variables for simulated conditions (N, pTrash, R2, lin_inter) 
ppsENETw <- idx2infoNew(ppsENETw)
ppsENETwo <- idx2infoNew(ppsENETwo)
ppsGBM <- idx2infoNew(ppsGBM)

# in ANOVA sample will be the ID; but all samples run from 1:100 in all simulated 
#   conditions that represent independent samples (between factors); identical sample
#   names suggest within factor in ANOVA! 
#   -> introduce additional ID variable to indicate independent samples
# independent observations for samples {1:100} in different simulated conditions
#   that all have the same sample numbers
ppsENETw$ID <- seq_len(dim(ppsENETw)[1])
ppsENETwo$ID <- seq_len(dim(ppsENETwo)[1])
ppsGBM$ID <- seq_len(dim(ppsGBM)[1])

# # check
# all(colnames(ppsENETw) == colnames(ppsENETwo))
# all(colnames(ppsGBM) == colnames(ppsENETw))

# merge ENET_w, ENET_wo and GBM data; concatenate data for all models
rSquaredTest <- rbind(ppsENETw, ppsENETwo, ppsGBM)
str(rSquaredTest)
rSquaredTest$testR2_diff <- as.numeric(rSquaredTest$R2) - as.numeric(rSquaredTest$Rsq_test)
rSquaredTest$testR2_relational <- as.numeric(rSquaredTest$Rsq_test)/as.numeric(rSquaredTest$R2)

#range(rSquaredTest$testR2_diff)
#range(rSquaredTest$testR2_relational)

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model")
rSquaredTest[col2fac] <- lapply(rSquaredTest[col2fac], factor)
# change variables to numeric
chr2num <- c("RMSE_train", "Rsq_train", "MAE_train", "RMSE_test", "Rsq_test", "MAE_test")
rSquaredTest[chr2num] <- lapply(rSquaredTest[chr2num], as.numeric)

################################################################################
# mixed ANOVA for R^2
################################################################################
# mixed ANOVA with ...
# ... id = sample but ID (independent samples between simulated conditions)
# ... dv = {R^2}
# ... between: 3 x 2 x 3 x 3 x 3
#   N (3) {100, 300, 1000}
#   pTrash (2) {10, 50}
#   rel (3) {0.6, 0.8, 1}
#   R2 (3) {0.2, 0.5, 0.8}
#   lin_inter (3) {0.2_0.8, 0.5_0.5, 0.8_0.2}
# ... within: 
#   model {Enetw; Enetwo, GBM}

# Type 3 sums of squaress (e.g., Maxwell and Delaney, 2004)
# contr.sum for effects-coding for the categorical variables
# outcome: generalized eta2 (Olejnik and Algina 2003)

# ANOVA: 
anovaTestR2 <- aov_ez(id = "ID", 
                  dv = "Rsq_test",
                  data = rSquaredTest, 
                  between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"),
                  within = "model")
# summary(anovaTestR2)
# nice(anovaTestR2)

eta2TestR2 <- eta_squared(
  anovaTestR2, # fitted model
  partial = FALSE, # not partial!
  generalized = TRUE, # generalized eta squared
  ci = 0.95,
  verbose = TRUE)

# sort generalized eta-squared results 
# which higher order interactions do we need to illustrate to report simulation results?
(eta2TestR2.ordered <- eta2TestR2[order(eta2TestR2$Eta2_generalized, decreasing = T),])

#### 2-way interaction
# R² x model: interactions with R² can be interesting
# lin vs. inter x model: 
eta2Thresh <- 0.1
plotEta2mixed <- eta2TestR2.ordered[eta2TestR2.ordered$Eta2_generalized > eta2Thresh,]

sortIdx <- order(plotEta2mixed$Eta2_generalized, decreasing = T)
plotEta2mixed$Parameter <- factor(plotEta2mixed$Parameter, 
                         levels = plotEta2mixed$Parameter[sortIdx])

# plot generalized eta^2 for mixed ANOVA with test R^2 as dependent variable and 
# ... model as within factor
# ... N, pTrash, rel, R2 and lin_inter as between factors
# comparison in test R2 between models with and without option to extract interaction effects is not fair
# -> next: plot generalized eta^2 for only ENETw and GBM which both can take interaction effects into account
(pEta2mixedFull <- ggplot(plotEta2mixed, aes(x = Parameter, y = Eta2_generalized)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  geom_text(aes(label=round(Eta2_generalized, 2)), 
            #angle = 90, hjust = 1.5, vjust=0.5, 
            angle = 0, vjust=1.5, 
            color="black", size=3.5)+
  ylab("generalisiertes eta^2") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.2, hjust=0.1)))

# ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_R2_thresh", eta2Thresh, ".png"),
#                 plot = pEta2mixedFull,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

eta2modelTable <- plotEta2mixed[plotEta2mixed$Eta2_generalized >= eta2Thresh,]

print(xtable::xtable(eta2modelTable, type = "latex"), 
      file = paste0(plotFolder, "/ANOVAresults/mixedANOVA_R2_thresh", eta2Thresh, ".tex"))

# save(rSquaredTest,
#      file = "results/finalResults/dependentMeasures/rSquaredData_eachSample.rda")

################################################################################
# calculate ANOVA for each model separately
################################################################################
modVec <- unique(rSquaredTest$model)
for (iModel in modVec) {
  # only between factor ANOVA (only within factor was model)
  anovaObjectName <- paste0("anovaRes", iModel)
  tmp <- aov_ez(id = "ID",
                dv = "Rsq_test", # absolut R2
                #dv = "testR2_relational", # relational R2
                data = rSquaredTest[rSquaredTest$model == iModel,],
                between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"))
  assign(anovaObjectName, tmp)

  tmpEta2 <- eta_squared(
    tmp, # fitted model
    partial = FALSE, # not partial!
    generalized = TRUE, # generalized eta squared
    ci = 0.95,
    verbose = TRUE)

  # sort generalized eta-squared results
  # which higher order interactions do we need to illustrate to report simulation results?
  etaObjectName <- paste0("eta2", iModel)
  assign(etaObjectName, tmpEta2[order(tmpEta2$Eta2_generalized, decreasing = T),])
}

summary(anovaResENETw)
summary(anovaResENETwo)
summary(anovaResGBM)

# merge datasets
eta2ENETw$model <- rep("ENETw", dim(eta2ENETw)[1])
eta2ENETwo$model <- rep("ENETwo", dim(eta2ENETwo)[1])
eta2GBM$model <- rep("GBM", dim(eta2GBM)[1])

eta2Thresh <- 0.1
##### eta2 table #####
eta2Table <- rbind(eta2ENETw, eta2ENETwo, eta2GBM) 
eta2Table <- tidyr::pivot_wider(eta2Table[,c("Parameter", "Eta2_generalized", "model")], 
                                names_from = model, 
                                values_from = Eta2_generalized)

eta2Table <- subset(eta2Table, eta2Table$ENETw >= eta2Thresh | 
                      eta2Table$ENETwo >= eta2Thresh |
                      eta2Table$GBM >= eta2Thresh)

eta2Table$sumEta2 <- apply(eta2Table[,c("ENETw", "ENETwo", "GBM")], 1, sum)
eta2Table$M <- apply(eta2Table[,c("ENETw", "ENETwo", "GBM")], 1, mean)
eta2Table <- dplyr::arrange(eta2Table, desc(M))
eta2Table <- eta2Table[,c("Parameter", "M", "GBM", "ENETw", "ENETwo")]

print(xtable::xtable(eta2Table, type = "latex"), 
      file = paste0(plotFolder, "/ANOVAresults/betweenANOVA_R2.tex"))

##### plot eta2 #####
# check how many parameters are still in depending on threshold
# dim(eta2ENETw[eta2ENETw$Eta2_generalized >= eta2Thresh,])
# dim(eta2ENETwo[eta2ENETwo$Eta2_generalized >= eta2Thresh,])
# dim(eta2GBM[eta2GBM$Eta2_generalized >= eta2Thresh,])

# eta2 <- rbind(eta2ENETw[eta2ENETw$Eta2_generalized >= eta2Thresh,], 
#               eta2ENETwo[eta2ENETwo$Eta2_generalized >= eta2Thresh,],
#               eta2GBM[eta2GBM$Eta2_generalized >= eta2Thresh,])

eta2 <- rbind(eta2ENETw, eta2ENETwo, eta2GBM)

# sort variables according to total eta2 across all models
eta2Sums <- aggregate(Eta2_generalized ~ Parameter, data = eta2, sum)
sortIdx <- order(eta2Sums$Eta2_generalized, decreasing = T)
eta2$Parameter <- factor(eta2$Parameter, 
                         levels = eta2Sums$Parameter[sortIdx])

eta2$model <- factor(eta2$model, 
                     levels = c("GBM", "ENETw", "ENETwo")) # obvious order of models

# as barplot with cut-off: gen. eta² = .1 
(pEta2Model <- ggplot(eta2[eta2$Eta2_generalized >= eta2Thresh,], 
                      aes(x = Parameter, y = Eta2_generalized,
                          group = model, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    geom_text(aes(label=round(Eta2_generalized, 2), group = model), 
              angle = 90, hjust = 1.5, vjust=0.5, 
              #angle = 0, vjust=1.5, 
              position = position_dodge(width = .9), 
              color="black", size=3.5)+
    ylab("generalisiertes eta^2") +
    scale_fill_manual(values = c("#990000", "#006600", "#009999")) +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))


(pEta2_stacked <- ggplot(eta2[eta2$Eta2_generalized >= eta2Thresh,], 
                      aes(x = model, y = Eta2_generalized,
                          group = Parameter, fill = Parameter)) +
    geom_col(position = position_stack(reverse = TRUE)) +
    ylab("generalisiertes eta^2") +
  geom_text(aes(label=round(Eta2_generalized, 2), group = Parameter), 
            angle = 0, vjust=1.75, 
            position = position_stack(reverse = TRUE), 
            color="black", size=3.5)+
    scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3',
                                 '#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd')) +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))


# ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/betweenANOVA_R2.png"),
#                 plot = pEta2Model,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/betweenANOVA_R2stacked.png"),
#                 plot = pEta2_stacked,
#                 width = 13.63,
#                 height = 12.07,
#                 units = "in")