# calculate ANOVA to check for effects of manipulated variables in simulation 
#   interactions between manipulations 
library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²
library(emmeans)
library(ggplot2)

# load parameters and helper functions 
source("setParameters.R")
source("analysisTools.R")

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

################################################################################
# ANOVA - R² in test sample
################################################################################
# to do:
#   - für test R2 ist sample immer dasselber aber Hyperparameter kommen von anderem training

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

# ANOVA: 
## mit ENET - ohne
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

##### Haupteffekte 
# R2: (Manipulationscheck) Effekt von simulated R2 auf R2 banal 
#     ! aber Interaktionen mit R2 können informativ sein
# lin_inter: je mehr in linearen Effekt, desto größer erklärtes R2 
#     vielleicht getrieben durch ENETwo, in dem keine Interaktionen gefunden werden können
# rel: reliablere Items führen zu höherem R2 (siehe Jacobucci)?
# N: Je mehr Daten, desto höher R2
# model: R2 unterscheidet sich je nach Modell 
#     wenig überraschend, weil Modelle ohne Interaktionen schlechter sein müssen (ohne ENETwo?)
# pTrash: kleinster Effekt! Anzahl der Trash-Variablen macht geringsten Unterschied! 
#####

#### 2-fach Interaktionen
# R² x model: Interaktionen mit R² können aufschlussreich sein
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

# ggplot2::ggsave(filename = paste0(plotFolder, "/mixedANOVA_R2.eps"),
#                 plot = pEta2mixedFull,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/mixedANOVA_R2.png"),
#                 plot = pEta2mixedFull,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")


# pEta2TestR2.ordered <- print_html(
#   eta2TestR2.ordered[eta2TestR2.ordered$Eta2_generalized >= 0.01,],
#   digits = 2)

# pEtaFormat <- format(
#   eta2TestR2.ordered,
#   digits = 2,
#   output = c("text", "markdown", "html"))

################################################################################
# mixed ANOVA for R^2 (but only models with interactions {ENETw & GBM})
################################################################################
## ohne ENET - ohne
anovaTestR2sub <- aov_ez(id = "ID", 
                      dv = "Rsq_test",
                      data = rSquaredTest[rSquaredTest$model != "ENETwo",], 
                      between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"),
                      within = "model")
# summary(anovaTestR2)
# nice(anovaTestR2)

eta2TestR2sub <- eta_squared(
  anovaTestR2sub, # fitted model
  partial = FALSE, # not partial!
  generalized = TRUE, # generalized eta squared
  ci = 0.95,
  verbose = TRUE)

# sort generalized eta-squared results 
# which higher order interactions do we need to illustrate to report simulation results?
(eta2TestR2sub.ordered <- eta2TestR2sub[order(eta2TestR2sub$Eta2_generalized, decreasing = T),])

plotEta2mixedSub <- eta2TestR2sub.ordered[eta2TestR2sub.ordered$Eta2_generalized > eta2Thresh,]

sortIdx <- order(plotEta2mixedSub$Eta2_generalized, decreasing = T)
plotEta2mixedSub$Parameter <- factor(plotEta2mixedSub$Parameter, 
                                  levels = plotEta2mixedSub$Parameter[sortIdx])

# plot generalized eta^2 for mixed ANOVA with test R^2 as dependent variable and 
# ... model as within factor
# ... N, pTrash, rel, R2 and lin_inter as between factors
# comparison in test R2 between models with and without option to extract interaction effects is not fair
# -> next: plot generalized eta^2 for only ENETw and GBM which both can take interaction effects into account
(pEta2mixedSub <- ggplot(plotEta2mixedSub, aes(x = Parameter, y = Eta2_generalized)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  geom_text(aes(label=round(Eta2_generalized, 2)), 
            #angle = 90, hjust = 1.5, vjust=0.5, 
            angle = 0, vjust=1.5, 
            color="black", size=3.5)+
  ylab("generalisiertes eta^2") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.2, hjust=0.1)))

# ggplot2::ggsave(filename = paste0(plotFolder, "/mixedANOVA_sub_R2.eps"),
#                 plot = pEta2mixedSub,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/mixedANOVA_sub_R2.png"),
#                 plot = pEta2mixedSub,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

# compare mixed ANOVA results with and without ENETwo as level in within factor 
plotEta2mixedSub$algo <- rep("sub", dim(plotEta2mixedSub)[1]) 
plotEta2mixed$algo <- rep("full", dim(plotEta2mixed)[1]) 

plotEta2mixedCompare <- rbind(plotEta2mixedSub, plotEta2mixed)
plotEta2mixedCompare$algo <- factor(plotEta2mixedCompare$algo,
                                    levels = c("full", "sub"))

eta2Sums <- aggregate(Eta2_generalized ~ Parameter,
                      data = plotEta2mixedCompare, sum)
sortIdx <- order(eta2Sums$Eta2_generalized, decreasing = T)
plotEta2mixedCompare$Parameter <- factor(plotEta2mixedCompare$Parameter, 
                                        levels = eta2Sums$Parameter[sortIdx])

(pEta2mixedCompare <- ggplot(plotEta2mixedCompare, aes(x = Parameter, y = Eta2_generalized,
                 group = algo, fill = algo)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  geom_text(aes(label=round(Eta2_generalized, 2), group = algo), 
            angle = 90, hjust = 1.5, vjust=0.5, 
            #angle = 0, vjust=1.5, 
            position = position_dodge(width = .9), 
            color="black", size=3.5)+
  ylab("generalisiertes eta^2") +
  scale_fill_manual(values = c("#990000", "#006600")))

# ggplot2::ggsave(filename = paste0(plotFolder, "/mixedANOVA_compare_R2.eps"),
#                 plot = pEta2mixedCompare,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/mixedANOVA_compare_R2.png"),
#                 plot = pEta2mixedCompare,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

# save(rSquaredTest, 
#      file = "results/finalResults/dependentMeasures/rSquared_fromANOVA.rda")
################################################################################
# calculate ANOVA for each model separately
################################################################################
modVec <- unique(rSquaredTest$model)
for (iModel in modVec) {
  # only between factor ANOVA (only within factor was model)
  anovaObjectName <- paste0("anovaRes", iModel)
  tmp <- aov_ez(id = "ID",
                dv = "Rsq_test",
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

# Datensätze zusammenführen
eta2ENETw$model <- rep("ENETw", dim(eta2ENETw)[1])
eta2ENETwo$model <- rep("ENETwo", dim(eta2ENETwo)[1])
eta2GBM$model <- rep("GBM", dim(eta2GBM)[1])

# check how many parameters are still in depending on threshold
eta2Thresh <- 0.1
# dim(eta2ENETw[eta2ENETw$Eta2_generalized >= eta2Thresh,])
# dim(eta2ENETwo[eta2ENETwo$Eta2_generalized >= eta2Thresh,])
# dim(eta2GBM[eta2GBM$Eta2_generalized >= eta2Thresh,])

eta2 <- rbind(eta2ENETw[eta2ENETw$Eta2_generalized >= eta2Thresh,], 
              eta2ENETwo[eta2ENETwo$Eta2_generalized >= eta2Thresh,],
              eta2GBM[eta2GBM$Eta2_generalized >= eta2Thresh,])

# sort variables according to total eta2 across all models
eta2Sums <- aggregate(Eta2_generalized ~ Parameter, data = eta2, sum)
sortIdx <- order(eta2Sums$Eta2_generalized, decreasing = T)
eta2$Parameter <- factor(eta2$Parameter, 
                         levels = eta2Sums$Parameter[sortIdx])

eta2$model <- factor(eta2$model, 
                     levels = c("GBM", "ENETw", "ENETwo")) # obvious order of models

# als Balkendiagramme mit Cut-Off bei gen. eta² von .1 plotten
#   y-Achse unterschiedliche Variablen; Farben unterschiedliche Modelle
(pEta2Model <- ggplot(eta2, aes(x = Parameter, y = Eta2_generalized,
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
  
# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_R2.eps"),
#                 plot = pEta2Model,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_R2.png"),
#                 plot = pEta2Model,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")
################################################################################
# simuliertes R² als zusätzlicher Split
################################################################################
modVec <- unique(rSquaredTest$model)
r2Vec <- unique(rSquaredTest$R2)
for (iModel in modVec) {
  for (iR2 in r2Vec) {
    # only between factor ANOVA (only within factor was model)
    anovaObjectName <- paste0("anovaRes", iModel, "_", iR2)
    tmp <- aov_ez(id = "ID",
                  dv = "Rsq_test",
                  data = rSquaredTest[rSquaredTest$model == iModel &
                                        rSquaredTest$R2 == iR2,],
                  between = c("N" , "pTrash" , "rel" , "lin_inter"))
    assign(anovaObjectName, tmp)
    
    tmpEta2 <- eta_squared(
      tmp, # fitted model
      partial = FALSE, # not partial!
      generalized = TRUE, # generalized eta squared
      ci = 0.95,
      verbose = TRUE)
    
    tmpEta2$model <- rep(iModel, dim(tmpEta2)[1])
    tmpEta2$R2 <- rep(iR2, dim(tmpEta2)[1])
    
    # sort generalized eta-squared results
    # which higher order interactions do we need to illustrate to report simulation results?
    etaObjectName <- paste0("eta2", iModel, "_", iR2)
    assign(etaObjectName, tmpEta2[order(tmpEta2$Eta2_generalized, decreasing = T),])
  }
}

# #
# eta2Thresh <- 0.1
# dim(eta2ENETw[eta2ENETw$Eta2_generalized >= eta2Thresh,])
# dim(eta2ENETwo[eta2ENETwo$Eta2_generalized >= eta2Thresh,])
# dim(eta2GBM[eta2GBM$Eta2_generalized >= eta2Thresh,])

eta2_ModelR2 <- rbind(eta2ENETw_0.2, eta2ENETw_0.5, eta2ENETw_0.8,
                      eta2ENETwo_0.2, eta2ENETwo_0.5, eta2ENETwo_0.8,
                      eta2GBM_0.2, eta2GBM_0.5, eta2GBM_0.8)

eta2Thresh <- 0.05
eta2_ModelR2 <- eta2_ModelR2[eta2_ModelR2$Eta2_generalized >= eta2Thresh,]

# sort variables according to total eta2 across all models
eta2Sums <- aggregate(Eta2_generalized ~ Parameter, data = eta2_ModelR2, sum)
sortIdx <- order(eta2Sums$Eta2_generalized, decreasing = T)
eta2_ModelR2$Parameter <- factor(eta2_ModelR2$Parameter, 
                                 levels = eta2Sums$Parameter[sortIdx])

eta2_ModelR2$model <- factor(eta2_ModelR2$model, 
                     levels = c("GBM", "ENETw", "ENETwo")) # obvious order of models
eta2_ModelR2$R2 <- factor(eta2_ModelR2$R2, 
                             levels = c("0.2", "0.5", "0.8")) # obvious order of models

(pEta2ModelR2 <- ggplot(eta2_ModelR2, aes(x = Parameter, y = Eta2_generalized,
                                group = model, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    geom_text(aes(label=round(Eta2_generalized, 2), group = model), 
              angle = 90, hjust = 1.5, vjust=0.5, 
              #angle = 0, vjust=1.5, 
              position = position_dodge(width = .9), 
              color="black", size=3.5)+
    ylab("generalisiertes eta^2") +
    scale_fill_manual(values = c("#990000", "#006600", "#009999")) +
    facet_wrap(~ R2) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.2, hjust=0.1)))

# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_modelR2_R2.eps"),
#                 plot = pEta2ModelR2,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/betweenANOVA_modelR2_R2.png"),
#                 plot = pEta2ModelR2,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

################################################################################
##### Voraussetzungsprüfungen #####
################################################################################
# # residual vs fitted plot (Punkte sollten sich vertikal unsystematisch und 
# # über die gesamte x-Achse gleich streuend um die 0 verteilen)
# plot(x = fitted(anovaTestR2), y = resid(anovaTestR2)) 
# 
# ## Heteroskedastizität prüfen
# # standardabweichungen pro Zelle deskriptiv prüfen
# allSDs <- aggregate(Rsq_test ~ N * pTrash * R2 * rel * lin_inter * model, 
#                     rSquaredTest,
#                     function(x) sd(x)) 
# 
# plot(allSDs$Rsq_test)
# allSDs[which(allSDs$Rsq_test > 0.07),]

################################################################################
# ANOVA - overfit (ENET vs. GBM)
################################################################################
rSquaredTest$Rsq_overfit <- rSquaredTest$Rsq_train - rSquaredTest$Rsq_test 

anovaTestoverfit <- aov_ez(id = "ID", 
                      dv = "Rsq_overfit",
                      data = rSquaredTest, 
                      between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"),
                      within = "model")
# summary(anovaTestR2)
# nice(anovaTestR2)

eta2Testoverfit <- eta_squared(
  anovaTestoverfit, # fitted model
  partial = FALSE, # not partial!
  generalized = TRUE, # generalized eta squared
  ci = 0.95,
  verbose = TRUE)

# sort generalized eta-squared results 
# which higher order interactions do we need to illustrate to report simulation results?
(eta2Testoverfit.ordered <- eta2Testoverfit[order(eta2Testoverfit$Eta2_generalized, 
                                                  decreasing = T),])
# -> depending on N and model the overfit varies
# plot overfit generalized eta^2
eta2Thresh = 0.05
plotEta2overfit <- eta2Testoverfit.ordered[eta2Testoverfit.ordered$Eta2_generalized > eta2Thresh,]
sortIdx <- order(plotEta2overfit$Eta2_generalized, decreasing = T)
plotEta2overfit$Parameter <- factor(plotEta2overfit$Parameter, 
                                  levels = plotEta2overfit$Parameter[sortIdx])

# plot generalized eta^2 for mixed ANOVA with test R^2 as dependent variable and 
# ... model as within factor
# ... N, pTrash, rel, R2 and lin_inter as between factors
# comparison in test R2 between models with and without option to extract interaction effects is not fair
# -> next: plot generalized eta^2 for only ENETw and GBM which both can take interaction effects into account
(pEta2overfit <- ggplot(plotEta2overfit, aes(x = Parameter, y = Eta2_generalized)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    geom_text(aes(label=round(Eta2_generalized, 2)), 
              #angle = 90, hjust = 1.5, vjust=0.5, 
              angle = 0, vjust=1.5, 
              color="black", size=3.5)+
    ylab("generalisiertes eta^2") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.2, hjust=0.1)))

# ggplot2::ggsave(filename = paste0(plotFolder, "/mixedANOVA_overfit.eps"),
#                 plot = pEta2overfit,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/mixedANOVA_overfit.png"),
#                 plot = pEta2overfit,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")
##### post-tests ##### 
# for model 

# post tests für N nicht sonderlich interessant: je größer N, desto kleiner overfit 
# post tests für model: welches Model zeigt den größten Overfit? 
#   most overfit in GBM > ENET - mit > ENET - ohne
# averaged over the levels of: N, pTrash, R2, rel, lin_inter
(postOverfitModel <- emmeans(anovaTestoverfit, specs = "model"))
pairs(postOverfitModel, adjust="holm")
(pPostOverfitN <- plot(postOverfitModel, comparisons = TRUE))

# ggplot2::ggsave(filename = paste0(plotFolder, "/mixedANOVA_overfit_postN.png"),
#                 plot = pPostOverfitN,
#                 width = 865,
#                 height = 410,
#                 units = "px")

# Interaktion ändert nichts an der Reihenfolge! 
# bedingt auf die entsprechende sample size (bei Kontrolle für sample size) rank order
#     of models in terms of overfit still applies (GBM > ENETw > ENETwo)
(postOverfitModelxN <- emmeans(anovaTestoverfit, ~ N * model))
(pPostOverfitModelxN <- plot(postOverfitModelxN, comparisons = TRUE))

# ggplot2::ggsave(filename = paste0(plotFolder, "/mixedANOVA_overfit_postModelxN.png"),
#                 plot = pPostOverfitModelxN,
#                 width = 865,
#                 height = 810,
#                 units = "px")

# ##### obvious post test #####
# # je mehr Trash Variablen, desto mehr overfit
# (postOverfitpTrash <- emmeans(anovaTestoverfit, ~pTrash))
# plot(postOverfitpTrash, comparisons = TRUE)
# 
# # je größer R², desto weniger overfit
# (postOverfitR2 <- emmeans(anovaTestoverfit, ~ R2))
# plot(postOverfitR2, comparisons = TRUE)
# 
# # weniger overfit, wenn Interaktionen mehr Effekt bekommen 
# #   weil Interaktionen eher runterregularisiert werden
# (postOverfit_linInter <- emmeans(anovaTestoverfit, ~ lin_inter))
# plot(postOverfit_linInter, comparisons = TRUE)

# # pairwise comparison for model + Bonferroni-Holm-correction for these 3 comparissons
# pairs(postOverfitModel, adjust="holm") 
