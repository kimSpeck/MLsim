# load parameters and helper functions 
source("utils/setParameters.R")
source("utils/analysisTools.R")

# ANOVA
library(afex) # für aov_ez()
library(effectsize) # für Berechnung von Effektstärken; generalisiertes eta²

# plot results
library(ggplot2)
library(ggh4x)
library(patchwork)
library(see)

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
allLinVars <- paste0("Var", seq_len(max(setParam$dgp$pTrash)))
# every beta coefficient != 0 that is ...
#   ... an interaction with simulated effects or is a main effect receives 0
#   ... an interaction without simulated effects and is not a main effect receives 1
pviENETw$interFP <- ifelse(pviENETw$pviRank %in% interLabels |
                               pviENETw$pviRank %in% allLinVars, 0, 1)

# sum up number of true positives and false negatives within each sample
interENETw <- aggregate(cbind(interTP, interFP) ~ sample + idxCondLabel + model + N_pTrash, 
                          data = pviENETw, sum)

# 3. calculate positive predictive value (TP / (TP + FP)) for every sample
interENETw$interPPV <- interENETw$interTP / (interENETw$interTP + interENETw$interFP) 

interENETw <- idx2infoNew(interENETw)

# 2. add measures... 
#   ... interTN =  interactions without simulated effect and which are not extracted 

# all possible interactions - true interactions (with simulated effects) - false positive interactions
# as this only counts interactions as true negative options
interENETw$interTN <- choose(as.numeric(interENETw$pTrash) + length(setParam$dgp$linEffects), 
                                 setParam$dgp$interDepth) - 
  length(setParam$dgp$interEffects) - interENETw$interFP

#   ... interFN = simulated interaction effect that is not extracted
interENETw$interFN <- length(setParam$dgp$interEffects) - interENETw$interTP

# 2b. add measures...
#   ... accuracy = (TP + TN) / (TP + TN + FP + FN)
interENETw$interACC <- (interENETw$interTP + interENETw$interTN) / 
  (interENETw$interTP + interENETw$interTN + interENETw$interFP + interENETw$interFN)

#   ... sensitivity = TP / (TP + FN)
interENETw$interSensitivity <- interENETw$interTP / (interENETw$interTP + interENETw$interFN)

#   ... specificity = TN / (TN + FP)
interENETw$interSpecificity <- interENETw$interTN / (interENETw$interTN + interENETw$interFP)

#   ... balanced accuracy = (sensitivity + specificity) / 2 
interENETw$interBalACC <- (interENETw$interSpecificity + interENETw$interSensitivity) / 2

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter")
interENETw[col2fac] <- lapply(interENETw[col2fac], factor)
str(interENETw)

# average sensitivity across all simulated conditions (for interaction effects in the ENETinter)
mean(interENETw$interSensitivity)

# # save data for sensitivity and specificity result plots
# save(interENETw, file = paste0(resFolder, "/interENETinterMeasures.rda"))
################################################################################
# run ANOVAs for ...
#     ... different dependent variables {interTP, interFP, ...}
################################################################################
# there are only between-sample factors for the ANOVA
# generate ID along data set
interENETw$ID <- seq_len(dim(interENETw)[1])

dvVec <- c("interTP", "interFP", "interPPV", 
           "interACC", "interSensitivity", 
           "interSpecificity", "interBalACC")
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

pSensitivity <- plotEta2_4eachModel(eta2interSensitivity, eta2Thresh, "#006600") + 
  ggtitle("Sensitivity") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pSpecificity <- plotEta2_4eachModel(eta2interSpecificity, eta2Thresh, "#006600") + 
  ggtitle("Specificity") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pBalACC <- plotEta2_4eachModel(eta2interBalACC, eta2Thresh, "#006600") + 
  ggtitle("BalACC") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pACC <- plotEta2_4eachModel(eta2interACC, eta2Thresh, "#006600") + 
  ggtitle("ACC") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

# pSensitivity + pSpecificity + pBalACC + pACC 

pTP <- plotEta2_4eachModel(eta2interTP, eta2Thresh, "#006600") + 
  ggtitle("TP") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pFP <- plotEta2_4eachModel(eta2interFP, eta2Thresh, "#006600") + 
  ggtitle("FP") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

pPPV <- plotEta2_4eachModel(eta2interPPV, eta2Thresh, "#006600") + 
  ggtitle("PPV") + # for the main title
  theme(plot.title = element_text(hjust = 0.5, size = 25))

# pTP + pFP + pPPV

pPviInter <- pSensitivity + pSpecificity + pBalACC + 
  pTP + pFP + pPPV + plot_layout(ncol = 3)

# ggplot2::ggsave(filename = paste0(plotFolder, "/detectInteractions/betweenANOVA_pviInterMeasures.png"),
#                 plot = pPviInter,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

# results for pvi and beta are exactly the same

