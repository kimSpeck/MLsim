# calculate ANOVA for R² to check for effects of manipulated variables in simulation 
#   interactions between manipulations 
library(afex) # for aov_ez()
library(effectsize) # for efect size calculation; generalized eta²
library(emmeans)
library(ggplot2)
library(tidyverse)
library(cowplot)

# load parameters and helper functions 
source("utils/setParameters.R")
source("utils/analysisTools.R")

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

resFolder <- "results/"

depMeasureFolder = paste0(resFolder, "dependentMeasures")
if (!file.exists(depMeasureFolder)){
  dir.create(depMeasureFolder)
}

# load results files
dgpVec <- c("inter", "pwlinear", "nonlinear3")
resFolderVec <- paste0("results/", dgpVec, "/dependentMeasures")

################################################################################
# ANOVA - R² in test sample
################################################################################
# load perform per sample data
listDir <- dir(resFolderVec)
dataList <- listDir[stringr::str_detect(listDir, "^performPerSample")]

dgps <- stringr::str_extract(dataList, "_[:alpha:]*3{0,1}_")
dgps <- stringr::str_sub(dgps, start = 2L, end = -2)

models <- stringr::str_extract(dataList, "_[:alpha:]*.rda$")
models <- stringr::str_sub(models, start = 2L, end = -5)

for (iData in seq_len(length(dataList))){
  objectName <- paste0("pps", models[iData], "_", dgps[iData])
  assign(objectName, loadRData(paste0(resFolder, "/", dgps[iData], 
                                      "/dependentMeasures/", dataList[iData])))
}

ppsList <- sapply(grep("pps", ls(), value = TRUE), get, simplify = FALSE)

ppsList <- lapply(ppsList, function(iPPS) {
  # pull data from nested list of all results (fullData)
  rbindSingleResults(iPPS)
})

ppsList <- lapply(ppsList, function(iPPS) {
  # get informative variables for simulated conditions (N, pTrash, R2, lin_inter)
  idx2infoNew(iPPS)
})

# concatenate data for all models
rSquaredTest <- do.call(rbind, ppsList)

# remove single objects since all information is in rSquaredTest, now
rmObjects <- ls(pattern = "pps")
rm(list = rmObjects)

rSquaredTest$testR2_relE <- (as.numeric(rSquaredTest$Rsq_test) - as.numeric(rSquaredTest$R2)) / as.numeric(rSquaredTest$R2) 

# change variables to factors
col2fac <- c("N", "pTrash" , "R2" , "rel" , "lin_inter", "model", "dgp")
rSquaredTest[col2fac] <- lapply(rSquaredTest[col2fac], factor)
# change variables to numeric
chr2num <- c("RMSE_train", "Rsq_train", "MAE_train", "RMSE_test", "Rsq_test", "MAE_test")
rSquaredTest[chr2num] <- lapply(rSquaredTest[chr2num], as.numeric)

################################################################################
# mixed ANOVA for R^2
################################################################################
# mixed ANOVA with ...
# ... id = sample but ID (independent samples between simulated conditions)
# ... dv = {rel R^2 error}
# ... between: 3 x 2 x 3 x 3 x 3
#   N (3) {100, 300, 1000}
#   pTrash (2) {10, 50}
#   rel (3) {0.6, 0.8, 1}
#   R2 (3) {0.2, 0.5, 0.8}
#   lin_inter (3) {0.2_0.8, 0.5_0.5, 0.8_0.2}
#   dgp (4) {inter, pwlinear, nonlinear3, nonlinear}
# ... within:
#   model {Enetw, Enetwo, GBM, RF}

# Type 3 sums of squaress (e.g., Maxwell and Delaney, 2004)
# contr.sum for effects-coding for the categorical variables
# outcome: generalized eta2 (Olejnik and Algina 2003)

rSquaredTest <- rSquaredTest %>%
  group_by(model) %>%
  mutate(ID = seq_len(n())) %>%
  ungroup()

# ANOVA:
anovaTestR2 <- aov_ez(id = "ID",
                      dv = "testR2_relE",
                      data = rSquaredTest,
                      between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter", "dgp"),
                      within = "model")

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
eta2Thresh <- 0.05
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

# ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_relErrorR2_thresh", eta2Thresh, ".png"),
#                 plot = pEta2mixedFull,
#                 width = 17.52,
#                 height = 10.76,
#                 units = "in")

eta2modelTable <- plotEta2mixed[plotEta2mixed$Eta2_generalized >= eta2Thresh,]

print(xtable::xtable(eta2modelTable, type = "latex"),
      file = paste0(plotFolder, "/ANOVAresults/mixedANOVA_relErrR2_thresh", eta2Thresh, ".tex"))

################################################################################
# expected marginal means (post hoc)
################################################################################
# estimated marginal means (EMMs), which are the predicted means adjusted for the 
#   effects of other factors in the model; only categorical variables and therefore
#   equivalent to means
# the more negative values == less proportion of R² found in the data

##### main effects as expected #####
(emRel <- emmeans(anovaTestR2, ~ rel))
# rel   emmean         SE     df lower.CL upper.CL
# 0.6 -0.62060 0.00013744 485514 -0.62087 -0.62033
# 0.8 -0.48041 0.00013744 485514 -0.48068 -0.48014
# 1   -0.33096 0.00013744 485514 -0.33123 -0.33069

(emN <- emmeans(anovaTestR2, ~ N))
# N      emmean         SE     df lower.CL upper.CL
# 100  -0.62958 0.00013744 485514 -0.62985 -0.62931
# 300  -0.44885 0.00013744 485514 -0.44912 -0.44858
# 1000 -0.35354 0.00013744 485514 -0.35381 -0.35327

(emR2 <- emmeans(anovaTestR2, ~ R2))
# R2    emmean        SE     df lower.CL upper.CL
# 0.2 -0.61509 0.00013744 485514 -0.61536 -0.61482
# 0.5 -0.44743 0.00013744 485514 -0.44770 -0.44716
# 0.8 -0.36945 0.00013744 485514 -0.36972 -0.36918

(emLin_inter <- emmeans(anovaTestR2, ~ lin_inter))
# lin_inter   emmean         SE     df lower.CL upper.CL
# 0.2_0.8   -0.58481 0.00013744 485514 -0.58508 -0.58454
# 0.5_0.5   -0.47576 0.00013744 485514 -0.47603 -0.47549
# 0.8_0.2   -0.37140 0.00013744 485514 -0.37167 -0.37113

##### interactions #####
(emDGPxModel <- emmeans(anovaTestR2, "model", by = "dgp"))
# dgp = inter:
#   model   emmean        SE     df lower.CL upper.CL
# ENETw  -0.3712 0.0002210 485514  -0.3716  -0.3708
# ENETwo -0.6225 0.0001707 485514  -0.6228  -0.6221
# GBM    -0.5446 0.0001903 485514  -0.5449  -0.5442
# RF     -0.4862 0.0001503 485514  -0.4865  -0.4859
# 
# dgp = nonlinear3:
#   model   emmean        SE     df lower.CL upper.CL
# ENETw  -0.5150 0.0002210 485514  -0.5154  -0.5146
# ENETwo -0.4236 0.0001707 485514  -0.4240  -0.4233
# GBM    -0.5089 0.0001903 485514  -0.5093  -0.5085
# RF     -0.4629 0.0001503 485514  -0.4632  -0.4626
# 
# dgp = pwlinear:
#   model   emmean        SE     df lower.CL upper.CL
# ENETw  -0.4840 0.0002210 485514  -0.4844  -0.4836
# ENETwo -0.3865 0.0001707 485514  -0.3868  -0.3862
# GBM    -0.4875 0.0001903 485514  -0.4879  -0.4872
# RF     -0.4350 0.0001503 485514  -0.4353  -0.4347

(emLinInterxModelxDGP <- emmeans(anovaTestR2, "lin_inter", by = c("model", "dgp")))
# model = ENETw, dgp = inter:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.4089 0.0003828 485514  -0.4096  -0.4081
# 0.5_0.5   -0.3806 0.0003828 485514  -0.3813  -0.3798
# 0.8_0.2   -0.3242 0.0003828 485514  -0.3249  -0.3234
# 
# model = ENETwo, dgp = inter:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.8794 0.0002956 485514  -0.8800  -0.8788
# 0.5_0.5   -0.6229 0.0002956 485514  -0.6235  -0.6223
# 0.8_0.2   -0.3651 0.0002956 485514  -0.3657  -0.3645
# 
# model = GBM, dgp = inter:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.6909 0.0003297 485514  -0.6916  -0.6903
# 0.5_0.5   -0.5328 0.0003297 485514  -0.5334  -0.5321
# 0.8_0.2   -0.4100 0.0003297 485514  -0.4106  -0.4093
# 
# model = RF, dgp = inter:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.6254 0.0002603 485514  -0.6259  -0.6249
# 0.5_0.5   -0.4759 0.0002603 485514  -0.4764  -0.4753
# 0.8_0.2   -0.3573 0.0002603 485514  -0.3578  -0.3567
# 
# model = ENETw, dgp = nonlinear3:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.6385 0.0003828 485514  -0.6393  -0.6378
# 0.5_0.5   -0.5188 0.0003828 485514  -0.5196  -0.5181
# 0.8_0.2   -0.3876 0.0003828 485514  -0.3884  -0.3869
# 
# model = ENETwo, dgp = nonlinear3:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.5438 0.0002956 485514  -0.5444  -0.5432
# 0.5_0.5   -0.4155 0.0002956 485514  -0.4160  -0.4149
# 0.8_0.2   -0.3116 0.0002956 485514  -0.3122  -0.3110
# 
# model = GBM, dgp = nonlinear3:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.5733 0.0003297 485514  -0.5739  -0.5726
# 0.5_0.5   -0.5143 0.0003297 485514  -0.5149  -0.5136
# 0.8_0.2   -0.4392 0.0003297 485514  -0.4399  -0.4386
# 
# model = RF, dgp = nonlinear3:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.5423 0.0002603 485514  -0.5428  -0.5418
# 0.5_0.5   -0.4614 0.0002603 485514  -0.4619  -0.4609
# 0.8_0.2   -0.3849 0.0002603 485514  -0.3854  -0.3844
# 
# model = ENETw, dgp = pwlinear:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.5878 0.0003828 485514  -0.5885  -0.5870
# 0.5_0.5   -0.4875 0.0003828 485514  -0.4882  -0.4867
# 0.8_0.2   -0.3768 0.0003828 485514  -0.3775  -0.3760
# 
# model = ENETwo, dgp = pwlinear:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.4842 0.0002956 485514  -0.4848  -0.4837
# 0.5_0.5   -0.3784 0.0002956 485514  -0.3789  -0.3778
# 0.8_0.2   -0.2969 0.0002956 485514  -0.2975  -0.2963
# 
# model = GBM, dgp = pwlinear:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.5423 0.0003297 485514  -0.5430  -0.5417
# 0.5_0.5   -0.4900 0.0003297 485514  -0.4906  -0.4893
# 0.8_0.2   -0.4303 0.0003297 485514  -0.4310  -0.4297
# 
# model = RF, dgp = pwlinear:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.5009 0.0002603 485514  -0.5014  -0.5004
# 0.5_0.5   -0.4314 0.0002603 485514  -0.4319  -0.4308
# 0.8_0.2   -0.3729 0.0002603 485514  -0.3734  -0.3724

(emNxR2 <- emmeans(anovaTestR2, "N", by = "R2"))
# R2 = 0.2:
#   N     emmean        SE     df lower.CL upper.CL
# 100  -0.8008 0.0002381 485514  -0.8012  -0.8003
# 300  -0.6049 0.0002381 485514  -0.6054  -0.6045
# 1000 -0.4396 0.0002381 485514  -0.4400  -0.4391
# 
# R2 = 0.5:
#   N     emmean        SE     df lower.CL upper.CL
# 100  -0.6055 0.0002381 485514  -0.6060  -0.6050
# 300  -0.3981 0.0002381 485514  -0.3986  -0.3977
# 1000 -0.3387 0.0002381 485514  -0.3391  -0.3382
# 
# R2 = 0.8:
#   N     emmean        SE     df lower.CL upper.CL
# 100  -0.4825 0.0002381 485514  -0.4829  -0.4820
# 300  -0.3435 0.0002381 485514  -0.3440  -0.3430
# 1000 -0.2824 0.0002381 485514  -0.2829  -0.2819

(emNxModel <- emmeans(anovaTestR2, "N", by = "model"))
# model = ENETw:
#   N      emmean         SE     df lower.CL upper.CL
# 100  -0.65398 0.00022102 485514 -0.65442 -0.65355
# 300  -0.41422 0.00022102 485514 -0.41465 -0.41379
# 1000 -0.30201 0.00022102 485514 -0.30245 -0.30158
# 
# model = ENETwo:
#   N      emmean         SE     df lower.CL upper.CL
# 100  -0.58416 0.00017069 485514 -0.58449 -0.58382
# 300  -0.44599 0.00017069 485514 -0.44632 -0.44565
# 1000 -0.40244 0.00017069 485514 -0.40277 -0.40210
# 
# model = GBM:
#   N      emmean         SE     df lower.CL upper.CL
# 100  -0.67920 0.00019035 485514 -0.67957 -0.67883
# 300  -0.49620 0.00019035 485514 -0.49658 -0.49583
# 1000 -0.36560 0.00019035 485514 -0.36598 -0.36523
# 
# model = RF:
#   N      emmean         SE     df lower.CL upper.CL
# 100  -0.60098 0.00015029 485514 -0.60127 -0.60069
# 300  -0.43898 0.00015029 485514 -0.43927 -0.43868
# 1000 -0.34412 0.00015029 485514 -0.34441 -0.34382

(emNxPtrash <- emmeans(anovaTestR2, "N", by = "pTrash"))
# pTrash = 10:
#   N     emmean        SE     df lower.CL upper.CL
# 100  -0.5754 0.0001944 485514  -0.5758  -0.5750
# 300  -0.4469 0.0001944 485514  -0.4473  -0.4465
# 1000 -0.3567 0.0001944 485514  -0.3571  -0.3563
# 
# pTrash = 50:
#   N     emmean        SE     df lower.CL upper.CL
# 100  -0.6838 0.0001944 485514  -0.6841  -0.6834
# 300  -0.4508 0.0001944 485514  -0.4512  -0.4504
# 1000 -0.3504 0.0001944 485514  -0.3508  -0.3500

# save expected marginal means into r object
df_emDGPxModel <- as.data.frame(emDGPxModel)
df_emLinInterxModelxDGP <- as.data.frame(emLinInterxModelxDGP)
df_emNxR2 <- as.data.frame(emNxR2)
df_emNxModel <- as.data.frame(emNxModel)
df_emNxPtrash <- as.data.frame(emNxPtrash)
rm(emDGPxModel, emLinInterxModelxDGP, emNxR2, emNxModel, emNxPtrash)
save(df_emDGPxModel, df_emLinInterxModelxDGP, 
     df_emNxR2, df_emNxModel, df_emNxPtrash, file = paste0(depMeasureFolder, "/mixedANOVA_postHocEMMEANS.rda"))


##### plot data #####
load(paste0(depMeasureFolder, "/mixedANOVA_postHocEMMEANS.rda"))

colRamp <- colorRampPalette(c("#6f9c3d", "#b8c36b", "#ffb366", "#ff8829", "#fe6b40"))(n = 200)
limit.bias <- 1


##### utility functions #####
themeFunction <- function(plot) {
  plot + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=15),
        legend.key.size = unit(2.5, 'cm'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(), # remove facet background
        #strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"), # only black frame
        strip.text = element_text(size = 15), # remove facet text
        axis.line = element_line(colour = "white"))
}

plotHeatmap <- function(data, xVar, yVar, xLabel = "", yLabel= "") {
  ggplot(data, 
         aes(x = xVar, y = yVar, fill = emmean)) + 
    geom_tile() +
    geom_text(aes(x = xVar, y = yVar, label = round(emmean, 2)), 
              color="black", size=rel(5)) +
    scale_fill_gradientn("relative R² error",
                         # colours = c(colRamp[length(colRamp):1], colRamp),
                         colours = c(colRamp[length(colRamp):1]),
                         values = scales::rescale(
                           # limited do +/- limit.bias
                           x = seq(from = -limit.bias,
                                   to = 0, length.out = 200),
                           from = c(-limit.bias, 0)),
                         limits = c(-limit.bias, 0)) +
    xlab(xLabel) + ylab(yLabel) + guides(fill = "none")
}

##### plot data #####
df_emDGPxModel$model <- factor(df_emDGPxModel$model, 
                               levels = c("ENETw", "RF", "GBM", "ENETwo"), 
                               labels = c("ENETint", "RF", "GBM", "ENETlin"))
df_emDGPxModel$dgp <- factor(df_emDGPxModel$dgp, 
                             levels = c("inter", "pwlinear", "nonlinear3"), 
                             labels = c("inter", "pw", "sw")) 

(plot4guide <- ggplot(df_emDGPxModel, 
                      aes(x = dgp, y = model, fill = emmean)) + 
    geom_tile() +
    geom_text(aes(x = dgp, y = model, label = round(emmean, 2)), 
              color="black", size=rel(5)) +
    scale_fill_gradientn("",
                         colours = colRamp[length(colRamp):1],
                         values = scales::rescale(
                           # limited do +/- limit.bias
                           x = seq(from = -limit.bias,
                                   to = 0, length.out = 200),
                           from = c(-limit.bias, 0)),
                         limits = c(-limit.bias, 0))) 

(plot4guide <- themeFunction(plot4guide))
# cowplot::get_legend(plot4guide)

(pDGPxModel <- plotHeatmap(df_emDGPxModel, df_emDGPxModel$dgp, df_emDGPxModel$model,
                           xLabel = "DGP", yLabel = "Model"))
(pDGPxModel <- themeFunction(pDGPxModel))


df_emLinInterxModelxDGP$model <- factor(df_emLinInterxModelxDGP$model, 
                               levels = c("ENETw", "RF", "GBM", "ENETwo"), 
                               labels = c("ENETint", "RF", "GBM", "ENETlin"))
df_emLinInterxModelxDGP$dgp <- factor(df_emLinInterxModelxDGP$dgp, 
                             levels = c("inter", "pwlinear", "nonlinear3"), 
                             labels = c("inter", "pw", "sw")) 
df_emLinInterxModelxDGP$lin_inter <- factor(df_emLinInterxModelxDGP$lin_inter,
                                            levels = c("0.8_0.2", "0.5_0.5", "0.2_0.8"),
                                            labels = c("80:20", "50:50", "20:80"))
(pLinInterxDGPxModel <- plotHeatmap(df_emLinInterxModelxDGP,
                                    df_emLinInterxModelxDGP$dgp, 
                                    interaction(df_emLinInterxModelxDGP$lin_inter, df_emLinInterxModelxDGP$model, sep = "x")))
(pLinInterxDGPxModel <- themeFunction(pLinInterxDGPxModel))


df_emNxR2$N <- factor(df_emNxR2$N, 
                      levels = c(100, 300, 1000))
(pNxR2 <- plotHeatmap(df_emNxR2, df_emNxR2$N, df_emNxR2$R2,  
                      xLabel = "", yLabel = "Simulated R²")) 
(pNxR2 <- themeFunction(pNxR2))


df_emNxPtrash$N <- factor(df_emNxPtrash$N, 
                         levels = c(100, 300, 1000))
(pNxPtrash <- plotHeatmap(df_emNxPtrash, df_emNxPtrash$N, df_emNxPtrash$pTrash, 
                         xLabel = "Sample Size", yLabel = "# Noise Variables")) 
(pNxPtrash <- themeFunction(pNxPtrash))


pColR <- cowplot::plot_grid(pNxR2, pNxPtrash,
                   labels = c("R² x N", "noise x N"), ncol = 1)
pColR_woLabel <- cowplot::plot_grid(pNxR2, pNxPtrash, 
                            labels = c("", ""), ncol = 1)
pEMMfull <- cowplot::plot_grid(pDGPxModel, pLinInterxDGPxModel, pColR, cowplot::get_legend(plot4guide), 
  labels = c("model x DGP", "[EC x model] x DGP", "", ""), ncol = 4, rel_widths = c(1.3, 1.8, 1, 0.8))

# ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_EMMfull.png"),
#                 plot = pEMMfull,
#                 width = 13.63,
#                 height = 8.14,
#                 units = "in")

legend <- cowplot::get_legend(plot4guide)
# pEMM <- cowplot::plot_grid(pDGPxModel, pColR, legend, 
#                            labels = c("model x DGP", "", "rel. R²\n error"), 
#                            # label_x = 0.1,
#                            # label_y = 0.95,
#                            ncol = 3, rel_widths = c(1, 1, 0.5))

pEMM_woLabel <- cowplot::plot_grid(pDGPxModel, pColR_woLabel, legend, 
                                   labels = c("", "", "rel. R²\n error"), 
                                   ncol = 3, rel_widths = c(1, 1, 0.5))


# ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_EMM.png"),
#                 plot = pEMM,
#                 width = 13.33,
#                 height = 6.38,
#                 units = "in")

ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_EMM_woLabel.png"),
                plot = pEMM_woLabel,
                width = 13.33,
                height = 6.38,
                units = "in")

# df_emNxModel$N <- factor(df_emNxModel$N, 
#                          levels = c(100, 300, 1000))
# df_emNxModel$model <- factor(df_emNxModel$model, 
#                              levels = c("ENETw", "RF", "GBM", "ENETwo"), 
#                              labels = c("ENETint", "RF", "GBM", "ENETlin"))
# (pNxModel <- plotHeatmap(df_emNxModel, df_emNxModel$N, df_emNxModel$model, 
#                          xLabel = "sample size")) 
# (pNxModel <- themeFunction(pNxModel))

################################################################################
# in ANOVA sample will be the ID; but all samples run from 1:100 in all simulated
#   conditions that represent independent samples (between factors); identical sample
#   names suggest within factor in ANOVA!
#   -> introduce additional ID variable to indicate independent samples
# independent observations for samples {1:100} in different simulated conditions
#   that all have the same sample numbers
rSquaredTest <- rSquaredTest %>%
  group_by(model, dgp) %>%
  mutate(ID = seq_len(n())) %>%
  ungroup()

# save(rSquaredTest,
#      file = paste0(depMeasureFolder, "/rSquaredData_eachSample.rda"))
################################################################################
# calculate ANOVA for each model and dgp separately
################################################################################
modVec <- unique(rSquaredTest$model)
dgpVec <- unique(rSquaredTest$dgp)

for (iDGP in dgpVec) {
  for (iModel in modVec) {
    # only between factor ANOVA (only within factor was model)
    # anovaObjectName <- paste0("anovaRes_", iDGP, "_", iModel)
    tmp <- aov_ez(id = "ID",
                  dv = "testR2_relE", 
                  data = rSquaredTest[rSquaredTest$model == iModel &
                                        rSquaredTest$dgp == iDGP,],
                  between = c("N" , "pTrash" , "R2" , "rel" , "lin_inter"))
    # assign(anovaObjectName, tmp)
    
    tmpEta2 <- eta_squared(
      tmp, # fitted model
      partial = FALSE, # not partial!
      generalized = TRUE, # generalized eta squared
      ci = 0.95,
      verbose = TRUE)
    
    # prepare merging of eta results
    tmpEta2[,"dgp"] <- iDGP
    tmpEta2[,"model"] <- iModel
    
    # sort generalized eta-squared results
    # which higher order interactions do we need to illustrate to report simulation results?
    etaObjectName <- paste0("eta2_", iDGP, "_", iModel)
    assign(etaObjectName, tmpEta2[order(tmpEta2$Eta2_generalized, decreasing = T),])
  
  }
}

##### eta2 table #####
# eta2res <- ls(pattern = "eta2") # check objects that match the pattern
rm(eta2Thresh, eta2modelTable, eta2TestR2, eta2TestR2.ordered, eta2res)
eta2List <- lapply(grep("eta2", ls(), value = TRUE), get)
eta2Table <- do.call(rbind, eta2List)

eta2Table <- tidyr::pivot_wider(eta2Table[,c("Parameter", "Eta2_generalized", "model", "dgp")], 
                                names_from = c(model, dgp),  
                                values_from = Eta2_generalized,
                                names_sep = "_")

# # do not remove variables below a certain threshold
eta2Thresh <- 0.1
eta2Table$thresh <- apply(eta2Table[,!(colnames(eta2Table) %in% "Parameter")], 1, function(r) any(r >= eta2Thresh))
eta2Table$sumEta2 <- apply(eta2Table[,!(colnames(eta2Table) %in% c("Parameter", "thresh"))], 1, sum)
eta2Table$M <- apply(eta2Table[,!(colnames(eta2Table) %in% c("Parameter", "thresh", "sumEta2"))], 1, mean)
eta2Table <- dplyr::arrange(eta2Table, desc(M))

colnames(eta2Table)
texTable <- eta2Table[eta2Table$thresh == TRUE,!(colnames(eta2Table) %in% c("thresh", "sumEta2", "M"))]
orderCols <- c("Parameter", "RF_inter", "GBM_inter", "ENETw_inter", "ENETwo_inter",
               "RF_pwlinear", "GBM_pwlinear", "ENETw_pwlinear", "ENETwo_pwlinear",
               "RF_nonlinear3", "GBM_nonlinear3", "ENETw_nonlinear3", "ENETwo_nonlinear3")
texTable <- texTable[, orderCols]
      
# print(xtable::xtable(texTable, type = "latex"),
#       file = paste0(plotFolder, "/ANOVAresults/betweenANOVA_relErrorR2_", eta2Thresh, ".tex"))

rm(eta2List, eta2Table)
