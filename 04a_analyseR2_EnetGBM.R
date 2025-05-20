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

# # save anova result object
# save(anovaTestR2, file = paste0(depMeasureFolder, "/mixedANOVA_aovResult_relErrorR2.rda"))
# load(paste0(depMeasureFolder, "/mixedANOVA_aovResult_relErrorR2.rda"))

eta2TestR2 <- eta_squared(
  anovaTestR2, # fitted model
  partial = FALSE, # not partial!
  generalized = TRUE, # generalized eta squared
  ci = 0.95,
  verbose = TRUE)

# # save eta² result object
# save(eta2TestR2, file = paste0(depMeasureFolder, "/mixedANOVA_eta2Result_relErrorR2.rda"))
load(paste0(depMeasureFolder, "/mixedANOVA_eta2Result_relErrorR2.rda"))

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

ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_relErrorR2_thresh", eta2Thresh, ".png"),
                plot = pEta2mixedFull,
                width = 17.52,
                height = 10.76,
                units = "in")

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
# rel   emmean        SE     df lower.CL upper.CL
# 0.6 -0.61855 0.0001368 485514 -0.61882 -0.61828
# 0.8 -0.47829 0.0001368 485514 -0.47855 -0.47802
# 1   -0.32881 0.0001368 485514 -0.32908 -0.32855

(emN <- emmeans(anovaTestR2, ~ N))
# N      emmean        SE     df lower.CL upper.CL
# 100  -0.62557 0.0001368 485514 -0.62584 -0.62530
# 300  -0.44708 0.0001368 485514 -0.44735 -0.44681
# 1000 -0.35300 0.0001368 485514 -0.35327 -0.35274

(emR2 <- emmeans(anovaTestR2, ~ R2))
# R2    emmean        SE     df lower.CL upper.CL
# 0.2 -0.61180 0.0001368 485514 -0.61207 -0.61153
# 0.5 -0.44535 0.0001368 485514 -0.44562 -0.44508
# 0.8 -0.36850 0.0001368 485514 -0.36877 -0.36824

(emLin_inter <- emmeans(anovaTestR2, ~ lin_inter))
# lin_inter   emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.58257 0.0001368 485514 -0.58284 -0.58230
# 0.5_0.5   -0.47358 0.0001368 485514 -0.47385 -0.47331
# 0.8_0.2   -0.36951 0.0001368 485514 -0.36977 -0.36924

##### interactions #####
(emDGPxModel <- emmeans(anovaTestR2, "model", by = "dgp"))
# dgp = inter:
#   model   emmean        SE     df lower.CL upper.CL
# ENETwo -0.6225 0.0001686 485514  -0.6228  -0.6221
# GBM    -0.5446 0.0001903 485514  -0.5449  -0.5442
# RF     -0.4862 0.0001503 485514  -0.4865  -0.4859
# ENETw  -0.3712 0.0002201 485514  -0.3716  -0.3708
# 
# dgp = nonlinear3:
#   model   emmean        SE     df lower.CL upper.CL
# GBM    -0.5089 0.0001903 485514  -0.5093  -0.5085
# ENETw  -0.5081 0.0002201 485514  -0.5085  -0.5077
# RF     -0.4629 0.0001503 485514  -0.4632  -0.4626
# ENETwo -0.4181 0.0001686 485514  -0.4184  -0.4177
#
# dgp = pwlinear:
#   model   emmean        SE     df lower.CL upper.CL
# GBM    -0.4875 0.0001903 485514  -0.4879  -0.4872
# ENETw  -0.4769 0.0002201 485514  -0.4773  -0.4765
# RF     -0.4350 0.0001503 485514  -0.4353  -0.4347
# ENETwo -0.3808 0.0001686 485514  -0.3811  -0.3805

(emLinInterxModelxDGP <- emmeans(anovaTestR2, "lin_inter", by = c("model", "dgp")))
# model = ENETw, dgp = inter:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.4089 0.0003813 485514  -0.4096  -0.4081
# 0.5_0.5   -0.3806 0.0003813 485514  -0.3813  -0.3798
# 0.8_0.2   -0.3242 0.0003813 485514  -0.3249  -0.3234
# 
# model = ENETwo, dgp = inter:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.8794 0.0002920 485514  -0.8800  -0.8788
# 0.5_0.5   -0.6229 0.0002920 485514  -0.6235  -0.6223
# 0.8_0.2   -0.3651 0.0002920 485514  -0.3657  -0.3645
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
# 0.2_0.8   -0.6312 0.0003813 485514  -0.6319  -0.6305
# 0.5_0.5   -0.5115 0.0003813 485514  -0.5123  -0.5108
# 0.8_0.2   -0.3816 0.0003813 485514  -0.3823  -0.3808
# 
# model = ENETwo, dgp = nonlinear3:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.5379 0.0002920 485514  -0.5385  -0.5373
# 0.5_0.5   -0.4099 0.0002920 485514  -0.4104  -0.4093
# 0.8_0.2   -0.3064 0.0002920 485514  -0.3070  -0.3059
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
# 0.2_0.8   -0.5803 0.0003813 485514  -0.5810  -0.5795
# 0.5_0.5   -0.4798 0.0003813 485514  -0.4806  -0.4791
# 0.8_0.2   -0.3705 0.0003813 485514  -0.3713  -0.3698
# 
# model = ENETwo, dgp = pwlinear:
#   lin_inter  emmean        SE     df lower.CL upper.CL
# 0.2_0.8   -0.4781 0.0002920 485514  -0.4787  -0.4775
# 0.5_0.5   -0.3726 0.0002920 485514  -0.3732  -0.3720
# 0.8_0.2   -0.2917 0.0002920 485514  -0.2922  -0.2911
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
# 100  -0.7956 0.0002369 485514  -0.7960  -0.7951
# 300  -0.6016 0.0002369 485514  -0.6021  -0.6011
# 1000 -0.4382 0.0002369 485514  -0.4387  -0.4378
# 
# R2 = 0.5:
#   N     emmean        SE     df lower.CL upper.CL
# 100  -0.6011 0.0002369 485514  -0.6015  -0.6006
# 300  -0.3965 0.0002369 485514  -0.3970  -0.3961
# 1000 -0.3384 0.0002369 485514  -0.3389  -0.3380
# 
# R2 = 0.8:
#   N     emmean        SE     df lower.CL upper.CL
# 100  -0.4801 0.0002369 485514  -0.4805  -0.4796
# 300  -0.3431 0.0002369 485514  -0.3436  -0.3426
# 1000 -0.2824 0.0002369 485514  -0.2828  -0.2819

(emNxModel <- emmeans(anovaTestR2, "N", by = "model"))
# model = ENETw:
#   N      emmean         SE     df lower.CL upper.CL
# 100  -0.64543 0.00022012 485514 -0.64586 -0.64500
# 300  -0.41001 0.00022012 485514 -0.41044 -0.40958
# 1000 -0.30076 0.00022012 485514 -0.30120 -0.30033
# 
# model = ENETwo:
#   N      emmean         SE     df lower.CL upper.CL
# 100  -0.57667 0.00016858 485514 -0.57700 -0.57633
# 300  -0.44312 0.00016858 485514 -0.44345 -0.44279
# 1000 -0.40153 0.00016858 485514 -0.40186 -0.40120
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
# 100  -0.5686 0.0001935 485514  -0.5690  -0.5682
# 300  -0.4442 0.0001935 485514  -0.4446  -0.4438
# 1000 -0.3558 0.0001935 485514  -0.3562  -0.3554
# 
# pTrash = 50:
#   N     emmean        SE     df lower.CL upper.CL
# 100  -0.6826 0.0001935 485514  -0.6829  -0.6822
# 300  -0.4500 0.0001935 485514  -0.4503  -0.4496
# 1000 -0.3502 0.0001935 485514  -0.3506  -0.3498

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

colRamp <- colorRampPalette(c('#55fa79','#FF9900', '#FF6666'))(n = 200)
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
                         colours = c(colRamp[length(colRamp):1], colRamp),
                         values = scales::rescale(
                           # limited do +/- limit.bias
                           x = seq(from = -limit.bias,
                                   to = limit.bias, length.out = 400),
                           from = c(-limit.bias, limit.bias)),
                         limits = c(-limit.bias, limit.bias)) +
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
                         colours = c(colRamp[length(colRamp):1], colRamp),
                         values = scales::rescale(
                           # limited do +/- limit.bias
                           x = seq(from = -limit.bias,
                                   to = limit.bias, length.out = 400),
                           from = c(-limit.bias, limit.bias)),
                         limits = c(-limit.bias, limit.bias))) 
(plot4guide <- themeFunction(plot4guide))
# cowplot::get_legend(plot4guide)

(pDGPxModel <- plotHeatmap(df_emDGPxModel, df_emDGPxModel$dgp, df_emDGPxModel$model))
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
                      xLabel = "", yLabel = "simulated R²")) 
(pNxR2 <- themeFunction(pNxR2))


df_emNxPtrash$N <- factor(df_emNxPtrash$N, 
                         levels = c(100, 300, 1000))
(pNxPtrash <- plotHeatmap(df_emNxPtrash, df_emNxPtrash$N, df_emNxPtrash$pTrash, 
                         xLabel = "sample size", yLabel = "# noise variables")) 
(pNxPtrash <- themeFunction(pNxPtrash))


pColR <- cowplot::plot_grid(pNxR2, pNxPtrash, 
                   labels = c("R² x N", "noise x N"), ncol = 1)
pEMMfull <- cowplot::plot_grid(pDGPxModel, pLinInterxDGPxModel, pColR, cowplot::get_legend(plot4guide), 
  labels = c("model x DGP", "[EC x model] x DGP", "", ""), ncol = 4, rel_widths = c(1.3, 1.8, 1, 0.8))

ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_EMMfull.png"),
                plot = pEMMfull,
                width = 13.63,
                height = 8.14,
                units = "in")

legend <- cowplot::get_legend(plot4guide)
pEMM <- cowplot::plot_grid(pDGPxModel, pColR, legend, 
                           labels = c("model x DGP", "", "rel. R²\n error"), ncol = 3, rel_widths = c(1, 1, 0.5))


ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/mixedANOVA_EMM.png"),
                plot = pEMM,
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
eta2List <- lapply(grep("eta2", ls(), value = TRUE), get)
eta2Table <- do.call(rbind, eta2List)

# eta2Thresh <- 0.1
eta2Table <- tidyr::pivot_wider(eta2Table[,c("Parameter", "Eta2_generalized", "model", "dgp")], 
                                names_from = c(model, dgp),  
                                values_from = Eta2_generalized,
                                names_sep = "_")

# # do not remove variables below a certain threshold
# eta2Table <- subset(eta2Table, eta2Table$ENETw >= eta2Thresh | 
#                       eta2Table$ENETwo >= eta2Thresh |
#                       #eta2Table$GBM >= eta2Thresh |
#                       eta2Table$RF >= eta2Thresh)


eta2Table$sumEta2 <- apply(eta2Table[,!(colnames(eta2Table) %in% "Parameter")], 1, sum)
eta2Table$M <- apply(eta2Table[,!(colnames(eta2Table) %in% c("Parameter", "sumEta2"))], 1, mean)
# eta2Table$sumEta2 <- apply(eta2Table[,c("ENETw", "ENETwo", "RF")], 1, sum)
# eta2Table$M <- apply(eta2Table[,c("ENETw", "ENETwo", "RF")], 1, mean)
eta2Table <- dplyr::arrange(eta2Table, desc(M))
# eta2Table <- eta2Table[,c("Parameter", "M", "RF", "GBM", "ENETw", "ENETwo")]
# eta2Table <- eta2Table[,c("Parameter", "M", "RF", "ENETw", "ENETwo")]

# print(xtable::xtable(eta2Table, type = "latex"),
#       file = paste0(plotFolder, "/ANOVAresults/betweenANOVA_R2.tex"))

# print(xtable::xtable(eta2Table, type = "latex"),
#       file = paste0(plotFolder, "/ANOVAresults/betweenANOVA_relErrorR2.tex"))

rm(eta2List, eta2Table)
##### plot eta2 #####
# check how many parameters are still in depending on threshold
# dim(eta2ENETw[eta2ENETw$Eta2_generalized >= eta2Thresh,])
# dim(eta2ENETwo[eta2ENETwo$Eta2_generalized >= eta2Thresh,])
# dim(eta2GBM[eta2GBM$Eta2_generalized >= eta2Thresh,])

# eta2 <- rbind(eta2ENETw[eta2ENETw$Eta2_generalized >= eta2Thresh,], 
#               eta2ENETwo[eta2ENETwo$Eta2_generalized >= eta2Thresh,],
#               eta2GBM[eta2GBM$Eta2_generalized >= eta2Thresh,])

eta2List <- lapply(grep("eta2", ls(), value = TRUE), get)
eta2 <- do.call(rbind, eta2List)

# sort variables according to total eta2 across all models
eta2Sums <- aggregate(Eta2_generalized ~ Parameter, data = eta2, sum)
sortIdx <- order(eta2Sums$Eta2_generalized, decreasing = T)
eta2$Parameter <- factor(eta2$Parameter, 
                         levels = eta2Sums$Parameter[sortIdx])

eta2$model <- factor(eta2$model, 
                     levels = c("RF", "GBM", "ENETw", "ENETwo")) # obvious order of models

eta2$dgp <- factor(eta2$dgp, 
                     levels = c("inter", "pwlinear", "nonlinear3", "nonlinear")) # obvious order of models

# as barplot with cut-off: gen. eta² = .1 
eta2Thresh <- 0.1
(pEta2Model <- ggplot(eta2[eta2$Eta2_generalized >= eta2Thresh,], 
# (pEta2Model <- ggplot(eta2, 
                      aes(x = Parameter, y = Eta2_generalized,
                          group = model, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    geom_text(aes(label=round(Eta2_generalized, 2), group = model), 
              angle = 90, hjust = 1.5, vjust=0.5, 
              #angle = 0, vjust=1.5, 
              position = position_dodge(width = .9), 
              color="black", size=3.5)+
    ylab("generalisiertes eta^2") +
    facet_wrap(~ dgp, ncol = 1) +
    scale_fill_manual(values = c("#FFA500", "#990000", "#006600", "#009999")) +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

(pEta2dgp <- ggplot(eta2[eta2$Eta2_generalized >= eta2Thresh,], 
                      # (pEta2Model <- ggplot(eta2, 
                      aes(x = Parameter, y = Eta2_generalized,
                          group = dgp, fill = dgp)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    geom_text(aes(label=round(Eta2_generalized, 2), group = dgp), 
              angle = 90, hjust = 1.5, vjust=0.5, 
              #angle = 0, vjust=1.5, 
              position = position_dodge(width = .9), 
              color="black", size=3.5)+
    ylab("generalisiertes eta^2") +
    facet_wrap(~ model, ncol = 1) +
    scale_fill_manual(values = c("#FFA500", "#990000", "#006600", "#009999")) +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/betweenANOVA_R2_Model.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/betweenANOVA_relErrorR2_Model.png"),
                plot = pEta2Model,
                width = 13.63,
                height = 23.89,
                units = "in")
ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/betweenANOVA_R2_DGP.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/betweenANOVA_relErrorR2_DGP.png"),
                plot = pEta2dgp,
                width = 13.63,
                height = 23.89,
                units = "in")

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
                                 '#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd', 
                                 "#ccebc5", "#ffed6f", "#f4cae4")) +
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# ggplot2::ggsave(filename = paste0(plotFolder, "/ANOVAresults/betweenANOVA_R2stacked.png"),
#                 plot = pEta2_stacked,
#                 width = 13.63,
#                 height = 12.07,
#                 units = "in")


