# plot DGP illustration
library(ggplot2)
library(rockchalk)
library(patchwork)

source("utils/simTools.R")

# path to plot folder
plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

################################################################################
##### plot utility functions #####
################################################################################

themeFunction <- function(plot) {
  plot +  
    xlim(c(-5, 5)) +
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
          legend.position = c(.15, .85), 
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.box = "horizontal")
}

################################################################################
##### illustrate linear effect #####
################################################################################
# load nonlinear data
# ... without measurement error
load("data/bigTestSamples/simDataN1e+06_pTrash10_rel1_nonlinear.rda")

# save subset for linear effect
df_lm <- data.frame(x = dataList[["X_int"]][1:1000, "Var1"], 
                 y = dataList[["yMat"]][1:1000, "R20.8lin_inter0.8_0.2"])


plotLM <- function(data, x, y) {
  ggplot(data, aes(x = x, y = y)) +
    geom_point(alpha = 0.4, colour = "blue", size = 2.5) +
    geom_smooth(data = data[data$x <= 0, ], 
                method = "lm", colour = "black") + 
    geom_smooth(data = data[data$x > 0, ], 
                method = "lm", colour = "black")
}
pLM_rel1 <- plotLM(df_lm, x, y)
pLM_rel1 <- themeFunction(pLM_rel1)

# load nonlinear data ... with measurement error
load("data/bigTestSamples/simDataN1e+06_pTrash10_rel0.6_nonlinear.rda")
df_lm_rel0.6 <- data.frame(x = dataList[["X_int"]][1:1000, "Var1"], 
                           y = dataList[["yMat"]][1:1000, "R20.8lin_inter0.8_0.2"])

pLM_rel0.6 <- plotLM(df_lm_rel0.6, x, y)
pLM_rel0.6 <- themeFunction(pLM_rel0.6)

################################################################################
##### illustrate non-linear effect #####
################################################################################
# ... without measurement error
load("data/bigTestSamples/simDataN1e+06_pTrash10_rel1_nonlinear.rda")

# save subset for nonlinear effect 
df_nl <- data.frame(x5 = dataList[["X_int"]][1:1000, "Var5"], 
                    x6 = dataList[["X_int"]][1:1000, "Var6"], 
                    y = dataList[["yMat"]][1:1000, "R20.8lin_inter0.2_0.8"])

df_nl$d5 <- createDummy(df_nl$x5, q = 0.5, effectCoding = T)
df_nl$d6 <- createDummy(df_nl$x6, q = 0.5, effectCoding = T)
df_nl$group <- paste0(df_nl$d5, "x", df_nl$d6)

plotNL <- function(data) {
  ggplot(data) +
    geom_point(aes(x = x5, y = y, group = group, colour = x6, shape = group), 
               alpha = 0.6, size = 2.5) +
    scale_shape_manual(values = c(1, 16, 2, 17)) +
    scale_colour_gradient2(
      low      = "#F44336",  # red
      mid      = "#9C27B0",  # purple
      high     = "#2196F3",  # blue
      midpoint = 0,          # numeric center of z
      space    = "Lab"
    ) +
    # scale_colour_continuous(type = "viridis") + # colorblind version
    geom_smooth(aes(x = d5, y = y, group = d6), method = "lm", colour = "black")
}

pNL_rel1 <- plotNL(df_nl)
pNL_rel1 <- themeFunction(pNL_rel1)

# ... with measurement error
load("data/bigTestSamples/simDataN1e+06_pTrash10_rel0.6_nonlinear.rda")
df_nl_rel0.6 <- data.frame(x5 = dataList[["X_int"]][1:1000, "Var5"], 
                    x6 = dataList[["X_int"]][1:1000, "Var6"], 
                    y = dataList[["yMat"]][1:1000, "R20.8lin_inter0.2_0.8"])

df_nl_rel0.6$d5 <- createDummy(df_nl_rel0.6$x5, q = 0.5, effectCoding = T)
df_nl_rel0.6$d6 <- createDummy(df_nl_rel0.6$x6, q = 0.5, effectCoding = T)
df_nl_rel0.6$group <- paste0(df_nl_rel0.6$d5, "x", df_nl_rel0.6$d6)

pNL_rel0.6 <- plotNL(df_nl_rel0.6)
pNL_rel0.6 <- themeFunction(pNL_rel0.6)

################################################################################
##### illustrate piecewise linear effect #####
################################################################################
# load piecewise linear data
# ... without measurement error
load("data/bigTestSamples/simDataN1e+06_pTrash10_rel1_pwlinear.rda")

df_pwl <- data.frame(x = dataList[["X_int"]][1:1000, "Var5"], 
                 y = dataList[["yMat"]][1:1000, "R20.8lin_inter0.2_0.8"])

plotPWL <- function(data) {
  ggplot(data, aes(x = x, y = y)) +
    geom_point(alpha = 0.4, size = 2.5, colour = "blue") +
    geom_smooth(data = data[data$x <= 0, ], 
                method = "lm", colour = "black") + 
    geom_smooth(data = data[data$x > 0, ], 
                method = "lm", colour = "black")
}

pPWL_rel1 <- plotPWL(df_pwl)
pPWL_rel1 <- themeFunction(pPWL_rel1)

# ... with measurement error
load("data/bigTestSamples/simDataN1e+06_pTrash10_rel0.6_pwlinear.rda")

df_pwl_rel0.6 <- data.frame(x = dataList[["X_int"]][1:1000, "Var5"], 
                            y = dataList[["yMat"]][1:1000, "R20.8lin_inter0.2_0.8"])

pPWL_rel0.6 <- plotPWL(df_pwl_rel0.6)
pPWL_rel0.6 <- themeFunction(pPWL_rel0.6)

################################################################################
##### illustrate interaction effects #####
################################################################################
# load linear data
# ... without measurement error
load("data/bigTestSamples/simDataN1e+06_pTrash10_rel1_inter.rda")

df_inter <- data.frame(x1 = dataList[["X_int"]][1:1000, "Var1"], 
                     x2 = dataList[["X_int"]][1:1000, "Var2"], 
                     y = dataList[["yMat"]][1:1000, "R20.8lin_inter0.2_0.8"])


plotInter <- function(data) {
  fit <- lm(y ~ x1 + x2 + x1:x2, data = data)
  
  b0 <- coef(fit)["(Intercept)"]              # Konstante
  b2.m <- coef(fit)["x2"]             # Reg.koeffizient des Moderators 
  b1 <- coef(fit)["x1"]
  b3 <- coef(fit)["x1:x2"]
  
  ggplot() +   
    geom_point(data = data[data$x2 <= -1, ], 
               aes(x = x1, y = y), alpha = 0.3, size = 2.5, colour = "#F44336") +                                 
    geom_point(data = data[data$x2 > -1 & data$x2 < 1, ], 
               aes(x = x1, y = y), alpha = 0.3, size = 2.5, colour = "#9C27B0") +                                 
    geom_point(data = data[data$x2 >= 1, ], 
               aes(x = x1, y = y), alpha = 0.3, size = 2.5, colour = "#2196F3") +
    geom_abline(aes(intercept = (b0 + b2.m * (0-1)),
                    slope = (b1 + b3 * (0-1)), colour = "#F44336"), linewidth = 1.5) + # MW - SD
    geom_abline(aes(intercept = (b0 + b2.m * 0),
                    slope = (b1 + b3 * 0), colour = "#9C27B0"), linewidth = 1.5) + # MW
    geom_abline(aes(intercept = (b0 + b2.m * (0+1)),
                    slope = (b1 + b3 * (0+1)), colour = "#2196F3"), linewidth = 1.5) +  # MW + SD
    scale_color_identity("",
                         breaks = c("#F44336", "#9C27B0", "#2196F3"),
                         labels = c("x2: low (M - SD)", "x2: mean (M)", "x2: high (M + SD)"),
                         guide = "legend")
}

pInter_rel1 <- plotInter(df_inter)
pInter_rel1 <- themeFunction(pInter_rel1)
  
# ... with measurement error
load("data/bigTestSamples/simDataN1e+06_pTrash10_rel0.6_inter.rda")

df_inter_rel0.6 <- data.frame(x1 = dataList[["X_int"]][1:1000, "Var1"], 
                       x2 = dataList[["X_int"]][1:1000, "Var2"], 
                       y = dataList[["yMat"]][1:1000, "R20.8lin_inter0.2_0.8"])

pInter_rel0.6 <- plotInter(df_inter_rel0.6)
pInter_rel0.6 <- themeFunction(pInter_rel0.6)

################################################################################
# in rows: simulated predictors without (top row) and with measurement error (bottom row)
# in columns: 
#     ... linear effects (without simulated interaction on the variable)
#     ... interaction effects shown as conditional regressions
#     ... nonlinear effects from the 2x2 factorial design with interaction between both dummies
#     ... piecewise linear effect 
# the data is simulated with R² = 0.8 and for the depicted effect (linear or other)
#   the condition with 80% of R² on the respectice effect is illustrated

pFull <- cowplot::plot_grid(
  pLM_rel1, pInter_rel1, pNL_rel1, pPWL_rel1,
  pLM_rel0.6, pInter_rel0.6, pNL_rel0.6, pPWL_rel0.6,
  labels = "AUTO", ncol = 4)

ggplot2::ggsave(filename = paste0(plotFolder, "/DGPoverview.png"),
                plot = pFull,
                width = 24.24,
                height = 13.28,
                units = "in")


