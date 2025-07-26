library(rpart)
library(partykit)
library(ggplot2)   # ← first
library(ggparty) 

# plot DGP illustration
library(rockchalk)
library(patchwork)

# for the tree plot
library(rpart.plot)
library(grid)
library(ggplotify)

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
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.position = c(.08, .92), 
          legend.title = element_text(size = 25),
          # legend.text = element_text(size = 20),
          legend.text = element_blank(),
          legend.direction = "horizontal",
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

##### line plot to only schematically illustrate DGP #####
# Simulate data (linear relationship)
set.seed(123)        # for reproducibility
n <- 50
x <- seq(0, 10, length.out = n)
beta0 <- 2           # intercept
beta1 <- 1           # slope
# Generate "true" y plus random noise
y <- beta0 + beta1 * x + rnorm(n, mean = 0, sd = 1)

# Store data in a data frame
df_linear <- data.frame(x = x, y = y)

# Plot
(linear <- ggplot(df_linear, aes(x = x, y = y)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black") +
  labs(x = "X", y = "Y") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "white", fill = "white"),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()))

################################################################################
##### illustrate non-linear effect (3 dummies, no interaction) #####
################################################################################
# ... without measurement error
load("data/bigTestSamples/simDataN1e+06_pTrash10_rel1_nonlinear3.rda")

# save subset for nonlinear effect 
df_nl3 <- data.frame(x5 = dataList[["X_int"]][1:1000, "Var5"], 
                    x6 = dataList[["X_int"]][1:1000, "Var6"], 
                    y = dataList[["yMat"]][1:1000, "R20.8lin_inter0.2_0.8"])

df_nl3$d5 <- createDummy(df_nl3$x5, q = 0.5, effectCoding = T)
df_nl3$d6 <- createDummy(df_nl3$x6, q = 0.5, effectCoding = T)
df_nl3$group <- paste0(df_nl3$d5, "x", df_nl3$d6)

plotNL <- function(data) {
  ggplot(data) +
    geom_point(aes(x = x5, y = y, group = group, colour = x6, shape = group), 
               alpha = 0.6, size = 2.5) +
    scale_shape_manual(values = c(1, 17, 16, 2)) +
    scale_colour_gradient2(low = "#fa9e1e", mid = "#9C27B0", high = "#25f7f0", 
                           midpoint = 0, space    = "Lab") +
    geom_point(aes(x = d5, y = y, group = group, shape = group), 
               alpha = 1, size = 5, stroke = 2, color = "black", stat = 'summary', fun = 'mean') +
    guides(shape = "none")
}

pNL3_rel1 <- plotNL(df_nl3)
pNL3_rel1 <- themeFunction(pNL3_rel1)

# ... with measurement error
load("data/bigTestSamples/simDataN1e+06_pTrash10_rel0.6_nonlinear3.rda")
df_nl3_rel0.6 <- data.frame(x5 = dataList[["X_int"]][1:1000, "Var5"], 
                           x6 = dataList[["X_int"]][1:1000, "Var6"], 
                           y = dataList[["yMat"]][1:1000, "R20.8lin_inter0.2_0.8"])

df_nl3_rel0.6$d5 <- createDummy(df_nl3_rel0.6$x5, q = 0.5, effectCoding = T)
df_nl3_rel0.6$d6 <- createDummy(df_nl3_rel0.6$x6, q = 0.5, effectCoding = T)
df_nl3_rel0.6$group <- paste0(df_nl3_rel0.6$d5, "x", df_nl3_rel0.6$d6)

pNL3_rel0.6 <- plotNL(df_nl3_rel0.6)
pNL3_rel0.6 <- themeFunction(pNL3_rel0.6) # + guides(color = "none")

##### line plot to only schematically illustrate DGP #####
# Define model parameters (no interaction): Y = β₀ + β₁*D1 + β₂*D2
beta0 <- 2
beta1 <- 1
beta2 <- 3

# Create all combinations of the two dummy variables
df_dum3 <- expand.grid(X5 = c(-1, 1), X6 = c(-1, 1))

# Calculate Y for each (D1, D2) pair
df_dum3$Y <- beta0 + beta1 * df_dum3$X5 + beta2 * df_dum3$X6

# Convert D1 and D2 to factors (for cleaner plotting)
df_dum3$X5 <- factor(df_dum3$X5, levels = c(-1, 1), labels = c("< 0", "> 0"))
df_dum3$X6 <- factor(df_dum3$X6, levels = c(-1, 1), labels = c("< 0", "> 0"))
# df_dum3$Y <- factor(df_dum3$Y, levels = c(-2, 0, 4, 6), 
#                     labels = c("E(-1,-1)", "E(-1,1)", "E(1,-1)", "E(1,1)"))

# Fit a small regression tree to illustrate the splitting logic
tree_fit <- rpart(
  Y ~ X5 + X6,
  data = df_dum3,
  method = "anova",
  control = rpart.control(
    minsplit = 2,     # minimum number of obs needed to attempt a split
    minbucket = 1,    # minimum number of obs in any terminal node
    cp = 0,           # cost-complexity parameter (lower means more splits)
    xval = 0          # turn off cross-validation
  )
)

# convert rpart object to party object for plotting
tree_party <- as.party(tree_fit)

terminal_ids <- nodeids(tree_party, terminal = TRUE)  # get terminal node IDs

# create node labels
leafLabels   <- c(
  "E(y|X5<0,X6<0)",
  "E(y|X5>0,X6<0)",
  "E(y|X5<0,X6>0)",
  "E(y|X5>0,X6>0)"
)

# create a named vector for leaf labels
leaf_map <- setNames(leafLabels, terminal_ids)

# plot the tree with custom labels

nonlinear3 <- ggparty(
  tree_party,
  add_vars = list(
    leaf_lab = function(data, node) {
      leaf_map[as.character(data$id)]
    }
  ),
  terminal_space = 0.1) + 
  ggparty::geom_edge() +
  ggparty::geom_edge_label(size = 6) +            # größere Split-Labels
  ggparty::geom_node_label(
    aes(label = splitvar),
    ids = "inner",
    size = 6                                       # größere Variablennamen
  ) +
  ggparty::geom_node_label(
    aes(label = leaf_lab),
    ids = "terminal",
    parse = TRUE,
    fill       = "white",
    colour     = "black",
    label.size = 0,
    size       = 5.5                                 # größere Blatt-Labels
  ) +
  coord_fixed(ratio = 0.55, clip = "off") +
  # theme_void() +
  theme(
    panel.background = element_rect(fill = "white", colour = "white", size = 0.5),
    plot.background  = element_rect(fill = "white", colour = "white"),
    plot.title   = element_blank(),
    plot.margin  = unit(c(0.5, 2, 0.5, 2), "cm")    # top, right, bottom, left
  )

# check plot
nonlinear3

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

##### line plot to only schematically illustrate DGP #####
# Set seed for reproducibility (optional)
set.seed(123)

# Simulate data
n <- 100                             # Number of points
x <- seq(-10, 10, length.out = n)      # X-values
cutpoint <- 0                       # Where the piecewise function changes
beta0 <- 12                           # Initial constant (intercept in first segment)
beta1 <- 2                           # Slope in second segment

# Piecewise linear function:
# If x < cutpoint: y = beta0
# If x >= cutpoint: y = beta0 + beta1*(x - cutpoint)
y <- ifelse(x < cutpoint,
            beta0,
            beta0 + beta1 * (x - cutpoint))

# Combine into a data frame
df_piecewise <- data.frame(x = x, y = y)

# Plot the piecewise linear function
(pwlinear <- ggplot(df_piecewise, aes(x = x, y = y)) +
    geom_line(color = "black", size = 1) +
    geom_vline(xintercept = cutpoint, linetype = "dashed", color = "red") +
    labs(x = "X", y = "Y") + 
    ylim(0, 42) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.x = element_blank(),
          strip.text.y = element_blank()))

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
                         # breaks = c("#F44336", "#9C27B0", "#2196F3"),
                         # labels = c("X2: low (M - SD)", "X2: mean (M)", "X2: high (M + SD)"),
                         breaks = c("#2196F3", "#9C27B0", "#F44336"),
                         labels = c("X2: high (M + SD)", "X2: mean (M)", "X2: low (M - SD)"),
                         guide = "legend") +
    guides(color = "none")
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

##### line plot to only schematically illustrate DGP #####
# Define model parameters: Y = β0 + β1*X + β2*Z + β3*(X*Z)
beta0 <- 0
beta1 <- 1
beta2 <- 1
beta3 <- 1

# Create a sequence for X
X1 <- seq(-10, 10, length.out = 100)

# Choose 3 distinct Z levels to illustrate different lines
X2 <- c(2, 6, 10)

# Build a data frame by computing Y for each combination (X, Z-level)
df_list <- lapply(X2, function(iX2) {
  data.frame(
    X1 = X1,
    X2 = paste("X2 =", iX2),  # label for plotting
    Y = beta0 + beta1*X1 + beta2*iX2 + beta3*(X1*iX2)
  )
})

df_combined <- do.call(rbind, df_list)
df_combined$X2 <- factor(df_combined$X2, levels = c("X2 = 10", "X2 = 6", "X2 = 2"))
# Plot multiple lines (one for each Z level)
(inter <- ggplot(df_combined, aes(x = X1, y = Y, color = X2)) +
  geom_line(size = 1) +
  labs(x = "X1", y = "Y") +
  scale_color_manual("",
                     values = c("#2196F3", "#9C27B0", "#F44336"),
                     labels = c("X2 = 10" = "X2: high (M + SD)",
                                "X2 = 6" = "X2: mean (M)",
                                "X2 = 2" = "X2: low (M - SD)"),
                       guide = "legend") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "white", fill = "white"),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position = c(.15, .85), 
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.box = "horizontal"))





################################################################################
# patchwork plot (with cowplot)
################################################################################
# in rows: simulated predictors without (top row) and with measurement error (bottom row)
# in columns: 
#     ... linear effects (without simulated interaction on the variable)
#     ... interaction effects shown as conditional regressions
#     ... nonlinear effects from the 2x2 factorial design with interaction between both dummies
#     ... piecewise linear effect 
# the data is simulated with R² = 0.8 and for the depicted effect (linear or other)
#   the condition with 80% of R² on the respectice effect is illustrated


(pFull <- cowplot::plot_grid(
  linear, inter, pwlinear, nonlinear3,
  pLM_rel1, pInter_rel1, pPWL_rel1, pNL3_rel1, 
  pLM_rel0.6, pInter_rel0.6, pPWL_rel0.6, pNL3_rel0.6,
  labels = c("Linear", "Interaction", "Piecewise", "Stepwise", 
             "rel = 1", "",  "", "",
             "rel = 0.6", "", "", ""), ncol = 4, 
  label_x = 0.5, hjust = 0.5, label_size = 20))

ggplot2::ggsave(filename = paste0(plotFolder, "/DGPoverview.png"),
                plot = pFull,
                width = 28.24,
                height = 13.28,
                units = "in",
                bg    = "white")

################################################################################

