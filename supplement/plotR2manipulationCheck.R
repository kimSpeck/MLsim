# plot check big data sample models
# linear vs. inter/nonlinear/piecewise linear
# different reliability measures

# "trueR2" = simulated R2 for all variables
# "R2_full" = fit data generating model 
# "R2lm_full" = fit data with original variables (no dummy variables or segmented process)
# "R2oracle_full" = fit only variables with simulated effect (only for pwlinear) 
# "trueLin" = simulated R2 only linear effects
# "R2_lin" = R2 when fitting only variables with linear effects
# "trueLinPC"    
# "trueOther" = simulated R2 only other effects     
# "R2_other" = R2 when fitting only variables with other effects (only interaction, only dummies, etc.)      
# "R2lm_other" = R2 when fitting only original variables with other effects 
# "trueOtherPC"   
# "pTrash" = 10       
# "rel" = 0.6, 0.8, 1          
# "data" = inter, nonlinear, pwlinear 
 
# originally: check data generating process!
# check effect of reliability on the linear models 
#   (separately from small sample error)

library(stringr)
library(ggplot2)
library(ggh4x)

# plotting colors across final paper plots
colValuesR2 <- c('#db4a07', '#850c0c', '#3c1518')
colValuesInter <- c('#050440', '#181ff2', '#0eb2e8')
colValuesLin <- c('#0eb2e8', '#181ff2', '#050440')

load("results/checkR2manipulation.rda")

condNames <- names(checkR2)
checkR2 <- do.call(rbind, checkR2)
checkR2 <- data.frame(cbind(checkR2, condNames = rep(condNames, each = 9)))

pattern <- "pTrash(\\d+)_rel(\\d+\\.?\\d*)_(inter|nonlinear|pwlinear)\\.rda$"

matches <- stringr::str_match(checkR2$condNames, pattern)

checkR2 <- cbind(checkR2, 
                 pTrash = matches[, 2],
                 rel = matches[, 3], 
                 data = matches[, 4])
checkR2$condNames <- NULL

checkR2[checkR2$data == "inter", "R2lm_full"] <- checkR2[checkR2$data == "inter", "R2_full"]
checkR2[checkR2$data == "inter", "R2lm_other"] <- checkR2[checkR2$data == "inter", "R2_other"]

str(checkR2)
chr2num <- c("trueR2", "R2_full", "R2lm_full", "R2oracle_full", "trueLin", "R2_lin", 
             "trueOther", "R2_other", "R2lm_other")
checkR2[, chr2num] <- sapply(checkR2[, chr2num], as.numeric)
checkR2$lin_inter <- paste0(checkR2$trueLinPC, "_", checkR2$trueOtherPC) 

checkR2$lin_inter <- factor(checkR2$lin_inter, 
                            levels = c("0.8_0.2", "0.5_0.5", "0.2_0.8"),
                            labels = c("80:20", "50:50", "20:80"))
  
plotCheckFunction <- function(data, y) {
  ggplot(data, aes(x = factor(rel), y = y, group = interaction(trueR2, data),
                    color = factor(trueR2), linetype = factor(data))) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(name = "data", values = c("solid", "dotted", "dashed")) +
  scale_color_manual(name = "R²", values = colValuesR2) +
  facet_grid2(~ lin_inter, 
              strip = strip_themed(
                background_x = list(element_rect(fill = alpha(colValuesLin[3], 0.4)),
                                    element_rect(fill = alpha(colValuesLin[2], 0.4)),
                                    element_rect(fill = alpha(colValuesLin[1], 0.4))))) + 
  geom_hline(yintercept = 0.2, col = colValuesR2[1],
             alpha = 0.6, linewidth = 1) +
  geom_hline(yintercept = 0.5, col = colValuesR2[2],
             alpha = 0.6, linewidth = 1) +
  geom_hline(yintercept = 0.8, col = colValuesR2[3],
             alpha = 0.6, linewidth = 1) +
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
        legend.position = "bottom", 
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.box = "horizontal")
}

# check Full R2
# here, I used the actual variables from the data generating process in the regression
#   - R2 matches trueR2 if reliability was 1
#   - irrespective of sample size, measurement error reduces performance dramatically
#   - for linear and piecewise linear DGPs it is identical (classic attenuation bias)
plotCheckFunction(checkR2, checkR2$R2_full)

# actual variables from the data generating process in the regression
# but only interaction, nonlinear and pwlinear effects
# --> check: R2 for rel = 1 is the same for every dgp!
plotCheckFunction(checkR2, checkR2$R2_other)

# here, I used the originally sampled continuous variables in the regression
#   - R2 matches trueR2 if reliability was 1 and we simulated linear vs. inter data
#   - for pwlinear data, the linear model still finds a cerain proportion of R2
#   - for nonlinear data, the linear model struggles to explain variance
# --> Warum ist performance für nonlinear data schlechter, wenn die linearen Zsmh. mehr bekommen?
plotCheckFunction(checkR2, checkR2$R2lm_full)

# only originally sampled continuous variables with simulated linear effects in the regression
#  - keine Unterschiede zwischen verschiedenen DGPs!
plotCheckFunction(checkR2, checkR2$R2_lin)

# only originally sampled continuous variables with other simulated effects in the regression
plotCheckFunction(checkR2, checkR2$R2lm_other)

