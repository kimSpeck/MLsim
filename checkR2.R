# check Rsquared in all data samples
# by simulating and adding measurement error to the predictor variables Rsquared increases
# presumably as a result of var(Rhat) = var(X %*% beta) increasing while var(y)
#     stays constant
# compare Rsquared with and without measurement error by plotting
#     ... bias between Rsquared with and without measurement error
#     ... bias between Rsquared without measurement error and target Rsquared
#     ... bias between Rsquared with measurement error and target Rsquared
library(ggplot2)

source("setParameters.R") # parameter values
source("simTools.R")

dataFolder <- "data"

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

condGrid <- expand.grid(N = setParam$dgp$N, 
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability)

# extract saved R2 from data and calculate R2 with ME from data
checkR2 <- lapply(seq_len(nrow(condGrid)), function(iSim) {
  
  # # test it
  # iSim <- 1
  
  # load data
  fileName <- paste0("simDataN", condGrid[iSim, "N"],
                     "_pTrash", condGrid[iSim, "pTrash"],
                     "_rel", condGrid[iSim,"reliability"], ".rda")
  load(paste0(dataFolder, "/", fileName))
  
  # get R^2 for X without measurement error
  R2 <- do.call(rbind, lapply(seq_along(data), function(subList) {
    data[[subList]][["R2"]]
  }))
  R2 <- cbind(R2, iSim)
  colnames(R2) <- c(setParam$dgp$condLabels, "idxCondGrid")
  
  # calculate R^2 for X with measurement error
  R2_wME <- do.call(rbind, lapply(seq_along(data), function(subList) {
    bMatrix <- data[[subList]][["trueB"]]
    sapply(seq_len(ncol(bMatrix)), function(x) {
      getR2(data[[subList]][["X_int"]], 
            bMatrix[,x], 
            setParam$dgp$sigmaE)
    })
  }))
  R2_wME <- cbind(R2_wME, iSim)
  colnames(R2_wME) <- c(setParam$dgp$condLabels, "idxCondGrid")
  
  list(R2 = R2, 
       R2_wME = R2_wME)
}) 

# extract data
R2 <- data.frame(do.call(rbind, lapply(seq_along(checkR2), function(subList) {
  checkR2[[subList]][["R2"]]
})))

R2_wME <- data.frame(do.call(rbind, lapply(seq_along(checkR2), function(subList) {
  checkR2[[subList]][["R2_wME"]]
})))

R2$type <- rep("R2", times = dim(R2)[1])
R2$sample <- rep(seq_len(setParam$dgp$nSamples), 
                 times = length(setParam$dgp$N) * length(setParam$dgp$pTrash) * 
                   length(setParam$dgp$reliability))
R2_wME$type <- rep("R2_wME", times = dim(R2_wME)[1])
R2_wME$sample <- rep(seq_len(setParam$dgp$nSamples),
                     times = length(setParam$dgp$N) * length(setParam$dgp$pTrash) * 
                       length(setParam$dgp$reliability))

R2 <- tidyr::pivot_longer(R2, !c(idxCondGrid, type, sample), 
                          names_to = "condLabel", values_to = "values")

R2_wME <- tidyr::pivot_longer(R2_wME, !c(idxCondGrid, type, sample), 
                              names_to = "condLabel", values_to = "values")
plotCheckR2 <- rbind(R2, R2_wME)


gridFull <- expand.grid(N = setParam$dgp$N,
                        pTrash = setParam$dgp$pTrash,
                        reliability = setParam$dgp$reliability)
gridFull <- tidyr::unite(gridFull, "condLabel", c(N, pTrash, reliability), 
                         sep = "_", remove = F)

idx <- match(plotCheckR2$idxCondGrid, row.names(gridFull))
plotCheckR2$condGrid <- gridFull$condLabel[idx]

plotCheckR2 <- tidyr::separate(plotCheckR2, condGrid, c("N", "pTrash", "rel"), sep = "_")
plotCheckR2$N <- factor(plotCheckR2$N, levels = setParam$dgp$N)
plotCheckR2$pTrash <- factor(plotCheckR2$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))

plotCheckR2$trueR2 <- stringr::str_sub(plotCheckR2$condLabel, start = 3L, end = 5L)
plotCheckR2$trueR2 <- as.numeric(plotCheckR2$trueR2)
plotCheckR2$lin_inter <- stringr::str_sub(plotCheckR2$condLabel, start = -7L, end = -1L)

plotCheckR2 <- tidyr::pivot_wider(plotCheckR2, names_from = "type", values_from = "values")
plotCheckR2$biasME <- plotCheckR2$R2_wME - plotCheckR2$R2 
plotCheckR2$biasTrue <- plotCheckR2$R2 - plotCheckR2$trueR2 
plotCheckR2$biasTrueME <- plotCheckR2$R2_wME - plotCheckR2$trueR2 

# aggregate data across all samples
plotBiasR2 <- aggregate(cbind(biasME, biasTrue, biasTrueME) ~ N + pTrash + rel + trueR2 + lin_inter, data = plotCheckR2, 
                        FUN = mean, na.rm = TRUE)

plotBiasR2$trueR2 <- factor(plotBiasR2$trueR2)

colValues <- c("green3", "darkcyan", "darkblue", "darkmagenta")
(pBiasME <- ggplot(plotBiasR2,
                   aes(x = interaction(pTrash, N, sep = " x "), y = biasME, 
                       group = trueR2, colour = trueR2)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(rel ~ lin_inter) +
    ylab(paste0("bias in R2")) +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle("R2 with ME - R2 without ME") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# # save plots
# ggplot2::ggsave(filename = paste0(plotFolder, "/", "plotBiasR2_ME", ".eps"),
#                 plot = pBiasME,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 14.39,
#                 height = 10.83,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/", "plotBiasR2_ME", ".png"),
#                 plot = pBiasME,
#                 width = 14.39,
#                 height = 10.83,
#                 units = "in")

# essentially almost no deviation from target R2 on average
# note the scale!
(pBiasTrue <- ggplot(plotBiasR2,
                     aes(x = interaction(pTrash, N, sep = " x "), y = biasTrue, 
                         group = trueR2, colour = trueR2)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(rel ~ lin_inter) +
    ylab(paste0("bias in R2")) +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle("R2 without ME - target R2") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# this is bias due to measurement error + bias due to deviation from target R2
# as deviation from target R2 in data without ME is quasi non existant, this plot 
#     is basically the same as biasME
(pBiasTrueME <- ggplot(plotBiasR2,
                       aes(x = interaction(pTrash, N, sep = " x "), y = biasTrueME, 
                           group = trueR2, colour = trueR2)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_grid(rel ~ lin_inter) +
    ylab(paste0("bias in R2")) +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle("R2 with ME - target R2") +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))
