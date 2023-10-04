library(ggplot2)

source("setParameters.R")

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

resFolder <- "results"
load(paste0(resFolder, "/estB_initialResults.rda"))

condGrid <- expand.grid(N = setParam$dgp$N, 
                        pTrash = setParam$dgp$pTrash)

condLabels <- sapply(seq_len(nrow(condGrid)), function(x) {
  paste0("N", condGrid[x, "N"], "_pTrash", condGrid[x, "pTrash"])})
names(res_estB) <- condLabels

linEffects <- sapply(seq_len(setParam$dgp$p), function(x) paste0("Var", x))
# choose variables for interaction that have no linear effects (R2 budget)
# interEffects <- c("Var1:Var2", "Var1:Var4", "Var2:Var3", "Var3:Var4")
interEffects <- c("Var5:Var6", "Var5:Var8", "Var6:Var7", "Var7:Var8")
nEffects <- length(c(linEffects, interEffects))


tmp_effects <- lapply(res_estB, function(x) {
  idxM <- which(rownames(x[["estB_M"]]) %in% c(linEffects, interEffects))
  idxSD <- which(rownames(x[["estB_SD"]]) %in% c(linEffects, interEffects))
  idxSE <- which(rownames(x[["estB_SE"]]) %in% c(linEffects, interEffects))
  subList <- list(estB_M = x[["estB_M"]][idxM,],
                  estB_SD = x[["estB_SD"]][idxSD,],
                  estB_SE = x[["estB_SE"]][idxSE,])
  subList <- do.call(cbind, subList)
  addVarName <- rep(c("estB_M", "estB_SD", "estB_SE"), each = length(setParam$dgp$condLabels))
  colnames(subList) <- paste0(addVarName, "__", colnames(subList))
  subList
})
names(tmp_effects) <- condLabels
estB_effects <- do.call(rbind, tmp_effects)
estB_effects <- cbind(N_pTrash = rep(names(tmp_effects), each = nEffects), 
                      data.frame(estB_effects))
estB_effects$effect <- rep(c(linEffects, interEffects), times = length(condLabels))

estB_effects <- tidyr::pivot_longer(estB_effects, cols = !c(N_pTrash, effect), names_to = "type", values_to = "values")
estB_effects <- tidyr::separate(estB_effects, N_pTrash, c("N", "pTrash"), sep = "_")
estB_effects <- tidyr::separate(estB_effects, type, c("measure", "type"), sep = "__")

estB_effects$N <- stringr::str_sub(estB_effects$N, start = 2L)
estB_effects$pTrash <- stringr::str_sub(estB_effects$pTrash, start = 7L)

estB_effects$R2 <- stringr::str_sub(estB_effects$type, start = 3L, end = 5L)
estB_effects$lin_inter <- stringr::str_sub(estB_effects$type, start = -7L, end = -1L)
# estB_effects <- tidyr::separate(estB_effects, lin_inter, c("lin", "inter"), sep = "_")

estB_effects$type <- NULL

estB_effects <- tidyr::pivot_wider(estB_effects, names_from = measure, values_from = values)

str(estB_effects)
estB_effects$N <- factor(estB_effects$N, levels = setParam$dgp$N)
estB_effects$pTrash <- factor(estB_effects$pTrash, levels = sort(setParam$dgp$pTrash, decreasing = T))

estB_effects$effectType <- ifelse(estB_effects$effect %in% 
                                    unique(estB_effects$effect)[stringr::str_detect(unique(estB_effects$effect), ":")], "inter", "lin")

trueB <- data.frame(R2 = setParam$dgp$Rsquared, setParam$dgp$trueEffects)
trueB <- tidyr::pivot_longer(trueB, cols = !R2, names_to = "lin", values_to = "trueB_lin")
trueB$lin <- stringr::str_sub(trueB$lin, start = 2L)
trueB$lin <- as.numeric(trueB$lin)
trueB$inter <- (1 - trueB$lin)

trueB2 <- data.frame(R2 = setParam$dgp$Rsquared, setParam$dgp$trueEffects)
trueB2 <- tidyr::pivot_longer(trueB2, cols = !R2, names_to = "inter", values_to = "trueB_inter")
trueB2$inter <- stringr::str_sub(trueB2$inter, start = 2L)

trueB <- merge(trueB, trueB2, by = c("R2", "inter"))
trueB <- tidyr::unite(trueB, "lin_inter", c(lin, inter), sep = "_")
colnames(trueB) <- c("R2", "lin_inter", "lin", "inter")
trueB <- tidyr::pivot_longer(trueB, cols = c(lin, inter), names_to = "effectType", values_to = "trueB")

trueB$R2 <- as.character(trueB$R2) # identical variable types to merge correctly
estB_effects <- merge(estB_effects, trueB, by = c("R2", "lin_inter", "effectType"))

# absolute bias biases interpretation of deviation in estimated coefficients due to different sizes of 
#     simulated coefficients depending on R2 budget size
estB_effects$absBiasB <- estB_effects$estB_M - estB_effects$trueB
estB_effects$relBiasB <- (estB_effects$estB_M - estB_effects$trueB) /  estB_effects$trueB


linInterTypes <- unique(estB_effects$lin_inter)
linIntLabels <- paste0("lin = ", 
                       stringr::str_sub(linInterTypes, start = 1L, end = 3L),
                       " + inter = ", stringr::str_sub(linInterTypes, start = -3L, end = -1L))
estB_effects$lin_inter <- factor(estB_effects$lin_inter, 
                                 levels = c("0.2_0.8", "0.5_0.5", "0.8_0.2"),
                                 labels = linIntLabels)

# to do:
#   - plot for all variables

colValues <- c("green3", "darkcyan", "darkblue", "darkmagenta")

plotRelBias <- function(data) {
  ggplot(data,
         aes(x = interaction(pTrash, N, sep = " x "), y = relBiasB, group = R2, colour = R2)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = relBiasB - estB_SE, ymax = relBiasB + estB_SE), width=.2) +
    scale_y_continuous(limits = c(-1, 0), breaks = seq(-1, 0, 0.2)) +
    scale_color_manual(values = colValues) +
    geom_hline(aes(yintercept = 0)) +
    facet_wrap(~lin_inter) +
    ylab("relative bias in est. \u03B2") +
    xlab("pTrash (decreasing) x N (increasing)") +
    ggtitle(paste0("relativ bias for ", iVar)) +
    theme(axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20))
}
  
for (iVar in c(linEffects, interEffects)) {
  subset <- estB_effects[which(estB_effects$effect == iVar), ]
  tmp_p <- plotRelBias(subset)
  
  plotName <- paste0("pRelBias_", sub(':','.', iVar))
  assign(plotName, tmp_p)
}

plotNames <- ls(pattern = "^pRelBias")
for (iPlot in plotNames) {
  pName <- paste0(iPlot, "_initialResults")
    
  ggplot2::ggsave(filename = paste0(plotFolder, "/", pName, ".eps"),
                    plot = eval(parse(text = iPlot)),
                    device = cairo_ps,
                    dpi = 300,
                    width = 11.3,
                    height = 8.22,
                    units = "in")
    
  ggplot2::ggsave(filename = paste0(plotFolder, "/", pName, ".png"),
                    plot = eval(parse(text = iPlot)),
                    width = 11.3,
                    height = 8.22,
                    units = "in") 
}

# # test plots
# # absolute bias
# ggplot(estB_effects[which(estB_effects$effect == "Var1"), ],
#        aes(x = interaction(pTrash, N), y = absBiasB, group = R2, colour = R2)) +
#   geom_point() +
#   geom_line() +
#   scale_color_manual(values = colValues) +
#   geom_hline(aes(yintercept = 0)) +
#   geom_errorbar(aes(ymin = absBiasB - estB_SE, ymax = absBiasB + estB_SE), width=.2) +
#   facet_wrap(~lin_inter) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# # averaged coeffients
# ggplot(estB_effects[which(estB_effects$effect == "Var1"), ],
#        aes(x = interaction(pTrash, N), y = estB_M, group = R2, colour = R2)) +
#   geom_point() +
#   geom_line() +
#   scale_color_manual(values = colValues) +
#   geom_errorbar(aes(ymin = estB_M - estB_SE, ymax = estB_M + estB_SE), width=.2) +
#   facet_wrap(~lin_inter) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

