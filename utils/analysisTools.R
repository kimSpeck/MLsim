# rather specific stuff to deal with simulation results 
# but needs to be done repeatedly -> functions to reduce amount of code in script
#     and to avoid copy-paste mistakes

loadRData <- function(fileName){
  ###
  # load rda file and assign it to any chosen object name instead of the original 
  #   object name
  ###
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

rowNames2col <- function(data, variableName = "rowNames"){
  ###
  # add rownames as another variable to the data frame/matrix 
  ## input:
  # data          - [data frame] data to add row names variable to
  # variableName  - [chr] variable name for new row names variable
  ## output:
  # data          - [data frame] data with row names variable
  ###
  # rownames to column in data frame
  data$rowNames <- rownames(data)
  # remove numbering (rownames need to be unique values, thus suffices .1, .2, ...)
  data$rowNames <- stringr::str_replace(data$rowNames, "\\.[:digit:]{1,}$", "")
  # name new column
  names(data)[names(data) == 'rowNames'] <- variableName
  return(data)
}

# rownames2col <- function(data, varName = "rownames") {
#   ###
#   # variables from rownames to own column to work woth variable information
#   ## input:
#   # data    - [data frame] specific result in one data frame 
#   # varName - [character string] name for new column with rowname information
#   ## output:
#   # data    - [data.frame] new column with rownames information
#   ###
#   # rownames to column in data frame
#   data$tmp <- rownames(data)
#   # remove numbering (rownames need to be unique values, thus suffices .1, .2, ...)
#   data$tmp <- stringr::str_replace(data$tmp, "\\.[:digit:]{1,}$", "")
#   # name new column
#   names(data)[names(data) == 'tmp'] <- varName
#   return(data)
# }

rbindResults <- function(data, result) {
  ###
  # rbind data from nested list + add simulated condition info (N_pTrash)
  # choose one result matrix by name out of different results lists
  ## input:
  # data    - [list] nested list with N_pTrash as lists
  # result  - [character] name of specific outcome
  ## output:
  #         - [data.frame] result in one data frame including index for N_Trash Labels
  ###
  tmpRes <- do.call(rbind, lapply(seq_along(data), function(subList) {
    cbind(data[[subList]][[result]], 
          idxN_pTrash = rep(subList, nrow(data[[subList]][[result]])))
  }))
  return(data.frame(tmpRes))
}

rbindSingleResults <- function(data) {
  ###
  # rbind data from list of matrices 
  ## input:
  # data    - [list] list with matrix for every N_pTrash_rel condition
  ## output:
  #         - [data.frame] result in one data frame 
  ###
  tmpRes <- do.call(rbind, lapply(seq_along(data), function(subList) {
    data[[subList]]
  }))
  return(data.frame(tmpRes))
}




idx2info <- function(data, cond_Np, type) {
  ###
  # replace idx values with informative labels (used to be matrices therefore idx **numbers**)
  ## input:
  # data    - [data frame] specific result in one data frame including index for labels
  ## output:
  # data    - [data.frame] numeric idx columns replaced by character type label columns
  ###
  data$idxCondLabel <- as.numeric(data$idxCondLabel)
  data$idxN_pTrash <- as.numeric(data$idxN_pTrash)
  
  data$condLabel <- setParam$dgp$condLabel[data$idxCondLabel]
  data$N_pTrash <- cond_Np[data$idxN_pTrash]
  
  if (type == "enet") {
    data <- tidyr::separate(data, N_pTrash, c("N", "pTrash", "rel", "factor", "indicators"), sep = "_")
    data$factor <- stringr::str_sub(data$factor, start = 2L)
    data$indicators <- stringr::str_sub(data$indicators, start = 4L)
  } else if (type == "gbm") {
    data <- tidyr::separate(data, N_pTrash, c("N", "pTrash", "rel"), sep = "_")  
  }
  
  data$N <- stringr::str_sub(data$N, start = 2L)
  data$pTrash <- stringr::str_sub(data$pTrash, start = 7L)
  data$rel <- stringr::str_sub(data$rel, start = 4L)
  
  
  data$R2 <- stringr::str_sub(data$condLabel, start = 3L, end = 5L)
  data$lin_inter <- stringr::str_sub(data$condLabel, start = -7L, end = -1L)
  
  data$condLabel <- NULL
  data$idxCondLabel <- NULL
  data$idxN_pTrash <- NULL
  return(data)
}

idx2infoNew <- function(data) {
  ###
  # replace idx values with informative labels (used to be matrices therefore idx **numbers**)
  ## input:
  # data    - [data frame] specific result in one data frame including index for labels
  ## output:
  # data    - [data.frame] numeric idx columns replaced by character type label columns
  ###
  data$idxCondLabel <- as.numeric(data$idxCondLabel)
  
  data$condLabel <- setParam$dgp$condLabel[data$idxCondLabel]
  
  data <- tidyr::separate(data, N_pTrash, c("N", "pTrash", "rel"), sep = "_")  
  
  data$N <- stringr::str_sub(data$N, start = 2L)
  data$pTrash <- stringr::str_sub(data$pTrash, start = 7L)
  data$rel <- stringr::str_sub(data$rel, start = 4L)
  
  data$R2 <- stringr::str_sub(data$condLabel, start = 3L, end = 5L)
  data$lin_inter <- stringr::str_sub(data$condLabel, start = -7L, end = -1L)
  
  data$condLabel <- NULL
  data$idxCondLabel <- NULL
  return(data)
}

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
