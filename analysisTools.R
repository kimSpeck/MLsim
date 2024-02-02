# rather specific stuff to deal with simulation results 
# but needs to be done repeatedly -> functions to reduce amount of code in script
#     and to avoid copy-paste mistakes
rbindResults <- function(data, result) {
  ###
  # rbind data from nested list + add simulated condition info (N_pTrash)
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

idx2info <- function(data) {
  ###
  # replace idx values with informative labels (used to be matrices therefore idx **numbers**)
  ## input:
  # data    - [data frame] specific result in one data frame including index for labels
  ## output:
  # data    - [data.frame] numeric idx columns replaced by character type label columns
  ###
  data$condLabel <- setParam$dgp$condLabel[data$idxCondLabel]
  data$N_pTrash <- condN_pTrash[data$idxN_pTrash]
  
  data <- tidyr::separate(data, N_pTrash, c("N", "pTrash", "rel", "factor", "indicators"), sep = "_")
  data$N <- stringr::str_sub(data$N, start = 2L)
  data$pTrash <- stringr::str_sub(data$pTrash, start = 7L)
  data$rel <- stringr::str_sub(data$rel, start = 4L)
  data$factor <- stringr::str_sub(data$factor, start = 2L)
  data$indicators <- stringr::str_sub(data$indicators, start = 4L)
  
  data$R2 <- stringr::str_sub(data$condLabel, start = 3L, end = 5L)
  data$lin_inter <- stringr::str_sub(data$condLabel, start = -7L, end = -1L)
  
  data$condLabel <- NULL
  data$idxCondLabel <- NULL
  data$idxN_pTrash <- NULL
  return(data)
}

