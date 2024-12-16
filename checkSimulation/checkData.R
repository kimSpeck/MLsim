load("data/simDataN100_pTrash10_rel0.8_dalek.rda")
dataDalek <- data

load("data/simDataN100_pTrash10_rel0.8_nerdholland.rda")
dataNerdy <- data

dataDalek[["1"]]$yMat - dataNerdy[["1"]]$yMat
dataDalek[["1"]]$X_int == dataNerdy[["1"]]$X_int

all(dataDalek[["1"]]$X_int == dataNerdy[["1"]]$X_int)

resGBM <- readRDS("results/resultsModelGBM_N100_pTrash10_rel0.6.rds")
resENETw <- readRDS("results/resultsModelENETw_N100_pTrash10_rel0.6.rds")
resENETwo <- readRDS("results/resultsModelENETwo_N100_pTrash10_rel0.6.rds")         
