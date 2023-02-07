
#### Parameter festlegen ####
# für Reproduzierbarkeit
set.seed(42)

# Zahl der VP
N <- 500

# Zahl der Prädiktoren
P <- 20

# Anteil der Regressionsgewichte != 0
p_features <- 0.2

# mittleres Regressionsgewicht und SD der Gewichte
mean_b <- 0.3
sd_b <- 0.2

# mittlere Korrelation der Prädiktoren und deren SD
mean_r <- 0
sd_r = 0.2

# Fehlervarianz
sigma_e <- 1
####

#### Daten simulieren
library(mvtnorm)
library(truncnorm)
library(lavaan)
library(glmnet)

# Prädiktormittelwerte
mean_x <- rep(0,P)

# zufällig gezogene Prädiktorkorrelationen
R_x <- lav_matrix_upper2full(rtruncnorm(P*(P-1)/2, mean = mean_r, sd = sd_r, a = -0.5, b = 0.5), diagonal = F)
diag(R_x) <- 1

# ziehe die Prädiktoren aus einer Normalverteilung
X <- rmvnorm(n = N, mean = mean_x, sigma = R_x)

pop_model <- paste0("~ (0 +", paste("X", 1:P, sep = "", collapse = "+"), ")^4",
                    paste("+ poly(X", 1:P, ",2)", sep = "", collapse = "+")
                    )
X_int <- model.matrix(as.formula(pop_model),data.frame(X))




# Regressionskoeffizienten ziehen
b <- rtruncnorm(n = ncol(X_int), a = -0.5, b = 0.5, mean = mean_b, sd = sd_b)
# 0 features auf 0 setzen
b[1:round((1-p_features)*length(b))] <- 0

# AV berechnen
y <- X_int %*% b + rnorm(n = N, mean = 0, sd = sigma_e)

# Data.frame erzeugen
dat <- data.frame(y,X)

# fitte eine regularisierte Regression
fit_cv <- cv.glmnet(X, y)
fit <- glmnet(X, y, lambda = fit_cv$lambda.1se)


b_est <- fit$beta
# wie viele falsch positive?
length(which(b == 0 & b_est != 0))/sum(b == 0)

b
b_est


###################### glm via caret
library(caret)
library(gbm)
preds <- names(dat)
res_train <- NULL
res_test <- NULL
res_vi <- data.frame(matrix(ncol = 1, nrow = ncol(dat)))

tictoc::tic()
for (i in 1:100) {
  set.seed(i)
  print(i)
  intrain <- createDataPartition(dat$y, p = .80, list = FALSE) #partitioning the data in training (80%) and testing data (20%)
  intrain <- as.vector(intrain)
  
  train <- dat[ intrain,] 
  test  <- dat[-intrain,]

  #model <-  train(y ~ .,
  #                train[,c(preds)],
  #                method = "glmnet", #elastic net regularized regression
  #                metric = "RMSE",
  #                preProc = "scale",
  #                trControl = trainControl(method="cv", number=10), tuneLength = 21) #model training using ten-fold cross-validation 

  #model <-  train(y ~ .,
  #                train[,c(preds)],
  #                method = "gbm", #gradient boosting machines
  #                metric = "RMSE",
  #                preProc = "scale",
  #                trControl = trainControl(method="cv", number=10))
  
  grid <- expand.grid(interaction.depth = c(1,2,3,4,5),   #hängt von simulation ab
                      n.minobsinnode = c(5,10,20,50),
                      n.trees = c(50,100,150,300,500,1000),
                      shrinkage = seq(.001, .201, .02))
  
  model <- train(y ~ .,
                 train[,c(preds)],
                 method = "gbm", #gradient boosting machines
                 metric = "RMSE",
                 trControl = trainControl(method="cv", number=10), #model training using ten-fold cross-validation
                 tuneGrid = grid) #setting the tuning grid
  
  
  predictions <- predict(model, train) 
  #returns explained variance (R?), Root Mean Squared Error (RMSE), and Mean Absolute Error (MAE) for model training:
  eval_train <- caret::postResample(pred = predictions, obs = train$y)  
  
  predictions <- predict(model, test) 
  #returns explained variance (R?), Root Mean Squared Error (RMSE), and Mean Absolute Error (MAE) for model testing:
  eval_test <- caret::postResample(pred = predictions, obs = test$y)

  vi <- t(varImp(model)$importance)  
  
  res_train <-  rbind(res_train, eval_train)
  res_vi <- cbind(res_vi, vi)
  res_test <- rbind(res_test, eval_test)
  
}
tictoc::toc()


