library(SuperLearner)

library(future)
library(npcausalML)
source("simRR.R")
SL.gam3 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 3
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam4 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 4
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam5 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 5
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
list_of_sieves_uni   <- list(
  "no_sieve" = NULL,
  fourier_basis$new(orders = c(1,0)),
  fourier_basis$new(orders = c(2,0)),
  fourier_basis$new(orders = c(3,0)),
  fourier_basis$new(orders = c(4,0))
)

lrnr_gam3 <- Lrnr_pkg_SuperLearner$new("SL.gam3" , family = binomial())
lrnr_gam4 <- Lrnr_pkg_SuperLearner$new("SL.gam4", family = binomial() )
lrnr_gam5 <- Lrnr_pkg_SuperLearner$new("SL.gam5", family = binomial() )

lrnr_gam3pois <- Lrnr_pkg_SuperLearner$new("SL.gam3" , family = poisson())
lrnr_gam4pois <- Lrnr_pkg_SuperLearner$new("SL.gam4", family = poisson() )
lrnr_gam5pois <- Lrnr_pkg_SuperLearner$new("SL.gam5", family = poisson() )

hard <- F
pos <- F

onesim <- function(n) {

  sieve_list <- list_of_sieves_uni

  data <- as.data.frame(sim.RR(n, hard, pos))
  W <- data[,grep("^W", colnames(data))]
  A <- data$A
  Y <- data$Y
  W <- data.frame(W=data$W)
  EY1Wtrue <- data$EY1W
  EY0Wtrue <- data$EY0W
  pA1Wtrue <- data$pA1W
  EYWtrue <- ifelse(A==1, EY1Wtrue, EY0Wtrue)

  LRR <- EY1Wtrue - EY0Wtrue



  # sieve method
  lrnr_Y <- make_learner(Pipeline, Lrnr_cv$new(Stack$new(
    Lrnr_xgboost$new(max_depth =4, objective = "count:poisson"),
    Lrnr_xgboost$new(max_depth =5, objective = "count:poisson"),
    Lrnr_xgboost$new(max_depth =6, objective = "count:poisson")))#Stack$new(
    #Lrnr_stratified$new(Lrnr_gam$new(), "A"))
    , Lrnr_cv_selector$new(loss_squared_error))


  lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(
    Stack$new(
      Lrnr_xgboost$new(max_depth =4),
      Lrnr_xgboost$new(max_depth =5),
      Lrnr_xgboost$new(max_depth =6)
    )
  ), Lrnr_cv_selector$new(loss_squared_error))


  data_train <-  data #as.data.frame(sim.LRR(n, hard, pos))

  initial_likelihood <- npcausalML:::estimate_initial_likelihood(W=data_train[,c("W"), drop = F], data_train$A, data_train$Y,  weights = rep(1,n), lrnr_A, lrnr_Y, folds = 10, outcome_type = "continuous")
  data1 <- data
  data0 <- data
  data1$A <- 1
  data0$A <- 0
  taskY <- sl3_Task$new(data, covariates = c("W", "A"), outcome = "Y", folds = origami::folds_vfold(n), outcome_type = "continuous")
  folds <- taskY$folds
  taskY0 <- sl3_Task$new(data0, covariates = c("W", "A"), outcome = "Y", folds = folds, outcome_type = "continuous")
  taskY1 <- sl3_Task$new(data1, covariates = c("W", "A"), outcome = "Y", folds = folds, outcome_type = "continuous")
  taskA <- sl3_Task$new(data, covariates = c("W"), outcome = "A", folds = folds)

  pA1W_est <- initial_likelihood$internal$sl3_Learner_pA1W_trained$predict(taskA)
  EY1W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY1)
  EY0W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY0)

  pA1W_est <- pmax(pA1W_est, 0.01)
  pA1W_est <- pmin(pA1W_est, 0.99)

  data$weightsIPW <- data$Y/ifelse(data$A==1,pA1W_est, 1 - pA1W_est)
  sl3_Task_IPW <- sl3_Task$new(data, covariates = "W", outcome = "A", weights = "weightsIPW")

  LRR_library <- list(    Lrnr_earth$new(family = binomial()), Lrnr_xgboost$new(max_depth = 3, verbosity = 0, objective = "binary:logistic"), Lrnr_xgboost$new(max_depth = 5, verbosity = 0, objective = "binary:logistic"),   Lrnr_xgboost$new(max_depth = 7, verbosity = 0, objective = "binary:logistic"),  lrnr_gam3, lrnr_gam5, Lrnr_glm$new(family = binomial())    )

  IPW_learner <- Stack$new(LRR_library)
  IPW_learner <- IPW_learner$train(sl3_Task_IPW)
  preds_IPW <- apply(IPW_learner$predict(sl3_Task_IPW), 2, qlogis)


  LRR_library_plugin <- list(   Lrnr_earth$new(family = poisson()), Lrnr_xgboost$new(max_depth = 3, verbosity = 0, objective = "count:poisson"), Lrnr_xgboost$new(max_depth = 5, verbosity = 0, objective = "count:poisson"),   Lrnr_xgboost$new(max_depth = 7, verbosity = 0, objective = "count:poisson"),   lrnr_gam3pois, lrnr_gam5pois, Lrnr_glm$new(family = poisson())    )

  subst_compare <- Stack$new(LRR_library_plugin)
  LRR_library_strat <- lapply(LRR_library_plugin , function (lrnr) {
    Lrnr_stratified$new(lrnr, "A")
  })

  subst_compare <- Stack$new(LRR_library_strat)
  subst_compare <- subst_compare$train(taskY)
  subst_EY1W <- pmax(subst_compare$predict(taskY1), 1e-5)
  subst_EY0W <- pmax(subst_compare$predict(taskY0), 1e-5)
  subst_LRR <- log(subst_EY1W/ subst_EY0W)




  fit_npcausalML <- npcausalML(LRR_library,
                               W= W, A = A, Y = Y, V = data.frame(W = W$W),
                               EY1W = EY1W_est, EY0W = EY0W_est,  pA1W = pA1W_est,
                               sl3_Learner_EYAW = NULL, sl3_Learner_pA1W = NULL, outcome_type = "continuous", list_of_sieves = sieve_list,
                               outcome_function_plugin = outcome_function_plugin_LRR, weight_function_plugin = weight_function_plugin_LRR,
                               outcome_function_IPW = outcome_function_IPW_LRR, weight_function_IPW = weight_function_IPW_LRR,
                               design_function_sieve_plugin = design_function_sieve_plugin_LRR,
                               weight_function_sieve_plugin = weight_function_sieve_plugin_LRR,
                               design_function_sieve_IPW = design_function_sieve_IPW_LRR, weight_function_sieve_IPW = weight_function_sieve_IPW_LRR,
                               family_risk_function = poisson(),
                               efficient_loss_function = efficient_loss_function_LRR,
                               use_sieve_selector = FALSE,
                               transform_function = qlogis,
                               cross_validate_ERM = T, folds = origami::folds_vfold(length(A), 5))


  preds <- predict(fit_npcausalML,  data.frame(W = W$W), F)


  # Compute least-squares risk of predictions using oracle loss function.
  risks_oracle <- as.vector(apply(preds, 2, function(theta) {
    mean((theta -  LRR)^2)
  })[grep("plugin", colnames(preds))])

  # Compute estimated cross-validated one-step risk of predictions
  cvrisksDR <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
    loss <- efficient_loss_function_LRR(W, theta, A, Y, EY1W_est,EY0W_est, pA1W_est )
    mean(loss)
  })[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])

  # Compute estimated cross-validated oracle one-step risk of predictions

  cvrisksDRoracle <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
    loss <- efficient_loss_function_LRR(W, theta, A, Y, EY1Wtrue,EY0Wtrue, pA1Wtrue )
    mean(loss)
  })[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])






  risk_subst<-  apply(subst_LRR, 2, function(pred) {
    mean((pred - LRR)^2)
  })

  risk_IPW <-  apply(preds_IPW, 2, function(pred) {
    mean((pred - LRR)^2)
  })


  list( risk_IPW = risk_IPW,  risk_subst = risk_subst,    sieve =data.frame(grep("plugin", colnames(fit_npcausalML$cv_predictions), value = T), cvrisksDRoracle, cvrisksDR, risks_oracle))
}

nsims <- 50
print(500)
simresults500 <- lapply(1:nsims, function(i){
  print(i)
  onesim(500)
})

save(simresults500, file = paste0("simsLRR", hard,pos, "n500_3"))


print(1000)
simresults1000 <- lapply(1:nsims, function(i){
  print(i)
 onesim(1000)
})

save(simresults1000, file = paste0("simsLRR", hard,pos, "n1000_3"))

print(2500)
simresults2500 <- lapply(1:nsims, function(i){
  print(i)
  onesim(2500)
})

save(simresults2500, file = paste0("simsLRR", hard,pos, "n2500_3"))

print(5000)
simresults5000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(5000)
})

save(simresults5000, file = paste0("simsLRR", hard,pos, "n5000_3"))


print(10000)
simresults10000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(10000)
})

save(simresults10000, file = paste0("simsLRR", hard,pos, "n10000_3"))
