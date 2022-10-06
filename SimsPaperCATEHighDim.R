library(SuperLearner)
library(npcausalML)
library(future)
source("simCATEHighDim.R")


hard <- F
pos <- F

onesim <- function(n) {

  sieve_list <- list_of_sieves_high_dim

  data <- as.data.frame(sim.CATEHighDim(n, hard, pos))
  W <- data[,grep("^W", colnames(data))]
  A <- data$A
  Y <- data$Y

  EY1Wtrue <- data$EY1W
  EY0Wtrue <- data$EY0W
  pA1Wtrue <- data$pA1W
  EYWtrue <- ifelse(A==1, EY1Wtrue, EY0Wtrue)

  CATE <- EY1Wtrue - EY0Wtrue



  # sieve method
  lrnr_Y <- make_learner(Pipeline, Lrnr_cv$new(Stack$new(
    Lrnr_glmnet$new(),
    Lrnr_gam$new(),

    Lrnr_xgboost$new(max_depth =7),
    Lrnr_xgboost$new(max_depth =4),
    Lrnr_xgboost$new(max_depth =5),
    Lrnr_xgboost$new(max_depth =6)))#Stack$new(
    #Lrnr_stratified$new(Lrnr_gam$new(), "A"))
    , Lrnr_cv_selector$new(loss_squared_error))


  lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(
    Stack$new(
      Lrnr_glmnet$new(),
      Lrnr_gam$new(),
      Lrnr_xgboost$new(max_depth =5)
    )
  ), Lrnr_cv_selector$new(loss_squared_error))


  data_train <-  data #as.data.frame(sim.CATE(n, hard, pos))
  covariates <- paste0("W", 1:20)
  initial_likelihood <- npcausalML:::estimate_initial_likelihood(W=data_train[,covariates], data_train$A, data_train$Y,  weights = rep(1,n), lrnr_A, lrnr_Y, folds = 10)
  data1 <- data
  data0 <- data
  data1$A <- 1
  data0$A <- 0
  taskY <- sl3_Task$new(data, covariates = c(covariates, "A"), outcome = "Y")
  taskY0 <- sl3_Task$new(data0, covariates = c(covariates, "A"), outcome = "Y")
  taskY1 <- sl3_Task$new(data1, covariates = c(covariates, "A"), outcome = "Y")
  taskA <- sl3_Task$new(data, covariates = c(covariates), outcome = "A")

  pA1W_est <- initial_likelihood$internal$sl3_Learner_pA1W_trained$predict(taskA)
  EY1W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY1)
  EY0W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY0)

  pA1W_est <- pmax(pA1W_est, 0.05)
  pA1W_est <- pmin(pA1W_est, 0.95)

  CATE_library <- list(  Lrnr_ranger$new(max.depth = 7), Lrnr_ranger$new(max.depth = 10),     Lrnr_xgboost$new(max_depth = 3, verbosity = 0), Lrnr_xgboost$new(max_depth = 5, verbosity = 0),   Lrnr_xgboost$new(max_depth = 7, verbosity = 0), Lrnr_gam$new(), Lrnr_glmnet$new()    )

  subst_compare <- Stack$new(CATE_library)
  CATE_library_strat <- lapply(CATE_library , function (lrnr) {
    Lrnr_stratified$new(lrnr, "A")
  })

  subst_compare <- Stack$new(CATE_library_strat)
  subst_compare <- subst_compare$train(taskY)
  subst_EY1W <-subst_compare$predict(taskY1)
  subst_EY0W <- subst_compare$predict(taskY0)
  subst_CATE <- subst_EY1W - subst_EY0W




  fit_npcausalML <- npcausalML(CATE_library,
                               W= W, A = A, Y = Y, V = W,
                               EY1W = EY1W_est, EY0W = EY0W_est,  pA1W = pA1W_est,
                               sl3_Learner_EYAW = NULL, sl3_Learner_pA1W = NULL, outcome_type = "continuous", list_of_sieves = sieve_list,
                               outcome_function_plugin = outcome_function_plugin_CATE, weight_function_plugin = weight_function_plugin_CATE,
                               outcome_function_IPW = outcome_function_IPW_CATE, weight_function_IPW = weight_function_IPW_CATE,
                               design_function_sieve_plugin = design_function_sieve_plugin_CATE,
                               weight_function_sieve_plugin = weight_function_sieve_plugin_CATE,
                               design_function_sieve_IPW = design_function_sieve_IPW_CATE, weight_function_sieve_IPW = weight_function_sieve_IPW_CATE, transform_function = function(x){x},
                               family_risk_function = gaussian(),
                               efficient_loss_function = efficient_loss_function_CATE,
                               use_sieve_selector = FALSE,
                               cross_validate_ERM = T, folds = origami::folds_vfold(length(A), 5))


  preds <- predict(fit_npcausalML, W, F)


  # Compute least-squares risk of predictions using oracle loss function.
  risks_oracle <- as.vector(apply(preds, 2, function(theta) {
    mean((theta -  CATE)^2)
  })[grep("plugin", colnames(preds))])

  # Compute estimated cross-validated one-step risk of predictions
  cvrisksDR <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
    loss <- efficient_loss_function_CATE(W, theta, A, Y, EY1W_est,EY0W_est, pA1W_est )
    mean(loss)
  })[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])

  # Compute estimated cross-validated oracle one-step risk of predictions

  cvrisksDRoracle <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
    loss <- efficient_loss_function_CATE(W, theta, A, Y, EY1Wtrue,EY0Wtrue, pA1Wtrue )
    mean(loss)
  })[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])

  CATEonestepbench <- DR_learner(CATE_library, W , A, Y, EY1W_est, EY0W_est, pA1W_est, NULL, NULL)
  CATEonestepbench <- apply(CATEonestepbench, 2, function(pred) {
    mean((pred - CATE)^2)
  })

  CATEonestepbenchoracle <- DR_learner(CATE_library, W , A, Y, EY1Wtrue, EY0Wtrue, pA1Wtrue, NULL, NULL)
  CATEonestepbenchoracle <- apply(CATEonestepbenchoracle, 2, function(pred) {
    mean((pred - CATE)^2)
  })


  risk_subst<-  apply(subst_CATE, 2, function(pred) {
    mean((pred - CATE)^2)
  })

  Y.hat <- EY1W_est * pA1W_est + EY0W_est * (1-pA1W_est)
  W.hat <- pA1W_est
  fit <- grf::causal_forest(X = W, Y  = Y, W = A, Y.hat = Y.hat, W.hat = W.hat)
  preds_cf <-  fit$predictions
  risk_cf <- mean((CATE - preds_cf)^2)


  list(risk_cf = risk_cf, risk_subst = risk_subst, CATEonestepbenchoracle =CATEonestepbenchoracle, CATEonestepbench = CATEonestepbench, sieve =data.frame(grep("plugin", colnames(fit_npcausalML$cv_predictions), value = T), cvrisksDRoracle, cvrisksDR, risks_oracle))
}

print(500)
simresults <- lapply(1:20, function(i){
  print(i)
  onesim(500)
})

save(simresults, file = paste0("simsCATEHighDim", hard,pos, "n500_2"))


print(1000)
simresults <- lapply(1:20, function(i){
  print(i)
  onesim(1000)
})

save(simresults, file = paste0("simsCATEHighDim", hard,pos, "n1000_2"))

print(2500)
simresults <- lapply(1:20, function(i){
  print(i)
  onesim(2500)
})

save(simresults, file = paste0("simsCATEHighDim", hard,pos, "n2500_2"))

print(5000)
simresults <- lapply(1:20, function(i){
  print(i)
  onesim(5000)
})

save(simresults, file = paste0("simsCATEHighDim", hard,pos, "n5000_2"))


print(10000)
simresults <- lapply(1:20, function(i){
  print(i)
  onesim(10000)
})

save(simresults, file = paste0("simsCATEHighDim", hard,pos, "n10000_2"))

onestepbenchoracle <- rowMeans(do.call(cbind, lapply(simresults, `[[`, 1)))
onestepbench  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, 2)))



lrnr_names <- simresults[[1]]$sieve[[1]]

cvrisksDRoracle <- rowMeans(do.call(cbind, lapply(simresults, function(item) {
  item$sieve$cvrisksDRoracle
})) )

cvrisksDR <- rowMeans(do.call(cbind, lapply(simresults, function(item) {
  item$sieve$cvrisksDR
})) )

risks_oracle <- rowMeans(do.call(cbind, lapply(simresults, function(item) {
  item$sieve$risks_oracle
})) )

dt <- data.table(lrnr_full = lrnr_names, cvrisksDRoracle,cvrisksDR,  risks_oracle)
dt$lrnr[ grep("gam", dt$lrnr_full)] <- "gam"
dt$lrnr[ grep("glm", dt$lrnr_full)] <- "glm"
dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
dt$degree <- as.numeric(stringr::str_match(dt$lrnr_full, "fourier_basis_([0-9]+)")[,2])
dt$degree[grep("no_sieve", dt$lrnr_full)] <- 0

library(ggplot2)

ggplot(dt, aes(x = degree, y = cvrisksDRoracle, color = lrnr, group = lrnr)) + geom_line() + geom_hline(aes(yintercept = ) )
ggplot(dt, aes(x = degree, y = cvrisksDR, color = lrnr, group = lrnr)) + geom_line()
ggplot(dt, aes(x = degree, y = risks_oracle, color = lrnr, group = lrnr)) + geom_line()  + geom_hline(yintercept = onestepbench, color = c("red", "green", "blue"))

