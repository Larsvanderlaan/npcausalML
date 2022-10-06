library(SuperLearner)
library(npcausalML)
library(future)
source("simCATE.R")
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
  fourier_basis$new(orders = c(1,0)),
  fourier_basis$new(orders = c(2,0)),
  fourier_basis$new(orders = c(3,0)),
  fourier_basis$new(orders = c(4,0))
)

lrnr_gam3 <- Lrnr_pkg_SuperLearner$new("SL.gam3" )
lrnr_gam4 <- Lrnr_pkg_SuperLearner$new("SL.gam4" )
lrnr_gam5 <- Lrnr_pkg_SuperLearner$new("SL.gam5" )


hard <- F
pos <- F

onesim <- function(n) {

  sieve_list <- list_of_sieves_uni

  data <- as.data.frame(sim.CATE(n, hard, pos))
  W <- data[,grep("^W", colnames(data))]
  A <- data$A
  Y <- data$Y
  W1 <- data$W1
  EY1Wtrue <- data$EY1W
  EY0Wtrue <- data$EY0W
  pA1Wtrue <- data$pA1W
  EYWtrue <- ifelse(A==1, EY1Wtrue, EY0Wtrue)

  CATE <- EY1Wtrue - EY0Wtrue



  # sieve method
  lrnr_Y <- make_learner(Pipeline, Lrnr_cv$new(Stack$new(
    Lrnr_xgboost$new(max_depth =4),
    Lrnr_xgboost$new(max_depth =5),
    Lrnr_xgboost$new(max_depth =6)))#Stack$new(
    #Lrnr_stratified$new(Lrnr_gam$new(), "A"))
    , Lrnr_cv_selector$new(loss_squared_error))


  lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(
    Stack$new(
      Lrnr_xgboost$new(max_depth =4),
      Lrnr_xgboost$new(max_depth =5),
      Lrnr_xgboost$new(max_depth =6)
    )
  ), Lrnr_cv_selector$new(loss_squared_error))


  data_train <-  data #as.data.frame(sim.CATE(n, hard, pos))

  initial_likelihood <- npcausalML:::estimate_initial_likelihood(W=data_train[,c("W1", "W2","W3")], data_train$A, data_train$Y,  weights = rep(1,n), lrnr_A, lrnr_Y, folds = 10)
  data1 <- data
  data0 <- data
  data1$A <- 1
  data0$A <- 0
  taskY <- sl3_Task$new(data, covariates = c("W1", "W2", "W3", "A"), outcome = "Y")
  taskY0 <- sl3_Task$new(data0, covariates = c("W1", "W2", "W3", "A"), outcome = "Y")
  taskY1 <- sl3_Task$new(data1, covariates = c("W1", "W2", "W3", "A"), outcome = "Y")
  taskA <- sl3_Task$new(data, covariates = c("W1", "W2", "W3"), outcome = "A")

  pA1W_est <- initial_likelihood$internal$sl3_Learner_pA1W_trained$predict(taskA)
  EY1W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY1)
  EY0W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY0)

  pA1W_est <- pmax(pA1W_est, 0.01)
  pA1W_est <- pmin(pA1W_est, 0.99)

  CATE_library <- list(    Lrnr_xgboost$new(max_depth = 1, verbosity = 0, nrounds = 10),    Lrnr_xgboost$new(max_depth = 2, verbosity = 0, nrounds = 10), Lrnr_xgboost$new(max_depth = 3, verbosity = 0, nrounds = 10), Lrnr_xgboost$new(max_depth = 4, verbosity = 0, nrounds = 10),   Lrnr_xgboost$new(max_depth = 5, verbosity = 0, nrounds = 10),   Lrnr_xgboost$new(max_depth = 6, verbosity = 0, nrounds = 10),   Lrnr_xgboost$new(max_depth = 7, verbosity = 0, nrounds = 10) )
  CATE_library <- list(Lrnr_sl$new(CATE_library, metalearner = Lrnr_cv_selector$new(loss_squared_error)))
  subst_compare <- Stack$new(CATE_library)
  #CATE_library_strat <- lapply(CATE_library , function (lrnr) {
   # Lrnr_stratified$new(lrnr, "A")
  #})

  subst_compare <- Stack$new(CATE_library)
  subst_compare <- subst_compare$train(taskY)
  subst_EY1W <-subst_compare$predict(taskY1)
  subst_EY0W <- subst_compare$predict(taskY0)
  subst_CATE <- subst_EY1W - subst_EY0W




  fit_npcausalML <- npcausalML(CATE_library,
                               W= W, A = A, Y = Y, V = data.frame(W1 = W$W1),
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


  preds <- predict(fit_npcausalML,  data.frame(W1 = W$W1), F)


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

  CATEonestepbench <- DR_learner(CATE_library, W1, A, Y, EY1W_est, EY0W_est, pA1W_est, NULL, NULL)
  CATEonestepbench <- apply(CATEonestepbench, 2, function(pred) {
    mean((pred - CATE)^2)
  })

  CATEonestepbenchoracle <- DR_learner(CATE_library, W1, A, Y, EY1Wtrue, EY0Wtrue, pA1Wtrue, NULL, NULL)
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

nsims <- 50
print(500)
simresults500 <- lapply(1:nsims, function(i){
  print(i)
  onesim(500)
})

save(simresults500, file = paste0("mainSimResults/","simsCATE", hard,pos, "n500_3"))


print(1000)
simresults1000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(1000)
})

save(simresults1000, file = paste0("mainSimResults/","simsCATE", hard,pos, "n1000_3"))

print(2500)
simresults2500 <- lapply(1:nsims, function(i){
  print(i)
  onesim(2500)
})

save(simresults2500, file = paste0("mainSimResults/","simsCATE", hard,pos, "n2500_3"))

print(5000)
simresults5000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(5000)
})

save(simresults5000, file = paste0("mainSimResults/", "simsCATE", hard,pos, "n5000_3"))


print(10000)
simresults10000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(10000)
})

save(simresults10000, file = paste0("mainSimResults/","simsCATE", hard,pos, "n10000_3"))
