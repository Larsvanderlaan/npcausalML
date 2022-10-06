library(SuperLearner)
library(npcausalML)
library(future)
source("simCATE.R")
SL.gam1 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 1
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam2 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 2
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
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

lrnr_gam1 <- Lrnr_pkg_SuperLearner$new("SL.gam1" , name = "Lrnr_gam_s1")
lrnr_gam2 <- Lrnr_pkg_SuperLearner$new("SL.gam2", name = "Lrnr_gam_s2")
lrnr_gam3 <- Lrnr_pkg_SuperLearner$new("SL.gam3", name = "Lrnr_gam_s3")
lrnr_gam4 <- Lrnr_pkg_SuperLearner$new("SL.gam4" , name = "Lrnr_gam_s4")
lrnr_gam5 <- Lrnr_pkg_SuperLearner$new("SL.gam5" , name = "Lrnr_gam_s5")


hard <- T
pos <- T

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

  lrnr_gam <- list(   lrnr_gam1, lrnr_gam2, lrnr_gam3, lrnr_gam4, lrnr_gam5 )
  lrnr_gam_sl <-  Lrnr_sl$new(lrnr_gam, metalearner = Lrnr_cv_selector$new(loss_squared_error))


  lrnr_lm <- list(Lrnr_earth$new(), Lrnr_glm$new())



  CATE_library <- c(lrnr_gam,   lrnr_lm  )

  CATE_library_subst <- c(CATE_library ,list(lrnr_gam_sl  ))


  subst_compare <- Stack$new(CATE_library_subst)
  subst_compare_trained <- subst_compare$train(taskY)
  subst_EY1W_trained <-subst_compare$train(taskY1[A==1]$next_in_chain(covariates = c("W1", "W2", "W3")))
  subst_EY0W_trained <- subst_compare$train(taskY0[A==0]$next_in_chain(covariates = c("W1", "W2", "W3")))

  subst_EY1W <-subst_EY1W_trained$predict(taskY1$next_in_chain(covariates = c("W1", "W2", "W3")))
  subst_EY0W <- subst_EY0W_trained$predict(taskY0$next_in_chain(covariates = c("W1", "W2", "W3")))
  subst_CATE <- subst_EY1W - subst_EY0W #subst_compare_trained$predict(taskY1) - subst_compare_trained$predict(taskY0)

  # apply(subst_EY1W -subst_EY0W , 2, function(p) {mean((p - CATE)^2)})
  # apply(subst_compare_trained$predict(taskY1) - subst_compare_trained$predict(taskY0), 2, function(p) {mean((p - CATE)^2)})

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
                               cross_validate_ERM = T, folds = origami::folds_vfold(length(A), 8))


  preds <- predict(fit_npcausalML,  data.frame(W1 = W$W1), F)


  # Compute least-squares risk of predictions using oracle loss function.
  risks_oracle <- as.vector(apply(preds, 2, function(theta) {
    mean((theta -  CATE)^2)
  })[grep("plugin", colnames(preds))])

  # Compute estimated cross-validated one-step risk of predictions
  cvrisksDR <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
    loss <- efficient_loss_function_CATE(W, theta, A, Y, EY1W_est,EY0W_est, pA1W_est )
    mean(loss)
  }))#[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])

  # Compute estimated cross-validated oracle one-step risk of predictions

  cvrisksDRoracle <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
    loss <- efficient_loss_function_CATE(W, theta, A, Y, EY1Wtrue,EY0Wtrue, pA1Wtrue )
    mean(loss)
  }))#[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])
  lrnrs_full <-  colnames(fit_npcausalML$cv_predictions)
  lrnrs <- gsub("[._]fourier.+", "", lrnrs_full)
  lrnrs <- gsub("[._]no_sieve.+", "", lrnrs)
  degree <- as.numeric(stringr::str_match(lrnrs_full, "fourier_basis_([0-9]+)")[,2])
  degree[grep("no_sieve", lrnrs_full)] <- 0



  tmp <- data.table(lrnrs_full, lrnrs , degree, risk = cvrisksDR, risks_oracle = risks_oracle, cvrisksDR = cvrisksDR, cvrisksDRoracle)
  gam_keep <- tmp[grep("gam", lrnrs_full),risks_oracle[which.min(risk)], by = degree]$V1
  names(gam_keep) <- paste0("Lrnr_gam_cv", "_fourier.basis_", 0:4, "_plugin")
  risks_oracle <- c(risks_oracle, gam_keep)


  gam_keep <- tmp[grep("gam", lrnrs_full),cvrisksDR[which.min(risk)], by = degree]$V1
  names(gam_keep) <- paste0("Lrnr_gam_cv", "_fourier.basis_", 0:4, "_plugin")
  cvrisksDR <- c(cvrisksDR, gam_keep)

  gam_keep <- tmp[grep("gam", lrnrs_full),cvrisksDRoracle[which.min(risk)], by = degree]$V1
  names(gam_keep) <- paste0("Lrnr_gam_cv", "_fourier.basis_", 0:4, "_plugin")
  cvrisksDRoracle <- c(cvrisksDRoracle, gam_keep)

  sieve_names <- c(colnames(fit_npcausalML$cv_predictions),   names(gam_keep))

  CATE_library <- c(CATE_library, list(lrnr_gam_sl))
  CATEonestepbench <- DR_learner(CATE_library, W1, A, Y, EY1W_est, EY0W_est, pA1W_est, NULL, NULL)
  CATEonestepbench <- apply(CATEonestepbench, 2, function(pred) {
    mean((pred - CATE)^2)
  })
  names(CATEonestepbench) <- c(unique(tmp$lrnrs) , "Lrnr_gam_cv")

  CATEonestepbenchoracle <- DR_learner(CATE_library, W1, A, Y, EY1Wtrue, EY0Wtrue, pA1Wtrue, NULL, NULL)
  CATEonestepbenchoracle <- apply(CATEonestepbenchoracle, 2, function(pred) {
    mean((pred - CATE)^2)
  })
  names(CATEonestepbenchoracle) <- c(unique(tmp$lrnrs), "Lrnr_gam_cv")



  risk_subst<-  apply(subst_CATE, 2, function(pred) {
    mean((pred - CATE)^2)
  })

  Y.hat <- EY1W_est * pA1W_est + EY0W_est * (1-pA1W_est)
  W.hat <- pA1W_est
  fit <- grf::causal_forest(X = W, Y  = Y, W = A, Y.hat = Y.hat, W.hat = W.hat)
  preds_cf <-  fit$predictions
  risk_cf <- mean((CATE - preds_cf)^2)


  list(risk_cf = risk_cf, risk_subst = risk_subst, CATEonestepbenchoracle =CATEonestepbenchoracle, CATEonestepbench = CATEonestepbench, sieve =data.frame(sieve_names, cvrisksDRoracle, cvrisksDR, risks_oracle))
}

nsims <- 10
print(500)
simresults500 <- lapply(1:nsims, function(i){
  print(i)
  onesim(500)
})

save(simresults500, file = paste0("mainSimResults/","simsCATE", hard,pos, "n500_gam"))


print(1000)
simresults1000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(1000)
})

save(simresults1000, file = paste0("mainSimResults/","simsCATE", hard,pos, "n1000_gam"))

print(2500)
simresults2500 <- lapply(1:nsims, function(i){
  print(i)
  onesim(2500)
})

save(simresults2500, file = paste0("mainSimResults/","simsCATE", hard,pos, "n2500_gam"))

print(5000)
simresults5000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(5000)
})

save(simresults5000, file = paste0("mainSimResults/", "simsCATE", hard,pos, "n5000_gam"))


print(10000)
simresults10000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(10000)
})

save(simresults10000, file = paste0("mainSimResults/","simsCATE", hard,pos, "n10000_gam"))
