library(SuperLearner)
library(npcausalML)
library(future)
library(sl3)
# plan(cluster, workers = 3)
source("./FinalSimulationCode/simCATE.R")
SL.gam3 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 3
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam4 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 4
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam5 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 5
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
list_of_sieves_uni   <- list(
  "no_sieve" = NULL,
  fourier_basis$new(orders = c(1,0)),
  fourier_basis$new(orders = c(2,0)),
  fourier_basis$new(orders = c(3,0)),
  fourier_basis$new(orders = c(1,1))
)

lrnr_gam3 <- Lrnr_pkg_SuperLearner$new("SL.gam3" )
lrnr_gam4 <- Lrnr_pkg_SuperLearner$new("SL.gam4" )
lrnr_gam5 <- Lrnr_pkg_SuperLearner$new("SL.gam5" )




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
    Lrnr_xgboost$new(max_depth =3),
    Lrnr_xgboost$new(max_depth =4),
    Lrnr_xgboost$new(max_depth =5),
    Lrnr_xgboost$new(max_depth =6)))#Stack$new(
    #Lrnr_stratified$new(Lrnr_gam$new(), "A"))
    , Lrnr_cv_selector$new(loss_squared_error))


  lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(
    Stack$new(
      Lrnr_xgboost$new(max_depth =3),
      Lrnr_xgboost$new(max_depth =4),
      Lrnr_xgboost$new(max_depth =5),
      Lrnr_xgboost$new(max_depth =6)
    )
  ), Lrnr_cv_selector$new(loss_squared_error))


  data_train <-  data #as.data.frame(sim.CATE(n, hard, pos))

  initial_likelihood <- npcausalML:::estimate_initial_likelihood(W=data_train[,c("W1", "W2","W3")], data_train$A, data_train$Y,  weights = rep(1,n), lrnr_A, lrnr_Y, folds = 10)
  print("inirtal lik")
  data1 <- data
  data0 <- data
  data1$A <- 1
  data0$A <- 0
  folds <- initial_likelihood$internal$folds
  taskY <- sl3_Task$new(data, covariates = c(grep("W", colnames(data_train), value = T), "A"), outcome = "Y", folds = folds)
  taskY0 <- sl3_Task$new(data0, covariates = c(grep("W", colnames(data_train), value = T), "A"), outcome = "Y", folds = folds)
  taskY1 <- sl3_Task$new(data1, covariates = c(grep("W", colnames(data_train), value = T), "A"), outcome = "Y", folds = folds)
  taskA <- sl3_Task$new(data, covariates = c(grep("W", colnames(data_train), value = T)), outcome = "A", folds = 5)

  pA1W_est <- initial_likelihood$internal$sl3_Learner_pA1W_trained$predict(taskA)
  EY1W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY1)
  EY0W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY0)

  pA1W_est <- pmax(pA1W_est, 0.01)
  pA1W_est <- pmin(pA1W_est, 0.99)

  lrnr_xg <- list(    Lrnr_xgboost$new(max_depth = 1, verbosity = 0, nrounds = 10),    Lrnr_xgboost$new(max_depth = 2, verbosity = 0, nrounds = 10), Lrnr_xgboost$new(max_depth = 3, verbosity = 0, nrounds = 10), Lrnr_xgboost$new(max_depth = 4, verbosity = 0, nrounds = 10),   Lrnr_xgboost$new(max_depth = 5, verbosity = 0, nrounds = 10) )
  lrnr_xg_sl <-  Lrnr_sl$new(lrnr_xg, metalearner = Lrnr_cv_selector$new(loss_squared_error))


  lrnr_rf <- list(Lrnr_ranger$new(
    max.depth = 3, name = "Lrnr_rf_3_xg"),
    Lrnr_ranger$new(
      max.depth = 5, name = "Lrnr_rf_5_xg"),
    Lrnr_ranger$new(
      max.depth = 7, name = "Lrnr_rf_7_xg"),
    Lrnr_ranger$new(
      max.depth = 9, name = "Lrnr_rf_9_xg"))

  lrnr_rf_sl <-  Lrnr_sl$new(lrnr_rf, metalearner = Lrnr_cv_selector$new(loss_squared_error) )


  CATE_library <- c(lrnr_rf,   lrnr_xg  )

  CATE_library_subst <- c(CATE_library ,list(lrnr_xg_sl , lrnr_rf_sl  ))


  subst_compare <- Stack$new(CATE_library_subst)
  #subst_compare_trained <- subst_compare$train(taskY)
  subst_EY1W_trained <-subst_compare$train(taskY1[A==1]$next_in_chain(covariates = c("W1", "W2", "W3")))
  subst_EY0W_trained <- subst_compare$train(taskY0[A==0]$next_in_chain(covariates = c("W1", "W2", "W3")))

  subst_EY1W <-subst_EY1W_trained$predict(taskY1$next_in_chain(covariates = c("W1", "W2", "W3")))
  subst_EY0W <- subst_EY0W_trained$predict(taskY0$next_in_chain(covariates = c("W1", "W2", "W3")))
  subst_CATE <- subst_EY1W - subst_EY0W #subst_compare_trained$predict(taskY1) - subst_compare_trained$predict(taskY0)

  # apply(subst_EY1W -subst_EY0W , 2, function(p) {mean((p - CATE)^2)})
  # apply(subst_compare_trained$predict(taskY1) - subst_compare_trained$predict(taskY0), 2, function(p) {mean((p - CATE)^2)})
  print("initial CATE")
  t <- Sys.time()
  # V = data.frame(W1=W$W1)
  fit_npcausalML <- EP_learn(CATE_library,V = as.data.frame(W), A = A, Y = Y, EY1W = EY1W_est  , EY0W = EY0W_est  , pA1W = pA1W_est, sieve_basis_generator_list = sieve_list ,EP_learner_spec = EP_learner_spec_CATE, cross_validate = TRUE, nfolds = 5)
  print(  Sys.time() - t)
print("npcausal")

  preds <- fit_npcausalML$full_predictions


  # Compute least-squares risk of predictions using oracle loss function.
  risks_oracle <- as.vector(apply(preds, 2, function(theta) {
    mean((theta -  CATE)^2)
  }))

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
  lrnrs <- gsub("[._]sieve_fourier.+", "", lrnrs_full)
  lrnrs <- gsub("_no_sieve", "", lrnrs)
  degree <- (stringr::str_match(lrnrs_full, "fourier_basis_([0-9_]+)")[,2])
  degree[grep("no_sieve", lrnrs_full)] <- "0"
  uniq_degrees <- sort(unique(degree))


  tmp <- data.table(lrnrs_full, lrnrs , degree, risk = cvrisksDR, risks_oracle = risks_oracle, cvrisksDR = cvrisksDR, cvrisksDRoracle)

  rf_keep <- tmp[grep("rf", lrnrs_full),risks_oracle[which.min(risk)], by = degree]$V1
  xg_keep <- tmp[grep("xgboost", lrnrs_full),risks_oracle[which.min(risk)], by = degree]$V1
  print(xg_keep)
  names(xg_keep) <- paste0("Lrnr_xgboost_cv", "_fourier.basis_", uniq_degrees, "_plugin")
  names(rf_keep) <- paste0("Lrnr_rf_cv", "_fourier.basis_", uniq_degrees, "_plugin")
  risks_oracle <- c(risks_oracle, xg_keep, rf_keep)
print("here")
  rf_keep <- tmp[grep("rf", lrnrs_full),cvrisksDR[which.min(risk)], by = degree]$V1
  xg_keep <- tmp[grep("xgboost", lrnrs_full),cvrisksDR[which.min(risk)], by = degree]$V1
  names(xg_keep) <- paste0("Lrnr_xgboost_cv", "_fourier.basis_", uniq_degrees, "_plugin")
  names(rf_keep) <- paste0("Lrnr_rf_cv", "_fourier.basis_", uniq_degrees, "_plugin")
  cvrisksDR <- c(cvrisksDR, xg_keep, rf_keep)
  print("here")
  rf_keep <- tmp[grep("rf", lrnrs_full),cvrisksDRoracle[which.min(risk)], by = degree]$V1
  xg_keep <- tmp[grep("xgboost", lrnrs_full),cvrisksDRoracle[which.min(risk)], by = degree]$V1
  names(xg_keep) <- paste0("Lrnr_xgboost_cv", "_fourier.basis_", uniq_degrees, "_plugin")
  names(rf_keep) <- paste0("Lrnr_rf_cv", "_fourier.basis_", uniq_degrees, "_plugin")
  cvrisksDRoracle <- c(cvrisksDRoracle, xg_keep, rf_keep)
  print("here")
  sieve_names <- c(colnames(fit_npcausalML$cv_predictions),   names(xg_keep),   names(rf_keep))

  CATE_library <- c(CATE_library, list(lrnr_rf_sl, lrnr_xg_sl))
  CATEonestepbench <- DR_learner(CATE_library, as.data.frame(W), A, Y, EY1W_est, EY0W_est, pA1W_est, NULL, NULL)
  CATEonestepbench <- apply(CATEonestepbench, 2, function(pred) {
    mean((pred - CATE)^2)
  })
  names(CATEonestepbench) <- c(unique(tmp$lrnrs) , "Lrnr_rf_cv",  "Lrnr_xgboost_cv")

  CATEonestepbenchoracle <- DR_learner(CATE_library, as.data.frame(W), A, Y, EY1Wtrue, EY0Wtrue, pA1Wtrue, NULL, NULL)
  CATEonestepbenchoracle <- apply(CATEonestepbenchoracle, 2, function(pred) {
    mean((pred - CATE)^2)
  })
  names(CATEonestepbenchoracle) <- c(unique(tmp$lrnrs), "Lrnr_rf_cv",  "Lrnr_xgboost_cv")



  risk_subst<-  apply(subst_CATE, 2, function(pred) {
    mean((pred - CATE)^2)
  })
  risk_subst_cv <- mean((EY1W_est - EY0W_est - CATE)^2)

  Y.hat <- EY1W_est * pA1W_est + EY0W_est * (1-pA1W_est)
  W.hat <- pA1W_est
  print("cf")
  fit <- grf::causal_forest(X = W, Y  = Y, W = A, Y.hat = Y.hat, W.hat = W.hat, num.trees = 500, tune.parameters = "all")
  preds_cf <-  fit$predictions
  risk_cf <- mean((CATE - preds_cf)^2)


  list(risk_subst_cv = risk_subst_cv, risk_cf = risk_cf, risk_subst = risk_subst, CATEonestepbenchoracle =CATEonestepbenchoracle, CATEonestepbench = CATEonestepbench, sieve =data.frame(sieve_names, cvrisksDRoracle, cvrisksDR, risks_oracle))
}

# n <- 5000
# hard <- pos <- F
# nsims<- 1000

hard <- hard == "TRUE"
pos <- pos == "TRUE"
n <- as.numeric(n)
if(n==5000) {
  nsims <- 1000
}
simresults <- lapply(1:nsims, function(i){try({
  print(i)
  onesim(n)
})
})
save(simresults, file = paste0("mainSimResults2/","simsCATE", hard,pos, "n", n, "_xgboost"))


