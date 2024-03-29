print("OK")
print(.libPaths())
library(SuperLearner)
library(npcausalML)
library(future)
print(getwd())
nsims<- 1000
source("./FinalSimulationCode/simCATE.R")
print("OK")
SL.gam1 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 1
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam2 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 2
  SuperLearner::SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
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

lrnr_gam1 <- Lrnr_pkg_SuperLearner$new("SL.gam1" , name = "Lrnr_gam_s1_x")
lrnr_gam2 <- Lrnr_pkg_SuperLearner$new("SL.gam2", name = "Lrnr_gam_s2_x")
lrnr_gam3 <- Lrnr_pkg_SuperLearner$new("SL.gam3", name = "Lrnr_gam_s3_x")
lrnr_gam4 <- Lrnr_pkg_SuperLearner$new("SL.gam4" , name = "Lrnr_gam_s4_x")
lrnr_gam5 <- Lrnr_pkg_SuperLearner$new("SL.gam5" , name = "Lrnr_gam_s5_x")


onesim <- function(n) {
  library(sl3)
  library(SuperLearner)
  library(npcausalML)
  library(future)
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

  initial_likelihood <- npcausalML:::estimate_initial_likelihood(W=data_train[,c("W1", "W2","W3")], data_train$A, data_train$Y,  weights = rep(1,n), lrnr_A, lrnr_Y, folds = 5)
  data1 <- data
  data0 <- data
  data1$A <- 1
  data0$A <- 0
  folds <- initial_likelihood$internal$folds
  taskY <- sl3_Task$new(data, covariates = c(grep("W", colnames(data_train), value = T), "A"), outcome = "Y", folds = folds)
  taskY0 <- sl3_Task$new(data0, covariates = c(grep("W", colnames(data_train), value = T), "A"), outcome = "Y", folds = folds)
  taskY1 <- sl3_Task$new(data1, covariates = c(grep("W", colnames(data_train), value = T), "A"), outcome = "Y", folds = folds)
  taskA <- sl3_Task$new(data, covariates = c(grep("W", colnames(data_train), value = T)), outcome = "A", folds = folds)

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
  subst_compare_trained <- delayed_learner_train(subst_compare, taskY)$compute()
  #subst_compare_trained <- subst_compare$train(taskY)
  subst_EY1W_trained <- delayed_learner_train(subst_compare, taskY1[A==1]$next_in_chain(covariates = c("W1", "W2", "W3")))$compute()
  subst_EY0W_trained <- delayed_learner_train(subst_compare, taskY1[A==0]$next_in_chain(covariates = c("W1", "W2", "W3")))$compute()
  #subst_EY1W_trained <-subst_compare$train(taskY1[A==1]$next_in_chain(covariates = c("W1", "W2", "W3")))
  #subst_EY0W_trained <- subst_compare$train(taskY0[A==0]$next_in_chain(covariates = c("W1", "W2", "W3")))

  subst_EY1W <-subst_EY1W_trained$predict(taskY1$next_in_chain(covariates = c("W1", "W2", "W3")))
  subst_EY0W <- subst_EY0W_trained$predict(taskY0$next_in_chain(covariates = c("W1", "W2", "W3")))
  subst_CATE <- subst_EY1W - subst_EY0W #subst_compare_trained$predict(taskY1) - subst_compare_trained$predict(taskY0)

  # apply(subst_EY1W -subst_EY0W , 2, function(p) {mean((p - CATE)^2)})
  # apply(subst_compare_trained$predict(taskY1) - subst_compare_trained$predict(taskY0), 2, function(p) {mean((p - CATE)^2)})
  t <- Sys.time()
  fit_npcausalML <- EP_learn(CATE_library,V = as.data.frame(W), A = A, Y = Y, EY1W = EY1W_est  , EY0W = EY0W_est  , pA1W = pA1W_est, sieve_basis_generator_list = sieve_list ,EP_learner_spec = EP_learner_spec_CATE, cross_validate = TRUE, nfolds = 5)
  print(  Sys.time() - t)

  t <- Sys.time()
  # fit_npcausalML <- npcausalML(CATE_library,
  #                              W= W, A = A, Y = Y, V = data.frame(W1 = W$W1),
  #                              EY1W = EY1W_est, EY0W = EY0W_est,  pA1W = pA1W_est,
  #                              sl3_Learner_EYAW = NULL, sl3_Learner_pA1W = NULL, outcome_type = "continuous", list_of_sieves = sieve_list,
  #                              outcome_function_plugin = outcome_function_plugin_CATE, weight_function_plugin = weight_function_plugin_CATE,
  #                              outcome_function_IPW = outcome_function_IPW_CATE, weight_function_IPW = weight_function_IPW_CATE,
  #                              design_function_sieve_plugin = design_function_sieve_plugin_CATE,
  #                              weight_function_sieve_plugin = weight_function_sieve_plugin_CATE,
  #                              design_function_sieve_IPW = design_function_sieve_IPW_CATE, weight_function_sieve_IPW = weight_function_sieve_IPW_CATE, transform_function = function(x){x},
  #                              family_risk_function = gaussian(),
  #                              efficient_loss_function = efficient_loss_function_CATE,
  #                              use_sieve_selector = FALSE,
  #                              cross_validate_ERM = T, folds = origami::folds_vfold(length(A), 5))
  #
  #
  # print(  Sys.time() - t)
  #preds <- predict(fit_npcausalML,  data.frame(W1 = W$W1), F)
  preds <- fit_npcausalML$full_predictions

  # Compute least-squares risk of predictions using oracle loss function.
  risks_oracle <- as.vector(apply(preds, 2, function(theta) {
    mean((theta -  CATE)^2)
  }) )
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
  gam_keep <- tmp[grep("gam", lrnrs_full),risks_oracle[which.min(risk)], by = degree]$V1
  names(gam_keep) <- paste0("Lrnr_gam_cv", "_fourier_basis_", uniq_degrees, "_plugin")
  risks_oracle <- c(risks_oracle, gam_keep)


  gam_keep <- tmp[grep("gam", lrnrs_full),cvrisksDR[which.min(risk)], by = degree]$V1
  names(gam_keep) <- paste0("Lrnr_gam_cv", "_fourier_basis_", uniq_degrees, "_plugin")
  cvrisksDR <- c(cvrisksDR, gam_keep)

  gam_keep <- tmp[grep("gam", lrnrs_full),cvrisksDRoracle[which.min(risk)], by = degree]$V1
  names(gam_keep) <- paste0("Lrnr_gam_cv", "_fourier_basis_", uniq_degrees, "_plugin")
  cvrisksDRoracle <- c(cvrisksDRoracle, gam_keep)

  sieve_names <- c(colnames(fit_npcausalML$cv_predictions),   names(gam_keep))

  CATE_library <- c(CATE_library, list(lrnr_gam_sl))
  CATEonestepbench <- DR_learner(CATE_library_subst, as.data.frame(W), A, Y, EY1W_est, EY0W_est, pA1W_est, NULL, NULL)
  CATEonestepbench <- apply(CATEonestepbench, 2, function(pred) {
    mean((pred - CATE)^2)
  })
  names(CATEonestepbench) <- c(unique(tmp$lrnrs) , "Lrnr_gam_cv")

  CATEonestepbenchoracle <- DR_learner(CATE_library, as.data.frame(W), A, Y, EY1Wtrue, EY0Wtrue, pA1Wtrue, NULL, NULL)
  CATEonestepbenchoracle <- apply(CATEonestepbenchoracle, 2, function(pred) {
    mean((pred - CATE)^2)
  })
  names(CATEonestepbenchoracle) <- c(unique(tmp$lrnrs), "Lrnr_gam_cv")



  risk_subst<-  apply(subst_CATE, 2, function(pred) {
    mean((pred - CATE)^2)
  })
  risk_subst_cv <- mean((EY1W_est - EY0W_est - CATE)^2)




  list(risk_subst_cv = risk_subst_cv,  risk_subst = risk_subst, CATEonestepbenchoracle =CATEonestepbenchoracle, CATEonestepbench = CATEonestepbench, sieve =data.frame(sieve_names, cvrisksDRoracle, cvrisksDR, risks_oracle))
}

hard <- hard == "TRUE"
pos <- pos == "TRUE"
n <- as.numeric(n)
simresults <- lapply(1:nsims, function(i){try({
  print(i)
  onesim(n)
})
})
save(simresults, file = paste0("mainSimResults2/","simsCATE", hard,pos, "n", n, "_gam"))

