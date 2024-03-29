library(SuperLearner)

library(future)
library(npcausalML)
source("simRR.R")

library(earth)
SL.gam1 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 1
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam2 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 3
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam3 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 5
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam4 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 7
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam5 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 9
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
list_of_sieves_uni   <- list(
  "no_sieve" = NULL,
  fourier_basis$new(orders = c(1,0)),
  fourier_basis$new(orders = c(2,0)),
  fourier_basis$new(orders = c(3,0)),
  fourier_basis$new(orders = c(4,0))
)

lrnr_gam1 <- Lrnr_pkg_SuperLearner$new("SL.gam1" , name = "Lrnr_gam_s1_x", family = binomial())
lrnr_gam2 <- Lrnr_pkg_SuperLearner$new("SL.gam2", name = "Lrnr_gam_s3_x", family = binomial())
lrnr_gam3 <- Lrnr_pkg_SuperLearner$new("SL.gam3", name = "Lrnr_gam_s5_x", family = binomial())
lrnr_gam4 <- Lrnr_pkg_SuperLearner$new("SL.gam4" , name = "Lrnr_gam_s7_x", family = binomial())
lrnr_gam5 <- Lrnr_pkg_SuperLearner$new("SL.gam5" , name = "Lrnr_gam_s9_x", family = binomial())


 hard <- T
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

  LRR <- log(EY1Wtrue/EY0Wtrue)



  # sieve method
  lrnr_Y <- make_learner(Pipeline, Lrnr_cv$new(Stack$new(
    Lrnr_xgboost$new(max_depth =2, verbosity = 0, nrounds = 10, objective = "binary:logistic"),
    Lrnr_xgboost$new(max_depth =3, verbosity = 0, nrounds = 10, objective = "binary:logistic"),
    Lrnr_xgboost$new(max_depth =4, verbosity = 0, nrounds = 10, objective = "binary:logistic"),
    Lrnr_xgboost$new(max_depth =5, verbosity = 0, nrounds = 10, objective = "binary:logistic"),
    Lrnr_xgboost$new(max_depth =6, verbosity = 0, nrounds = 10, objective = "binary:logistic")
))
    , Lrnr_cv_selector$new(loss_squared_error))


  lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(
    Stack$new(


      Lrnr_xgboost$new(max_depth =2, verbosity = 0, nrounds = 10),
      Lrnr_xgboost$new(max_depth =3, verbosity = 0, nrounds = 10),
      Lrnr_xgboost$new(max_depth =4, verbosity = 0, nrounds = 10),
      Lrnr_xgboost$new(max_depth =5, verbosity = 0, nrounds = 10),
      Lrnr_xgboost$new(max_depth =6, verbosity = 0, nrounds = 10)
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


  lrnr_gam <- list(   lrnr_gam1, lrnr_gam2, lrnr_gam3, lrnr_gam4, lrnr_gam5 )
  lrnr_gam_sl <-  Lrnr_sl$new(lrnr_gam, metalearner = Lrnr_cv_selector$new(loss_squared_error))

  lrnr_lm <- list(   Lrnr_glm$new(family = binomial()))

  LRR_library <- c(lrnr_gam,  lrnr_lm  )


  IPW_learner <- Stack$new(c(LRR_library, list(lrnr_gam_sl)))
  IPW_learner <- IPW_learner$train(sl3_Task_IPW)
  preds_IPW <- apply(IPW_learner$predict(sl3_Task_IPW), 2, qlogis)

  lrnr_gam_plugin <- list(   lrnr_gam1, lrnr_gam2, lrnr_gam3, lrnr_gam4, lrnr_gam5 )
  lrnr_gam_sl_plugin <-  Lrnr_sl$new(lrnr_gam_plugin, metalearner = Lrnr_cv_selector$new(loss_squared_error))
  lrnr_lm_plugin <- list(  Lrnr_glm$new(family = binomial()))

  LRR_library_plugin <- c(lrnr_gam_plugin, lrnr_lm_plugin, list(lrnr_gam_sl_plugin))

  subst_compare <-Stack$new(c(LRR_library, list(lrnr_gam_sl)))
  subst_EY1W_trained <-subst_compare$train(taskY1[A==1]$next_in_chain(covariates = c("W")))
  subst_EY0W_trained <- subst_compare$train(taskY0[A==0]$next_in_chain(covariates = c("W")))

  subst_EY1W <-subst_EY1W_trained$predict(taskY1$next_in_chain(covariates = c("W")))
  subst_EY0W <- subst_EY0W_trained$predict(taskY0$next_in_chain(covariates = c("W")))


  subst_EY1W <- pmax(subst_EY1W, 1e-5)
  subst_EY0W <- pmax(subst_EY0W, 1e-5)
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
                               family_risk_function = binomial(),
                               efficient_loss_function = efficient_loss_function_LRR,
                               use_sieve_selector = FALSE,
                               transform_function =qlogis,
                               cross_validate_ERM = T, folds = origami::folds_vfold(length(A), 8))


  preds <- predict(fit_npcausalML,  data.frame(W = W$W), F)


  # Compute least-squares risk of predictions using oracle loss function.
  risks_oracle <- as.vector(apply(preds, 2, function(theta) {
    mean((theta -  LRR)^2)
  })[grep("plugin", colnames(preds))])

  # Compute estimated cross-validated one-step risk of predictions
  cvrisksDR <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {

    loss <-  efficient_loss_function_LRR(W, theta, A, Y, EY1W_est,EY0W_est, pA1W_est )
    colMeans(loss)

  }))
  names(cvrisksDR) <- colnames(fit_npcausalML$cv_predictions)

  # Compute estimated cross-validated oracle one-step risk of predictions

  cvrisksDRoracle <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
    loss <- efficient_loss_function_LRR(W, theta, A, Y, EY1Wtrue,EY0Wtrue, pA1Wtrue )
    colMeans(loss)
  }))

  lrnrs_full <-  colnames(fit_npcausalML$cv_predictions)
  lrnrs <- gsub("[._]fourier.+", "", lrnrs_full)
  lrnrs <- gsub("[._]no_sieve.+", "", lrnrs)
  degree <- as.numeric(stringr::str_match(lrnrs_full, "fourier_basis_([0-9]+)")[,2])
  degree[grep("no_sieve", lrnrs_full)] <- 0


  tmp <- data.table(lrnrs_full, lrnrs , degree, risk = cvrisksDR, risks_oracle = risks_oracle, cvrisksDR = cvrisksDR, cvrisksDRoracle)
  gam_keep <- tmp[grep("gam", lrnrs_full),risks_oracle[which.min(risk)], by = degree]$V1
  names(gam_keep) <- paste0("Lrnr_gam_cv", "_fourier.basis_", 0:4, "_plugin")
  risks_oracle <- c(risks_oracle, gam_keep)
  sieve_names <- c(colnames(fit_npcausalML$cv_predictions),   names(gam_keep))


  gam_keep <- tmp[grep("gam", lrnrs_full),cvrisksDR[which.min(risk)], by = degree]$V1
  names(gam_keep) <- paste0("Lrnr_gam_cv", "_fourier.basis_", 0:4, "_plugin")
  cvrisksDR <- c(cvrisksDR, gam_keep)

  gam_keep <- tmp[grep("gam", lrnrs_full),cvrisksDRoracle[which.min(risk)], by = degree]$V1
  names(gam_keep) <- paste0("Lrnr_gam_cv", "_fourier.basis_", 0:4, "_plugin")
  cvrisksDRoracle <- c(cvrisksDRoracle, gam_keep)






  risk_subst<-  apply(subst_LRR, 2, function(pred) {
    mean((pred - LRR)^2)
  })

  risk_IPW <-  apply(preds_IPW, 2, function(pred) {
    mean((pred - LRR)^2)
  })


  list( risk_IPW = risk_IPW,  risk_subst = risk_subst,    sieve =data.frame(sieve_names, cvrisksDRoracle, cvrisksDR, risks_oracle))
}

nsims <- 25
print(250)
simresults250 <- lapply(1:nsims, function(i){
  print(i)
  onesim(250)
})

save(simresults250, file = paste0("mainSimResults/","simsLRR_GAM_", hard,pos, "n250"))


print(500)
simresults500 <- lapply(1:nsims, function(i){
  print(i)
  onesim(500)
})

save(simresults500, file = paste0("mainSimResults/","simsLRR_GAM_", hard,pos, "n500"))

print(1000)
simresults1000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(1000)
})

save(simresults1000, file = paste0("mainSimResults/","simsLRR_GAM_", hard,pos, "n1000"))

print(2500)
simresults2500 <- lapply(1:nsims, function(i){
  print(i)
  onesim(2500)
})

save(simresults2500, file = paste0("mainSimResults/","simsLRR_GAM_", hard,pos, "n2500"))


print(5000)
simresults5000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(5000)
})

save(simresults5000, file = paste0("mainSimResults/","simsLRR_GAM_", hard,pos, paste0("n5000")))
