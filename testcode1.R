n <- 2500
sieve_list <- list_of_sieves_uni
library(sl3)
list_cateonestep <- list()
list_cateplugin_adaptive <- list()
list_cateplugin_adaptivetmle <- list()
list_cateplugin_oracle <- list()
list_cateonestep_oracle <- list()
list_catesubst <- list()
list_causalforest <- list()
for(i in 1:nsims) {
  print(i)
  data <- sim.CATE(n, F, T)
  W <- data$W
  W <- as.matrix(data.frame(W = W))
  A <- data$A
  Y <- data$Y
  EY1Wtrue <- data$EY1W
  EY0Wtrue <- data$EY0W
  pA1Wtrue <- data$pA1W
  EYWtrue <- ifelse(A==1, EY1Wtrue, EY0Wtrue)

  CATE <- EY1Wtrue - EY0Wtrue



  # sieve method
  lrnr_Y <- make_learner(Pipeline, Lrnr_cv$new(
    Stack$new(
      Lrnr_hal9001$new(smoothness_orders = 0, max_degree = 2, num_knots = c(10,10))
    )), Lrnr_cv_selector$new(loss_squared_error))


  lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(
    Stack$new(
      Lrnr_hal9001$new(smoothness_orders = 0, max_degree = 2, num_knots = c(10,10))
    )
  ), Lrnr_cv_selector$new(loss_squared_error))

  initial_likelihood <- estimate_initial_likelihood(W, A, Y, weights = rep(1,n), lrnr_A, lrnr_Y, folds = 10)
  EY1W_est <- initial_likelihood$EY1W
  EY0W_est <- initial_likelihood$EY0W
  pA1W_est <- initial_likelihood$pA1W
  EYW_est <- ifelse(A==1, EY1W_est, EY0W_est)

  CATE_library <- list(Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 1, num_knots = 10, lambda = 1e-10, fit_control = list(cv_select = F)))
  fit_npcausalML <- npcausalML(CATE_library,
                               W= as.matrix(data.table(W)), A = A, Y = Y, V = as.matrix(as.data.table(W)),
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
                               cross_validate_ERM = TRUE)


   preds <- predict(fit_npcausalML, data.frame(W1=W$W1), F)

   mean((EY1W_est - EY0W_est - CATE)^2)
   mean((preds[,1] - CATE)^2)
   Ypseudo <- EY1W_est - EY0W_est
   task_tmp <- sl3_Task$new(data.table(W1, Ypseudo), covariates = "W1", outcome = "Ypseudo")
   plugin <- Lrnr_gam$new()$train(task_tmp)$predict(task_tmp)
   mean((plugin - CATE)^2)


   risks_oracle <- as.vector(apply(preds, 2, function(theta) {
     mean((theta -  CATE)^2)
   })[grep("plugin", colnames(preds))])
   as.vector(risks_oracle)
   plot(0:10, risks_oracle, type = "l")


   risks_oracle2 <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
     mean(theta ^2-  2*theta*CATE)
   })[-grep("IPW", colnames(preds))])
   as.vector(risks_oracle2)
   plot(0:10, risks_oracle2, type = "l")





   risksDR <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
     loss <- efficient_loss_function_CATE(W, theta, A, Y, EY1W_est,EY0W_est, pA1W_est )
     mean(loss)
   })[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])

   plot(0:10, risksDR, type = "l")

   risksDR2 <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
     loss <- efficient_loss_function_CATE(W, theta, A, Y, EY1Wtrue,EY0Wtrue, pA1Wtrue )

     mean(loss)
   })[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])


   plot(0:10, risksDR2, type = "l")


   CATEonestep <- DR_learner(CATE_library, W1, A, Y, EY1W_est, EY0W_est, pA1W_est, NULL, NULL)
   c(risks_oracle, mean((CATEonestep[[1]] - CATE)^2))

   risksDR <- as.vector(apply(cbind(preds, CATEonestep), 2, function(theta) {
     loss <- efficient_loss_function_CATE(W, theta, A, Y, EY1W_est,EY0W_est, pA1W_est )
     mean(loss)
   })[-grep("IPW", colnames(preds))])

   plot(0:11, risksDR, type = "l")




   data_val <-  sim.CATE(n, F, F)
   Wval <- data_val$W
   Wval <- as.matrix(data.frame(W = Wval))
   Aval <- data_val$A
   Yval <- data_val$Y
   EY1Wtrue_val <- data_val$EY1W
   EY0Wtrue_val <- data_val$EY0W
   pA1Wtrue_val <- data_val$pA1W
   CATE_val <- EY1Wtrue_val - EY0Wtrue_val
   preds_val <- predict(fit_npcausalML, Wval, F)
   risksDRoracle <- as.vector(apply(preds_val, 2, function(theta) {
     loss <- efficient_loss_function_CATE(Wval, theta, Aval, Yval, EY1Wtrue_val,EY0Wtrue_val, pA1Wtrue_val )
     mean(loss)
   })[grep("plugin", colnames(preds_val))])
   plot(0:10, risksDRoracle, type = "l")



  CATEonestep <- DR_learner(CATE_library, W, A, Y, EY1W_est, EY0W_est, pA1W_est, NULL, NULL)
  CATEonesteporacle <-DR_learner(CATE_library, W, A, Y, EY1Wtrue, EY0Wtrue, pA1Wtrue, NULL, NULL)

  risks_onestep <- apply(CATEonestep, 2, function(est) {
    mean((CATE - est)^2)
  })

  risks_onestep_oracle <- apply(CATEonesteporacle, 2, function(est) {
    mean((CATE - est)^2)
  })

  risk_plugin_adaptive <- apply(preds_adaptive, 2, function(est) {
    mean((CATE - est)^2)
  })
  risk_plugin_adaptivetmle <- apply(preds_adaptivetmle, 2, function(est) {
    mean((CATE - est)^2)
  })
  risk_plugin_oracle <- apply(preds_oracle, 2, function(est) {
    mean((CATE - est)^2)
  })



  risk_subst<-  mean((CATE - (EY1W_est  - EY0W_est))^2)

  Y.hat <- EY1W_est * pA1W_est + EY0W_est * (1-pA1W_est)
  W.hat <- pA1W_est
  fit <- grf::causal_forest(X = W, Y  = Y, W = A, Y.hat = Y.hat, W.hat = W.hat)
  preds_cf <-  fit$predictions
  risk_cf <- mean((CATE - preds_cf)^2)

  list_causalforest[[i]] <- risk_cf
  list_cateonestep[[i]] <- risks_onestep
  list_cateonestep_oracle[[i]] <- risks_onestep_oracle
  list_cateplugin_oracle[[i]] <- risk_plugin_oracle
  list_cateplugin_adaptive[[i]] <- risk_plugin_adaptive
  list_cateplugin_adaptivetmle[[i]] <- risk_plugin_adaptivetmle
  list_catesubst[[i]] <- risk_subst
}
cateonestep <- rowMeans(do.call(cbind, list_cateonestep))
cateonestep_oracle <- rowMeans(do.call(cbind, list_cateonestep_oracle))
cateplugin_oracle <- rowMeans(do.call(cbind, list_cateplugin_oracle))
cateplugin_adaptive <- rowMeans(do.call(cbind, list_cateplugin_adaptive))
cateplugin_adaptivetmle <- rowMeans(do.call(cbind, list_cateplugin_adaptivetmle))

causalforest <- rowMeans(do.call(cbind, list_causalforest))
subst <- mean(unlist(list_catesubst))
