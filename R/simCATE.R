

sim.CATE <- function(n, hard = TRUE, positivity = TRUE, ...) {
  if(!positivity & !hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    CATE <- 1 + W
    EY0W <- W + sin(5*W) + 1/(W + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 0.5)
  }

  ## hardCATE
  if(!positivity & hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    CATE <- 1 + W + sin(5*W)
    EY0W <- W + sin(5*W) + 1/(W + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 0.5)
  }
  ##  positivity
  ## easy CATE
  if(positivity & !hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(3*W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    CATE <- 1 + W
    EY0W <- W + sin(5*W) + 1/(W + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 0.5)
  }
  ## hardCATE
  if(positivity & hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(3*W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(pA1W)
    CATE <- 1 + W + sin(5*W)
    EY0W <- W + sin(5*W) + 1/(W + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 0.5)

  }
  return(data.table(W, A, Y, EY0W, EY1W, pA1W))
}


library(sl3)
library(data.table)
lrnr_stack <- list(
  Lrnr_glm$new(),
  Lrnr_earth$new(),
  Lrnr_gam$new(),
  Lrnr_rpart$new(),
  Lrnr_ranger$new(),
  Lrnr_xgboost$new(max_depth = 3),
  Lrnr_xgboost$new(max_depth = 4),
  Lrnr_xgboost$new(max_depth = 5)
)








list_of_sieves <- list_of_sieves <- list(
  NULL,
  fourier_basis$new(orders = c(1,0)),
  fourier_basis$new(orders = c(2,0)),
  fourier_basis$new(orders = c(1,1)),
  fourier_basis$new(orders = c(2,1)),
  fourier_basis$new(orders = c(3,2))
)

list_of_sieves_uni <- list_of_sieves <- list(
  NULL,
  fourier_basis$new(orders = c(1,0)),
  fourier_basis$new(orders = c(2,0)),
  fourier_basis$new(orders = c(3,0)),
  fourier_basis$new(orders = c(4,0)),
  fourier_basis$new(orders = c(5,0))
)

out <- estCATE(250,  lrnr_stack, 1, sim.CATE, list_of_sieves_uni, positivity = FALSE, hard = FALSE)


estCATE <- function(n, CATE_library, nsims, sim_function, sieve_list, ...) {


  library(sl3)
  list_cateonestep <- list()
  list_cateplugin <- list()
  list_cateonestep_oracle <- list()
  list_catesubst <- list()
  list_causalforest <- list()
  for(i in 1:nsims) {
    data <- sim_function(n, ...)
    W <- data$W
    W <- as.matrix(data.frame(W = W))
    A <- data$A
    Y <- data$Y
    EY1Wtrue <- data$EY1W
    EY0Wtrue <- data$EY0W
    pA1Wtrue <- data$pA1W
    CATE <- EY1Wtrue - EY0Wtrue

    # sieve method
    lrnr_Y <- make_learner(Pipeline, Lrnr_cv$new(
      Stack$new(
        Lrnr_stratified$new(Lrnr_glm$new(), variable_stratify = "A"),
        Lrnr_stratified$new(Lrnr_gam$new(), variable_stratify = "A"),
        Lrnr_earth$new(),
        Lrnr_xgboost$new(max_depth = 3),
        Lrnr_xgboost$new(max_depth = 4),
        Lrnr_xgboost$new(max_depth = 5)
      )), Lrnr_cv_selector$new(loss_squared_error))


    lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(
      Stack$new(
        Lrnr_gam$new(),
        Lrnr_earth$new(),
        Lrnr_xgboost$new(max_depth = 3),
        Lrnr_xgboost$new(max_depth = 4),
        Lrnr_xgboost$new(max_depth = 5)
      )
    ), Lrnr_cv_selector$new(loss_squared_error))


    fit_npcausalML <- npcausalML(CATE_library,
                                 W= as.matrix(data.table(W)), A = A, Y = Y, V = as.matrix(as.data.table(W)), sl3_Learner_EYAW = lrnr_Y, sl3_Learner_pA1W = lrnr_A, outcome_type = "continuous", list_of_sieves = sieve_list,cross_validate = FALSE,
                                 outcome_function_plugin = outcome_function_plugin_CATE, weight_function_plugin = weight_function_plugin_CATE,
                                 outcome_function_IPW = outcome_function_IPW_CATE, weight_function_IPW = weight_function_IPW_CATE,
                                 design_function_sieve_plugin = design_function_sieve_plugin_CATE,
                                 weight_function_sieve_plugin = weight_function_sieve_plugin_CATE,
                                 design_function_sieve_IPW = design_function_sieve_IPW_CATE, weight_function_sieve_IPW = weight_function_sieve_IPW_CATE, transform_function = function(x){x},
                                 family_risk_function = gaussian(),
                                 efficient_loss_function = efficient_loss_function_CATE, use_sieve_selector = TRUE)

    preds <- predict(fit_npcausalML, W)
    EY1W_est <- fit_npcausalML$EY1W
    EY0W_est <- fit_npcausalML$EY0W
    pA1W_est <- fit_npcausalML$pA1W





    CATEonestep <- DR_learner(CATE_library, W, A, Y, EY1W_est, EY0W_est, pA1W_est, NULL, NULL)
    CATEonesteporacle <-DR_learner(CATE_library, W, A, Y, EY1Wtrue, EY0Wtrue, pA1Wtrue, NULL, NULL)

    risks_onestep <- apply(CATEonestep, 2, function(est) {
      mean((CATE - est)^2)
    })

    risks_onestep_oracle <- apply(CATEonesteporacle, 2, function(est) {
      mean((CATE - est)^2)
    })

    risk_plugin <- apply(preds, 2, function(est) {
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
    list_cateplugin[[i]] <- risk_plugin
    list_catesubst[[i]] <- risk_subst
  }
  cateonestep <- rowMeans(do.call(cbind, list_cateonestep))
  cateplugin <- rowMeans(do.call(cbind, list_cateplugin))
  subst <- mean(unlist(list_subst))
  initial <- mean(unlist(list_initial))
  print(initial)
  return(data.table(initial = initial, subst = subst, cateonestep, cateplugin, lrnr = names(cateplugin)))




}





DR_learner <- function(CATE_library, W, A, Y, EY1W, EY0W, pA1W, lrnr_A, lrnr_Y) {
  if(missing(EY1W)) {
    data <-data.table(W,A,Y)
    covariates <- colnames(W)
    task_A <- sl3_Task$new(data, covariates = covariates, outcome = "A")
    task_Y <- sl3_Task$new(data, covariates = c(covariates, "A"), outcome = "Y")
    data1 <- data
    data1$A <- 1
    task_Y1 <- sl3_Task$new(data, covariates = c(covariates, "A"), outcome = "Y")
    data0 <- data
    data0$A <- 0
    task_Y0 <- sl3_Task$new(data, covariates = c(covariates, "A"), outcome = "Y")
    lrnr_Y <- lrnr_Y$train(task_Y)
    lrnr_A <- lrnr_A$train(task_A)
    pA1W <- lrnr_A$predict(task_A)
    EY1W <- lrnr_Y$predict(task_Y1)
    EY0W <- lrnr_Y$predict(task_Y0)
    pA1W <- pmax(pmin(pA1, 0.98), 0.02)
  }
  EY <- ifelse(A==1, EY1W, EY0W)
  pA0 <- 1- pA1W
  Ytilde <- EY1W - EY0W + (A/pA1W - (1-A)/pA0)*(Y - EY)

  data <- data.table(W, Ytilde)
  task_onestep <- sl3_Task$new(data, covariates = c("W"), outcome = "Ytilde")
  lrnr <- Stack$new(CATE_library)
  lrnr_1 <- lrnr$train(task_onestep)
  CATEonestep <- lrnr_1$predict(task_onestep)
  return(CATEonestep)
}

