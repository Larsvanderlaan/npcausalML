

lrnr_stack_LRR <- list(
  Lrnr_glm$new(family = binomial()),
  Lrnr_earth$new(family = binomial()),
  Lrnr_gam$new(family = binomial())
  ,Lrnr_xgboost$new(max_depth = 3, nrounds = 1, num_parallel_tree = 1),
  Lrnr_xgboost$new(max_depth = 3, nrounds = 1,
    'colsample_bynode'= 0.8,
    'learning_rate'= 1,
    'max_depth'= 5,
    'num_parallel_tree'= 100,
    'objective'= 'reg:logistic',
    'subsample'= 0.8
  ),
  Lrnr_xgboost$new(max_depth = 3, 'objective'= 'reg:logistic'),
  Lrnr_xgboost$new(max_depth = 4, 'objective'= 'reg:logistic'),
  Lrnr_xgboost$new(max_depth = 5, 'objective'= 'reg:logistic')
)



# Sims CATE
## no positivity
## easy RR
sim.RR <- function(n, hard = TRUE, positivity = TRUE) {
  if(!positivity & !hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    LRR <- -1 + W
    EY0W <- plogis(-1+W + sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rpois(n, EY0W * exp(A* LRR))
  }


  ## hard RR
  if(!positivity & hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    LRR <- -1 + W + sin(5*W)
    EY0W <- plogis(-1+W + sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rpois(n, EY0W * exp(A* LRR))
  }
  ##  positivity
  ## easy RR
  if(positivity & !hard) {

    W <- runif(n, -1 , 1)
    pA1W <- plogis(3*W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    LRR <- -1 + W
    EY0W <- plogis(-1+W + sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rpois(n, EY0W * exp(A* LRR))
  }
  ## hard RR
  if(positivity & hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(3*W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    LRR <- -1 + W + sin(5*W)
    EY0W <- plogis(-1+W + sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rpois(n, EY0W * exp(A* LRR))
  }
  return(data.table(W, A, Y, pA1W, EY1W, EY0W))
}





estRR <- function(n, LRR_library, nsims, sim_function, sieve_list, ...) {


  library(sl3)
  list_RRIPW <- list()
  list_RRplugin <- list()
  list_RRIPW_oracle <- list()
  list_RRsubst <- list()
  list_causalforest <- list()
  for(i in 1:nsims) {
    try({
    print(i)
    data <- sim_function(n, ...)
    W <- data$W
    W <- as.matrix(data.frame(W = W))
    A <- data$A
    Y <- data$Y
    EY1Wtrue <- data$EY1W
    EY0Wtrue <- data$EY0W
    pA1Wtrue <- data$pA1W
    LRR <- log(EY1Wtrue/ EY0Wtrue)

    # sieve method
    lrnr_Y <- make_learner(Pipeline, Lrnr_cv$new(
      Stack$new(
        Lrnr_stratified$new(Lrnr_glm$new(family = poisson()), variable_stratify = "A"),
        Lrnr_stratified$new(Lrnr_gam$new(family = poisson()), variable_stratify = "A"),
        #Lrnr_earth$new(family = poisson()),
        Lrnr_xgboost$new(max_depth = 3, objective = "count:poisson"),
        Lrnr_xgboost$new(max_depth = 4, objective = "count:poisson"),
        Lrnr_xgboost$new(max_depth = 5, objective = "count:poisson")
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


    fit_npcausalML <- npcausalML(LRR_library,
                                 W= as.matrix(data.table(W)), A = A, Y = Y, V = as.matrix(as.data.table(W)), sl3_Learner_EYAW = lrnr_Y, sl3_Learner_pA1W = lrnr_A, outcome_type = "continuous", list_of_sieves = sieve_list,cross_validate = FALSE,
                                 outcome_function_plugin = outcome_function_plugin_LRR, weight_function_plugin = weight_function_plugin_LRR,
                                 outcome_function_IPW = outcome_function_IPW_LRR, weight_function_IPW = weight_function_IPW_LRR,
                                 design_function_sieve_plugin = design_function_sieve_plugin_LRR,
                                 weight_function_sieve_plugin = weight_function_sieve_plugin_LRR,
                                 design_function_sieve_IPW = design_function_sieve_IPW_LRR, weight_function_sieve_IPW = weight_function_sieve_IPW_LRR, transform_function = function(x){x},
                                 family_risk_function = binomial(),
                                 efficient_loss_function = efficient_loss_function_LRR, use_sieve_selector = TRUE)

    preds <- predict(fit_npcausalML, W)

    EY1W_est <- fit_npcausalML$EY1W
    EY0W_est <- fit_npcausalML$EY0W
    pA1W_est <- fit_npcausalML$pA1W





    LRRIPW <- IPW_learner(LRR_library, W, A, Y, EY1W_est, EY0W_est, pA1W_est, NULL, NULL)
    LRRIPWoracle <-IPW_learner(LRR_library, W, A, Y, EY1Wtrue, EY0Wtrue, pA1Wtrue, NULL, NULL)

    risks_IPW <- apply(LRRIPW, 2, function(est) {
      mean((LRR - est)^2)
    })

    risks_IPW_oracle <- apply(LRRIPWoracle, 2, function(est) {
      mean((LRR - est)^2)
    })

    risk_plugin <- apply(preds, 2, function(est) {
      mean((LRR - est)^2)
    })



    risk_subst<-  mean((LRR - (log(EY1W_est)  - log(EY0W_est)))^2)



    list_RRIPW[[i]] <- risks_IPW
    list_RRIPW_oracle[[i]] <- risks_IPW_oracle
    list_RRplugin[[i]] <- risk_plugin
    list_RRsubst[[i]] <- risk_subst
    })
  }
  RRIPW <- rowMeans(do.call(cbind, list_RRIPW))
  RRIPW_oracle <- rowMeans(do.call(cbind, list_RRIPW_oracle))
  RRplugin <- rowMeans(do.call(cbind, list_RRplugin))
  subst <- mean(unlist(list_RRsubst))

  return(data.table(subst,  RRIPW_oracle, RRIPW, RRplugin, lrnr = names(RRplugin)))




}




IPW_learner <- function(LRR_library, W, A, Y, EY1W, EY0W, pA1W, lrnr_A, lrnr_Y) {
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
  pA0W <- 1- pA1W
  Ytilde <- A
  weights = Y/ifelse(A==1, pA1W, pA0W)
  data <- data.table(W, Ytilde, weights = weights)
  task_IPW <- sl3_Task$new(data, covariates = c("W"), outcome = "Ytilde", weights = "weights")
  lrnr <- Stack$new(LRR_library)
  lrnr_1 <- lrnr$train(task_IPW)
  CATEIPW <- lrnr_1$predict(task_IPW)
  return(CATEIPW)
}



