

sim.CATE <- function(n, hard = TRUE, positivity = TRUE, ...) {
  if(!positivity & !hard) {
    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- plogis((W1 + W2 + W3)/3)
    quantile(pA1W)
    A <- rbinom(n, 1 ,  pA1W)
    CATE <- 1 + W1
    EY0W <- 0.5*(W1 + W2 + W3) + sin(5*W1) + sin(5*W2) + sin(5*W3) + 1/(W1 + 1.2) + 1/(W2 + 1.2) + 1/(W3 + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 2)
  }

  ## hardCATE
  if(!positivity & hard) {


    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- plogis((W1 + W2 + W3)/3)
    quantile(pA1W)
    A <- rbinom(n, 1 ,  pA1W)
    CATE <- 1 + W1 + sin(5*W1)
    EY0W <- 0.5*(W1 + W2 + W3) + sin(5*W1) + sin(5*W2) + sin(5*W3) + 1/(W1 + 1.2) + 1/(W2 + 1.2) + 1/(W3 + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 2)


  }
  ##  positivity
  ## easy CATE
  if(positivity & !hard) {
    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- plogis((W1 + W2 + W3))
    quantile(pA1W)
    A <- rbinom(n, 1 ,  pA1W)
    CATE <- 1 + W1
    EY0W <- 0.5*(W1 + W2 + W3) + sin(5*W1) + sin(5*W2) + sin(5*W3) + 1/(W1 + 1.2) + 1/(W2 + 1.2) + 1/(W3 + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 2)
  }
  ## hardCATE
  if(positivity & hard) {
    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- plogis((W1 + W2 + W3))
    quantile(pA1W)
    A <- rbinom(n, 1 ,  pA1W)
    CATE <- 1 + W1 + sin(5*W1)
    EY0W <- 0.5*(W1 + W2 + W3) + sin(5*W1) + sin(5*W2) + sin(5*W3) + 1/(W1 + 1.2) + 1/(W2 + 1.2) + 1/(W3 + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 2)

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

list_of_sieves_uni   <- list(
  NULL,
  fourier_basis$new(orders = c(1,0)),
  fourier_basis$new(orders = c(2,0)),
  fourier_basis$new(orders = c(3,0)),
  fourier_basis$new(orders = c(4,0)),
  fourier_basis$new(orders = c(5,0))
)

#out <- estCATE(250,  lrnr_stack, 1, sim.CATE, list_of_sieves_uni, positivity = FALSE, hard = FALSE)


estCATE <- function(n, CATE_library, nsims, sim_function, sieve_list, ...) {


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
        Lrnr_hal9001$new(smoothness_orders = 0, max_degree = 2, num_knots = c(30,30))
      )), Lrnr_cv_selector$new(loss_squared_error))


    lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(
      Stack$new(
        Lrnr_hal9001$new(smoothness_orders = 0, max_degree = 2, num_knots = c(30,30))
      )
    ), Lrnr_cv_selector$new(loss_squared_error))


    fit_npcausalML <- npcausalML(CATE_library,
                                 W= as.matrix(data.table(W)), A = A, Y = Y, V = as.matrix(as.data.table(W)), sl3_Learner_EYAW = lrnr_Y, sl3_Learner_pA1W = lrnr_A, outcome_type = "continuous", list_of_sieves = sieve_list,
                                 outcome_function_plugin = outcome_function_plugin_CATE, weight_function_plugin = weight_function_plugin_CATE,
                                 outcome_function_IPW = outcome_function_IPW_CATE, weight_function_IPW = weight_function_IPW_CATE,
                                 design_function_sieve_plugin = design_function_sieve_plugin_CATE,
                                 weight_function_sieve_plugin = weight_function_sieve_plugin_CATE,
                                 design_function_sieve_IPW = design_function_sieve_IPW_CATE, weight_function_sieve_IPW = weight_function_sieve_IPW_CATE, transform_function = function(x){x},
                                 family_risk_function = gaussian(),
                                 efficient_loss_function = efficient_loss_function_CATE, use_sieve_selector = FALSE,
                                 cross_validate_ERM = FALSE)


    fit_npcausalML_adaptive <- fit_npcausalML
    fit_npcausalML_adaptivetmle <- fit_npcausalML
    fit_npcausalML_oracle <- fit_npcausalML


    #fit_npcausalML_adaptive$learners <- subset_best_sieve(fit_npcausalML$learners, fit_npcausalML$learner_names, A, Y, fit_npcausalML$EY1W, fit_npcausalML$EY0W, fit_npcausalML$pA1W, rep(1, n), efficient_loss_function_CATE, V = W)


    oracle_risk <- function(V, theta, A , Y, EY1W , EY0W , pA1W, oracle = FALSE){
       theta^2 - 2*theta*CATE
    }
    print("oracle")
    fit_npcausalML_oracle$learners <- subset_best_sieve(fit_npcausalML$learners, fit_npcausalML$learner_names, A, Y,
                                             fit_npcausalML$EY1W, fit_npcausalML$EY0W, fit_npcausalML$pA1W, rep(1,n),
                                             oracle_risk, V = W)
    print("ADAPTIVEtmle")
    fit_npcausalML_adaptivetmle$learners <- subset_best_sieve(fit_npcausalML$learners, fit_npcausalML$learner_names, A, Y,
                                                          fit_npcausalML$EY1W, fit_npcausalML$EY0W, fit_npcausalML$pA1W, rep(1,nrow(data)),
                                                        function(V, theta, A , Y, EY1W , EY0W , pA1W, oracle = FALSE){
                                                          pA <-  ifelse(A==1,pA1W, 1- pA1W)
                                                          EY <- ifelse(A==1, EY1W, EY0W)
                                                          CATE <- EY1W - EY0W
                                                          beta <- coef(lm(Y ~ H - 1, offset = EY, data = list(Y = Y, H = theta*(A - (1-A))), weights = 1/pA))
                                                          EY1W <- EY1W +  theta*beta
                                                          EY0W <- EY0W - theta*beta
                                                          CATE <- EY1W - EY0W
                                                          EY <- ifelse(A==1, EY1W, EY0W)
                                                          #print("score")
                                                          #print(mean(theta*(1/pA)*(A - (1-A))*(Y - EY)))
                                                          #loss <- theta^2 - theta*(CATE + (1/pA)*(A - (1-A))*(Y - EY))
                                                          loss <- (theta - (CATE))^2
                                                        }, V = data$W)

    print("ADAPTIVE")
    fit_npcausalML_adaptive$learners <- subset_best_sieve(fit_npcausalML$learners, fit_npcausalML$learner_names, A, Y,
                                                          fit_npcausalML$EY1W, fit_npcausalML$EY0W, fit_npcausalML$pA1W, rep(1,nrow(data)),
                                                          function(V, theta, A , Y, EY1W , EY0W , pA1W, oracle = FALSE){
                                                            pA <-  ifelse(A==1,pA1W, 1- pA1W)
                                                            EY <- ifelse(A==1, EY1W, EY0W)
                                                            CATE <- EY1W - EY0W

                                                            loss <- theta^2 - 2*theta*(CATE + (1/pA)*(A - (1-A))*(Y - EY))

                                                          }, V = data$W)
    # print("ORACLEDR")
    # data_val <- sim_function(n, ...)
    # preds_all <- predict(fit_npcausalML, data.frame(W=data_val$W), subset_cv = FALSE)
    #
    # risks <- apply(preds_all, 2, efficient_risk_function, A = data_val$A, Y = data_val$Y,
    #       EY1W = data_val$EY1W, EY0W = data_val$EY0W, pA1W = data_val$pA1W, weights = rep(1, nrow(data_val)),
    #      efficient_loss_function = efficient_loss_function_CATE)
    # print(as.vector(risks))
    # print(which.min(as.vector(risks)))
    #
    # print("leakage")
    # preds_all <- predict(fit_npcausalML, data.frame(W=data$W), subset_cv = FALSE)
    #
    # risks <- apply(preds_all, 2, efficient_risk_function, A = data$A, Y = data$Y,
    #                EY1W = data$EY1W, EY0W = data$EY0W, pA1W = data$pA1W, weights = rep(1, n),
    #                efficient_loss_function = efficient_loss_function_CATE)
    # print(as.vector(risks))
    # print(which.min(as.vector(risks)))



    preds_adaptive <- as.matrix(predict(fit_npcausalML_adaptive, W))
    preds_adaptivetmle <- as.matrix(predict(fit_npcausalML_adaptivetmle, W))

    preds_oracle <- as.matrix(predict(fit_npcausalML_oracle, W))


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

  return(data.table(subst, causalforest,  cateonestep_oracle, cateonestep, cateplugin_adaptive, cateplugin_adaptivetmle, cateplugin_oracle, lrnr = names(cateplugin_adaptive)))




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

