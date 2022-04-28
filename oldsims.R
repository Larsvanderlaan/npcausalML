

simpleparametricCATE <- function(n, CATE_library, nsims) {

  library(sl3)
  list_cateonestep <- list()
  list_cateplugin <- list()
  list_subst <- list()
  list_initial <- list()
  for(i in 1:nsims) {
    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)# rbinom(n, 1 , plogis(W1))
    W <- cbind(W1, W2)
    A <- rbinom(n, 1 , plogis(2*(W1 + W2 )))
    Y <- rnorm(n, W1 + W2 + A*(1 + W1 + W2), 1)
    CATE <- (1 + W1 + W2)
    quantile(plogis(2*(W1 + W2 )))
    data <- data.table(W1,W2,A,Y)

    task_A <- sl3_Task$new(data, covariates = c("W1", "W2"), outcome = "A")
    task_Y <- sl3_Task$new(data, covariates = c("W1", "W2", "A"), outcome = "Y")
    data1 <- data
    data1$A <- 1
    task_Y1 <- sl3_Task$new(data, covariates = c("W1", "W2", "A"), outcome = "Y")
    data0 <- data
    data0$A <- 0
    task_Y0 <- sl3_Task$new(data, covariates = c("W1", "W2", "A"), outcome = "Y")
    lrnr_Y <- Lrnr_glm$new(formula = ~ W1 + W2 + A *(W1 + W2))
    lrnr_A <- Lrnr_glm$new()
    lrnr_Y <- lrnr_Y$train(task_Y)
    lrnr_A <- lrnr_A$train(task_A)
    pA1 <- lrnr_A$predict(task_A)
    EY1 <- lrnr_Y$predict(task_Y1)
    EY0 <- lrnr_Y$predict(task_Y0)
    EY <- ifelse(A==1, EY1, EY0)
    pA1 <- pmax(pmin(pA1, 0.98), 0.02)
    pA0 <- 1- pA1
    Ytilde <- EY1 - EY0 + (A/pA1 - (1-A)/pA0)*(Y - EY)
    quantile(Ytilde)
    beta <- coef(glm(Y~X.W1 + X.W2 - 1, data= data.frame(Y = Y, X = (A - (1-A))*W), offset = EY, weights = ifelse(A==1, 1/pA1, 1/pA0)))

    EY1Star <- EY1 +(1 - (1-1))*W %*% beta
    EY0Star <- EY0 +(0 - (1-0))*W %*% beta
    Ytilde_star =  EY1Star - EY0Star
    data <- data.frame(W1, W2, Ytilde = Ytilde, Ytilde_star =  EY1Star - EY0Star)
    task_onestep <- sl3_Task$new(data, covariates = c("W1", "W2"), outcome = "Ytilde")
    task_plugin <- sl3_Task$new(data, covariates = c("W1", "W2"), outcome = "Ytilde_star")

    lrnr <- CATE_library
    lrnr_1 <- lrnr$train(task_onestep)
    CATEonestep <- lrnr_1$predict(task_onestep)
    cateonestep <- apply(CATEonestep, 2, function(est) {
      mean((CATE - est)^2)
    })

    lrnr_1 <- lrnr$train(task_plugin)
    CATEplugin <- lrnr_1$predict(task_plugin)
    cateplugin <- apply(CATEplugin, 2, function(est) {
      mean((CATE - est)^2)
    })

    list_cateonestep[[i]] <- cateonestep
    list_cateplugin[[i]] <- cateplugin
    list_subst[[i]] <- mean((CATE - (EY1Star - EY0Star))^2)
    list_initial[[i]] <- mean((CATE - (EY1 - EY0))^2)
  }
  cateonestep <- rowMeans(do.call(cbind, list_cateonestep))
  cateplugin <- rowMeans(do.call(cbind, list_cateplugin))
  subst <- mean(unlist(list_subst))
  initial <- mean(unlist(list_initial))
  print(initial)
  return(data.table(initial = initial, subst = subst, cateonestep, cateplugin, lrnr = names(cateplugin)))




}







HardMeanSimpleCATE <- function(n, CATE_library, nsims) {

  library(sl3)
  list_cateonestep <- list()
  list_cateplugin <- list()
  list_subst <- list()
  list_initial <- list()
  for(i in 1:nsims) {
    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)# rbinom(n, 1 , plogis(W1))
    W <- cbind(W1, W2)
    A <- rbinom(n, 1 , plogis((sin(5*W1) + cos(5*W2) )))
    Y <- rnorm(n, sin(5*W1) + cos(5*W2) + A*(1 + W1 + W2), 1)
    CATE <- (1 + W1 + W2)
    quantile(plogis(2*(W1 + W2 )))
    data <- data.table(W1,W2,A,Y)

    task_A <- sl3_Task$new(data, covariates = c("W1", "W2"), outcome = "A")
    task_Y <- sl3_Task$new(data, covariates = c("W1", "W2", "A"), outcome = "Y")
    data1 <- data
    data1$A <- 1
    task_Y1 <- sl3_Task$new(data, covariates = c("W1", "W2", "A"), outcome = "Y")
    data0 <- data
    data0$A <- 0
    task_Y0 <- sl3_Task$new(data, covariates = c("W1", "W2", "A"), outcome = "Y")
    lrnr_Y <- make_learner(Pipeline, Lrnr_cv$new(
      Stack$new(
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


    lrnr_Y <- lrnr_Y$train(task_Y)
    lrnr_A <- lrnr_A$train(task_A)
    pA1 <- lrnr_A$predict(task_A)
    EY1 <- lrnr_Y$predict(task_Y1)
    EY0 <- lrnr_Y$predict(task_Y0)
    EY <- ifelse(A==1, EY1, EY0)
    pA1 <- pmax(pmin(pA1, 0.98), 0.02)
    pA0 <- 1- pA1
    Ytilde <- EY1 - EY0 + (A/pA1 - (1-A)/pA0)*(Y - EY)
    quantile(Ytilde)
    beta <- coef(glm(Y~X.W1 + X.W2 - 1, data= data.frame(Y = Y, X = (A - (1-A))*W), offset = EY, weights = ifelse(A==1, 1/pA1, 1/pA0)))

    EY1Star <- EY1 +(1 - (1-1))*W %*% beta
    EY0Star <- EY0 +(0 - (1-0))*W %*% beta
    Ytilde_star =  EY1Star - EY0Star
    data <- data.frame(W1, W2, Ytilde = Ytilde, Ytilde_star =  EY1Star - EY0Star)
    task_onestep <- sl3_Task$new(data, covariates = c("W1", "W2"), outcome = "Ytilde")
    task_plugin <- sl3_Task$new(data, covariates = c("W1", "W2"), outcome = "Ytilde_star")

    lrnr <- CATE_library
    lrnr_1 <- lrnr$train(task_onestep)
    CATEonestep <- lrnr_1$predict(task_onestep)
    cateonestep <- apply(CATEonestep, 2, function(est) {
      mean((CATE - est)^2)
    })

    lrnr_1 <- lrnr$train(task_plugin)
    CATEplugin <- lrnr_1$predict(task_plugin)
    cateplugin <- apply(CATEplugin, 2, function(est) {
      mean((CATE - est)^2)
    })

    list_cateonestep[[i]] <- cateonestep
    list_cateplugin[[i]] <- cateplugin
    list_subst[[i]] <- mean((CATE - (EY1Star - EY0Star))^2)
    list_initial[[i]] <- mean((CATE - (EY1 - EY0))^2)
  }
  cateonestep <- rowMeans(do.call(cbind, list_cateonestep))
  cateplugin <- rowMeans(do.call(cbind, list_cateplugin))
  subst <- mean(unlist(list_subst))
  initial <- mean(unlist(list_initial))
  print(initial)
  return(data.table(initial = initial, subst = subst, cateonestep, cateplugin, lrnr = names(cateplugin)))




}




