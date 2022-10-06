

sim.CATE <- function(n, hard = TRUE, positivity = TRUE, randomized = F,  ...) {
  if(randomized) {
    W1 <- runif(n, -1 , 1)
    W2 <- runif(n, -1 , 1)
    W3 <- runif(n, -1 , 1)
    W <- data.table(W1, W2, W3)
    pA1W <- 0.5
    A <- rbinom(n, 1 ,  pA1W)
    CATE <- 1 + W1
    EY0W <- 0.5*(W1 + W2 + W3) + sin(5*W1) + sin(5*W2) + sin(5*W3) + 1/(W1 + 1.2) + 1/(W2 + 1.2) + 1/(W3 + 1.2)
    EY1W <- EY0W + CATE
    Y <- rnorm(n, EY0W + A*CATE, 2)
    return(data.table(W, A, Y, EY0W, EY1W, pA1W))

  }
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

  data <- data.table(W, Y=Ytilde)
  task_onestep <- sl3_Task$new(data, covariates = colnames(W), outcome = "Y")
  lrnr <- Stack$new(CATE_library)
  lrnr_1 <- lrnr$train(task_onestep)
  CATEonestep <- lrnr_1$predict(task_onestep)
  return(CATEonestep)
}

