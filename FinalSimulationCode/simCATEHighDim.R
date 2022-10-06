
library(tmvtnorm)



sim.CATEHighDim <- function(n, hard = TRUE, positivity = TRUE, ...) {
  d <- 20
  Sigma <- outer(1:d, 1:d, Vectorize(function(i,j) {
    if(i==j) {return(1)}
    else {
      return(0.4)
    }
  }))

  W <- rtmvnorm(n, mean = rep(0, d),
                sigma = Sigma,
                lower=rep(-2, d),
                upper=rep( 2, d), algorithm = "gibbs")/2


  colnames(W) <- paste0("W", 1:d)
  W1 <- W[,1]
  W5 <- W[,5]
  W9 <- W[,9]
  W11 <- W[,11]
  W19 <- W[,19]
  W10 <- W[,10]
  W15 <- W[,15]
  if(positivity){
    pA1W <- plogis((W1 + W5 + W9 + W11 + W19)/5  )
  } else {
    pA1W <- plogis((W1 + W5 + W9 + W11 + W19) /1.3  )

  }
  quantile(pA1W)
  A <- rbinom(n, 1 ,  pA1W)
  if(!hard) {
    CATE <- 1 + (W1 + W5 + W9 + W15 + W10)/5

  } else {
    CATE <- 1 + (sin(4*W1) + sin(4*W5) + cos(4*W9) + 1.5*(W15^2 - W10^2))/5

  }

  EY0W <- (cos(4*W1) + cos(4*W5) + sin(4*W9) + 1/(1.5+W15) + 1/(1.5+ W10))/5 +
    (sin(5*W1) + sin(5*W5) + sin(5*W9))/2

  EY1W <- EY0W + CATE
  Y <- rnorm(n, EY0W + A*CATE, 1)

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








list_of_sieves_high_dim <-   list(
  NULL,
  fourier_basis$new(orders = c(1,0)),
  fourier_basis$new(orders = c(2,0)),
  fourier_basis$new(orders = c(3,0))
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

  data <- data.table(W, Ytilde)
  task_onestep <- sl3_Task$new(data, covariates = colnames(W), outcome = "Ytilde")
  lrnr <- Stack$new(CATE_library)
  lrnr_1 <- lrnr$train(task_onestep)
  CATEonestep <- lrnr_1$predict(task_onestep)
  return(CATEonestep)
}

