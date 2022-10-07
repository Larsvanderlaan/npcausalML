# efficient_loss_function <- function(theta, A , Y, EY1W , EY0W , pA1W){
#   LRR <- theta
#   EY <- ifelse(A==1, EY1W, EY0W)
#   plugin_risk <- (EY0W + EY1W) * log(1 + exp(LRR)) - EY1W * LRR
#   score_comp <- (A/pA1W)*(log(1 + exp(LRR)) - LRR)*(Y - EY) + ((1-A)/(1-pA1W))*(log(1 + exp(LRR)) - LRR)*(Y - EY)
#   plugin_risk + score_comp
# }
#
#
# n <- 35
# W <- runif(n, -1 , 1)
# A <- rbinom(n, 1 , plogis(2*W))
# Y <- rpois(n,exp( W + A*(W)))
# LRR <- W
#
#
# data <- data.table(W,A,Y)
# task_A <- sl3_Task$new(data, covariates = c("W"), outcome = "A")
# task_Y <- sl3_Task$new(data, covariates = c("W", "A"), outcome = "Y", outcome_type = "continuous")
# data1 <- data
# data1$A <- 1
# task_Y1 <- sl3_Task$new(data, covariates = c("W", "A"), outcome = "Y")
# data0 <- data
# data0$A <- 0
# task_Y0 <- sl3_Task$new(data, covariates = c("W", "A"), outcome = "Y")
# lrnr_Y <- Lrnr_cv$new( Lrnr_xgboost$new(max_depth = 5, objective = "count:poisson"))
# lrnr_A <- Lrnr_cv$new(Lrnr_glm$new())
# lrnr_Y <- lrnr_Y$train(task_Y)
# lrnr_A <- lrnr_A$train(task_A)
# pA1 <- lrnr_A$predict(task_A)
# EY1 <- lrnr_Y$predict(task_Y1)
# EY0 <- lrnr_Y$predict(task_Y0)
# EY <- ifelse(A==1, EY1, EY0)
# pA1 <- pmax(pmin(pA1, 0.98), 0.02)
# pA0 <- 1 - pA1
#
#
# betas <- seq(-10, 10, length = 25)
# losses <- do.call(cbind, lapply(betas, function(beta){
#    efficient_loss_function(beta * W, A , Y, EY1 , EY0 , pA1)[order(W)]
# }))
#
# plot(betas, colMeans(losses))
#
# #print(persp(betas,W[order(W)],t(losses), xlab = "beta", ylab = "Obs", zlab = "loss", theta = 45,
#          #   phi = 20, col = "lightblue", shade = 0.1))
#
#
