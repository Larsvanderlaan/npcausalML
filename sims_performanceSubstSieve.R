library(SuperLearner)
library(npcausalML)
library(future)
source("simCATE.R")
SL.gam3 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 3
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam4 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 4
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
SL.gam5 <- function(Y, X, newX, family, obsWeights, cts.num = 4,...) {
  deg.gam <- 5
  SL.gam(Y, X, newX, family, obsWeights, deg.gam, cts.num,... )
}
list_of_sieves_uni   <- list(
  "no_sieve" = NULL,
  fourier_basis$new(orders = c(1,0)),
  fourier_basis$new(orders = c(2,0)),
  fourier_basis$new(orders = c(3,0)),
  fourier_basis$new(orders = c(4,0))
)

lrnr_gam3 <- Lrnr_pkg_SuperLearner$new("SL.gam3" )
lrnr_gam4 <- Lrnr_pkg_SuperLearner$new("SL.gam4" )
lrnr_gam5 <- Lrnr_pkg_SuperLearner$new("SL.gam5" )


hard <- F
pos <- F

onesim <- function(n) {

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
    Lrnr_xgboost$new(max_depth =4),
    Lrnr_xgboost$new(max_depth =5),
    Lrnr_xgboost$new(max_depth =6)))#Stack$new(
    #Lrnr_stratified$new(Lrnr_gam$new(), "A"))
    , Lrnr_cv_selector$new(loss_squared_error))


  lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(
    Stack$new(
      Lrnr_xgboost$new(max_depth =4),
      Lrnr_xgboost$new(max_depth =5),
      Lrnr_xgboost$new(max_depth =6)
    )
  ), Lrnr_cv_selector$new(loss_squared_error))


  data_train <-  data #as.data.frame(sim.CATE(n, hard, pos))

  initial_likelihood <- npcausalML:::estimate_initial_likelihood(W=data_train[,c("W1", "W2","W3")], data_train$A, data_train$Y,  weights = rep(1,n), lrnr_A, lrnr_Y, folds = 10)
  data1 <- data
  data0 <- data
  data1$A <- 1
  data0$A <- 0
  taskY <- sl3_Task$new(data, covariates = c("W1", "W2", "W3", "A"), outcome = "Y")
  taskY0 <- sl3_Task$new(data0, covariates = c("W1", "W2", "W3", "A"), outcome = "Y")
  taskY1 <- sl3_Task$new(data1, covariates = c("W1", "W2", "W3", "A"), outcome = "Y")
  taskA <- sl3_Task$new(data, covariates = c("W1", "W2", "W3"), outcome = "A")

  pA1W_est <- initial_likelihood$internal$sl3_Learner_pA1W_trained$predict(taskA)
  EY1W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY1)
  EY0W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY0)

  pA1W_est <- pmax(pA1W_est, 0.01)
  pA1W_est <- pmin(pA1W_est, 0.99)

  CATE_library <- list(  Lrnr_xgboost$new(max_depth =4),
                         Lrnr_xgboost$new(max_depth =5),
                         Lrnr_xgboost$new(max_depth =6) )

  subst_compare <- Stack$new(CATE_library)
  CATE_library_strat <- lapply(CATE_library , function (lrnr) {
    Lrnr_stratified$new(lrnr, "A")
  })

  subst_compare <- Stack$new(CATE_library)
  subst_compare <- subst_compare$train(taskY)
  subst_EY1W <-subst_compare$predict(taskY1)
  subst_EY0W <- subst_compare$predict(taskY0)
  subst_CATE <- subst_EY1W - subst_EY0W






  est_sieve_nuisances <- npcausalML:::compute_plugin_and_IPW_sieve_nuisances(data.frame(W1=data$W1), A, Y, EY1W_est, EY0W_est, pA1W_est, rep(1,n), list_of_sieves_uni[[4]], design_function_sieve_plugin_CATE, weight_function_sieve_plugin_CATE, design_function_sieve_IPW_CATE, weight_function_sieve_IPW_CATE, family_for_targeting = binomial(), debug = FALSE)
  mean((est_sieve_nuisances$EY1W_star - est_sieve_nuisances$EY0W_star - CATE)^2)
  mean((EY1W_est - EY0W_est - CATE)^2)


  risk_subst<-  apply(subst_CATE, 2, function(pred) {
    mean((pred - CATE)^2)
  })

  risk_subst_cv <- mean((EY1W_est - EY0W_est - CATE)^2)



  list(risk_subst_cv = risk_subst_cv, risk_subst = risk_subst)
}

nsims <- 50
print(500)
#simresults500 <- lapply(1:nsims, function(i){
# print(i)
#onesim(500)
#})

#save(simresults500, file = paste0("simsCATE", hard,pos, "n500_3"))


print(1000)
simresults1000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(1000)
})

save(simresults1000, file = paste0("simsCATE", hard,pos, "n1000_3"))

print(2500)
simresults2500 <- lapply(1:nsims, function(i){
  print(i)
  onesim(2500)
})

save(simresults2500, file = paste0("simsCATE", hard,pos, "n2500_3"))

print(5000)
simresults5000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(5000)
})

save(simresults5000, file = paste0("simsCATE", hard,pos, "n5000_3"))


print(10000)
simresults10000 <- lapply(1:nsims, function(i){
  print(i)
  onesim(10000)
})

save(simresults10000, file = paste0("simsCATE", hard,pos, "n10000_3"))
