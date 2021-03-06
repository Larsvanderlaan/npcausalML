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
hard <- T
pos <- T

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
      Lrnr_gam$new()
    )
  ), Lrnr_cv_selector$new(loss_squared_error))


  data_train <-  as.data.frame(sim.CATE(n, hard, pos))

  initial_likelihood <- npcausalML:::estimate_initial_likelihood(W=data_train[,c("W1", "W2","W3")], data_train$A, data_train$Y,  weights = rep(1,n), lrnr_A, lrnr_Y, folds = 10)
  data1 <- data
  data0 <- data
  data1$A <- 1
  data0$A <- 0
  taskY0 <- sl3_Task$new(data0, covariates = c("W1", "W2", "W3", "A"), outcome = "Y")
  taskY1 <- sl3_Task$new(data1, covariates = c("W1", "W2", "W3", "A"), outcome = "Y")
  taskA <- sl3_Task$new(data, covariates = c("W1", "W2", "W3"), outcome = "A")

  pA1W_est <- initial_likelihood$internal$sl3_Learner_pA1W_trained$predict(taskA)
  EY1W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY1)
  EY0W_est <- initial_likelihood$internal$sl3_Learner_EYAW_trained$predict(taskY0)
  pA1W_est <- pmax(pA1W_est, 0.05)
  pA1W_est <- pmin(pA1W_est, 0.95)


  lrnr_gam3 <- Lrnr_pkg_SuperLearner$new("SL.gam3" )
  lrnr_gam4 <- Lrnr_pkg_SuperLearner$new("SL.gam4" )
  lrnr_gam5 <- Lrnr_pkg_SuperLearner$new("SL.gam5" )



  CATE_library <- list(Lrnr_svm$new(), Lrnr_ranger$new(max.depth = 8), Lrnr_ranger$new(max.depth = 10), Lrnr_ranger$new(max.depth = 15), Lrnr_ranger$new(max.depth = 25), Lrnr_earth$new(), Lrnr_xgboost$new(max_depth = 3, verbosity = 0), Lrnr_xgboost$new(max_depth = 5, verbosity = 0),   Lrnr_xgboost$new(max_depth = 7, verbosity = 0), Lrnr_rpart$new(), lrnr_gam3, lrnr_gam4 , lrnr_gam5, Lrnr_glm$new()    )



  fit_npcausalML <- npcausalML(CATE_library,
                               W= W, A = A, Y = Y, V = data.frame(W1 = W$W1),
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
                               cross_validate_ERM = T, folds = origami::folds_vfold(length(A), 5))


  preds <- predict(fit_npcausalML,  data.frame(W1 = W$W1), F)


  # Compute least-squares risk of predictions using oracle loss function.
  risks_oracle <- as.vector(apply(preds, 2, function(theta) {
    mean((theta -  CATE)^2)
  })[grep("plugin", colnames(preds))])

  # Compute estimated cross-validated one-step risk of predictions
  cvrisksDR <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
    loss <- efficient_loss_function_CATE(W, theta, A, Y, EY1W_est,EY0W_est, pA1W_est )
    mean(loss)
  })[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])

  # Compute estimated cross-validated oracle one-step risk of predictions

  cvrisksDRoracle <- as.vector(apply(fit_npcausalML$cv_predictions, 2, function(theta) {
    loss <- efficient_loss_function_CATE(W, theta, A, Y, EY1Wtrue,EY0Wtrue, pA1Wtrue )
    mean(loss)
  })[-grep("IPW", colnames(fit_npcausalML$cv_predictions))])


  #CATE_library <- list(Lrnr_ranger$new(max.depth = 8), Lrnr_ranger$new(max.depth = 10), Lrnr_ranger$new(max.depth = 15), Lrnr_ranger$new(max.depth = 25))
  CATEonestepbench <- DR_learner(CATE_library, W1, A, Y, EY1W_est, EY0W_est, pA1W_est, NULL, NULL)
  CATEonestepbench <- apply(CATEonestepbench, 2, function(pred) {
    mean((pred - CATE)^2)
  })
  #CATEonestepbench

  CATEonestepbenchoracle <- DR_learner(CATE_library, W1, A, Y, EY1Wtrue, EY0Wtrue, pA1Wtrue, NULL, NULL)
  CATEonestepbenchoracle <- apply(CATEonestepbenchoracle, 2, function(pred) {
    mean((pred - CATE)^2)
  })


  risk_subst<-  mean((CATE - (EY1W_est  - EY0W_est))^2)

  Y.hat <- EY1W_est * pA1W_est + EY0W_est * (1-pA1W_est)
  W.hat <- pA1W_est
  fit <- grf::causal_forest(X = W, Y  = Y, W = A, Y.hat = Y.hat, W.hat = W.hat)
  preds_cf <-  fit$predictions
  risk_cf <- mean((CATE - preds_cf)^2)


  list(risk_cf = risk_cf, risk_subst = risk_subst, CATEonestepbenchoracle =CATEonestepbenchoracle, CATEonestepbench = CATEonestepbench, sieve =data.frame(grep("plugin", colnames(fit_npcausalML$cv_predictions), value = T), cvrisksDRoracle, cvrisksDR, risks_oracle))
}
print(500)
simresults <- lapply(1:20, function(i){
  print(i)
  onesim(500)
})

save(simresults, file = "simsCATETTn500")

print(1000)
simresults <- lapply(1:20, function(i){
  print(i)
  onesim(1000)
})

save(simresults, file = "simsCATETTn1000")

print(2500)
simresults <- lapply(1:20, function(i){
  print(i)
  onesim(2500)
})

save(simresults, file = "simsCATETTn2500")

print(5000)
simresults <- lapply(1:20, function(i){
  print(i)
  onesim(5000)
})

save(simresults, file = "simsCATETTn5000")


print(10000)
simresults <- lapply(1:20, function(i){
  print(i)
  onesim(10000)
})

save(simresults, file = "simsCATETTn10000")


load("simsCATEFFn10000_2")
onestepbenchoracle <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbenchoracle")))
onestepbench  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbench")))

causalforestrisks <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_cf")))
substrisks  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst")))



lrnr_names <- simresults[[1]]$sieve[[1]]

cvrisksDRoracle <- rowMeans(do.call(cbind, lapply(simresults, function(item) {
  item$sieve$cvrisksDRoracle
})) )

cvrisksDR <- rowMeans(do.call(cbind, lapply(simresults, function(item) {
  item$sieve$cvrisksDR
})) )

risks_oracle <- rowMeans(do.call(cbind, lapply(simresults, function(item) {
  item$sieve$risks_oracle
})) )

dt <- data.table(lrnr_full = lrnr_names, cvrisksDRoracle,cvrisksDR,  risks_oracle)
dt$degree <- as.numeric(stringr::str_match(dt$lrnr_full, "fourier_basis_([0-9]+)")[,2])
dt$degree[grep("no_sieve", dt$lrnr_full)] <- 0
#tmp <- data.table(risks_oracle = onestepbench, lrnr_full = names(onestepbench), degree = "onestep")
#dt <- rbind(dt, tmp, fill = T)
#tmp <- data.table(risks_oracle = onestepbenchoracle, lrnr_full = names(onestepbench), degree = "onesteporacle")
#dt <- rbind(dt, tmp, fill = T)

dt$lrnr[ grep("gam", dt$lrnr_full)] <- "gam"
dt$lrnr[ grep("glm", dt$lrnr_full)] <- "glm"
dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
dt$lrnr[ grep("rpart", dt$lrnr_full)] <- "rpart"
dt$lrnr[ grep("xgboost_20_1_5", dt$lrnr_full)] <- "xgboost_5"
dt$lrnr[ grep("xgboost_20_1_3", dt$lrnr_full)] <- "xgboost_3"
dt$lrnr[ grep("xgboost_1", dt$lrnr_full)] <- "randomforest"
dt$lrnr[ grep("xgboost_1", dt$lrnr_full)] <- "randomforest"
dt$lrnr[ grep("SL.gam", dt$lrnr_full)] <- "gam5"
dt$lrnr[ grep("Lrnr_gam", dt$lrnr_full)] <- "gamcv"
dt$lrnr[ grep("svm", dt$lrnr_full)] <- "svm"

dt <- dt[dt$lrnr != "randomforest",]

library(ggplot2)

ggplot(dt, aes(x = degree, y = cvrisksDRoracle, color = lrnr, group = lrnr)) + geom_line()
ggplot(dt, aes(x = degree, y = cvrisksDR, color = lrnr, group = lrnr)) + geom_line()
ggplot(dt, aes(x = degree, y = risks_oracle, color = lrnr, group = lrnr)) + geom_line() # + geom_hline(yintercept = onestepbench )


ns <- c( 1000, 2500, 5000, 10000)
sims_list <- lapply(ns, function(n) {
  load(paste0("simsCATEFTn", n,"_2"))
  onestepbenchoracle <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbenchoracle")))
  onestepbench  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbench")))

  causalforestrisks <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_cf")))
  substrisks  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst")))



  lrnr_names <- simresults[[1]]$sieve[[1]]

  cvrisksDRoracle <- rowMeans(do.call(cbind, lapply(simresults, function(item) {
    item$sieve$cvrisksDRoracle
  })) )

  cvrisksDR <- rowMeans(do.call(cbind, lapply(simresults, function(item) {
    item$sieve$cvrisksDR
  })) )

  risks_oracle <- rowMeans(do.call(cbind, lapply(simresults, function(item) {
    item$sieve$risks_oracle
  })) )

  dt <- data.table(lrnr_full = lrnr_names, cvrisksDRoracle,cvrisksDR,  risks_oracle)
  dt$degree <- as.numeric(stringr::str_match(dt$lrnr_full, "fourier_basis_([0-9]+)")[,2])
  dt$degree[grep("no_sieve", dt$lrnr_full)] <- 0
  tmp <- data.table(risks_oracle = onestepbench, lrnr_full = names(onestepbench), degree = "DR")
  dt <- rbind(dt, tmp, fill = T)
  tmp <- data.table(risks_oracle = onestepbenchoracle, lrnr_full = names(onestepbench), degree = "DRoracle")
  dt <- rbind(dt, tmp, fill = T)
  tmp <- data.table(risks_oracle = causalforestrisks, lrnr = "causalforest", degree = "causalforest")
  dt <- rbind(dt, tmp, fill = T)
  tmp <- data.table(risks_oracle = substrisks, lrnr = "subst", degree = "subst")
  dt <- rbind(dt, tmp, fill = T)

  dt$lrnr[ grep("gam3", dt$lrnr_full)] <- "gam3"
  dt$lrnr[ grep("gam4", dt$lrnr_full)] <- "gam4"
  dt$lrnr[ grep("gam5", dt$lrnr_full)] <- "gam5"
  dt$lrnr[ grep("glm", dt$lrnr_full)] <- "glm"
  dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
  dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
  dt$lrnr[ grep("rpart", dt$lrnr_full)] <- "rpart"
  dt$lrnr[ grep("xgboost_20_1_5", dt$lrnr_full)] <- "xgboost_5"
  dt$lrnr[ grep("xgboost_20_1_3", dt$lrnr_full)] <- "xgboost_3"
  dt$lrnr[ grep("xgboost_1", dt$lrnr_full)] <- "randomforest"
  dt$lrnr[ grep("svm", dt$lrnr_full)] <- "svm"
  dt$lrnr[ grep("ranger", dt$lrnr_full)] <- "ranger"
  print(dt$lrnr[37:42])
  dt$lrnr[37:42] <- "gam3"
  dt$lrnr[43:48] <- "gam4"
  dt$lrnr[49:54] <- "gam5"
  dt[!is.na(as.numeric(dt$degree)), risks_best := min(risks_oracle), by = c("lrnr")]
  dt[is.na(as.numeric(dt$degree)), risks_best := risks_oracle, by = c("lrnr")]
  dt$type[!is.na(as.numeric(dt$degree))] <- "sieve"
  dt$type[is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

  dt2 <- dt[,c("lrnr", "risks_best", "type"), with = F]
  dt2 <- unique(dt2)
  dt2$n <- n
  return(dt2)
})


dt <- rbindlist(sims_list)
dt <- dt[dt$lrnr != "subst"]
dt <- dt[dt$lrnr != "causalforest"]

plt <- ggplot(dt[!(dt$lrnr %in% c("glm", "earth", "gam3", "gam4", "gam5")),], aes(x = n, y = risks_best, group = type, color = type, linetype = type)) + geom_line() +
  facet_wrap(~lrnr, scales = "free") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ylab("MSE")
ggsave("performancePlotSieveForest2.pdf")
plt
