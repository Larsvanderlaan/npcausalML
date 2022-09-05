

#' Estimates the log relative risk function using a user-supplied binomial sl3 Learner object using either the IPW or plugin empirical risk function.
#' @param V A matrix of observations of a subset of the covariates `W` for which to estimate the (possibly semi-marginalized) log relative risk (LRR).
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' @param EY1W A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0W A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1W A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
#' @param  sl3_LRR_Learner_binomial A \code{sl3_Learner} object that minimizes the binomial/logistic risk function.
#' This function will automatically add a `family = binomial()` parameter to the internal params of the inputted \code{sl3_Learner}.
#' If the learner predicts the predictions at the probability scale (i.e. expit transform of link predictor) then set the argument `logit_transform = TRUE` (default). (Almost all Learners do this by default)
#' If the learner predicts the link predictor then set the argument `logit_transform = TRUE`.
#' Note to users familiar with \code{sl3}: if `learning_method = "plugin"`, the outcome_type of the sl3_Task is `quasibinomial`. Therefore, a family object should be passed to the learner to ensure the outcome type is correct.
#' @param learning_method A string being either "plugin" or "IPW". Whether the LRR should be estimated by minimizing the plugin or IPW empirical risk function.
#' @param Vpred A matrix of covariates observations at which to predict the LRR. By default, \code{Vpred} equals \code{W}.
#' @param logit_transform An internal argument specifying whether the predictions of \code{sl3_LRR_Learner_binomial} should be logit-transformed.
#' This argument is needed since the LRR predictions correspond with the logit-scale predictor and not probability-scale predictions of the binomial learner.
#' For most \code{sl3_Learner}s, the default `logit_transform = TRUE` is necessary for this method to work correctly.
#' @export
estimate_using_ERM <- function(V, A, Y,  EY1W, EY0W, pA1W, weights, family_risk_function, sl3_Learner,  outcome_function_plugin, weight_function_plugin, outcome_function_IPW, weight_function_IPW , learning_method = c("plugin", "IPW"), Vpred = V, transform_function = function(x){x}) {

  learning_method <- match.arg(learning_method)
  data <- as.data.table(V)

  covariates <- colnames(data)
  if(learning_method == "plugin") {

    #pseudo_outcome <- EY1W / (EY1W + EY0W)
    pseudo_outcome <- outcome_function_plugin(A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
    pseudo_weights <- weights * weight_function_plugin(A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
    #pseudo_weights <- weights * (EY1W + EY0W)
    data$pseudo_outcome <- pseudo_outcome
    data$pseudo_weights <- pseudo_weights
    params <- sl3_Learner$params
    params$family <- family_risk_function
    sl3_Learner <- sl3_Learner$clone()$reparameterize(params)
    outcome_type <- NULL

    task_ERM <- sl3_Task$new(data, covariates = covariates, outcome = "pseudo_outcome", weights = "pseudo_weights", outcome_type = outcome_type)
    if(task_ERM$outcome_type$type %in% c("constant", "categorical")) {
      task_ERM <- sl3_Task$new(data, covariates = covariates, outcome = "pseudo_outcome", weights = "pseudo_weights", outcome_type = "continuous")
    }



  }
  else if(learning_method == "IPW") {
    sl3_Learner <- Lrnr_mean$new()
    pseudo_outcome <- outcome_function_IPW(A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
    pseudo_weights <- weights * weight_function_IPW(A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
    data$pseudo_outcome <- pseudo_outcome
    data$pseudo_weights <- pseudo_weights


    params <- sl3_Learner$params
    params$family <- family_risk_function
    sl3_Learner <- sl3_Learner$clone()$reparameterize(params)
    task_ERM <- sl3_Task$new(data, covariates = covariates, outcome = "pseudo_outcome", weights = "pseudo_weights")
    if(task_ERM$outcome_type$type %in% c("constant", "categorical")) {
      task_ERM <- sl3_Task$new(data, covariates = covariates, outcome = "pseudo_outcome", weights = "pseudo_weights", outcome_type = "continuous")
    }

  }

  task_ERM_pred <- sl3_Task$new(as.data.table(Vpred), covariates = covariates, outcome = c(), outcome_type = "continuous")

  sl3_Learner_trained <- sl3_Learner$train(task_ERM[task_ERM$weights!=0])
  ERM <- sl3_Learner_trained$predict(task_ERM)
  ERM_pred <- sl3_Learner_trained$predict(task_ERM_pred)
  if(!is.null(transform_function)) {
    ERM <- transform_function(ERM)
    ERM_pred <- transform_function(ERM_pred)
  }

  output <- list(ERM_train = as.matrix(ERM), ERM_pred = as.matrix(ERM_pred), ERM_learner = sl3_Learner_trained)
  return(output)
}

#' The double-robust one-step efficient empirical risk function for the log relative risk.
#' Note this risk function is non-convex.
#' @param ERM A vector or matrix of log relative risk (LRR) estimates whose risk is to be evaluated using the one-step efficient double-robust risk function.
#' @param W A matrix of covariate observations
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' @param EY1W A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0W A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1W A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
#' @param debug ...
#' @param return_loss Boolean for whether to return loss function values or the risk value (i.e. average of the losses)
efficient_risk_function <- function(theta, A, Y, EY1W, EY0W, pA1W, weights, efficient_loss_function, debug = FALSE, return_loss = FALSE, V = NULL, oracle = FALSE) {
  LRR <- as.matrix(theta)
  if(!(nrow(LRR) == length(A) && nrow(LRR) == length(EY1W))) {
    stop("Input lengths dont match")
  }
  #EY <- ifelse(A==1, EY1W, EY0W)
  #plugin_risk <- (EY0W + EY1W) * log(1 + exp(LRR)) - EY1W * LRR
  #score_comp <- (A/pA1W)*(log(1 + exp(LRR)) - LRR)*(Y - EY) + ((1-A)/(1-pA1W))*(log(1 + exp(LRR)) - LRR)*(Y - EY)
  DR_loss <- weights*apply(as.matrix(LRR), 2, efficient_loss_function, V = V, A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
  if(debug){
    #print(colMeans(weights * score_comp))
  }
  #DR_loss <- weights * (plugin_risk + score_comp)
  if(return_loss) {
    return(DR_loss)
  } else {
    return(colMeans(DR_loss))
  }
}
