
#' @export
npcausalML <- function(learners, W, A, Y, V = W, weights, Vpred = V, EY1W, EY0W, pA1W, sl3_Learner_EYAW, sl3_Learner_pA1W, outcome_function_plugin, weight_function_plugin, outcome_function_IPW, weight_function_IPW, design_function_sieve_plugin, weight_function_sieve_plugin, design_function_sieve_IPW, weight_function_sieve_IPW, outcome_type = c("binomial", "continuous", "nonnegative"), family_risk_function, family_for_targeting, transform_function, efficient_loss_function, list_of_sieves, cross_validate_ERM = FALSE, use_sieve_selector = TRUE, folds = origami::make_folds(n=length(A))) {
  if(!is.list(learners)) {
    learners <- list(learners)
  }
  list_of_learners <- learners
  learner_names <- names(list_of_learners)
  if(is.null(learner_names) || length(learner_names)!= length(list_of_learners)) {
    learner_names <- sapply(list_of_learners, `[[`, "name")
    names(list_of_learners) <- learner_names
  }
  W <- as.matrix(W)
  A <- as.vector(A)
  Y <- as.vector(Y)
  if(missing(weights)) {
    weights <- rep(1, length(A))
  }
  weights <- as.vector(weights)

  try({
    likelihood <- estimate_initial_likelihood(W, A, Y, weights, sl3_Learner_EYAW = sl3_Learner_EYAW, sl3_Learner_pA1W = sl3_Learner_pA1W, folds = folds, outcome_type = outcome_type)
    EY1W <- likelihood$EY1W
    EY0W <- likelihood$EY0W
    pA1W <- likelihood$pA1W
  })
  EY1W <- as.vector(EY1W)
  EY0W <- as.vector(EY0W)
  pA1W <- as.vector(pA1W)

  if(missing(list_of_sieves)) {
    if(ncol(W))
      list_of_sieves <- list(
        NULL,
        fourier_basis$new(orders = c(1,0)),
        fourier_basis$new(orders = c(2,0)),
        fourier_basis$new(orders = c(3,0)),
        fourier_basis$new(orders = c(1,1)),
        fourier_basis$new(orders = c(2,1)),
        ourier_basis$new(orders = c(3,1))
      )
  }


  full_fit_ERM <- train_learners(V, A, Y, EY1W, EY0W, pA1W, weights, family_risk_function, outcome_function_plugin, weight_function_plugin,  outcome_function_IPW, weight_function_IPW, transform_function, design_function_sieve_plugin, weight_function_sieve_plugin, design_function_sieve_IPW, weight_function_sieve_IPW, family_for_targeting,  list_of_learners, list_of_sieves, Vpred = V)
  if(use_sieve_selector) {
    all_ERM_full_best <- subset_best_sieve(full_fit_ERM, learner_names, A, Y, EY1W, EY0W, pA1W, weights, efficient_loss_function)
  } else {
    all_ERM_full_best <- full_fit_ERM
  }
   all_ERM_full_predictions <- do.call(cbind, lapply(all_ERM_full_best, `[[`, "ERM_pred"))
  output_list <- list(learners = all_ERM_full_best)
  if(cross_validate_ERM) {
    if(use_sieve_selector) {
      print("use_sieve_selector is overwritten to TRUE when cross_validate_ERM = TRUE")
      use_sieve_selector <- TRUE
    }
    learners_all_folds <- lapply(folds, train_learners_using_fold, V, A, Y, EY1W, EY0W, pA1W, weights, family_risk_function, outcome_function_plugin, weight_function_plugin,  outcome_function_IPW, weight_function_IPW, transform_function, design_function_sieve_plugin, weight_function_sieve_plugin, design_function_sieve_IPW, weight_function_sieve_IPW, family_for_targeting,  list_of_learners, list_of_sieves, Vpred = V )
    names(learners_all_folds) <- seq_along(folds)
    learners_all_folds <- unlist(learners_all_folds)
    learner_sieve_names <- names(learners_all_folds)
    learners_all_folds_delayed <- bundle_delayed(learners_all_folds)
    learners_all_folds <- learners_all_folds_delayed$compute()
    learners_best_sieve_all_folds <- subset_best_sieve_all_folds(folds, learners_all_folds, learner_names, A, Y, EY1W, EY0W, pA1W, weights, efficient_loss_function = efficient_loss_function)
    cv_predictions <- cv_predict_learner(folds, learners_best_sieve_all_folds)
    print(data.table(cv_predictions))
    best_learner_index_cv <- which.min(efficient_risk_function(cv_predictions, A , Y, EY1W, EY0W, pA1W, weights, efficient_loss_function))
    print(dim(cv_predictions))
    all_ERM_full_predictions <- all_ERM_full_predictions[,best_learner_index_cv]
    output_list$cv_index <- best_learner_index_cv
  }

  return(output_list)




}

#' @export
predict <- function(output, Wpred) {

  learners <- output$learners
  task <- sl3_Task$new(as.data.table(Wpred), covariates = colnames(Wpred), outcome = c())
  pred_list <- lapply(seq_along(learners), function(index) {
    learner <- learners[[index]]
    name <-  names(learners)[index]
    lrnrs <- learner[["ERM_learner"]]
    if(!is.list(lrnrs)) {
      lrnrs <- list(lrnrs)
    }
    out <- do.call(cbind, lapply(seq_along(lrnrs), function(index) {
      lrnr <- lrnrs[[index]]
      preds <- as.matrix(lrnr$predict(task))[,index]
    }))
    colnames(out) <- paste0(name, "_", seq_along(colnames(out)))
    out

  })
  preds <- do.call(cbind, pred_list)
  cv_index <- output$cv_index
  if(!is.null(cv_index)) {
    preds <- preds[, cv_index]
  }

  return(preds)
}


fit_npcausalML
