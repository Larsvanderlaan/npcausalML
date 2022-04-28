
#' Train the LRR learners for all sieves and base learners on the training set specified by a fold object.
#' @param fold \code{fold} object from the \link{origami} package specifying the fold-specific training set at which to estimate the LRR.
#' @param V A matrix of observations of a subset of the covariates `W` for which to estimate the (possibly semi-marginalized) log relative risk (LRR).
#'  This should include all observations including both the training and validation observations.
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' This should include all observations including both the training and validation observations.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' This should include all observations including both the training and validation observations.
#' @param EY1W A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0W A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1W A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
#' @param list_of_learners A list of untrained \code{sl3_Learner} objects to be used to estimate the log relative risk LRR using the function \link{estimate_LRR_using_ERM}.
#' @param list_of_sieves A list of basis_generator objects specifying the sieve. See, for example, \code{fourier_basis} for an example and template.
train_learners_using_fold <- function(fold, V, A, Y, EY1W, EY0W, pA1W, weights, family_risk_function, outcome_function_plugin, weight_function_plugin,  outcome_function_IPW, weight_function_IPW, transform_function, design_function_sieve_plugin, weight_function_sieve_plugin, design_function_sieve_IPW, weight_function_sieve_IPW, family_for_targeting,  list_of_learners, list_of_sieves, Vpred = V) {
  list_of_learners <- list_of_learners
  index <- origami::training(fold = fold)
  index_val <- origami::validation(fold = fold)
  Vfull <- V
  V <- Vfull[index, ,drop = F]
  Vval <- Vfull[index_val,  ,drop = F]
  A <- A[index]
  Y <- Y[index]
  EY1W <- EY1W[index]
  EY0W <- EY0W[index]
  pA1W <- pA1W[index]
  weights <- weights[index]

  all_learners_delayed <- train_learners(V, A, Y, EY1W, EY0W, pA1W, weights, family_risk_function, outcome_function_plugin, weight_function_plugin,  outcome_function_IPW, weight_function_IPW, transform_function, design_function_sieve_plugin, weight_function_sieve_plugin, design_function_sieve_IPW, weight_function_sieve_IPW, family_for_targeting,  list_of_learners, list_of_sieves, Vpred = Vval, compute = FALSE)
  return(all_learners_delayed)
}

#' For each LRR learner and fold, the best sieve is chosen by minimizing the double-robust one-step efficient empirical risk function over the choices of sieve-spaces.
#' @param folds A \code{folds} object (list of \code{fold} objects) from the \link{origami} package specifying the fold-specific training set at which to estimate the LRR.
#' @param trained_learner_list A unnested list of trained LRR learners as
#' @param learner_names Names of LRR learners. Should be `learner_names <- sapply(list_of_learners, `[[`, "name)`.
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' This should include all observations including both the training and validation observations.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' This should include all observations including both the training and validation observations.
#' @param EY1W A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0W A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1W A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
subset_best_sieve_all_folds <- function(folds, trained_learner_list, learner_names, A, Y, EY1W, EY0W, pA1W, weights, efficient_loss_function) {
  trained_learner_list <- trained_learner_list
  output <- lapply(seq_along(folds), function(fold_number) {
    fold <- folds[[fold_number]]
    training_index <- origami::training(fold = fold)
    keep <- which(stringr::str_detect(names(trained_learner_list), paste0("^", fold_number, "\\.", "+")))
    learners <- trained_learner_list[keep]
    learners <- subset_best_sieve(learners, learner_names, A[training_index], Y[training_index], EY1W[training_index], EY0W[training_index], pA1W[training_index], weights[training_index], efficient_loss_function)

    return(learners)

  })
  names(output) <- seq_along(folds)
  return(output)
}



cv_predict_learner <- function(folds, learners_all_folds) {
  cv_fun <- function(fold) {
    fold_number <- fold_index()
    index <- validation()
    v <- origami::fold_index(fold = fold)
    list(index = index,
         fold_index = rep(fold_index(), length(index)), predictions=as.data.table(do.call(cbind, lapply(learners_all_folds[[v]] , `[[`, "ERM_pred"))))
  }
  comb_ctrl <- list(combiners = list(
    index = combiner_c, fold_index = combiner_c,
    predictions = function(x) rbindlist(x, fill = TRUE)
  ))
  cv_predictions <- origami::cross_validate(cv_fun, folds, .combine_control = comb_ctrl)
  cv_predictions <- as.matrix(cv_predictions$predictions[order(cv_predictions$index),] )
}
