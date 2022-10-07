# Outline
# Cross-validation external
# Base function
# Construct pseudo-outcomes



#' @export
compute_EP_nuisances <- function(V, A, Y, EY1W, EY0W, pA1W, weights, basis_generator, design_transform, weight_transform, debug = FALSE) {
  EY <- ifelse(A==1, EY1W, EY0W)
  # If no sieve generator given then return original nuisance estimators
  if(is.null(basis_generator)) {
    return(list(pA1W_star = pA1W, EY1W_star = EY1W, EY0W_star = EY0W, sieve = "no_sieve"))
  }


  # Get design matrices and weights for targeting of nuisance estimators
  X <- basis_generator$set(V)$eval(V)
  X_pseudo_plugin <- design_transform(X,A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
  X1_pseudo_plugin <- design_transform(X,A = rep(1, length(A)), Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
  X0_pseudo_plugin <- design_transform(X,A = rep(0, length(A)), Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
  weights_pseudo_plugin <- weights * weight_transform(A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)

  # Prepare targeting using the logistic model
  ## Rescaling and centering of outcome and outcome regression estimators [0,1]
  upper_bound <- max(c(Y, EY1W, EY0W))
  lower_bound <- min(c(Y, EY1W, EY0W))
  if(lower_bound < 0 || upper_bound > 1) {
    lower_bound <- min(0.9 * lower_bound, 1.1 * lower_bound)
    upper_bound <- max(0.9 * upper_bound, 1.1 * upper_bound)
  } else {
    lower_bound <- 0
    upper_bound <- 1
  }
  Y_scaled <- (Y - lower_bound) / (upper_bound - lower_bound)
  EY1W_scaled <- (EY1W - lower_bound) / (upper_bound - lower_bound)
  EY0W_scaled <- (EY0W - lower_bound) / (upper_bound - lower_bound)
  EY_scaled <- (EY - lower_bound) / (upper_bound - lower_bound)
  ## Fit sieve-adjusted mu-estimator
  suppressWarnings(sieve_fit_plugin <- glm.fit(X_pseudo_plugin, Y_scaled, weights = weights_pseudo_plugin, offset = qlogis(EY_scaled), family = binomial(), intercept = F))
  beta_plugin <- coef(sieve_fit_plugin)
  beta_plugin[is.na(beta_plugin)] <- 0
  # Get unscaled and uncentered sieve-adjusted predictions
  EY1W_star <- lower_bound + (upper_bound - lower_bound) * as.vector(plogis(qlogis(EY1W_scaled) + X1_pseudo_plugin %*% beta_plugin))
  EY0W_star <- lower_bound + (upper_bound - lower_bound) * as.vector(plogis(qlogis(EY0W_scaled) + X0_pseudo_plugin %*% beta_plugin))

  output <- list(pA1W_star = pA1W, EY1W_star = EY1W_star, EY0W_star = EY0W_star, sieve = basis_generator$name)


}

#EP_learner_spec <- list(outcome_type, sieve_design_transform, sieve_weight_transform, EP_family, EP_outcome_transform, EP_weight_transform)
#' @export
get_EP_learner_task_helper <- function(V, A, Y, EY1W, EY0W, pA1W, weights, sieve_basis_generator, EP_learner_spec, folds) {

  EP_nuisances <- compute_EP_nuisances(V, A, Y, EY1W, EY0W, pA1W, weights, sieve_basis_generator, EP_learner_spec$sieve_design_transform, EP_learner_spec$sieve_weight_transform )
  EY1W_star <- EP_nuisances$EY1W_star
  EY0W_star <- EP_nuisances$EY0W_star

  pseudo_outcome <- EP_learner_spec$EP_outcome_transform(A = A, Y = Y, EY1W = EY1W_star, EY0W = EY0W_star, pA1W = pA1W)
  pseudo_weights <- weights * EP_learner_spec$EP_weight_transform(A = A, Y = Y, EY1W = EY1W_star, EY0W = EY0W_star, pA1W = pA1W)

  data <- as.data.table(V)
  data$pseudo_outcome <- pseudo_outcome
  data$pseudo_weights <- pseudo_weights
  task <- sl3_Task$new(data, covariates = colnames(V), outcome = "pseudo_outcome", weights = "pseudo_weights", outcome_type = EP_learner_spec$outcome_type, folds = folds)
  return(task)

}
#' @export
get_EP_learner_task <- function(V, A, Y, EY1W, EY0W, pA1W, weights, sieve_basis_generator, EP_learner_spec, folds = NULL) {
  task_full <-  get_EP_learner_task_helper(V, A, Y, EY1W, EY0W, pA1W, weights, sieve_basis_generator, EP_learner_spec, folds = folds)
  if(is.null(folds)) {
    return(task_full)
  }
  EP_learner_tasks_all_folds <- lapply(folds, function(fold) {
    train <- training()
    val <- validation()
    n <- length(train) + length(val)
    sieve_weight_transform_fold <- function(...) {
      (1:n %in% train) * EP_learner_spec$sieve_weight_transform(...)
    }
    EP_learner_spec_fold <- EP_learner_spec
    EP_learner_spec_fold$sieve_weight_transform <- sieve_weight_transform_fold
    get_EP_learner_task_helper(V, A, Y, EY1W, EY0W, pA1W, weights , sieve_basis_generator, EP_learner_spec_fold, folds = folds)
  })
  task_generator <- function(task, fold_number) {
    if(fold_number == "full") {
      return(task)
    } else {
      return(EP_learner_tasks_all_folds[[fold_number]])
    }
  }
  return(sl3_revere_Task$new(task_generator, task_full ))
}


#
#' @export
EP_learn <- function(EP_learner, V, A, Y, EY1W, EY0W, pA1W, weights = rep(1, length(A)), sieve_basis_generator_list, EP_learner_spec, train_learners = TRUE, cross_validate = TRUE, nfolds = 10){
  if(cross_validate) {
    folds = origami::folds_vfold(length(A), nfolds)
  } else {
    folds <- NULL
  }
  all_tasks <- lapply(sieve_basis_generator_list, function(sieve_basis_generator) {
    task <- get_EP_learner_task(V, A, Y, EY1W, EY0W, pA1W, weights, sieve_basis_generator, EP_learner_spec, folds = folds)
    return(task)
  })
  if(class(EP_learner) == "Stack"){
    EP_learner <- EP_learner$learners
  } else if(!is.list(EP_learner)) {
    EP_learner <- list(EP_learner)
  }
  all_names <- c()
  EP_learner <- lapply(EP_learner, function(lrnr) {
    params <- lrnr$params
    params$family <- EP_learner_spec$EP_family
    all_names <<- c(all_names,  lrnr$name[1])
    lrnr <- lrnr$clone()$reparameterize(params)

  })
  all_names <- as.vector(sapply( sapply(sieve_basis_generator_list, `[[`, "name" ), function(name){
    if(is.null(name)) {
      return(paste0(all_names, "_no_sieve"))
    }
    paste0(all_names, "_sieve_",name)
  }))


  EP_learner <- Stack$new(EP_learner)

  if(cross_validate) {
    EP_learner <- Lrnr_cv$new(EP_learner, full_fit = TRUE, folds = folds)
  }


   all_EP_learners <- lapply(seq_along(sieve_basis_generator_list), function(i) {
    lrnr <- EP_learner$clone()
     name <-  paste0(sapply(EP_learner$learners, `[[`, "name"), "_sieve_", sieve_basis_generator_list[[i]]$name)
     lrnr$.__enclos_env__$private$.name <-name
    return(delayed_learner_train(lrnr, all_tasks[[i]]))
  })



  all_EP_learners <- unlist(all_EP_learners)
  all_EP_learners <- bundle_delayed(all_EP_learners)
  all_EP_learners <- all_EP_learners$compute()
  #EP_learner <- Stack$new(all_EP_learners)
  output <- list(EP_learners_trained = all_EP_learners, all_tasks = all_tasks )
  class(output) <- "EP_learn"

  if(cross_validate) {
    cv_predictions <- do.call(cbind,lapply(all_EP_learners, function(EP_learner) {
      task <- all_tasks[[1]]
      EP_learner$predict_fold(task, "validation")
    }))


    colnames(cv_predictions) <- all_names

    output$cv_predictions <- cv_predictions
    output$cv_risks <- efficient_risk_function( cv_predictions, A , Y, EY1W, EY0W, pA1W, weights, EP_learner_spec$efficient_loss_function)

  }
  full_predictions <- do.call(cbind,lapply(all_EP_learners, function(EP_learner) {
    task <- all_tasks[[1]]
    EP_learner$predict_fold(task, "full")
  }))
  colnames(full_predictions) <- all_names
  output$full_predictions <- full_predictions

  return(output)

}
#' @export
predict.EP_learn <- function(output, V) {
  task <-sl3_Task$new(as.data.table(V), covariates  = colnames(V), outcome = c(), outcome_type = "continuous")

}


