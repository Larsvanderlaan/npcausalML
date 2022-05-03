


#' For a given base LRR learner, train the LRR learners on all sieve-based plugin and IPW empirical risk functions
#' @param sl3_learner_binomial See argument \code{sl3_learner_binomial} of \link{estimate_LRR_using_ERM} for details.
#' @param V A matrix of observations of a subset of the covariates `W` for which to estimate the (possibly semi-marginalized) log relative risk (LRR).
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' @param EY1W A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0W A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1W A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
#' @param list_of_sieve_nuisances A list of sieve nuisance estimates as returned by the function \link{compute_plugin_and_IPW_sieve_nuisances}.
#' @param Vpred A matrix of covariates observations at which to predict the LRR. By default, \code{Vpred} equals \code{W}.
#' @param compute Whether to `compute` the list of delayed trained learners. See the package \code{delayed} for details.
train_learner_all_sieves <- function(sl3_Learner, V, A, Y, weights, family_risk_function, outcome_function_plugin, weight_function_plugin,  outcome_function_IPW, weight_function_IPW, transform_function, list_of_sieve_nuisances, Vpred = V, compute = FALSE) {
  list_of_sieve <- lapply(list_of_sieve_nuisances, function(sieve_nuisances) {
    EY1W_star <- sieve_nuisances$EY1W_star
    EY0W_star <- sieve_nuisances$EY0W_star
    pA1W_star <- sieve_nuisances$pA1W_star
    #print(data.table(sieve = EY1W_star))

    delayed_plugin_ERM <- delayed_fun(estimate_using_ERM)(V, A, Y,  EY1W_star, EY0W_star, pA1W_star, weights, family_risk_function, sl3_Learner,  outcome_function_plugin, weight_function_plugin, outcome_function_IPW, weight_function_IPW , learning_method = c("plugin"), Vpred = Vpred, transform_function = transform_function)


    delayed_IPW_ERM <- delayed_fun(estimate_using_ERM)(V, A, Y,  EY1W_star, EY0W_star, pA1W_star, weights, family_risk_function, sl3_Learner,  outcome_function_plugin, weight_function_plugin, outcome_function_IPW, weight_function_IPW , learning_method = c("IPW"), Vpred = Vpred, transform_function = transform_function)

    output <- (list( plugin = delayed_plugin_ERM, IPW = delayed_IPW_ERM))
  })
  sieve_names <- sapply(list_of_sieve_nuisances, `[[`, "sieve")
  names(list_of_sieve) <- sieve_names
  if(compute) {
    list_of_sieve <- bundle_delayed(unlist(list_of_sieve))
    list_of_sieve <- list_of_sieve$compute()
  }
  return(list_of_sieve)
}

#' Train the LRR learners for all sieves and base learners on the full data.
#' @param V A matrix of covariate observations.
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' @param EY1W A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0W A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1W A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
#' @param list_of_sieve_nuisances A list of sieve nuisance estimates as returned by the function \link{compute_plugin_and_IPW_sieve_nuisances}.
#' @param list_of_learners A list of untrained \code{sl3_Learner} objects to be used to estimate the log relative risk LRR using the function \link{estimate_LRR_using_ERM}.
#' @param list_of_sieves A list of basis_generator objects specifying the sieve. See, for example, \code{fourier_basis} for an example and template.
#' @param Vpred A matrix of covariates observations at which to predict the LRR. By default, \code{Vpred} equals \code{W}.
#' @param compute Whether to `compute` the list of delayed trained learners. See the package \code{delayed} for details.
train_learners <- function(V, A, Y, EY1W, EY0W, pA1W, weights, family_risk_function, outcome_function_plugin, weight_function_plugin,  outcome_function_IPW, weight_function_IPW, transform_function, design_function_sieve_plugin, weight_function_sieve_plugin, design_function_sieve_IPW, weight_function_sieve_IPW, family_for_targeting,  list_of_learners, list_of_sieves, Vpred = V, compute = TRUE) {
    #print("train_learners")

  list_of_sieve_nuisances <- lapply(list_of_sieves, function(sieve){
    compute_plugin_and_IPW_sieve_nuisances(basis_generator = sieve, V = V, A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W, weights = weights, design_function_sieve_plugin = design_function_sieve_plugin, weight_function_sieve_plugin = weight_function_sieve_plugin, design_function_sieve_IPW = design_function_sieve_IPW, weight_function_sieve_IPW = weight_function_sieve_IPW, family_for_targeting = family_for_targeting)})
  all_learners_delayed <- lapply(list_of_learners, train_learner_all_sieves, list_of_sieve_nuisances = list_of_sieve_nuisances, V = V, A = A, Y = Y, weights = weights, family_risk_function = family_risk_function, outcome_function_plugin = outcome_function_plugin, weight_function_plugin = weight_function_plugin,  outcome_function_IPW = outcome_function_IPW, weight_function_IPW = weight_function_IPW, transform_function = transform_function, Vpred = Vpred)
  learner_names <- names(list_of_learners)
  names(all_learners_delayed) <- paste0(learner_names)
  if(compute) {
    all_learners_delayed <- bundle_delayed(unlist(all_learners_delayed))
    suppressWarnings(suppressMessages(all_learners_delayed <- all_learners_delayed$compute()))
  }
  return(all_learners_delayed)
}

#' For each LRR learner, the best sieve is chosen by minimizing the double-robust one-step efficient empirical risk function over the choices of sieve-spaces.
#' @param trained_learner_list A unnested list of trained LRR learners as
#' @param learner_names Names of LRR learners. Should be `learner_names <- sapply(list_of_learners, `[[`, "name)`.
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' @param EY1W A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0W A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1W A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
subset_best_sieve <- function(trained_learner_list, learner_names, A, Y, EY1W, EY0W, pA1W, weights, efficient_loss_function, V = NULL) {
  learners <- trained_learner_list


  learners <- lapply(learner_names, function(learner_name) {
    print(learner_name)
    keep <- which(stringr::str_detect(names(learners), quotemeta(learner_name)))
    sieve_learners <- learners[keep]

    # Get LRR predictions on fold-specific training set for all sieves
    # Single learner

    all_ERM_training <- lapply(sieve_learners, `[[`, "ERM_train")

    ntrain <- nrow(all_ERM_training[[1]])
    tmp <- do.call(rbind, all_ERM_training)
    list_of_sieve_training <- lapply(seq_len(ncol(tmp)), function(j) {
      matrix(tmp[,j], nrow = ntrain, byrow = F)
    })

    all_ERM_pred <- lapply(sieve_learners, `[[`, "ERM_pred")
    npred <- nrow(all_ERM_pred[[1]])
    tmp <- do.call(rbind, all_ERM_pred)
    list_of_sieve_pred <- lapply(seq_len(ncol(tmp)), function(j) {
      matrix(tmp[,j], nrow = npred, byrow = F)
    })

    all_best_index <- rep(NA, length(list_of_sieve_training))
    all_best_ERM_training <- as.list( rep(NA, length(list_of_sieve_training)))
    all_best_ERM_pred <- as.list( rep(NA, length(list_of_sieve_training)))



    lapply(seq_along(list_of_sieve_training), function(j) {
      sieve_ERM_training <- list_of_sieve_training[[j]]

      all_training_risks <-  efficient_risk_function(sieve_ERM_training, A, Y, EY1W, EY0W, pA1W, weights, efficient_loss_function = efficient_loss_function, V = V)

      print(all_training_risks)

      best_index <- which.min(all_training_risks)
      print(best_index)
      all_best_ERM_training[[j]] <<- sieve_ERM_training[,best_index]

      all_best_ERM_pred[[j]] <<- list_of_sieve_pred[[j]][,best_index]
      all_best_index[j] <<-  best_index
      return(NULL)
    })

    all_best_ERM_training <- do.call(cbind, all_best_ERM_training)
    colnames(all_best_ERM_training) <- colnames(all_ERM_training[[1]])

    all_best_ERM_pred <- do.call(cbind, all_best_ERM_pred)
    colnames(all_best_ERM_pred) <- colnames(all_ERM_pred[[1]])

    best_sieve_learners <- lapply(all_best_index, function(best_index) {
      #print("best_index")
      #print(all_best_index)
      #print(names(sieve_learners[[best_index]]))
      sieve_learners[[best_index]]$ERM_learner
    })
    return(list(ERM_train = all_best_ERM_training, ERM_pred = all_best_ERM_pred, ERM_learner = best_sieve_learners, choice = names(sieve_learners)[all_best_index]))


  })
  names(learners) <- learner_names
  learners
}


