

#'
#' Machine learning of the covariate-adjusted conditional relative risk function.
#'
#' This function allows one to data-adaptively estimate the conditonal relative risk function using machine learning.
#' Given two sets of covariates `V` and `X`, a binary treatment `A`, and a binary outcome variable`Y`. This function estimates `E[E[Y=1|A=1, V, X] |V]/ E[E[Y=1|A=0, V, X] |V]`. When `X` is omitted, this reduces to `E[Y=1|A=1, V]/ E[Y=1|A=0, X]``
#'
#'
#' @param V A vector or matrix of voariates for which to estimate the adjusted conditional relative risk function.
#' @param X A vector or matrix of variables to adjust for/marginalize over when computing the conditional relative risk function. If NULL then the unadjusted conditional relative function is computed using \cpde{V}
#' @param A A binary vector representing treatment. A = 1 corresponds with recieving treatment and A = 0 corresponds with no treatment. The relative risk is computing comparing as the quotient of the A=1 risk and A=0 risk.
#' @param Y A binary vector representing the outcome variable.
#' @param weights A vector of weight values to incorporate in all estimation procedures.
#' @param library_RR An \code{sl3} \code{Learner} or \code{Stack}  object representing the machine-learning library for estimating the relative risk function. See the R github package \code{\link[sl3]} at \url{https://github.com/tlverse/sl3/tree/master/R}
#' @param library_A An \code{sl3} \code{Learner} or \code{Stack}  object representing the machine-learning library for estimating the distribution of \code{A} nuisance parameter, P(A=1|X, V).
#' @param library_Y An \code{sl3} \code{Learner} or \code{Stack}  object representing the machine-learning library for estimating the distribution of \code{Y} nuisance parameter, P(Y=1|A, X, V).
#' @param library_V An optional \code{sl3} \code{Learner} or \code{Stack}  object representing the machine-learning library for estimating the sequential regression E[E[Y|A=a, V, X]|X]. This need only be supplied when \code{V} is one-dimensional and one wishes to call \code{\link{npRR_inference}}.
#' @param cv_nuisance A boolean for whether or not to use cross-validation/cross-validated predictions for the nuisance parameter estimation. This must be true if \code{library_A} or \code{library_V} are \code{sl3} \code{Stack} objects.
#' If \code{library_A} or \code{library_V} are already cross-validated (Super)learners then it is computationally beneficial to set this parameter tp \code{FALSE}.
#' @param cv_RR A boolean for whether to return the cross-validated-selected predictions for the relative risk obtained by cross-validating \code{library_RR}.
#' @param basis_list Internal-use variable. list of sieve basis objects to use to convexify the double-robust RR-risk function. See \code{fourier_basis} object as an example.
#' @param likelihood A \code{tmle3} fitted likelihood object (as returned by this function) to use as estimates for the nuisance parameters.
npRR <- function(V, X = NULL, A, Y, weights = NULL, library_RR = "autoML", library_A = "autoML", library_Y = "autoML", library_V =  "autoML", cv_nuisance = TRUE, cv_RR = TRUE, basis_list = make_fourier_basis_list(V), likelihood = NULL) {


  V <- as.matrix(V)
  X <- cbind(V,as.matrix(X))
  duplicated.columns <- duplicated(t(X))
  X <- X[, !duplicated.columns]
  A <- as.vector(A)
  Y <- as.vector(Y)

  # if(is.character(library_V) && library_V == "autoML" && ncol(V)== ncol(X) && ncol(X) > 1) {
  #   library_V <- NULL
  # }
  if(ncol(V)>1){
    library_V <- NULL
  }
  # if(is.list(library_RR)) {
  #   for(i in seq_along(library_RR)) {
  #     lrnr <- library_RR[[i]]
  #     if(!("LRR" %in% lrnr$properties)) {
  #       if(!(all(c("weights", "binomial") %in% lrnr$properties))) {
  #         warning(paste0("RR Learner ", lrnr$name, " does not have `weights` or `binomial` in properties. Make sure the learner supports these."))
  #       }
  #       library_RR[[i]] <- binomial_to_LRR_learner(lrnr)
  #     }
  #   }
  # }


  if(library_RR == "autoML") {
    library_RR <- default_library_RR
  }

  if(is.character(library_A) && library_A == "autoML") {
    library_A <-  default_library_nuisance
  }
  if(is.character(library_Y) && library_Y == "autoML") {
    library_Y <-  default_library_nuisance
  }
  if(ncol(V) > 10) {
    screen.glmnet <- TRUE
  }
  if(is.character(library_V) && library_V == "autoML") {
    library_V <-  default_library_V
    if(screen.glmnet) {
      library_V <- make_learner(Pipeline, Lrnr_LRR_glmnet.screener$new(), library_V)
    }
  } else if (is.character(library_V) &&library_V == "NULL") {
    library_V <- NULL
  }
  task <- make_task(V, X, A, Y, weights = weights, folds = 10)
  print("Fitting likelihoods")
  if(is.null(likelihood)) {
    likelihood <- make_likelihood(task, library_A ,library_Y, library_V, cv = cv_nuisance )
  }
  print("Done Fitting likelihoods")
  genr <- make_generator(likelihood)
  task_RR <- genr(task, "validation")

  base_lrnr_lrr <- generate_plugin_learners(task_RR, basis_list, library_RR, task, likelihood, select_sieve = TRUE, screen.glmnet = screen.glmnet, cv = cv_RR)
  lrnr_lrr <- base_lrnr_lrr$train(task_RR)
  predictions <- exp(lrnr_lrr$predict(task_RR))
  if(cv_RR) {
    reduced_lrnr <- lrnr_lrr$fit_object$full_fit$fit_object$learner_fits$Stack$fit_object$learner_fits[[which(out$learner_LRR_trained$fit_object$cv_meta_fit$coefficients==1)]]
    cv_predictions <- exp(lrnr_lrr$predict_fold(task_RR, "validation"))
  } else {
    reduced_lrnr <- lrnr_lrr
    cv_predictions <- NULL
  }
  output <- list(RR = predictions, tmle_task = task, task_RR = task_RR, likelihood = likelihood, learner_LRR_trained = lrnr_lrr, learner_LRR_trained_squashed = reduced_lrnr, RR_cv = cv_predictions)
  class(output) <- "npRR"
  return(output)
}


#' Function that provides efficient estimates and inference for a data-adaptively kernel smoothed approximation of the covariate-adjusted conditional relative risk function. This method only works when the argument \code{V} provided to \code{\link{npRR}}is univariate.
#'
#' @param output An output object returned by \code{\link{npRR}}.
#' @param kernel_function The kernel function to use for kernel smoothing. By default, the guassian kernel is smooth.
#' Must be a function of the form `f(x,y,sigma)` where `x` is a vector of values to evaluate the kernel function at, `y` is the centering value of the kernel, and `sigma` is a real-valued smoothness parameter/bandwidth of the kernel.
#' @param npoints The number of points at which to estimate the kernel-smoothed relative risk function (equally spaced on the quantile scale).
#' @param points A vector of points at which to estimate the kernel smoothed RR estimate. This will override \code{npoints} if supplied.
#' @param sigmas_grid A vector of bandwidth values of the kernel to search through. Lepsky's method is used to data-adaptively select an optimal bandwidth.
#' @param min_samples The minimum number of samples to be contained in the support of the kernel. This is used to filter out values of \code{sigma_grid} that lead to very narrow kernels.
#' @param max_samples The maximum number of samples to be contained in the support of the kernel. This is used to filter out values of \code{sigma_grid} that lead to very wide kernels.
#' @param VarOverBias The ratio of the variance and bias of the kernel-smoothed relative risk estimator at a point relative to the true relative risk function.
#' High values of \code{VarOverBias} lead to undersmoothed, more variable and less biased estimates of the relative risk function. Under suitable conditions, undersmoothing the kernel can lead to the kernel-smoothed relative risk estimates and inference being correct for the true unsmoothed relative risk function.
#' We stress that the estimates and inference are for the kernel smoothed relative risk function, and not necessarily for the true relative risk function.
#'
#' @return A list containing:
#' 1. A matrix that contains point-estimates and 95`%` pointwise confidence intervals for the kernel-smoothed relative risk function for a number of points. The bandwidth `sigma` used for the kernel smooth at each point is also given.
#' 2. A targeted/updated (TMLE) estimator for the relative risk evaluated at \code{V}
#' 3. A matrix containing the estimated efficient influence function for the kernel-smoothed relative risk function point evaluations. This can be used to construct simultaneous confidence bands.
#'
npRR_inference <- function(output, kernel_function = function(x,y, sigma){exp(-(y-x)^2/(2*sigma^2)) /  sqrt(2*3.14*sigma^2)}, npoints = 35, points = NULL, sigmas_grid = NULL, min_samples = 100, max_samples = Inf, VarOverBias = 3) {
  #npoints <- max(npoints, 20)
  V <- as.matrix(output$task_RR$X)
  if(!is.null(output$RR_cv)) {
    LRR <- log(output$RR_cv)
  } else {
    LRR <- log(output$RR)
  }
  output_inference <- compute_TMLE(LRR, output$task_RR, output$tmle_task, output$likelihood, points = points, npoints = npoints, sigmas_grid = sigmas_grid, min_samples= min_samples, max_samples = max_samples, VarOverBias = VarOverBias, kernel_function = kernel_function)
  df <- as.data.frame(output_inference$estimates)
  plt <- ggplot(df, aes(x = V, y = RR)) + geom_line() + geom_ribbon(aes(ymin=lower_CI, ymax=upper_CI), alpha=0.2) + xlab("Value of V") + ylab("Relative Risk")
  output_inference$plot <- plt
  return(output_inference)
}

#' Predict relative risk values at new observations.
#'
#' Given an npRR object returned by \code{\link{npRR}}, this functions returns the relative risk values for given values of `V`.
#' @param object A \code{\link{npRR}} object
#' @param newV A vector or matrix values at which to evaluate the relative risk function. The order of columns must match the argument \code{V} passed to \code{\link{npRR}} exactly.
#' @return Returns a vector or matrix of predicted RR values.
predict.npRR <- function(object, newV) {
  newV <- as.matrix(newV)
  Anew <- rep(1, nrow(newV))
  Ynew <- rep(1, nrow(newV))
  X <- as.matrix(object$tmle_task$get_tmle_node("X"))
  Xnew <- apply(X, 2, function(x) {
    if(length(x) < nrow(newV)) {
      x <- suppressWarnings(x + rep(0, nrow(newV)))
    } else {
      x <- x[1:nrow(newV)]
    }
    x <- as.matrix(x)
    return(x)
  })
  if(!is.matrix(Xnew)) {
    Xnew <- t(as.matrix(Xnew))
  }

  colnames(Xnew) <- colnames(X)


  task <- make_task(newV, Xnew, Anew, Ynew, folds = 10)

  likelihood <- object$likelihood
  genr <- make_generator(likelihood)
  task_RR <- genr(task, "validation")
  output <- as.matrix(exp(object$learner_LRR_trained_squashed$predict(task_RR)))
  if(ncol(output)==1){
    output <- as.vector(output)
  }
  return(output)
}



