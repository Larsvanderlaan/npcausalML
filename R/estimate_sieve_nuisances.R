



#' Computes the data-adaptive sieve-based nuisance estimators of E[Y|A,W] and P(A=1|W)
#' that are constructed to be used with the IPW and plugin empirical risk functions.
#' @param V A matrix of observations of a subset of the covariates `W` for which to estimate the (possibly semi-marginalized) log relative risk (LRR).
#' The sieve basis will only be generated for these covariates.
#' @param A A binary vector specifying the treatment assignment. The values should be in {0,1}.
#' @param Y A numeric vector of binary or nonnegative observations of the outcome variable.
#' @param EY1W A numeric vector containing initial cross-fitted estimates of E[Y|A=1,W] for all observations.
#' @param EY0W A numeric vector containing initial cross-fitted estimates of E[Y|A=0,W] for all observations.
#' @param pA1W A numeric vector containing initial cross-fitted estimates of P(A=1|W) for all observations.
#' @param weights A numeric vector of observation weights. If no special weighting desired, supply a vector of 1's.
#' @param basis_generator A basis_generator object (see \code{fourier_basis} and \code{bspline_basis}) ...
#' that specifies a sieve basis transformation of the design matrix \code{W}.
#' @param family An R \code{family} object that specifies the link and risk function used in the data-adaptive sieve update step.
#' `family` should be a \code{binomial} object if `Y` is binary and a \code{poisson} object if `Y` is non-binary and nonnegative (e.g. a count).
#' @param debug ...
compute_plugin_and_IPW_sieve_nuisances <- function(V, A, Y, EY1W, EY0W, pA1W, weights, basis_generator, design_function_sieve_plugin, weight_function_sieve_plugin, design_function_sieve_IPW, weight_function_sieve_IPW, family_for_targeting = binomial(), debug = FALSE) {
  #print("compute_plugin_and_IPW_sieve_nuisances")
  if(all(Y %in% c(0,1))) {
    family_for_targeting <- binomial()
  } else if(all(Y >=0)) {
    family_for_targeting <- poisson()
  } else {
    family_for_targeting <- gaussian()
  }
  if(is.null(basis_generator)) {
    return(list(pA1W_star = pA1W, EY1W_star = EY1W, EY0W_star = EY0W, sieve = "no_sieve"))
  }
  # Compute sieve-transformed design matrix
  basis_generator <- basis_generator$clone()

  X <- basis_generator$set(V)$eval(V)
  # Compute data-adaptive sieve
  X_pseudo_plugin <- design_function_sieve_plugin(X,A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
  X1_pseudo_plugin <- design_function_sieve_plugin(X,A = rep(1, length(A)), Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
  X0_pseudo_plugin <- design_function_sieve_plugin(X,A = rep(0, length(A)), Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)

  weights_pseudo_plugin <- weights * weight_function_sieve_plugin(A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
  #X_plugin <- cbind(A*X, (1-A)*X)
  EY <- ifelse(A==1, EY1W, EY0W)
  #pA0 <- 1 - pA1W
  #pA <- ifelse(A==1, pA1W, pA0)
  #print("plugin")

  upper_bound <- max(c(Y, EY1W, EY0W))
  lower_bound <- min(c(Y, EY1W, EY0W))
  lower_bound <- min(0.9 * lower_bound, 1.1 * lower_bound)
  upper_bound <- max(0.9 * upper_bound, 1.1 * upper_bound)

  Y_scaled <- (Y - lower_bound) / (upper_bound - lower_bound)
  EY1W_scaled <- (EY1W - lower_bound) / (upper_bound - lower_bound)
  EY0W_scaled <- (EY0W - lower_bound) / (upper_bound - lower_bound)
  EY_scaled <- (EY - lower_bound) / (upper_bound - lower_bound)

  # print(quantile(EY_scaled))
  # print(quantile(EY0W_scaled))
  # print(quantile(EY1W_scaled))
  # print(quantile(Y_scaled))
  # print(any(is.na(Y_scaled) | is.infinite(Y_scaled)))
  # print(any(is.na(X_pseudo_plugin) | is.infinite(X_pseudo_plugin)))
  # print(any(is.na(weights_pseudo_plugin) | is.infinite(weights_pseudo_plugin)))
  # print(any(is.na(qlogis(EY_scaled)) | is.infinite(qlogis(EY_scaled))))

  suppressWarnings(sieve_fit_plugin <- glm.fit(X_pseudo_plugin, Y_scaled, weights = weights_pseudo_plugin, offset = qlogis(EY_scaled), family = binomial(), intercept = F))

  beta_plugin <- coef(sieve_fit_plugin)
  beta_plugin[is.na(beta_plugin)] <- 0


  EY1W_star <- lower_bound + (upper_bound - lower_bound) * as.vector(plogis(qlogis(EY1W_scaled) + X1_pseudo_plugin %*% beta_plugin))
  EY0W_star <- lower_bound + (upper_bound - lower_bound) * as.vector(plogis(qlogis(EY0W_scaled) + X0_pseudo_plugin %*% beta_plugin))
  #print("sieve")

  # print(quantile(EY1W_star))
  # print(quantile(EY0W_star))
  # print( range ((upper_bound - lower_bound) * as.vector(plogis(qlogis(EY1W_scaled) + X1_pseudo_plugin %*% beta_plugin))))
  #
  #   print(c(lower_bound, upper_bound))
  #print(quantile(as.vector(plogis(qlogis(EY1W_scaled) + X1_pseudo_plugin %*% beta_plugin))))
  #print(quantile(as.vector(plogis(qlogis(EY0W_scaled) + X0_pseudo_plugin %*% beta_plugin))))

  if(F) {
    EYstar <- ifelse(A==1, EY1W_star, EY0W_star )
    print("Sieve scores plugin")
    print(dim(X_pseudo_plugin))
    print(quantile(abs(colMeans(weights_pseudo_plugin*X_pseudo_plugin*(Y - EYstar)))))

  }

  #X_IPW <- cbind(EY1W/pA1W * X, EY0W/pA0 * X)
  X_pseudo_IPW <- design_function_sieve_IPW(X, A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)
  weights_pseudo_IPW <- weights * weight_function_sieve_IPW(A = A, Y = Y, EY1W = EY1W, EY0W = EY0W, pA1W = pA1W)

  suppressWarnings(sieve_fit_IPW <- glm.fit(X_pseudo_IPW, A, weights = weights_pseudo_IPW, offset = qlogis(pA1W), family = binomial(), intercept = F))
  beta_IPW <- coef(sieve_fit_IPW)
  beta_IPW[is.na(beta_IPW)] <- 0
  pA1W_star <- as.vector(plogis(qlogis(pA1W) + X_pseudo_IPW %*% beta_IPW))
  pA1W_star <- pmax(pmin(pA1W_star, 0.995), 0.005)
  if(debug) {
    print("Sieve scores IPW")
    print(colMeans(weights_pseudo_IPW*X_pseudo_IPW*(A - pA1W_star)))
  }

  output <- list(pA1W_star = pA1W_star, EY1W_star = EY1W_star, EY0W_star = EY0W_star, sieve = basis_generator$name)
}








