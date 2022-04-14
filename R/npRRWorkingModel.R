



#' Function to compute initial estimates of nuisance functions.
#' @param Formula_LRR A formula specifying the working model for the log relative risk (or msm for the log relative risk)
#' @param W A column-named matrix of baseline variables.
#' @param A A binary vector with values in {0,1} encoding the treatment assignment.
#' @param Y A numeric vector storing the outcome values.
#' @param weights A numeric vector of weights
#' @param EY1W An optional numeric vector of initial estimates of `E[Y | A = 1, W]` for each individual.
#' @param EY0W An optional numeric vector of initial estimates of `E[Y | A = 0, W]` for each individual.
#' @param pA1W An optional numeric vector of initial estimates of `P(A = 1} W)` for each individual.
#' @param sl3_Learner_pA1W A \code{sl3_Learner} object from the \code{tlverse/sl3} R github package that specifies the machine-learning algorithm for learning the propensity score `P(A = 1 | W)`
#' @param sl3_Learner_EYAW A \code{sl3_Learner} object from the \code{tlverse/sl3} R github package that specifies the machine-learning algorithm for learning the outcome conditional mean `E[Y | A, W]`. NOTE: the treatment arms are pooled in the regression. See the preprocessing sl3_Learner \code{Lrnr_stratified} if you wish to stratify the estimation by treatment.
#' @param folds A number representing the number of folds to use in cross-fitting or a fold object from the package \code{tlverse/origami}. This parameter will be passed to internal \code{sl3_Task} objects that are fed to the code{sl3_Learner}s.
#' @param outcome_type Internal use only.
#' @export
npRRWorkingModel <- function(formula_LRR, W, A, Y, weights, EY1W, EY0W, pA1W, sl3_Learner_EYAW, sl3_Learner_pA1W, folds = 10, outcome_type = NULL) {
  try({
    likelihood <- estimate_initial_likelihood(W, A, Y, weights, sl3_Learner_EYAW = sl3_Learner_EYAW, sl3_Learner_pA1W = sl3_Learner_pA1W, folds = folds, outcome_type = outcome_type)
    EY1W <- likelihood$EY1
    EY0W <- likelihood$EY0
    pA1W <- likelihood$pA1
    })
  pA1W <- as.vector(pA1W)
  EY0W <- as.vector(EY0W)
  EY1W <- as.vector(EY1W)

  n <- length(A)
  V <- model.matrix(formula_LRR, as.data.frame(W))
  pA0W <- 1 - pA1W
  pAW <- ifelse(A==1, pA1W, pA0W)
  EYAW <- ifelse(A==1, EY1W, EY0W)



  # Compute initial EIF for variance estimation


  beta <- suppressWarnings(coef(glm.fit(V, EY1W, offset = log(EY0W), family = poisson(), weights = weights)))
  RR_beta <- as.vector(exp(V %*% beta))
  H <- V * (A / pA1W - (1 - A) * RR_beta * (1 / pA0W))

  scale <- apply(V, 2, function(v) {
    apply(weights * V * (v) * RR_beta * EY0W, 2, mean)
  })
  scaleinv <- solve(scale)

  EIF_Y_initial <- weights * (H %*% scaleinv) * as.vector(Y - EYAW)
  EIF_WA_initial <- apply(V, 2, function(v) {
    weights * (v * (RR_beta * EY0W - EY1W) - mean(weights * v * (RR_beta * EY0W - EY1W)))
  }) %*% scaleinv
  ses <- sqrt(diag(var(EIF_Y_initial + EIF_WA_initial)))

  max_iter <- 100
  max_eps <- 0.01
  for(i in seq_len(max_iter)) {
    beta <- suppressWarnings(coef(glm.fit(V, EY1W, offset = log(EY0W), family = poisson(), weights = weights)))
    RR_beta <- as.vector(exp(V %*% beta))

    H <- V * (A  - (1 - A) * RR_beta )
    H1 <- V
    H0 <- - V  *  RR_beta
    EIF_Y <- weights/pAW * (H %*% scaleinv) * as.vector(Y - EYAW)
    print(max(abs(colMeans(EIF_Y))))
    if(all(abs(colMeans(EIF_Y)) <= 0.5 * ses / sqrt(n) / log(n))) {
      #print(colMeans(EIF_Y))
      print("Converged.")
      break
    }
    dir <- colMeans(EIF_Y)
    dir <- dir / norm(dir, type = "2")
    H_1d <- H %*% dir
    # Log-link submodel
    apply_sub_model <- function(pred, H, eps) {
      return(as.vector(exp(log(pred) + H * eps)))
    }
    # Weighted poisson risk function
    risk_function <- function(eps) {
      EYAW_eps <- apply_sub_model(EYAW, H_1d, eps)
      risk <- mean(weights/pAW * (EYAW_eps - Y * log(EYAW_eps)))
      return(risk)
    }

    optim_fit <- optim(par = list(epsilon = c(rep(0, ncol(V)))), fn =  risk_function, method = "Brent", lower = 0, upper = max_eps)
    epsilon <- optim_fit$par
    EY1W <- apply_sub_model(EY1W, H1 %*% dir , epsilon)
    EY0W <- apply_sub_model(EY0W, H0 %*% dir , epsilon)
    EYAW <- ifelse(A==1, EY1W, EY0W)

  }

  tmle_scores <- max(abs(colMeans(EIF_Y)))
  beta_tmle <- suppressWarnings(coef(glm.fit(V, EY1W, offset = log(EY0W), family = poisson(), weights = weights)))
  EIF <- EIF_Y_initial + EIF_WA_initial
  Zscore <- abs(sqrt(n) * beta_tmle / ses)
  pvalue <- signif(2 * (1 - pnorm(Zscore)), 5)
  coefs <- data.frame(coef = beta_tmle, se = ses/sqrt(n),
                      ci_left = beta - 1.96*ses/sqrt(n),  ci_right = beta + 1.96*ses/sqrt(n),
                      Z_score = Zscore, p_value = pvalue
  )
  output <- list(coefficients = coefs, EIF = EIF, tmle_scores = tmle_scores )
}
