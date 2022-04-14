compute_TMLE <- function(initial_LRR, task_RR, tmle_task, likelihood, points = NULL, npoints = 25, min_samples = 250, max_samples = Inf, sigmas_grid = NULL, VarOverBias = 4, kernel_function = NULL) {

  task <- tmle_task
  if(is.null(points)) {
    points <- quantile(unlist(task_RR$X), seq(0.05, 0.95, length.out = npoints))
  }
  if(ncol(task_RR$X) !=1) {
    stop("Inference only available for one dimensional `V`.")
  }
  ys <- points
  V <- unlist(task_RR$X)
  task <- tmle_task
  preds <- initial_LRR
  sigmas <- c()
  previous_best <- NULL
  if(!is.null(sigmas_grid)) {
    grid <- sigmas_grid
  } else {
    grid <- sd(V) *  sort(c(0.00001, 0.00005, 0.0002, 0.0005, 0.001, 0.002, 0.003,  c(seq(0.0001, 2.5, length = 100), 2.6, 2.7, 2.8, 3 )))
  }

  if(length(grid) > 1) {
    for(y in ys) {
      sigma <- get_h(V, y, preds, task_RR, task, likelihood, grid = grid,  previous_best = previous_best, min_samples = min_samples, max_samples = max_samples, VarOverBias = VarOverBias, kernel_function = kernel_function)
      previous_best <- sigma
      sigmas <- c(sigmas, sigma)
    }
  } else {
    sigmas <- rep(sigmas_grid, length(ys))
  }

  kernel <- function(x) {
    fun <- function(i) {
      y <- c(ys)[i]
      sigma <- c(sigmas)[i]
      v <- kernel_function(x,y, sigma) #exp(-(y-x)^2/(2*sigma^2)) /  sqrt(2*3.14*sigma^2)
      v
    }
    sapply(1:length(ys), fun)
  }
  c <- kernel(V)

  Q1V <- task_RR$get_data(,"Q1V")[[1]]
  Q0V <- task_RR$get_data(,"Q0V")[[1]]
  Q1 <- task_RR$get_data(,"Q1")[[1]]
  Q0 <- task_RR$get_data(,"Q0")[[1]]
  g1 <- task_RR$get_data(,"g1")[[1]]
  A <- task_RR$get_data(,"A")[[1]]
  Y <- task_RR$get_data(,"Y")[[1]]
  Delta <- task_RR$get_data(,"Delta")[[1]]
  G <- task_RR$get_data(,"G")[[1]]

  update <- preds
  RR <- exp(update)
  clever_covariates <-   (Q1V + Q0V)/(Q1V*Q0V) * kernel(V)
  EIC <- (Delta/G) * A/(g1) *(-1/(1+ RR))* (Y-Q1) + (Delta/G) * (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1V) + (RR/(1+ RR))*(Q0V) - (1/(1+ RR))*(Q1 - Q1V) + (RR/(1+ RR))*(Q0 - Q0V)
  weights <- apply(as.vector(EIC)*clever_covariates,2,sd)
  eff_loss <- make_eff_loss(task, likelihood)

  max_eps <- 5e-3
  for(i in 1:50) {
    RR <- exp(update)
    clever_covariates <-   (Q1V + Q0V)/(Q1V*Q0V) * kernel(V)
    EIC <- (Delta/G) * A/(g1) *(-1/(1+ RR))* (Y-Q1) + (Delta/G) * (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1V) + (RR/(1+ RR))*(Q0V) - (1/(1+ RR))*(Q1 - Q1V) + (RR/(1+ RR))*(Q0 - Q0V)
    dir <- colMeans(as.vector(EIC)*clever_covariates) / weights
    dir <- dir/sqrt(mean(dir^2) )
    clever_covariates_one_dim <- clever_covariates %*% (dir / colMeans(kernel(V)))

    risk <- function(epsilon) {

      loss <- (eff_loss(update  + epsilon*clever_covariates_one_dim ))
      loss <- loss
      risk <- mean(loss)
    }

    optim_fit <- optim(
      par = list(epsilon = 0 ), fn = risk,
      lower = -max_eps, upper = max_eps,
      method = "Brent")

    epsilon <- optim_fit$par
    update <- update + epsilon* (clever_covariates%*% dir)

    RR <- exp(update)
    kern <- kernel(V)

    score1 <- (Delta/G) * A/(g1) *(-1/(1+ RR))* (Y-Q1) + (Delta/G) * (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1V) + (RR/(1+ RR))*(Q0V) - (1/(1+ RR))*(Q1 - Q1V) + (RR/(1+ RR))*(Q0 - Q0V)

    score1a <- clever_covariates*as.vector(score1)
    print("SCORE")

    if(i >=5 && sqrt(mean((colMeans(score1a)/ weights)^2)) <= 1/sqrt(n)/log(n)){
      break
    }
    if(abs(optim_fit$par) < 1e-8) {
      max_eps <- max_eps/3
    }
  }

  est_X <- kern * as.vector(update)
  psi <- colMeans(est_X)
  EIC <- score1a + t(t(est_X) - psi)
  EIC <- t(t(EIC) / colMeans(kern)) + t( -(psi/colMeans(kern)^2)*(t(kern) - colMeans(kern)))
  psi <- psi/colMeans(kern)
  se <- apply(EIC,2,sd)/sqrt(n)
  radius <- 1.96*se
  upper <- exp(psi + radius)
  lower <- exp(psi - radius)
  psi <- exp(psi)
  ests <- cbind(points, psi, lower, upper, sigmas, log(psi),  se)
  colnames(ests) <- c("V", "RR", "lower_CI", "upper_CI", "kernel_sigmas", "LRR", "se_LRR/sqrt(n)")
  output <- list(estimates = ests, RR_targeted = exp(update), EIC = EIC)
  return(output)
}

get_h <- function(V, y, preds, task_RR, tmle_task, likelihood, grid, previous_best = NULL, min_samples = 250, max_samples = Inf, VarOverBias = 4, kernel_function = NULL) {


  task <- tmle_task
  eff_loss <- make_eff_loss(task, likelihood)
  grid <- rev(sort(grid))
  min_nonzero <- min(length(V),min(min_samples, max(length(V) / 10, 30)))
  max_nonzero <- max_samples
  kernel_support <- function(sigma) {
    v <- exp(-(y-V)^2/(2*sigma^2)) /  sqrt(2*3.14*sigma^2)
    #num_support <- sum(kernel_function(V,y, sigma) > 1e-6)
     num_support <- sum(abs(y-V) <= 1.75*sigma)
    return(num_support)
  }
  supports <- sapply(grid, kernel_support)
  keep <-supports >= min_nonzero

  tmp <- grid[which.min(abs(supports-min(250, length(V)/5)))]

  grid <- grid[keep]

  kernel <- function(x,y) {
    v <- exp(-(y-x)^2/(2*tmp^2)) /  sqrt(2*3.14*tmp^2)
    v
  }
  left <- mean(kernel(V,y  - tmp) * preds) / mean(kernel(V,y - tmp) )
  center <- mean(kernel(V,y ) * preds) / mean(kernel(V,y ) )
  right <- mean(kernel(V,y + tmp ) * preds) / mean(kernel(V,y + tmp) )
  decreasing <- mean(c(right - center, center - left)) <= 0
  previous_best <- which.min(abs(previous_best - grid))

  Q1V <- task_RR$get_data(,"Q1V")[[1]]
  Q0V <- task_RR$get_data(,"Q0V")[[1]]
  Q1 <- task_RR$get_data(,"Q1")[[1]]
  Q0 <- task_RR$get_data(,"Q0")[[1]]
  g1 <- task_RR$get_data(,"g1")[[1]]
  A <- task_RR$get_data(,"A")[[1]]
  Y <- task_RR$get_data(,"Y")[[1]]
  Delta <- task_RR$get_data(,"Delta")[[1]]
  G <- task_RR$get_data(,"G")[[1]]
  upper_CIs <- c()
  lower_CIs <- c()
  ests <- c()
  u <- max(V)
  l <- min(V)

  for(sigma in grid) {
    update <- preds
    kernel <- function(x) {
      fun <- function(y) {
        v <- kernel_function(x,y, sigma) #exp(-(y-x)^2/(2*sigma^2)) /  sqrt(2*3.14*sigma^2)
        v
      }
      return(as.matrix(fun(y)))

    }

    for(i in 1:1) {
      clever_covariates <-   (Q1V + Q0V)/(Q1V*Q0V) * kernel(V)

      risk <- function(epsilon) {
        # RR <- exp(preds + epsilon*1)
        # score <- A/(g1) *(-1/(1+ RR))* (Y-Q1) + (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1) + (RR/(1+ RR))*(Q0)
        # score <- abs(mean(score * clever_covariates))
        loss <- (eff_loss(update  + epsilon*clever_covariates ))
        loss <- loss

        risk <- mean(loss)
      }


      max_eps <- 0.05
      optim_fit <- optim(
        par = list(epsilon = 0 ), fn = risk,
        lower = -max_eps, upper = max_eps,
        method = "Brent")

      epsilon <- optim_fit$par

      update <- update + epsilon* (clever_covariates)
    }
    RR <- exp(update)
    kern <- kernel(V)
    score1 <- (Delta/G) * A/(g1) *(-1/(1+ RR))* (Y-Q1) + (Delta/G) * (1-A)/((1-g1))*(RR/(1+ RR))*(Y - Q0) - (1/(1+ RR))*(Q1V) + (RR/(1+ RR))*(Q0V) - (1/(1+ RR))*(Q1 - Q1V) + (RR/(1+ RR))*(Q0 - Q0V)
    score1a <- clever_covariates*as.vector(score1)




    est_X <- kern * as.vector(update)
    psi <- colMeans(est_X)
    EIC <- score1a + t(t(est_X) - psi)
    psi <- psi/ mean(kern)
    EIC <- EIC/mean(kern) - (psi/mean(kern)^2)*(kern - mean(kern))
    radius <- VarOverBias*apply(EIC,2,sd)/sqrt(n)
    upper <- psi + radius
    lower <- psi - radius
    upper_CIs <- c(upper_CIs, upper)
    lower_CIs <- c(lower_CIs, lower)
    ests <- c(ests, psi)

  }
  # Estimates increasing
  len <- round(length(ests)*0.8)
  start <- round(len/2) + 1
  increasing <- mean(diff(ests)[10:len]) >= 0
  decreasing <- mean(diff(ests)[10:len]) <= 0

  slopes <- sign(diff(ests))[-1]
 print(slopes)
  local_optima <- slopes*(diff(sign(diff(lower_CIs))) + diff(sign(diff(upper_CIs))))
   print(local_optima)



  best_index <- which(local_optima == -2) + 1
  best_index <- max(best_index)
  # if(increasing) {
  #   best_index <- which(diff(sign(diff(lower_CIs)))==-2)+1
  #   best_index <- min(best_index)
  #
  # } else if(decreasing) {
  #
  #   best_index <- which(diff(sign(diff(-upper_CIs)))==-2)+1
  #   best_index <- min(best_index)
  #
  # } else {
  #   # print("FLAT")
  #   # best_index <- which(diff(sign(diff(-(upper_CIs - lower_CIs))))==-2)+1
  #   # best_index <- min(best_index)
  # }
  if(is.infinite(best_index)) {


    if(is.null(previous_best)) {
      best_index <- round(length(grid)*0.9)
    } else {
      best_index <- min(previous_best + 1  , length(grid))
    }

  }
  if(abs(mean(diff(ests)[10:len])) < 1e-3) {
    best_index <- min(best_index, previous_best + 15)
    best_index <- max(best_index, previous_best - 15)

  }
  best_index <- min(best_index, previous_best + 15)
  best_index <- max(best_index, previous_best - 15)
  print(best_index)

  return(grid[best_index])
}
