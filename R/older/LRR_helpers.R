

#' Generates a default \codE{basis_list} using the fourier basis.
#' @param V A vector or matrix of values of `V`. See the corresponding argument in \code{npRR}.
#' @param The maximum size of the list.
#' @returns A list of \code{fourier_basis} objects.
make_fourier_basis_list <- function(V, size = Inf) {
  V <- as.matrix(V)
  if(ncol(V) == 1) {
    basis1 <- fourier_basis$new(orders = c(1,0,0), max_degrees = c(1,2,3))
    basis2 <- fourier_basis$new(orders = c(2,0,0), max_degrees = c(1,2,3))
    basis3 <- fourier_basis$new(orders = c(3,0,0), max_degrees = c(1,2,3))
    basis4 <- fourier_basis$new(orders = c(4,0,0), max_degrees = c(1,2,3))
    basis5 <- fourier_basis$new(orders = c(5,0,0), max_degrees = c(1,2,3))
    basis6 <- fourier_basis$new(orders = c(7,0,0), max_degrees = c(1,2,3))
    basis7 <- fourier_basis$new(orders = c(9,0,0), max_degrees = c(1,2,3))
    basis8 <- fourier_basis$new(orders = c(12,0,0), max_degrees = c(1,2,3))

    basis_list <- list("k=0" = NULL,"k=1"= basis1,
                       "k=2" = basis2,
                       "k=3" =basis3,
                       "k=4" =basis4,
                       "k=5" = basis5,
                       "k=7"= basis6,
                       "k=9" = basis7,
                       "k=12" = basis8)
  } else {
    basis1 <- fourier_basis$new(orders = c(1,0,0), max_degrees = c(1,2,3))
    basis2 <- fourier_basis$new(orders = c(1,1,0), max_degrees = c(1,2,3))
    basis2a <- fourier_basis$new(orders = c(2,0,0), max_degrees = c(1,2,3))
    basis3 <- fourier_basis$new(orders = c(3,1,0), max_degrees = c(1,2,3))
    basis4 <- fourier_basis$new(orders = c(4,0,0), max_degrees = c(1,2,3))
    basis5 <- fourier_basis$new(orders = c(4,1,0), max_degrees = c(1,2,3))
    basis6 <- fourier_basis$new(orders = c(5,1,0), max_degrees = c(1,2,3))
    basis7 <- fourier_basis$new(orders = c(6,2,0), max_degrees = c(1,2,3))

    basis_list <- list("k=0" = NULL,"k=1,0"= basis1,
                       "k=2,0" = basis2,
                       basis2a,
                       "k=1,1" =basis3,
                       "k=2,1" =basis4,
                       "k=3,1" =basis5,
                       "k=4,1" =basis6,
                       "k=3,2" =basis7)
    if(ncol(V) > 2) {
      basis_list$`k=3,2,1` <- fourier_basis$new(orders = c(5,2,1), max_degrees = c(1,2,3))
    }
  }
  basis_list <- basis_list[1:min(length(basis_list), size)]
  return(basis_list)

}


generate_plugin_learners <- function(task_RR, basis_list, learner_list, task = NULL, likelihood = NULL, screen.glmnet = F, select_sieve = F, cv = F, include_subst = F) {
  plugin_chainers <- list()
  if(ncol(task_RR$X) <=1) {
    screen.glmnet <- FALSE
  }
  for(name in names(basis_list)) {
    basis <- basis_list[[name]]
    lrr_plugin <- LRR_plugin_task_generator$new(sieve_basis = basis, name = name)
    if(screen.glmnet) {
      glmnet_screener <- Lrnr_LRR_glmnet.screener$new()
      lrr_plugin <- make_learner(Pipeline, lrr_plugin, glmnet_screener)

    }
    lrr_plugin <- lrr_plugin$train(task_RR)
    plugin_chainers[[name]] <- lrr_plugin
  }
  lrr_learners_all <- list()
  for(lrnr in learner_list) {
    lrr_learners <- list()
    for(chainer in plugin_chainers) {
      lrr_learners <- c(lrr_learners, list(Pipeline$new(chainer, lrnr)))
    }

    if(select_sieve) {
      lrr_learners <- make_learner(Stack, lrr_learners)
      eff_loss <- make_eff_loss(task, likelihood)
      lrr_learners <- make_learner(Pipeline, lrr_learners, Lrnr_selector$new(eff_loss))
      lrr_learners <- list(lrr_learners)
    }
    lrr_learners_all <- c(lrr_learners_all, lrr_learners)
  }
  if(include_subst) {
    lrr_learners_all <- c(lrr_learners_all, Lrnr_LRR_subst$new())
  }
  if(length(lrr_learners_all) > 1) {
    lrr_learners_all <- make_learner(Stack, lrr_learners_all)
  } else {
    lrr_learners_all <- lrr_learners_all[[1]]

  }
  if(cv) {
    eff_loss <- make_eff_loss(task, likelihood)
    lrr_learners_all <- Lrnr_sl$new(lrr_learners_all, metalearner = Lrnr_cv_selector$new(eff_loss))
    #lrr_learners_all <- make_learner(Pipeline, Lrnr_cv$new(lrr_learners_all), Lrnr_cv_selector$new(eff_loss))
  }
  lrr_learners_all
}

