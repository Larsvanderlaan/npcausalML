
#' @import data.table
#' @import tmle3

make_task <- function(V, X = V, A, Y, Delta = NULL, weights = NULL,...) {

  if(is.vector(V)) {
    V <- as.matrix(V)
    colnames(V) <- "V"
  }
  if(is.vector(X)) {
    X <- as.matrix(X)
    colnames(X) <- "X"
  }
  colnames(V) <- paste0("V", 1:ncol(V))
  colnames(X) <- paste0("X", 1:ncol(X))

  if(is.null(weights)) {
    weights <- rep(1, length(Y))
  }
  data <- data.table(V, X, A = A, Y = Y, Delta = Delta, weights = weights)

  data <- data[,!duplicated(names(data)), with = F]

  npsem <- list(V = tmle3::define_node("V", colnames(V), c()),
                X = tmle3::define_node("X", colnames(X), c()),
                A = tmle3::define_node("A", "A", c("X")),
                Y = tmle3::define_node("Y", "Y", c("A", "X")),
                RR = tmle3::define_node("RR", c(), c("V")))
  if(!is.null(data$Delta)) {
    censoring_node <- tmle3::define_node("Delta", "Delta", c("A", "X"))
    npsem$Y$censoring_node <- censoring_node
    npsem$Delta <- censoring_node
  }
  task <- tmle3_Task$new(data, npsem, weights = "weights",...)
  return(task)
}

# Make likelihood
#' @import sl3
make_likelihood <- function(task, lrnr_A, lrnr_Y, lrnr_V = NULL, lrnr_Delta = NULL, cv = TRUE) {
  if(task$npsem$Y$variable_type$type == "binomial") {
    loss_Y <- loss_loglik_binomial
  } else {
    loss_Y <- loss_squared_error
  }
  if(cv) {
    lrnr_A <- make_learner(Pipeline, Lrnr_cv$new(lrnr_A), Lrnr_cv_selector$new(loss_loglik_binomial))
    lrnr_Y <- make_learner(Pipeline, Lrnr_cv$new(lrnr_Y), Lrnr_cv_selector$new(loss_Y))
    if(!is.null(lrnr_V)) {
      lrnr_V <- make_learner(Pipeline, Lrnr_cv$new(lrnr_V), Lrnr_cv_selector$new(loss_squared_error))
    }
  }
  factor_list <- list(LF_fit$new("A", lrnr_A, type = "density"),
                      LF_fit$new("Y", lrnr_Y, type = "mean", bound = c(0,1) ))
  likelihood <- Likelihood$new(factor_list)
  likelihood <- likelihood$train(task)

  if(!is.null(lrnr_V)) {
    lf1 <- LF_derived$new("Q1V",lrnr_V, likelihood, generator_Q1V, type = "mean", bound = c(0,1)  )
    lf0 <- LF_derived$new("Q0V", lrnr_V$clone(), likelihood, generator_Q0V, type = "mean",bound = c(0,1)   )
    likelihood$add_factors(list(lf0, lf1))
  }
  if(!is.null(task$npsem$Delta)) {
    lf_delta <- LF_fit$new("Delta", lrnr_Delta, type = "mean")
    likelihood$add_factors(list(lf_delta))
  }
  return(likelihood)
}



#### Sequential regression task generators
generator_Q1V <- function(tmle_task, likelihood) {
  cf_task <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A= rep(1, tmle_task$nrow)))

  Q1 <- likelihood$get_likelihood(cf_task, "Y", fold_number = "validation")
  task <- tmle_task$get_regression_task("RR")
  column_names <- task$add_columns(data.table(Q1 = Q1))
  task <- task$next_in_chain(column_names = column_names, outcome = "Q1")
  task <- sl3_Task$new(task$internal_data, nodes = task$nodes, column_names = task$column_names, folds = task$folds, row_index = task$row_index)

  return(task)
}

generator_Q0V <- function(tmle_task, likelihood) {
  cf_task <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A= rep(0, tmle_task$nrow)))

  Q0 <- likelihood$get_likelihood(cf_task, "Y", fold_number = "validation")
  task <- tmle_task$get_regression_task("RR")
  column_names <- task$add_columns(data.table(Q0 = Q0))
  task <- task$next_in_chain(column_names = column_names, outcome = "Q0")
  task <- sl3_Task$new(task$internal_data, nodes = task$nodes,  column_names = task$column_names, folds = task$folds, row_index = task$row_index)

  return(task)
}


# Generates revere
make_revere <- function(task, likelihood ) {
  genf <- make_generator(likelihood )
  return(sl3_revere_Task$new(genf, task))
}




make_generator <- function(likelihood) {




  gen_task_general <- function(tmle_task, fold_number) {
    task <- tmle_task$get_regression_task("RR")
    X <- task$X
    Y <- tmle_task$get_tmle_node("Y")
    A <- tmle_task$get_tmle_node("A")
    g <- likelihood$get_likelihood(tmle_task, "A", fold_number)
    if(!is.null(tmle_task$npsem$Delta)) {
      G <- likelihood$get_likelihood(tmle_task, "Delta", fold_number)
      Delta <- tmle_task$get_tmle_node("Delta")
    } else {
      G <- 1
      Delta <- 1
    }
    g <- tmle3::bound(g, 0.005)
    Q <- likelihood$get_likelihood(tmle_task, "Y", fold_number)
    cf_task1 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A= rep(1, length(g))))
    cf_task0 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A=rep(0, length(g))))
    Q1 <- bound(likelihood$get_likelihood(cf_task1, "Y", fold_number), 0.001)
    Q0 <- bound(likelihood$get_likelihood(cf_task0, "Y", fold_number),0.001)
    if(!is.null(likelihood$factor_list$Q1V)) {
      Q1V <- likelihood$get_likelihood(tmle_task, "Q1V", fold_number)
      Q0V <- likelihood$get_likelihood(tmle_task, "Q0V", fold_number)
    } else {
      Q1V <- -1
      Q0V <- -1
    }
    # weightsIPW <- Y/g * task$weights
    # weightsplugin <- (Q1 + Q0)
    # YIPW <- A
    # Yplugin <- Q1/weightsplugin
    covariates <- colnames(X)
    outcome <- c("Y")
    weights <- c(tmle_task$nodes$weights)

    data <- cbind(data.table(g1 = ifelse(A==1, g, 1-g), Q = Q, g = g, ginv = 1/g,  gGinv = 1/g/ G,  Y = Y, A = A, Q1 = Q1, Q0 = Q0, Q0V=Q0V, Q1V = Q1V, G = G, Delta = Delta), X)
    data <- cbind(task$get_data(,weights), data)
    #data$RR <- data$Q1 / data$Q0
    data$Qg1 <- data$Q1 / data$g1
    data$Qg0 <- data$Q0 / (1-data$g1)
    new_task <- sl3_Task$new(data, covariates = covariates, outcome = outcome, weights = weights, folds = task$folds, outcome_type = variable_type("binomial"))
    return(new_task)
  }


  return(gen_task_general)


}



make_eff_loss <- function(tmle_task, likelihood) {
  tmle_task <- tmle_task
  likelihood <- likelihood
  Y <- tmle_task$get_tmle_node("Y")
  A <- tmle_task$get_tmle_node("A")
  if(!is.null(tmle_task$npsem$Delta)) {
    Delta <- tmle_task$get_tmle_node("Delta")
    G <-  likelihood$get_likelihood(tmle_task, "Delta", "validation")
  } else{
    Delta <- rep(1, length(Y))
    G <- rep(1, length(Y))
  }
  g <- likelihood$get_likelihood(tmle_task, "A", "validation")
  lik <- likelihood

  cf_task1 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A= rep(1, length(g))))
  cf_task0 <- tmle_task$generate_counterfactual_task(uuid::UUIDgenerate(), data.table::data.table(A=rep(0, length(g))))
  Q1 <- lik$get_likelihood(cf_task1, "Y", "validation")
  Q0 <- lik$get_likelihood(cf_task0, "Y", "validation")
  Q <- lik$get_likelihood(tmle_task, "Y", "validation")

  if(!is.null(likelihood$factor_list$Q0V)) {
    Q1V <-  lik$get_likelihood(tmle_task, "Q1V", "validation")
    Q0V <-  lik$get_likelihood(tmle_task, "Q0V", "validation")

  } else {
    Q1V <- Q1
    Q0V <- Q0
  }
  efficient_loss = function(preds, Y, row_index = NULL) {
    Y <- tmle_task$get_tmle_node("Y")
    LRR <- preds

    if(!is.null(row_index)) {
      Q0V <- Q0V[row_index]
      Q1V <- Q1V[row_index]
      Q0 <- Q0[row_index]
      Q1 <- Q1[row_index]
      g <- g[row_index]
      A <- A[row_index]
      Y <- Y[row_index]
      Delta <- Delta[row_index]
      G <- G[row_index]
    }
    le <- log(1 + exp(LRR))
    term1 <- (Q0V)*le + (le - LRR)*Q1V
    term2 <- (le - LRR)*(Q1 - Q1V) + le*(Q0 - Q0V)
    term3 <- (Delta/G)*(A/g)*(le - LRR)*(Y - Q1) + (Delta/G)*((1-A)/g)*( le)*(Y - Q0)


    loss <- term1 + term2 + term3


    return( (loss))
  }

}

