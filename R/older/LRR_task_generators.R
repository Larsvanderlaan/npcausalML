


Lrnr_weight_helper <- R6Class(
  classname = "Lrnr_weight_helper", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function( lrnr, name="", ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
      private$.name <- paste0(self$params$name, "_", self$name)
    },

    print = function() {
      print(self$name)
    }
  ),

  private = list(
    .properties = c("continuous", "binomial", "categorical", "weights", "offset"),

    .train = function(task) {
      lrnr <- self$params$lrnr
      task <- task[abs(task$weights)>1e-4]
      lrnr <- lrnr$train(task)
      fit_object <- list(lrnr = lrnr)

      return(fit_object)
    },

    .predict = function(task = NULL) {
      lrnr <- self$fit_object$lrnr

      return(lrnr$predict(task))

    }
  )
)




LRR_plugin_task_generator <- R6Class(
  classname = "LRR_plugin_task_generator", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(sieve_basis = NULL,  name = "plugin", sequential = F, ...) {
      if(is.null(sieve_basis)) {
        sieve_learner <- NULL
        sieve_learner_Q1V <- NULL
        sieve_learner_Q0V <- NULL
      } else {
        sieve_learner <- Lrnr_adaptive_sieve$new(basis_generator = sieve_basis, stratify_by = "A", mult_by = "gGinv")
        sieve_learner_Q1V <- Lrnr_adaptive_sieve$new(basis_generator = sieve_basis)
        sieve_learner_Q0V <- Lrnr_adaptive_sieve$new(basis_generator = sieve_basis)

      }
      params <- list(sequential=sequential,name = name, sieve_basis = sieve_basis, Q_var = "Q", Q1_var = "Q1", Q0_var = "Q0", treatment_var = "A", outcome_var = "Y", g1_var = "g1", ginv_var = "gGinv")
      params$sieve_learner <- sieve_learner
      params$sieve_learner_Q1V <- sieve_learner_Q1V
      params$sieve_learner_Q0V <- sieve_learner_Q0V
      super$initialize(params = params, ...)
      private$.name <- paste0(self$params$name, "_", self$name)

    },

    print = function() {
      print( self$name)
    }
  ),

  private = list(
    .properties = c("continuous", "binomial", "categorical", "weights", "offset"),

    .train = function(task) {

      trt <- self$params$treatment_var
      outcome <- self$params$outcome_var
      Q_var <- self$params$Q_var
      Q0_var <- self$params$Q0_var
      Q1_var <- self$params$Q1_var
      lrnr <- self$params$sieve_learner
      if(is.null(lrnr)) {
        return(list())
      }
      task_q <- task$next_in_chain(offset = Q_var, covariates = union(task$nodes$covariates, trt), outcome = outcome)
      task_q_sub <- task_q[task$get_data(,"Delta")[[1]]==1]
      lrnr <- lrnr$train(task_q_sub)
      fit_object <- list()
      fit_object$lrnr <- lrnr
      if(!self$params$sequential || all(task$get_data(, "Q1V")[[1]]==-1)){
        return(fit_object)
      }
      g1 <- task$get_data(,self$params$g1_var)[[1]]
      g0 <- 1-g1
      G <- task$get_data(,"G")[[1]]
      cf_data1 <- data.table(rep(1, task_q$nrow))
      names(cf_data1) <- trt
      cf_data1[["gGinv"]] <-  1/g1/G
      cf_data0 <- data.table(rep(0, task_q$nrow))
      names(cf_data0) <- trt
      cf_data0[["gGinv"]] <-  1/g0/G
      column_map1 <- task_q$add_columns(cf_data1)
      column_map0 <- task_q$add_columns(cf_data0)
      task_q1 <- task_q$next_in_chain(column_names = column_map1, offset = Q1_var)
      task_q0 <- task_q$next_in_chain(column_names = column_map0, offset = Q0_var)
      Qnew0 <- lrnr$predict(task_q0)
      Qnew1 <- lrnr$predict(task_q1)

      column_map <- task_q$add_columns(data.table(Qnew1 = Qnew1, Qnew0=Qnew0))
      task_qV1 <- task_q$next_in_chain(column_names = column_map, outcome = "Qnew1", offset = "Q1V")
      task_qV0 <- task_q$next_in_chain(column_names = column_map, outcome = "Qnew0", offset = "Q0V")
      task_qV0 <- sl3_Task$new(task_qV0$internal_data, nodes = task_qV0$nodes, outcome_type = variable_type("binomial"), column_names = task_qV0$column_names, folds = task_qV0$folds, row_index = task_qV0$row_index)
      task_qV1 <- sl3_Task$new(task_qV1$internal_data, nodes = task_qV1$nodes, outcome_type = variable_type("binomial"), column_names = task_qV1$column_names, folds = task_qV1$folds, row_index = task_qV1$row_index)

      lrnr1 <- self$params$sieve_learner_Q1V
      lrnr0 <- self$params$sieve_learner_Q0V
      lrnr1 <- lrnr1$train(task_qV1)
      lrnr0 <- lrnr0$train(task_qV0)
      fit_object$lrnr1 <- lrnr1
      fit_object$lrnr0 <- lrnr0

      return(fit_object)
    },

    .predict = function(task = NULL) {
      stop("Nothing to predict")
    },
    .chain = function(task = NULL) {
      print("PLUGIN")
      trt <- self$params$treatment_var
      outcome <- self$params$outcome_var
      Q_var <- self$params$Q_var
      Q0_var <- self$params$Q0_var
      Q1_var <- self$params$Q1_var

      lrnr <- self$fit_object$lrnr

      if(!is.null(lrnr)) {
        g1 <- task$get_data(,self$params$g1_var)[[1]]
        g0 <- 1-g1
        G <- task$get_data(,"G")[[1]]
        task_q <- task$next_in_chain(offset = Q_var, covariates = union(task$nodes$covariates, trt), outcome = outcome)
        Qnew <- lrnr$predict(task_q)
        cf_data1 <- data.table(rep(1, task_q$nrow))
        names(cf_data1) <- trt
        cf_data1[[self$params$ginv_var]] <-  1/g1/G
        cf_data0 <- data.table(rep(0, task_q$nrow))
        names(cf_data0) <- trt
        cf_data0[[self$params$ginv_var]] <-  1/g0/G
        column_map1 <- task_q$add_columns(cf_data1)
        column_map0 <- task_q$add_columns(cf_data0)
        task_q1 <- task_q$next_in_chain(column_names = column_map1, offset = Q1_var)
        task_q0 <- task_q$next_in_chain(column_names = column_map0, offset = Q0_var)


        Qnew0 <- lrnr$predict(task_q0)
        Qnew1 <- lrnr$predict(task_q1)

      } else {

        Qnew <- task$get_data(,Q_var)
        Qnew0 <-  unlist(task$get_data(,Q0_var))
        Qnew1 <-  unlist(task$get_data(,Q1_var))
      }

      if(self$params$sequential && !is.null(self$fit_object$lrnr1)){

        column_map <- task_q$add_columns(data.table(Qnew1 = Qnew1, Qnew0=Qnew0))
        task_qV1 <- task_q$next_in_chain(column_names = column_map, outcome = "Qnew1", offset = "Q1V")
        task_qV0 <- task_q$next_in_chain(column_names = column_map, outcome = "Qnew0", offset = "Q0V")

        lrnr1 <- self$fit_object$lrnr1
        lrnr0 <- self$fit_object$lrnr0
        Q1Vnew <- lrnr1$predict(task_qV1)
        Q0Vnew <- lrnr0$predict(task_qV0)
        weights <- (Q0Vnew + Q1Vnew)*task$weights
        weights <- weights/mean(weights)
        outcome_val <- Q1Vnew / (Q0Vnew + Q1Vnew)
        new_data <- data.table(outcome_val, weights)
        colnames(new_data) <- c("outcome", "weights")
        new_data$Q0Vnew <- Q0Vnew
        new_data$Q1Vnew <- Q1Vnew


      } else {
        weights <- (Qnew0 + Qnew1)*task$weights
        weights <- weights/mean(weights)
        outcome_val <- Qnew1 / (Qnew0 + Qnew1)
        new_data <- data.table(outcome_val, weights)
        colnames(new_data) <- c("outcome", "weights")
        new_data$Qnew0 <- Qnew0
        new_data$Qnew1 <- Qnew1
        new_data$Qnew <- Qnew
      }



      covariates_plugin <- setdiff(task$covariates, trt)
      column_names <- task$add_columns(new_data)
      task <- task$next_in_chain(column_names = column_names, weights = "weights", outcome = "outcome", covariates = covariates_plugin)
      task <- sl3_Task$new(task$internal_data, nodes = task$nodes, outcome_type = variable_type("binomial"), column_names = task$column_names, folds = task$folds, row_index = task$row_index)

      return(task)
    }
  )
)

LRR_IPW_task_generator <- R6Class(
  classname = "LRR_IPW_task_generator", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(sieve_basis, name = "IPW", ...) {
      if(is.null(sieve_basis)) {
        sieve_learner <- NULL
      } else {
        sieve_learner <-  Lrnr_adaptive_sieve$new(basis_generator = sieve_basis, mult_by = c("Qg1", "Qg0"))
      }
      params <- list(name = name, g1_var = "g1", treatment_var = "A", outcome_var = "Y", sieve_learner = sieve_learner)

      super$initialize(params = params, ...)
      private$.name <- paste0(self$params$name, "_", self$name)
    },

    print = function() {
      print(self$name)
    }
  ),

  private = list(
    .properties = c("continuous", "binomial", "categorical", "weights", "offset"),

    .train = function(task) {
      trt <- self$params$treatment_var
      outcome <- self$params$outcome_var
      g1_var <- self$params$g1_var
      lrnr <- self$params$sieve_learner
      if(is.null(lrnr)) {
        return(list())
      }
      task_g <- task$next_in_chain(offset = g1_var, covariates = setdiff(task$nodes$covariates, trt), outcome = trt)
      lrnr <- lrnr$train(task_g)
      fit_object <- list()
      fit_object$lrnr <- lrnr
      return(fit_object)
    },

    .predict = function(task = NULL) {
      stop("Nothing to predict")
    },
    .chain = function(task = NULL) {
      print("IPW")
      trt <- self$params$treatment_var
      outcome <- self$params$outcome_var
      g1_var <- self$params$g1_var
      lrnr <- self$fit_object$lrnr
      if(!is.null(lrnr)) {
        task_g <- task$next_in_chain(offset = g1_var, covariates = setdiff(task$nodes$covariates, trt), outcome = trt)
        g1new <- lrnr$predict(task_g)
      } else {
        g1new <- task$get_data(,g1_var)[[1]]
      }
      #g1new <- tmle3::bound(g1new, 0.125)
      A <- task$get_data()[[trt]]
      gnew <- ifelse(A==1, g1new, 1 - g1new)
      covariates_IPW <- setdiff(task$covariates, trt)
      weights <-  task$get_data(,outcome)[[1]]/gnew * task$weights
      weights <- weights /mean(weights[weights!=0])
      new_data <- data.table(g1new, weights, gnew)
      colnames(new_data) <- c(g1_var, "weights", "gnew")
      column_names <- task$add_columns(new_data)

      task <- task$next_in_chain(column_names = column_names, weights = "weights", outcome = trt, covariates = covariates_IPW)
      task <- sl3_Task$new(task$internal_data, nodes = task$nodes, outcome_type = variable_type("binomial"), column_names = task$column_names, folds = task$folds, row_index = task$row_index)

      return(task)
    }
  )
)



#
# Lrnr_LRR_plugin_chainer <- R6Class(
#   classname = "Lrnr_LRR_plugin_chainer", inherit = Lrnr_base,
#   portable = TRUE, class = TRUE,
#   public = list(
#     initialize = function(Q_var = "Q", Q1_var = "ER1", Q0_var = "ER0", treatment_var = "A", outcome_var = "R", sieve_learner = Lrnr_fourier$new(fourier_basis(2,1)),  name = "plugin", ...) {
#       params <- args_to_list()
#       super$initialize(params = params, ...)
#       private$.name <- paste0(self$params$name, "_", self$name)
#
#     },
#
#     print = function() {
#       print( self$name)
#     }
#   ),
#
#   private = list(
#     .properties = c("continuous", "binomial", "categorical", "weights", "offset"),
#
#     .train = function(task) {
#       trt <- self$params$treatment_var
#       outcome <- self$params$outcome_var
#       Q_var <- self$params$Q_var
#       lrnr <- self$params$sieve_learner
#       if(is.null(lrnr)) {
#         return(list())
#       }
#       task_q <- task$next_in_chain(offset = Q_var, covariates = union(task$nodes$covariates, trt), outcome = outcome)
#       lrnr <- lrnr$train(task_q)
#       fit_object <- list()
#       fit_object$lrnr <- lrnr
#       return(fit_object)
#     },
#
#     .predict = function(task = NULL) {
#       stop("Nothing to predict")
#     },
#     .chain = function(task = NULL) {
#       trt <- self$params$treatment_var
#       outcome <- self$params$outcome_var
#       Q_var <- self$params$Q_var
#       Q0_var <- self$params$Q0_var
#       Q1_var <- self$params$Q1_var
#
#       lrnr <- self$fit_object$lrnr
#
#       if(!is.null(lrnr)) {
#         task_q <- task$next_in_chain(offset = Q_var, covariates = union(task$nodes$covariates, trt), outcome = outcome)
#         cf_data1 <- data.table(rep(1, task_q$nrow))
#         names(cf_data1) <- trt
#         cf_data0 <- data.table(rep(0, task_q$nrow))
#         names(cf_data0) <- trt
#         column_map1 <- task_q$add_columns(cf_data1)
#         column_map0 <- task_q$add_columns(cf_data0)
#         task_q1 <- task_q$next_in_chain(column_names = column_map1, offset = Q1_var)
#         task_q0 <- task_q$next_in_chain(column_names = column_map0, offset = Q0_var)
#
#         Qnew <- lrnr$predict(task_q)
#         Qnew0 <- lrnr$predict(task_q0)
#         Qnew1 <- lrnr$predict(task_q1)
#       } else {
#         Qnew <- task$get_data(,Q_var)
#         Qnew0 <-  task$get_data(,Q0_var)
#         Qnew1 <-  task$get_data(,Q1_var)
#       }
#
#
#       weights <- (Qnew0 + Qnew1)*task$weights
#       covariates_plugin <- setdiff(task$covariates, trt)
#       outcome_val <- Qnew1 / (Qnew0 + Qnew1)
#       new_data <- data.table(outcome_val, weights)
#       colnames(new_data) <- c("outcome", "weights")
#       column_names <- task$add_columns(new_data)
#       task <- task$next_in_chain(column_names = column_names, weights = "weights", outcome = "outcome", covariates = covariates_plugin)
#       task <- sl3_Task$new(task$internal_data, nodes = task$nodes, outcome_type = variable_type("binomial"), column_names = task$column_names, folds = task$folds, row_index = task$row_index)
#       return(task)
#     }
#   )
# )
#




Lrnr_chainer_link <- R6Class(
  classname = "Lrnr_chainer_link", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(...) {
      params <- args_to_list()
      super$initialize(params = params, ...)

    },

    print = function() {
      print( self$name)
    }
  ),

  private = list(
    .properties = c("continuous", "binomial", "categorical", "weights", "offset"),

    .train = function(task) {

      return(list())
    },

    .predict = function(task = NULL) {

     return(stats::qlogis(as.matrix(task$X)))
    }
)
)

















