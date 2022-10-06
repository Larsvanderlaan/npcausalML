#' Earth: Multivariate Adaptive Regression Splines
#'
#' This learner provides fitting procedures for building regression models thru
#' the spline regression techniques described in

Lrnr_earth <- R6Class(
  classname = "Lrnr_earth", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(degree = 2, penalty = 3, pmethod = "backward",
                          nfold = 0, ncross = 1, minspan = 0, endspan = 0,
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),
    .train = function(task) {
      args <- self$params
      outcome_type <- self$get_outcome_type(task)
      args$x <- task$X
      args$y <- outcome_type$format(task$Y)

      if (task$has_node("weights")) {
        args$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset
      }
      if(!is.null(args$family)) {
        glm <- list(family = args$family)
      } else if (outcome_type$type == "continuous") {
        glm <- list(family = stats::gaussian)
      } else if (outcome_type$type == "binomial") {
        glm <- list(family = stats::binomial)
      } else {
        stop("Unsupported outcome type for Lrnr_earth.")
      }
      args$glm <- glm
      earth_fun <- utils::getS3method("earth", "default",
                                      envir = getNamespace("earth")
      )

      # incorporate arguments defined by non-default earth function class
      default_args <- names(formals(earth_fun))
      extra_args <- names(formals(
        utils::getS3method("earth", "fit", envir = getNamespace("earth"))
      ))
      extra_args <- extra_args[!(extra_args %in% default_args)]

      fit_object <- sl3:::call_with_args(earth_fun, args, other_valid = extra_args)
      return(fit_object)
    },
    .predict = function(task) {
      preds <- stats::predict(
        object = private$.fit_object, newdata = task$X,
        type = "response"
      )
      return(preds)
    },
    .required_packages = c("earth")
  )
)
