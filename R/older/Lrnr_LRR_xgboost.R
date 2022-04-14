
#' @export
Lrnr_LRR_xgboost <- R6Class(
  classname = "Lrnr_LRR_xgboost", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(nrounds = 20, nthread = 1, num_parallel_tree = 1, max_depth = 6, eta = 0.3, gamma = 0, random_forest = F,
                          method = NULL,...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("LRR"),

    .train = function(task) {

      method <- self$params$method
      X <- task$X_intercept
      Y <- task$Y
      weights <- task$weights
      #weights <- weights * task$weights
      params <- list(verbose = 0, nrounds = self$params$nrounds, nthread = self$params$nthread, label = Y, data = as.matrix(X), weight = weights, max_depth = self$params$max_depth,
                     eta = self$params$eta, num_parallel_tree = self$params$num_parallel_tree, gamma = self$params$gamma, objective = "reg:logistic" )
      if(self$params$random_forest) {
        params$subsample <- 0.6
        params$colsample_bynode <- 0.8
        params$nrounds = 1
        params$eta = 1
        params$num_parallel_tree <- max(self$params$num_parallel_tree, 500)

      }
      fit_object <- do.call(xgboost::xgboost, params)
      return(fit_object)
    },
    .predict = function(task = NULL) {
      X <- task$X_intercept
      fit_obj <- self$fit_object
      predictions <- as.vector(predict(fit_obj, as.matrix(X), outputmargin = T))


      return(predictions)
    },
    .required_packages = c("xgboost")
  )
)
