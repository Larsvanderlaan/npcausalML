Lrnr_LRR_glmnet.screener <- R6Class(
  classname = "Lrnr_LRR_glmnet.screener", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(
      ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("LRR"),

    .train = function(task) {

      keep <- task$weights!=0
      task <- task[task$weights!=0]
      method <- self$params$method
      X <- task$X
      Y <- (task$Y)
      weights <- task$weights
      Y <- as.vector(Y)
      X <- as.matrix(X)

      family <-   binomial()
      fit_object <- cv.glmnet(as.matrix(X), Y, family =  family, weights = weights, intercept = T)
      coefs <- coef(fit_object, s = "lambda.min")


      selected <- task$nodes$covariates
      selected <- selected[which(abs(coefs[-1]) > 1e-9)]
      fit_object <- list(coef = as.vector(coefs), selected = selected)
      return(fit_object)
    },
    .predict = function(task = NULL) {

      task <- task$next_in_chain(covariates = self$fit_object$selected)


      return(task$X)
    },
    .chain = function(task = NULL) {
      return(task$next_in_chain(covariates = self$fit_object$selected))
    },
    .required_packages = c( "glmnet")
  )
)
