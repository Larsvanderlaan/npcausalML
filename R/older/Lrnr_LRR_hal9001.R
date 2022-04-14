
#' A learner that can be used to estimate the relative risk function. Valid for \code{library_RR}.
#' @param max_degree The maximum degree of interactions to generat
#' @param smoothness_orders The smoothness order of the fit.
#' @param num_knots The number of knots to generate.
#'
#' See the R package \code{\link[hal9001]} for a more detailed description of these parameters.
#' @export
Lrnr_LRR_hal9001 <- R6Class(
  classname = "Lrnr_LRR_hal9001", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 2,
                          smoothness_orders =1,
                          num_knots = c(100, 50),

                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("LRR"),

    .train = function(task) {
      params <- self$params
      method <- self$params$method
      X <- task$X
      family <-  binomial()
      weights <- task$weights
      Y <- task$Y
      basis_list <- hal9001::enumerate_basis(as.matrix(task$X),  max_degree = self$params$max_degree, smoothness_orders = self$params$smoothness_orders, num_knots = self$params$num_knots)
      x_basis <- as.matrix(hal9001::make_design_matrix(as.matrix(task$X),basis_list))
      print("HAL dim")
      print(dim(x_basis))
      fit <- NULL
      fit <- glmnet::cv.glmnet( x_basis, Y, standardize = F, family = binomial(), weights = weights, intercept = F)
      fit <- list(coef = coef(fit, s = "lambda.min")[-1])
      return(fit_object = list(fit = fit, basis_list = basis_list))
    },
    .predict = function(task = NULL) {
      fit_obj <- self$fit_object
      fit <- fit_obj$fit
      #preds <- predict(fit, new_data = as.matrix(task$X), type = "link")
       coef <- as.vector(fit$coef)
       coef[is.na(coef)] <- 0
       X <- as.matrix(task$X)
       x_basis <- as.matrix(hal9001::make_design_matrix(as.matrix(task$X),fit_obj$basis_list))

       preds <- x_basis %*% coef

      return(preds)
    },

    .required_packages = c("hal9001")
  )
)

