Lrnr_LRR_subst <- R6Class(
  classname = "Lrnr_LRR_subst", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c( "LRR"),

    .train = function(task) {
      fit_object <- list()
      return(fit_object)
    },
    .predict = function(task = NULL) {
      fit_obj <- self$fit_object
      if(all(task$get_data(,"Q1V")[[1]] == -1) ){
        ER1 <- bound(task$get_data(,"Q1V"), c(0.001))
        ER0 <- bound(task$get_data(,"Q0V"), c(0.001))
      } else {
        ER1 <- bound(task$get_data(,"Q1"), c(0.001))
        ER0 <- bound(task$get_data(,"Q0"), c(0.001))
      }

      return(log(ER1/ER0))
    }
  )
)
