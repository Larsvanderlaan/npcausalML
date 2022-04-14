
#' An \code{sl3} learner object that screens variables using \code{glmnet}.
glmnet_screener <- Lrnr_pkg_SuperLearner_screener$new("screen.glmnet")

#' A vector of names of valid \code{\link[SuperLearner]} learners that can be used in \code{Library_RR} for estimating the relative risk.
#' @export
valid_learners <- c("SL.xgboost", "SL.earth", "SL.gam", "SL.bayesglm", "SL.glm", "SL.glm.interaction", "SL.glmnet")
#' Prints the names of valid \code{\link[SuperLearner]} learners for relative risk estimation.
#' @export
valid_SL_learners <- function() {
  print(valid_learners)
}
#' Converts a name of a \code{\link[SuperLearner]} learner and additional parameters to an \code{sl3} learner that can be used in that can be used in \code{libraryRR}.
#' @param name A character with the name of a \code{\link[SuperLearner]}learner.
#' @param ... Additional arguments to be passed to the learner.
#' @export
SL_learner_to_LRR_learner <- function(name, ...) {
  if(F && !(name %in% valid_learners)) {
    stop(paste0(name, " is not a valid SL learner. Please call valid_SL_learners() to see valid learners,"))
  }
  lrnr <- make_learner(Lrnr_pkg_SuperLearner, name, family = binomial(), outcome_type = variable_type("binomial"), ...)
  binomial_to_LRR_learner(lrnr)
}
#' Converts a binomial learner that predicts `expit` transformed predictions to a relative risk learner that can be used in \code{libraryRR}
#' @param learner A binomial \code{sl3} learner.
#' @export
binomial_to_LRR_learner <- function(learner) {
  make_learner(Pipeline,learner , Lrnr_chainer_link$new())
}



# Default superlearner library for the nuisance parameters associated with A and Y
#' @export
default_library_nuisance <- make_learner(Pipeline, glmnet_screener, Stack$new(list(
  Lrnr_xgboost$new(max_depth = 3),
  Lrnr_xgboost$new(max_depth = 4),
  Lrnr_xgboost$new(max_depth = 5),
  Lrnr_xgboost$new(max_depth = 6),
  #Lrnr_ranger$new(),
  Lrnr_glmnet$new(),
  Lrnr_gam$new(),
  Lrnr_earth$new()
)))

# Default superlearner library for the sequential regression nuisance parameter associated with V
#' @export
default_library_V <-  Stack$new(list(
  make_pkg_SL_learner("SL.xgboost" , max_depth = 3, family = binomial(), outcome_type = variable_type("binomial")),
  make_pkg_SL_learner("SL.xgboost" , max_depth = 4, family = binomial(), outcome_type = variable_type("binomial")),
  make_pkg_SL_learner("SL.xgboost" , max_depth = 5, family = binomial(), outcome_type = variable_type("binomial")),
  make_pkg_SL_learner("SL.xgboost" , max_depth = 6, family = binomial(), outcome_type = variable_type("binomial")),
  make_pkg_SL_learner("SL.gam" , deg.gam = 3, family = binomial(), outcome_type = variable_type("binomial")),
  make_pkg_SL_learner("SL.gam" , deg.gam = 5, family = binomial(), outcome_type = variable_type("binomial")),
  make_pkg_SL_learner("SL.earth", degree = 2, pmethod = "forward", family = binomial(), outcome_type = variable_type("binomial")),
  make_pkg_SL_learner("SL.glm" , family = binomial(), outcome_type = variable_type("binomial")),
  make_pkg_SL_learner("SL.glm.interaction" , family = binomial(), outcome_type = variable_type("binomial"))
))


# Default superlearner library for relative risk function.
#' @export
default_library_RR <- list(
  "xgboost_3" = Lrnr_LRR_xgboost$new(max_depth = 3),
  "xgboost_4" =  Lrnr_LRR_xgboost$new( max_depth = 4),
  "xgboost_5" =  Lrnr_LRR_xgboost$new( max_depth = 5),
  "xgboost_6" =  Lrnr_LRR_xgboost$new( max_depth = 6),
  #"randomForest_6" =  Lrnr_LRR_xgboost$new(random_forest = TRUE, max_depth = 6),
  #"randomForest_8" =  Lrnr_LRR_xgboost$new(random_forest = TRUE, max_depth = 8),
  "gam_3" = SL_learner_to_LRR_learner("SL.gam", deg.gam = 3),
  "gam_5" = SL_learner_to_LRR_learner("SL.gam", deg.gam = 5),
  "earth_1" <- SL_learner_to_LRR_learner("SL.earth", degree = 1, pmethod = "forward"),
  "earth_2" <- SL_learner_to_LRR_learner("SL.earth", degree = 2, pmethod = "forward"),
  #"bayesglm" = SL_learner_to_LRR_learner("SL.bayesglm"),
  "glm" = SL_learner_to_LRR_learner("SL.glm"),
  "glm.interaction" = SL_learner_to_LRR_learner("SL.glm.interaction"),
  #"hal9001" = Lrnr_LRR_hal9001$new(max_degree = 2, smoothness_orders = 1, num_knots = c(30, 5)),
  Lrnr_LRR_subst$new()
  #"hal9001_2" = Lrnr_LRR_hal9001$new(max_degree = 2, smoothness_orders = 1, num_knots = c(40, 10)),
)
