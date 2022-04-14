

Lrnr_adaptive_sieve <- R6Class(
  classname = "Lrnr_adaptive_sieve", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(basis_generator = fourier_basis$new(),   mult_by = NULL, stratify_by = NULL, covariates_to_add = NULL,  ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(

    .train = function(task) {

      train <- function(task) {
        print("Fitting sieve")
        basis_generator <- self$params$basis_generator
        basis_generator <- basis_generator$clone()
        outcome_type <- task$outcome_type$glm_family(return_object = T)
        linkinv_fun <- outcome_type$linkinv
        link_fun <- outcome_type$linkfun


        if(task$has_node("weights")) {
          weights <- task$weights
        } else {
          weights <- NULL
        }
        if(task$has_node("offset")) {
          offset <- task$offset_transformed(link_fun, for_prediction = T)
        } else {
          offset <- NULL
        }

        covariates <- c(task$nodes$covariates, self$params$covariates_to_add)
        X <- as.matrix(task$get_data(, covariates))

        basis_generator <- basis_generator$train(X)
        x_basis <- basis_generator$predict(X)


        if(!is.null(self$params$mult_by)) {
          list_of_x <- list()
          for(column in c(self$params$mult_by)) {
            w <- task$get_data(,column)[[1]]


            list_of_x[[column]] <- w * x_basis
            colnames(list_of_x[[column]]) <- paste0(colnames(x_basis), "_", column)
          }
          x_basis <- do.call(cbind, list_of_x)

        }
        Y <- task$Y


        print("size: ")
        print(dim(x_basis))

        fit <- glm.fit( as.matrix(x_basis), Y, family = outcome_type, weights = weights, offset = offset, intercept = F)

        #coefs <- fit$coef
        coefs <- coef(fit)
        coefs[is.na(coefs)] <- 0
        if(length(coefs)!= ncol(x_basis)) {
          coefs <- coefs[-1]
        }

        fit$coef <- coefs
        fit$linkinv_fun <- linkinv_fun
        fit$link_fun <- link_fun
        fit_object = list(has_offset = !is.null(offset), basis_generator = basis_generator, fit = fit )
      }

      stratify <- self$params$stratify_by

      if(!is.null(stratify)) {
        covariates <- setdiff(task$nodes$covariates, stratify)
        A <- as.factor(task$get_data(,stratify)[[1]])
        levels <- levels(A)
        fit_objects <- list()
        for(level in levels) {
          new_task <- task[which(A==level)]$next_in_chain(covariates = covariates)

          fit_objects[[level]] <- train(new_task)
        }
      } else {
        levels <- NULL
        fit_objects <- train(task)
      }

      return(list(fit_objects = fit_objects, levels = levels))
    },
    .predict = function(task) {

      fit_objects <- self$fit_object
      levels <- fit_objects$levels
      fit_objects <- fit_objects$fit_objects
      predict_once <- function(task, fit_object) {

        fit <- fit_object$fit
        has_offset <- fit_object$has_offset
        basis_generator <- fit_object$basis_generator
        has_offset <- fit_object$has_offset

        covariates <- c(task$nodes$covariates, self$params$covariates_to_add)
        X <- as.matrix(task$get_data(, covariates))


        x_basis <- basis_generator$predict(X)





        if(!is.null(self$params$mult_by)) {
          list_of_x <- list()
          for(column in c(self$params$mult_by)) {
            w <- task$get_data(,column)[[1]]

            list_of_x[[column]] <- w * x_basis
            colnames(list_of_x[[column]]) <- paste0(colnames(x_basis), "_", column)
          }
          x_basis <- do.call(cbind, list_of_x)

        }

        coefs <- fit$coef
        coefs[is.na(coefs)] <- 0




        link <- as.matrix(x_basis) %*% as.vector(coefs)

        if(has_offset) {
          link <- link + task$offset_transformed(fit$link_fun, for_prediction = T)
        }


        pred <- fit$linkinv_fun(link)


        return(as.vector(pred))
      }
      stratify <- self$params$stratify_by
      if(!is.null(stratify)) {
        covariates <- setdiff(task$nodes$covariates, stratify)
        A <- as.factor(task$get_data(,stratify)[[1]])

        predictions <- rep(0, task$nrow)
        for(level in levels) {

          fit_obj <- fit_objects[[level]]
          index <- which(A==level)
          if(length(index) == 0){
            next
          }
          new_task <- task[index]$next_in_chain(covariates = covariates)
          predictions[index] <- predict_once(new_task, fit_obj)
        }
      } else {
        predictions <- predict_once(task, fit_objects)
      }
      return(predictions)

    },
    .required_packages = c("glmnet", "speedglm", "fda")
  ))


