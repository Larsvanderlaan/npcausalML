#' Object that computes bspline basis matrix
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @export
#'
#'
#' @return
#'
#' @format \code{\link{R6Class}} object.
#'
#' @family Basis
#'
#
fourier_basis <- R6Class(
  classname = "fourier_basis",
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function( nknots = c(3,1,0), max_degrees = c(1,2,3), smoothness = 3,  ...) {
      # learner is an already initialized learner
      params <- list(
        max_degrees = max_degrees, nknots = nknots, smoothness = smoothness,
        ...
      )
      private$.params <- params

    },
    train = function(X) {
      nbasis <- max(self$params$nknots)
      norder <- self$params$smoothness
      remove <- NULL
      X <- as.matrix(X)
      basis_list <- list()
      var_names <- c()

      for(name in colnames(X)) {
        sd_val <- sd(X[,name])
        if(sd_val < 1e-9) {
          next
        }
        var_names <- c(var_names, name)
        a <- min(X[,name])
        b <- max(X[,name])
        width <- abs(b-a)
        a <- a - width/8
        b <- b + width/8
        rangeval <- c(a, b)
        basis <- fda::create.bspline.basis(rangeval = rangeval, nbasis = nbasis, norder = norder)
        basis_list[[name]] <- basis
      }
      private$.basis_list <- basis_list
      data <- self$predict(X)
      remove <- which(abs(apply(data, 2, sd))<1e-6)
      private$.remove <- remove
      return(self)
    },
    predict = function(X) {
      remove <- private$.remove
      max_degrees <- self$params$max_degrees
      orders <- self$params$orders

      X <- as.matrix(X)
      basis_list <- private$.basis_list
      var_names <- names(basis_list)
      nbasis_per_var <- ncol(as.matrix(fda::eval.basis(X[,var_names[[1]]], basis_list[[var_names[[1]]]])))
      indices <- 1:nbasis_per_var
      basis_names <- lapply(var_names, function(var) { paste0(var, indices)})
      names(basis_names) <- var_names
      main_terms <- do.call(cbind, lapply(var_names, function(name) {
        basis_eval <- data.table::as.data.table(fda::eval.basis(X[,name], basis_list[[name]]))
        colnames(basis_eval) <- basis_names[[name]]
        return(basis_eval)
      }))
      formula_of_degree <- function(i) {
        mdegree <- max_degrees[i]
        order <- orders[i]
        if(order == 0){
          return(NULL)
        }
        terms_list <- c()
        for(deg in  mdegree) {
          terms <- combn(var_names, deg, function(...){
            basis_names_reduced <- lapply(basis_names[...], function( names ) {
              return(names[1:(2*order)])
            })
            fun <- function(x,y){
              vals <- sort(c(x,y))
              paste0(vals[1],"*",vals[2])
            }

            fun <- Vectorize(fun)
            as.vector(unlist(reduce(basis_names_reduced, outer,  FUN = fun)))
          })
          terms_list <- c(terms_list, as.vector(unlist(terms)))
        }
        terms_list

      }
      formula <- paste0("~ ", paste0(sort(unique(as.vector(unlist(sapply(1:length(orders), formula_of_degree))))), collapse = " + "), " -1")
      formula <- as.formula(formula)
      x_basis <- model.frame(formula, data = as.data.frame(main_terms))
      print(as.data.frame(x_basis))
      x_basis <- model.matrix(formula, data = x_basis)
      x_basis <- as.matrix(x_basis)
      if(length(remove) > 0) {
        x_basis <- x_basis[, -remove]
      }
      return(x_basis)

    }
  ),
  active = list(
    name = function() {
      "fourier_basis"
    },
    params = function(){
      private$.params
    }
  ),
  private = list(
    .required_packages = NULL,
    .basis_list = NULL,
    .remove = c(),
    .params = list()
  )
)
