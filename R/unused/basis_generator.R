hal_basis <- function(nbasis = 10, max_degree = 2, smoothness_orders = 0, include_zero_order = F, include_lower_order = F, ...) {
  max_degree <- max_degree
  bins <- nbasis
  smoothness_orders <- smoothness_orders
  include_zero_order <- include_zero_order
  include_lower_order <- include_lower_order
  generator <- function(X) {
    X <-as.matrix(X)
    smoothness_orders <- rep(smoothness_orders[[1]], ncol(X))
    basis_list <- hal9001fast::enumerate_basis(X, max_degree, bins = rep(bins[1], ncol(X)), order_map = smoothness_orders, include_zero_order = include_zero_order, include_lower_order = include_lower_order )
    gen <- function(X) {
      X <-as.matrix(X)
      x_basis <- hal9001fast::make_design_matrix(X, basis_list)
      return(x_basis)
    }
    return(gen)
  }
  return(generator)
}
class(hal_basis) <- "basis_config"

fourier_basis <- function(nbasis = 5, max_degree = 2, unpenalized = NULL, ...) {
  nbasis <- nbasis
  max_degree <- max_degree
  generator <- function(X) {
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
      basis <- fda::create.fourier.basis(rangeval = rangeval, nbasis = nbasis, dropind = 1)
      basis_list[[name]] <- basis

    }
    num_orig <- length(var_names)

    gen <- function(X) {
      k <- NULL
      data <- do.call(cbind, lapply(var_names, function(name) {
        out <- as.data.table(fda::eval.basis(X[,name], basis_list[[name]]))
        if(all(colnames(out) == c("sin", "cos"))) {
          colnames(out) <- c("sin1", "cos1")
        }
        k <<- ncol(out)
        return(out)
      }))



      var_key <- sapply(1:num_orig, function(k) {rep(paste0("X_",k), nbasis-1)})

      colnames(data) <- sapply(seq_along(colnames(data)), function(k) {
        paste0(var_key[k],"_", colnames(data)[k])
      })

      #data <- matrix(as.vector(as.matrix(data)), ncol = num_orig)
      #colnames(data) <- var_names
      if(max_degree == 1) {
        form <- formula(paste0("~. - 1" ))

      } else {
        form <- formula(paste0("~.^", max_degree," - 1" ))

      }
      data <- model.frame(form, data = as.data.frame(data))
      data <- model.matrix(form, data = data)
      data <- as.matrix(data)
      ns <- colnames(data)
      data <- cbind(rep(1, nrow(data)), data)
      colnames(data) <- c("inter", ns)
      var_names <- paste0("X_", 1:num_orig)
      index <- unlist(sapply(var_names, function(name) {
        grep(paste0(name,".+", name), colnames(data))
      }))
      index1 <- unlist(sapply(var_names, function(name) {
        grep(paste0(name,".+", name, ".+", name), colnames(data))
      }))
      index <- union(index, index1)
      if(length(index)>0){
        data <- data[, -index, drop = F]
      }
      if(!is.null(remove) & length(remove) >0){
        data <- data[,-remove]
      }
      #data$intercept <- 1
      #var_names <- colnames(data)
      #data <- do.call(cbind, unlist(apply(data, 2, function(v) {list(matrix(v, ncol = k))}), recursive = F))
      #data <- as.matrix(data)
      #new_names <- sapply(var_names, function(a) paste0(a, "_", 1:k))
      #colnames(data) <- new_names
      return(data)
    }
    data <- gen(as.matrix(X))
     remove <- which(abs(apply(data, 2, sd))<1e-6)
    keep <- setdiff(1:ncol(data), remove)
     remove2 <- findCorrelation(cor(data[,keep]), cutoff = 0.9)
    if(length(remove) >0) {
      remove <- remove[-1]
    }
    if(length(remove2)>0) {
      remove <- union( keep[remove2], remove)
    }
    return(gen)
  }
  return(generator)

}

class(fourier_basis) <- "basis_config"


bspline_basis <- function(nbasis = 5,  norder = 4, max_degree = 2, unpenalized = NULL, ...) {
  nbasis <- nbasis
  max_degree <- max_degree
  norder <- norder
  generator <- function(X) {
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

    num_orig <- length(var_names)
    remove <- NULL
    gen <- function(X) {
      k <- NULL
      data <- do.call(cbind, lapply(var_names, function(name) {
        out <- as.data.table(fda::eval.basis(X[,name], basis_list[[name]]))


        k <<- ncol(out)
        return(out)
      }))



      var_key <- sapply(1:num_orig, function(k) {rep(paste0("X_",k), nbasis)})

      colnames(data) <- sapply(seq_along(colnames(data)), function(k) {
        paste0(var_key[k],"_", colnames(data)[k])
      })

      #data <- matrix(as.vector(as.matrix(data)), ncol = num_orig)
      #colnames(data) <- var_names
      if(max_degree == 1) {
        form <- formula(paste0("~.", "-1"))

      } else {
        form <- formula(paste0("~.^", max_degree, "-1"))

      }

      data <- model.frame(form, data = as.data.frame(data))
      data <- model.matrix(form, data = data)

      data <- as.matrix(data)
      var_names <- paste0("X_", 1:num_orig)
      index <- unlist(sapply(var_names, function(name) {
        grep(paste0(name,".+", name), colnames(data))
      }))
      index1 <- unlist(sapply(var_names, function(name) {
        grep(paste0(name,".+", name, ".+", name), colnames(data))
      }))
      index <- union(index, index1)

      if(length(index)>0){
        data <- data[, -index, drop = F]
      }
      if(!is.null(remove) & length(remove) >0){
        data <- data[,-remove]
      }
      data$intercept <- 1
      #var_names <- colnames(data)
      #data <- do.call(cbind, unlist(apply(data, 2, function(v) {list(matrix(v, ncol = k))}), recursive = F))
      #data <- as.matrix(data)
      #new_names <- sapply(var_names, function(a) paste0(a, "_", 1:k))
      #colnames(data) <- new_names
      return(data)
    }

    data <- gen(as.matrix(X))
    remove <- which(abs(apply(data, 2, sd))<1e-6)
    keep <- setdiff(1:ncol(data), remove)
    remove2 <- findCorrelation(cor(data[,keep]), cutoff = 0.9)
    if(length(remove) >0) {
      remove <- remove[-1]
    }
    if(length(remove2)>0) {
      remove <- union( keep[remove2], remove)
    }
    return(gen)
  }
  return(generator)

}

# bspline_basis <- function(nbasis = 5, norder = 4, max_degree = 1, unpenalized = NULL, ...) {
#   nbasis <- nbasis
#   max_degree <- max_degree
#   norder <- norder
#   generator <- function(X) {
#     X <- as.matrix(X)
#     basis_list <- list()
#     var_names <- c()
#     print(colnames(X))
#     for(name in colnames(X)) {
#       sd_val <- sd(X[,name])
#       if(sd_val < 1e-9) {
#         next
#       }
#       var_names <- c(var_names, name)
#       rangeval <- c(min(X[,name]) - sd_val/2, max(X[,name]) +sd_val/2)
#       basis <- fda::create.bspline.basis(rangeval = rangeval, nbasis = nbasis, norder = norder)
#       basis_list[[name]] <- basis
#
#     }
#     num_orig <- length(var_names)
#
#     gen <- function(X) {
#       k <- NULL
#       data <- do.call(cbind, lapply(var_names, function(name) {
#         out <- as.data.table(fda::eval.basis(X[,name], basis_list[[name]]))
#         k <<- ncol(out)
#         return(out)
#       }))
#       data <- matrix(as.vector(as.matrix(data)), ncol = num_orig)
#       colnames(data) <- var_names
#       if(max_degree == 1) {
#         form <- formula(paste0("~."))
#
#       } else {
#         form <- formula(paste0("~.^", max_degree))
#
#       }
#       data <- model.frame(form, data = as.data.frame(data))
#       data <- model.matrix(form, data = data)
#
#       data <- as.matrix(data[,-1, drop = F])
#       var_names <- colnames(data)
#       data <- do.call(cbind, unlist(apply(data, 2, function(v) {list(matrix(v, ncol = k))}), recursive = F))
#       data <- as.matrix(data)
#       new_names <- sapply(var_names, function(a) paste0(a, "_", 1:k))
#       colnames(data) <- new_names
#       return(data)
#     }
#
#     return(gen)
#   }
#   return(generator)
#
# }
class(bspline_basis) <- "basis_config"

#
#
# poly_basis <- function(order = 1, max_degree = 1, unpenalized = NULL, ...) {
#   max_degree <- max_degree
#   exponents <- 1:order
#   generator <- function(X) {
#     X <- as.matrix(X)
#     basis_list <- list()
#     var_names <- c()
#     for(name in colnames(X)) {
#       sd_val <- sd(X[,name])
#       if(sd_val < 1e-9) {
#         next
#       }
#       var_names <- c(var_names, name)
#       rangeval <- c(min(X[,name]) - sd_val/2, max(X[,name]) +sd_val/2)
#       basis <- fda::create.monomial.basis(rangeval = rangeval, exponents = exponents)
#       basis_list[[name]] <- basis
#
#     }
#     num_orig <- length(var_names)
#
#     gen <- function(X) {
#       k <- NULL
#       data <- do.call(cbind, lapply(var_names, function(name) {
#         out <- as.data.table(fda::eval.basis(X[,name], basis_list[[name]]))
#         k <<- ncol(out)
#         return(out)
#       }))
#       data <- matrix(as.vector(as.matrix(data)), ncol = num_orig)
#       colnames(data) <- var_names
#       if(max_degree == 1) {
#         form <- formula(paste0("~."))
#
#       } else {
#         form <- formula(paste0("~.^", max_degree))
#
#       }
#
#       data <- model.frame(form, data = as.data.frame(data))
#       data <- model.matrix(form, data = data)
#
#       data <- as.matrix(data[,-1, drop = F])
#       var_names <- colnames(data)
#       data <- do.call(cbind, unlist(apply(data, 2, function(v) {list(matrix(v, ncol = k))}), recursive = F))
#       data <- as.matrix(data)
#       new_names <- sapply(var_names, function(a) paste0(a, "_", 1:k))
#       colnames(data) <- new_names
#       return(data)
#     }
#
#     return(gen)
#   }
#   return(generator)
#
# }
# class(poly_basis) <- "basis_config"
#



merge_basis <- function(gen1, gen2) {
  gen1 <- gen1
  gen2 <- gen2
  generator <- function(X) {
    X <- as.matrix(X)
    gen1 <- gen1(X)
    gen2 <- gen2(X)
    remove <- NULL
    gen <- function(X) {
      x1 <- gen1(X)
      x2 <- gen2(X)

      x <- cbind(x1, x2)

      #x <- x[,!duplicated(colnames(x))]
      colnames(x) <- make.unique(colnames(x))
      if(!is.null(remove) & length(remove) >0) {
        x <- x[,-remove]
      }
      return(as.matrix(x))
    }
    data <- gen(as.matrix(X))
    remove <- c()
    remove <- which(abs(apply(data, 2, sd))<1e-6)
    keep <- setdiff(1:ncol(data), remove)
     remove2 <- findCorrelation(cor(data[,keep]), cutoff = 0.9)
    if(length(remove) >0) {
      remove <- remove[-1]
    }
    if(length(remove2)>0) {
      remove <- union( keep[remove2], remove)
    }
    return(gen)
  }
}


