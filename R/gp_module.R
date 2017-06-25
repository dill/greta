# Gaussian process module

# check availability of GPflow and gpflowr, error nicely if they aren't
# installed, and put the module in the calling environment if they are
check_gpflowr <- function () {

  gpflowr_available <- requireNamespace('gpflowr', quietly = TRUE)

  if (gpflowr_available)
    gpflow_available <- gpflowr::gpflow_available()

  if (!gpflowr_available || !gpflow_available) {

    stop ('the GPflow python package and the gpflowr R package must be installed to plot greta models',
          call. = FALSE)

  } else {

    assign('gpflow',
           gpflowr::gpflow,
           envir = parent.frame())

  }

}

# create a greta kernel function (to create ops)
greta_kernel <- function (kernel_name, parameters, gpflow_method, input_dim) {

  # check GPflow is available and get the kernel method
  check_gpflowr()

  gpflow_method <- gpflow$kernels[[gpflow_method]]

  kernel_name <- paste(kernel_name, "kernel function")

  parameters <- lapply(parameters, as.greta_array)

  kernel <- list(name = kernel_name,
                 parameters = parameters,
                 gpflow_method = gpflow_method,
                 input_dim = input_dim)

  # check and get the dimension of a target matrix
  get_dim <- function (x, name = 'X') {

    x_dim <- dim(x)

    if (length(x_dim) != 2) {
      stop (name, "must be a 2D greta array",
            call. = FALSE)
    }

    x_dim

  }

  # return a function here, acting on either one or two datasets
  kernel_function <- function (X, X_prime = NULL) {

    X <- as.greta_array(X)

    if (is.null(X_prime)) {

      op_data_list <- list(operation = 'self-covariance matrix',
                           X = X)
      tf_op <- 'tf_self_K'

      dimfun <- function (elem_list) {
        X_dim <- get_dim(elem_list[[1]], 'X')
        rep(X_dim[1], 2)
      }

    } else {

      X_prime <- as.greta_array(X_prime)

      op_data_list <- list(operation = 'covariance matrix',
                           X = X,
                           X_prime = X_prime)
      tf_op <- 'tf_K'

      dimfun <- function (elem_list) {

        X_dim <- get_dim(elem_list[[1]], 'X')
        X_prime_dim <- get_dim(elem_list[[2]], 'X_prime')

        if (X_dim[2] != X_prime_dim[2]) {
          stop ('number of columns of X and X_prime do not match',
                call. = FALSE)
        }

        c(X_prime_dim[1], X_dim[1])
      }

    }

    # kernel parameters (as greta arrays) are getting fetched here anyway so
    # just need method to fetch/assign parameters across more complex kernels

    args <- c(op_data_list,
              kernel$parameters,
              list(dimfun = dimfun,
                   operation_args = list(greta_kernel = kernel),
                   tf_operation = tf_op))

    do.call(op, args)

  }

  # give it a class and return
  class(kernel_function) <- c('greta_kernel_function', class(kernel_function))
  kernel_function

}

print.greta_kernel_function <- function (x, ...)
  cat(kern <- environment(x)$kernel$name, "\n")

# overload addition and multiplication of greta kernels
# - grab their kernel objects and combine them; then return a kernel function with that information

# function to create gpflow kernel from a greta kernel; called when compiling
# the tf graph
compile_gpflow_kernel <- function (greta_kernel, tf_parameters) {

  # take kernel object for each sub-kernel, assign the relevant parts of the
  # parameters
  # for now just do flat version

  # get gpflow version
  dim <- as.integer(greta_kernel$input_dim)
  gpflow_kernel <- greta_kernel$gpflow_method(dim)

  # put tensor in the gpflow kernel object
  parameter_names <- names(greta_kernel$parameters)
  for (i in seq_along(tf_parameters)) {
    name <- parameter_names[i]
    gpflow_kernel[[name]] <- tf_parameters[[i]]
  }

  gpflow_kernel

}

# create gpflow kernel and evaluate with tensors
tf_K <- function (X, X_prime, ..., greta_kernel) {

  # evaluate with tensors
  tf_parameters <- list(...)
  gpflow_kernel <- compile_gpflow_kernel(greta_kernel, tf_parameters)
  gpflow_kernel$K(X, X_prime)

}

tf_self_K <- function (X, ..., greta_kernel) {

  # evaluate with tensors
  tf_parameters <- list(...)
  gpflow_kernel <- compile_gpflow_kernel(greta_kernel, tf_parameters)
  res <- gpflow_kernel$K(X)

  # create a cholesky factor representation of this
  chol_res <- chol(res)
  res$node$representations$cholesky_factor = chol_res$node

  res

}

# make kernel constructors

# return an R function to evaluate the kernel with these parameters
bias <- function (variance, dim = 1) {
  greta_kernel("bias",
               parameters = list(variance = variance),
               gpflow_method = "Bias",
               input_dim = dim)
}

# return an R function to evaluate the kernel with these parameters
rbf <- function (lengthscales, variance) {
  greta_kernel("radial basis",
               parameters = list(lengthscales = t(lengthscales),
                                 variance = variance),
               gpflow_method = 'RBF',
               input_dim = length(lengthscales))
}

#' @name gp-module
#' @aliases gp
#' @title methods for Gaussian process modelling
#'
#' @description A module providing composable kernel functions for use in
#'   Gaussian process models. Currently the only provided methods are the kernel
#'   constructor functions under gp$kernels. Each of these returns a
#'   \emph{function} which can be executed on greta arrays to compute the
#'   covariance matrix between points in the space of the Gaussian process. See
#'   the example for a demonstration.
#'
#' @param variance (scalar) the marginal variance of a Gaussian process prior
#'   under this kernel
#' @param lengthscales (column vector) the decay distance along each dimensions
#'   of the Gaussian process
#' @param dim the dimension of the Gaussian process (number of columns on which
#'   it acts)
#'
#' @section Usage:
#' \preformatted{
#'   gp$kernels$bias(variance, dim = 1)
#'   gp$kernels$rbf(lengthscales, variance)
#' }
#'
#' @examples
#' # create a radial basis function kernel on two dimensions
#' k1 <- gp$kernels$rbf(lengthscales = c(0.1, 0.2), variance = 0.6)
#'
#' # evaluate it on a greta array to get the variance-covariance matrix
#' x <- greta_array(rnorm(8), dim = c(4, 2))
#' k1(x)
#'
#' # non-symmetric covariance between two sets of points
#' x2 <- greta_array(rnorm(10), dim = c(5, 2))
#' k1(x, x2)
#'
#'
NULL

#' @export
gp <- list(kernels = list(bias = bias,
                          rbf = rbf))
