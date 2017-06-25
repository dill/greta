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
greta_kernel <- function (kernel_name, gpflow_name, parameters, arguments) {

  # check GPflow is available and get the kernel method
  check_gpflowr()

  gpflow_method <- gpflow$kernels[[gpflow_name]]

  kernel_name <- paste(kernel_name, "kernel function")

  parameters <- lapply(parameters, as.greta_array)

  kernel <- list(name = kernel_name,
                 parameters = parameters,
                 gpflow_method = gpflow_method,
                 arguments = arguments)

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

#' @export
print.greta_kernel_function <- function (x, ...)
  cat(environment(x)$kernel$name, "\n")

# overload addition and multiplication of greta kernels
# - grab their kernel objects and combine them; then return a kernel function with that information

# function to create gpflow kernel from a greta kernel; called when compiling
# the tf graph
compile_gpflow_kernel <- function (greta_kernel, tf_parameters) {

  # take kernel object for each sub-kernel, assign the relevant parts of the
  # parameters
  # for now just do flat version

  # get gpflow version
  gpflow_kernel <- do.call(greta_kernel$gpflow_method,
                           greta_kernel$arguments)

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

bias_kernel <- function (variance, dim = 1) {
  greta_kernel("bias",
               gpflow_name = "Bias",
               parameters = list(variance = variance),
               arguments = list(input_dim = as.integer(dim)))
}

white_kernel <- function (variance, dim = 1) {
  greta_kernel("white",
               gpflow_name = "White",
               parameters = list(variance = variance),
               arguments = list(input_dim = as.integer(dim)))
}

linear_kernel <- function (variances) {
  greta_kernel("linear",
               gpflow_name = 'Linear',
               parameters = list(variance = variances),
               arguments = list(input_dim = length(variances),
                                ARD = TRUE))
}

rbf_kernel <- function (lengthscales, variance) {
  greta_kernel("radial basis",
               gpflow_name = 'RBF',
               parameters = list(lengthscales = t(lengthscales),
                                 variance = variance),
               arguments = list(input_dim = length(lengthscales),
                                ARD = TRUE))
}

exponential_kernel <- function (lengthscales, variance) {
  greta_kernel("exponential",
               gpflow_name = 'Exponential',
               parameters = list(lengthscales = t(lengthscales),
                                 variance = variance),
               arguments = list(input_dim = length(lengthscales),
                                ARD = TRUE))
}

matern12_kernel <- function (lengthscales, variance) {
  greta_kernel("Matern 1/2",
               gpflow_name = 'Matern12',
               parameters = list(lengthscales = t(lengthscales),
                                 variance = variance),
               arguments = list(input_dim = length(lengthscales),
                                ARD = TRUE))
}

matern32_kernel <- function (lengthscales, variance) {
  greta_kernel("Matern 3/2",
               gpflow_name = 'Matern32',
               parameters = list(lengthscales = t(lengthscales),
                                 variance = variance),
               arguments = list(input_dim = length(lengthscales),
                                ARD = TRUE))
}

matern52_kernel <- function (lengthscales, variance) {
  greta_kernel("Matern 5/2",
               gpflow_name = 'Matern52',
               parameters = list(lengthscales = t(lengthscales),
                                 variance = variance),
               arguments = list(input_dim = length(lengthscales),
                                ARD = TRUE))
}

periodic_kernel <- function (period, lengthscale, variance, dim = 1) {
  greta_kernel("periodic",
               gpflow_name = 'PeriodicKernel',
               parameters = list(period = period,
                                 lengthscales = lengthscale,
                                 variance = variance),
               arguments = list(input_dim = as.integer(dim)))
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
#' @details The kernels are imported from the GPflow python package, using the
#'   gpflowr R package. Both of those need to be installed before you can use
#'   these methods. See the \href{gpflow.readthedocs.io}{GPflow website} for
#'   details of the kernels implemented.
#'
#' @param variance,variances (scalar/vector) the variance of a Gaussian process
#'   prior in all dimensions (\code{variance}) or in each dimensions
#'   (\code{variances}).
#' @param lengthscale,lengthscales (scalar/vector) the correlation decay
#'   distance along all dimensions (\code{lengthscale}) or each dimension
#'   ((\code{lengthscales})) of the Gaussian process
#' @param period (scalar) the period of the Gaussian process
#' @param dim (scalar integer, not a greta array) the dimension of the Gaussian
#'   process (number of columns on which it acts)
#'
#' @section Usage:
#' \preformatted{
#'   gp$kernels$bias(variance, dim = 1)
#'   gp$kernels$white(variance, dim = 1)
#'   gp$kernels$linear(variances)
#'   gp$kernels$rbf(lengthscales, variance)
#'   gp$kernels$exponential(lengthscales, variance)
#'   gp$kernels$matern12(lengthscales, variance)
#'   gp$kernels$matern32(lengthscales, variance)
#'   gp$kernels$matern52(lengthscales, variance)
#'   gp$kernels$periodic(period, lengthscale, variance)
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
gp <- list(kernels = list(bias = bias_kernel,
                          white = white_kernel,
                          rbf = rbf_kernel,
                          exponential = exponential_kernel,
                          linear = linear_kernel,
                          matern12 = matern12_kernel,
                          matern32 = matern32_kernel,
                          matern52 = matern52_kernel,
                          periodic = periodic_kernel))
