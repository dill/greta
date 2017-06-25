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
greta_kernel <- function (kernel_name,
                          gpflow_name,
                          parameters,
                          components = NULL,
                          arguments = list()) {

  # check GPflow is available and get the kernel method
  check_gpflowr()

  gpflow_method <- gpflow$kernels[[gpflow_name]]

  kernel_name <- paste(kernel_name, "kernel function")

  parameters <- lapply(parameters, as.greta_array)

  kernel <- list(name = kernel_name,
                 parameters = parameters,
                 gpflow_method = gpflow_method,
                 components = components,
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

#' @export
is.greta_kernel_function <- function (x)
  inherits(x, "greta_kernel_function")

# combine greta kernel function objects
combine_greta_kernel_function <- function(a, b, combine = c('additive', 'multiplicative')) {

  combine <- match.arg(combine)

  if (!is.greta_kernel_function(a) && !is.greta_kernel_function(b)) {
    stop ("can only combine a greta kernel function",
          "with another greta kernel function",
          call. = FALSE)
  }

  kernel_a <- environment(a)$kernel
  kernel_b <- environment(b)$kernel

  gpflow_name <- switch(combine,
                        additive = 'Add',
                        multiplicative = 'Prod')

  greta_kernel(kernel_name = combine,
               gpflow_name = gpflow_name,
               parameters = c(kernel_a$parameters, kernel_b$parameters),
               components = list(kernel_a, kernel_b))

  # need to get names of the components. Options:

  # 1. work out the names of the members in the top-level object at this point,
  # and store them in the greta_kernel object. Need to account for naming of
  # duplicates as e.g. "bias_1", and handle the different names between greta
  # and gpflowr

  # 2. write own composition class on the R side, then build each gpflow
  # function one at a time, replacing parameters before combining to the next
  # level.

  # 2 seems more feasible - just need more code in compile_gpflow_kernel. Would
  # need to expose all (greta array) component parameters via $parameters in the
  # combination functions, to have it built properly

  # the order of parameters is more obvious than the naming.

  # so the combination functions have a parameters member which concatenates
  # those of their children. It must also have two component functions (kernel_a
  # and kernel_b). In compile_greta_kernel, these should be checked for and
  # iterated to build up the full tree, copying over the tensors at each stage

  # write a simple recursive function which, if it is a combination, does the
  # children first.

  # how to ensure the order is correct? Need to return a counter up in the
  # recursive function, always start with kernel_a in a combination function,
  # and use the next tensors in order

}

#' @export
`+.greta_kernel_function` <- function (e1, e2)
  combine_greta_kernel_function(e1, e2, 'additive')

#' @export
`*.greta_kernel_function` <- function (e1, e2)
  combine_greta_kernel_function(e1, e2, 'multiplicative')

# overload addition and multiplication of greta kernels
# - grab their kernel objects and combine them; then return a kernel function with that information

# recursively iterate through nested greta kernels, creating corresponding
# gpflow kernels and replacing their parameters with tensors
recurse_kernel <- function (greta_kernel, tf_parameters, counter) {

  # if it's compound, recursively call this function on the components then
  # combine them
  if (!is.null(greta_kernel$components)) {

    a <- recurse_kernel(greta_kernel$components[[1]],
                        tf_parameters,
                        counter)

    b <- recurse_kernel(greta_kernel$components[[2]],
                        tf_parameters,
                        counter)

    gpflow_kernel <- greta_kernel$gpflow_method(list(a, b))

  } else {

    # get gpflow version of the basis kernel
    gpflow_kernel <- do.call(greta_kernel$gpflow_method,
                             greta_kernel$arguments)

    # find the relevant tensors
    n_param <- length(greta_kernel$parameters)
    previous <- counter$count
    counter$count <- counter$count + n_param
    idx <- previous + seq_len(n_param)
    tf_parameters <- tf_parameters[idx]

    # put tensors in the gpflow kernel object
    parameter_names <- names(greta_kernel$parameters)
    for (i in seq_along(tf_parameters)) {
      name <- parameter_names[i]
      gpflow_kernel[[name]] <- tf_parameters[[i]]
    }

  }

  gpflow_kernel

}

# function to create gpflow kernel from a greta kernel; called when compiling
# the tf graph
compile_gpflow_kernel <- function (greta_kernel, tf_parameters) {

  # take kernel object for each sub-kernel, assign the relevant parts of the
  # parameters
  # for now just do flat version

  counter <- new.env()
  counter$count <- 0
  gpflow_kernel <- recurse_kernel(greta_kernel, tf_parameters, counter)

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

  # # create a cholesky factor representation of this
  # chol_res <- chol(res)
  # res$node$representations$cholesky_factor = chol_res$node

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
#'   The \code{+} and \code{*} operators can be used to combine kernel functions
#'   from the existing basis kernel functions, as demonstrated in the example.
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
#' # create a bias kernel, with the variance as a variable
#' k2 <- gp$kernels$bias(variance = lognormal(0, 1))
#'
#' # combine two kernels and evaluate
#' K <- k1 + k2
#' K(x, x2)
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
