context('greta_array class')

test_that('print and summary work', {

  source('helpers.R')

  ga_data <- as_data(matrix(1:9, nrow = 3))
  ga_stochastic <- normal(0, 1)
  ga_operation <- ga_data * ga_stochastic

  # data arrays
  # print method
  expected_output <- "greta array (data)\n\n     [,1] [,2] [,3]\n[1,]    1    4    7\n[2,]    2    5    8\n[3,]    3    6    9"
  result <- evaluate_promise(ga_data, print = TRUE)
  expect_identical(result$output, expected_output)

  # summary method
  expected_output <- "'data' greta array with 9 elements (3 x 3)  \n\n       V1            V2            V3     \n Min.   :1.0   Min.   :4.0   Min.   :7.0  \n 1st Qu.:1.5   1st Qu.:4.5   1st Qu.:7.5  \n Median :2.0   Median :5.0   Median :8.0  \n Mean   :2.0   Mean   :5.0   Mean   :8.0  \n 3rd Qu.:2.5   3rd Qu.:5.5   3rd Qu.:8.5  \n Max.   :3.0   Max.   :6.0   Max.   :9.0  "
  result <- evaluate_promise(summary(ga_data), print = TRUE)
  expect_identical(result$output, expected_output)

  # stochastic arrays
  # print method
  expected_output <- "greta array (stochastic)\n\n     [,1]\n[1,]   ? "
  result <- evaluate_promise(ga_stochastic, print = TRUE)
  expect_identical(result$output, expected_output)

  # summary method
  expected_output <- "'stochastic' greta array with 1 element following a normal distribution \n\n  (values currently unknown)"
  result <- evaluate_promise(summary(ga_stochastic), print = TRUE)
  expect_identical(result$output, expected_output)

  # operation arrays
  # print method
  expected_output <- "greta array (operation)\n\n     [,1] [,2] [,3]\n[1,]   ?    ?    ? \n[2,]   ?    ?    ? \n[3,]   ?    ?    ? "
  result <- evaluate_promise(ga_operation, print = TRUE)
  expect_identical(result$output, expected_output)

  # summary method
  expected_output <- "'operation' greta array with 9 elements (3 x 3)  \n\n  (values currently unknown)"
  result <- evaluate_promise(summary(ga_operation), print = TRUE)
  expect_identical(result$output, expected_output)

})


test_that('length and dim work', {

  source('helpers.R')

  ga_data <- as_data(matrix(1:9, nrow = 3))
  ga_stochastic <- normal(0, 1, dim = c(3, 3))
  ga_operation <- ga_data * ga_stochastic

  # length
  expect_identical(length(ga_data), 9L)
  expect_identical(length(ga_stochastic), 9L)
  expect_identical(length(ga_operation), 9L)

  # dim
  expect_identical(dim(ga_data), c(3L, 3L))
  expect_identical(dim(ga_stochastic), c(3L, 3L))
  expect_identical(dim(ga_operation), c(3L, 3L))

})

test_that('head and tail work', {

  source('helpers.R')

  a <- randn(10, 1)
  b <- randn(10, 4)
  c <- randn(10, 3, 3)

  check_op(head, a)
  check_op(tail, a)

  check_op(head, b)
  check_op(tail, b)

  check_op(head, c)
  check_op(tail, c)

})


