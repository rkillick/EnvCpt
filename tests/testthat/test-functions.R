# testing functions, aim to get 100% test coverage on exported code
# cpt.reg should ideally be tested when moved to changepoint package.

context("function tests")

# plot.aic.envcpt
set.seed(98135)
x=rnorm(50)
out=envcpt(x)
test_that("plotAIC xlim main", expect_silent(plot.aic.envcpt(out,xlim=c(0,10),main="AIC Test")))


# envcpt
test_that("no messages",expect_message(envcpt(x,verbose=FALSE),NA))

x[24]=NA
test_that("NA data",expect_error(envcpt(x),"data has missing values, this function cannot handle missing values"))

x=LETTERS
test_that("non-numeric",expect_error(envcpt(x),"data must be a numeric vector"))




# aic.envcpt
test_that("AIC not list",expect_error(aic.envcpt(rnorm(100)),"object argument must be a list"))
test_that("AIC not matrix",expect_error(aic.envcpt(list(summary=rnorm(100))),"first element in the object list must be a matrix."))
test_that("AIC not numeric",expect_error(aic.envcpt(list(summary=matrix(LETTERS,nrow=2))),"First two rows in matrix in first element of object list must be numeric"))




