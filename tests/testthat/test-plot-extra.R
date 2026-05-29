# Additional coverage tests for plot.envcpt argument handling
library(testthat)
library(EnvCpt)

test_that("plot.envcpt errors when x is not envcpt", {
  expect_error(
    EnvCpt:::plot.envcpt(list()),
    "x must be an object with class envcpt"
  )
})

# These plotting branches are run locally / outside CRAN
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  set.seed(99)
  x <- c(rnorm(40, 0, 1), rnorm(40, 3, 1))
  out <- envcpt(x)
# The function accepts the British spelling "colours", although it currently warns
  test_that("plot.envcpt accepts colours alias but warns", {
    expect_warning(plot(out, type = "aic", colours = rep("black", 12)))
  })
  
  test_that("plot.envcpt errors when colors vector is too short", {
    expect_error(
      plot(out, type = "aic", colors = rep("black", 3)),
      "colors must be a vector of length 12"
    )
  })
  
  test_that("plot.envcpt errors when a color is invalid", {
    bad_cols <- c(rep("black", 11), "not_a_colour")
    expect_error(
      plot(out, type = "aic", colors = bad_cols),
      "Atleast one of your colours is not resolvable by col2rgb."
    )
  })
  
  test_that("plot.envcpt errors on unsupported plot type", {
    expect_error(
      plot(out, type = "weird"),
      "type supplied can only be 'aic', 'bic' or 'fit'"
    )
  })
}
