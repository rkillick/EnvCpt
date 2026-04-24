# New tests to increase code coverage for EnvCpt
# Author: Pratik Bangerwa
# Target: CptReg.R (0%), envcpt.R (0%), diagnostics.R (1.24%)

library(testthat)
library(EnvCpt)
library(changepoint)
set.seed(42)

# ── Shared test data ───────────────────────────────────────────────────────────
n <- 200
ar1_data <- as.numeric(c(
  arima.sim(list(ar = 0.8),  n = n/2, sd = 1),
  arima.sim(list(ar = -0.6), n = n/2, sd = 1)
))

ar2_data <- as.numeric(c(
  arima.sim(list(ar = c(0.5, 0.3)),  n = n/2, sd = 1),
  arima.sim(list(ar = c(-0.4, 0.2)), n = n/2, sd = 1)
))

trend_data <- c(0.1*(1:100) + rnorm(100, 0, 0.3),
                0.5*(1:100) + rnorm(100, 0, 0.3))

make_ar1_mat <- function(x) {
  n <- length(x)
  cbind(x[-1], rep(1, n-1), x[-n])
}

make_ar2_mat <- function(x) {
  n <- length(x)
  cbind(x[-c(1:2)], rep(1, n-2), x[2:(n-1)], x[1:(n-2)])
}

# ═══════════════════════════════════════════════════════════════════════════════
# 1. cpt.reg — INPUT VALIDATION
# ═══════════════════════════════════════════════════════════════════════════════
context("cpt.reg input validation")

test_that("cpt.reg errors on non-array data", {
  expect_error(EnvCpt:::cpt.reg(as.data.frame(make_ar1_mat(ar1_data))),
               "Argument 'data' must be a numerical matrix/array")
})

test_that("cpt.reg errors on character data", {
  expect_error(EnvCpt:::cpt.reg(matrix(rep("a", 100), nrow=50)),
               "Argument 'data' must be a numerical matrix/array")
})

test_that("cpt.reg errors on invalid penalty type", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, penalty = 123),
               "Argument 'penalty' is invalid")
})

test_that("cpt.reg errors on multiple penalty values", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, penalty = c("MBIC", "AIC")),
               "Argument 'penalty' is invalid")
})

test_that("cpt.reg errors on invalid method type", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, method = 123),
               "Argument 'method' is invalid")
})

test_that("cpt.reg errors on unrecognised method", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, method = "CROPS"),
               "Invalid method, must be AMOC or PELT")
})

test_that("cpt.reg errors on invalid dist type", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, dist = 123),
               "Argument 'dist' is invalid")
})

test_that("cpt.reg warns on unsupported dist and converts to Normal", {
  mat <- make_ar1_mat(ar1_data)
  expect_warning(EnvCpt:::cpt.reg(mat, dist = "Poisson"),
                 "is not supported. Converted to dist='Normal'")
})

test_that("cpt.reg errors on invalid class argument", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, class = "yes"),
               "Argument 'class' is invalid")
})

test_that("cpt.reg errors on multiple class values", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, class = c(TRUE, FALSE)),
               "Argument 'class' is invalid")
})

test_that("cpt.reg errors on invalid param.estimates", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, param.estimates = "yes"),
               "Argument 'param.estimates' is invalid")
})

test_that("cpt.reg errors on non-numeric minseglen", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, minseglen = "three"),
               "Argument 'minseglen' is invalid")
})

test_that("cpt.reg errors on negative minseglen", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, minseglen = -1),
               "must be positive integer")
})

test_that("cpt.reg errors on non-integer minseglen", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, minseglen = 1.5),
               "must be positive integer")
})

test_that("cpt.reg errors on invalid tol type", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, tol = "small"),
               "Argument 'tol' is invalid")
})

test_that("cpt.reg errors on multiple tol values", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, tol = c(1e-7, 1e-5)),
               "Argument 'tol' is invalid")
})

test_that("cpt.reg errors on negative tol", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, tol = -1e-7),
               "must be positive")
})

test_that("cpt.reg warns when minseglen too small", {
  mat <- make_ar1_mat(ar1_data)
  expect_warning(EnvCpt:::cpt.reg(mat, method = "PELT", minseglen = 1),
                 "minseglen is too small")
})

test_that("cpt.reg errors when minseglen too large", {
  mat <- make_ar1_mat(ar1_data)
  expect_error(EnvCpt:::cpt.reg(mat, method = "PELT", minseglen = 100),
               "Minimum segment length is too large")
})

# ═══════════════════════════════════════════════════════════════════════════════
# 2. cpt.reg — FUNCTIONALITY
# ═══════════════════════════════════════════════════════════════════════════════
context("cpt.reg functionality")

test_that("cpt.reg PELT AR1 returns cpt.reg S4 object", {
  mat <- make_ar1_mat(ar1_data)
  result <- EnvCpt:::cpt.reg(mat, method = "PELT", minseglen = 3)
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg PELT AR1 detects changepoint near truth", {
  mat <- make_ar1_mat(ar1_data)
  result <- EnvCpt:::cpt.reg(mat, method = "PELT", minseglen = 3)
  expect_true(any(abs(cpts(result) - 100) <= 10))
})

test_that("cpt.reg AMOC AR1 returns cpt.reg S4 object", {
  mat <- make_ar1_mat(ar1_data)
  result <- EnvCpt:::cpt.reg(mat, method = "AMOC", minseglen = 3)
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg PELT AR2 returns cpt.reg S4 object", {
  mat <- make_ar2_mat(ar2_data)
  result <- EnvCpt:::cpt.reg(mat, method = "PELT", minseglen = 4)
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg with class=FALSE returns non-class output", {
  mat <- make_ar1_mat(ar1_data)
  result <- EnvCpt:::cpt.reg(mat, method = "PELT", 
                             minseglen = 3, class = FALSE)
  expect_false(isS4(result))
})

test_that("cpt.reg with param.estimates=FALSE returns cpt.reg", {
  mat <- make_ar1_mat(ar1_data)
  result <- EnvCpt:::cpt.reg(mat, method = "PELT", 
                             minseglen = 3, param.estimates = FALSE)
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg with AIC penalty works", {
  mat <- make_ar1_mat(ar1_data)
  result <- EnvCpt:::cpt.reg(mat, method = "PELT", 
                             penalty = "AIC", minseglen = 3)
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg with BIC penalty works", {
  mat <- make_ar1_mat(ar1_data)
  result <- EnvCpt:::cpt.reg(mat, method = "PELT", 
                             penalty = "BIC", minseglen = 3)
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg with shape=-1 (RSS cost) works", {
  mat <- make_ar1_mat(ar1_data)
  result <- EnvCpt:::cpt.reg(mat, method = "PELT", 
                             minseglen = 3, shape = -1)
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg with shape>0 (fixed variance) works", {
  mat <- make_ar1_mat(ar1_data)
  result <- EnvCpt:::cpt.reg(mat, method = "PELT", 
                             minseglen = 3, shape = 1)
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg with trend data works", {
  n_t <- length(trend_data)
  mat <- cbind(trend_data, rep(1, n_t), 1:n_t)
  result <- EnvCpt:::cpt.reg(mat, method = "PELT", minseglen = 5)
  expect_s4_class(result, "cpt.reg")
})

test_that("cpt.reg returns list for 3D array input", {
  mat1 <- make_ar1_mat(ar1_data)
  mat2 <- make_ar1_mat(ar1_data)
  arr  <- array(c(mat1, mat2), 
                dim = c(2, nrow(mat1), ncol(mat1)))
  result <- suppressWarnings(EnvCpt:::cpt.reg(arr, method = "PELT", minseglen = 3))
  expect_type(result, "list")
  expect_length(result, 2)
})

# ═══════════════════════════════════════════════════════════════════════════════
# 3. check_data
# ═══════════════════════════════════════════════════════════════════════════════
context("check_data validation")

test_that("check_data errors on non-array", {
  df <- as.data.frame(make_ar1_mat(ar1_data))
  expect_error(EnvCpt:::check_data(df),
               "Argument 'data' must be a numerical matrix")
})

test_that("check_data errors on 1D array", {
  expect_error(EnvCpt:::check_data(array(1:10, dim = 10)),
               "Argument 'data' must be a numerical matrix")
})

test_that("check_data errors when no regressors", {
  mat <- as.array(matrix(1:10, ncol = 1))
  expect_error(EnvCpt:::check_data(mat),
               "no regressors found")
})

test_that("check_data errors when more regressors than observations", {
  # 2 rows, 5 columns (1 response + 1 intercept + 3 regressors)
  mat <- as.array(matrix(c(1,2, 1,1, 3,4, 5,6, 7,8), nrow = 2))
  expect_error(EnvCpt:::check_data(mat))
})

test_that("check_data warns and adds intercept when missing", {
  mat <- as.array(cbind(ar1_data[-1], ar1_data[-length(ar1_data)]))
  expect_warning(EnvCpt:::check_data(mat),
                 "Missing intercept regressor")
})

test_that("check_data returns matrix unchanged when valid", {
  mat <- as.array(make_ar1_mat(ar1_data))
  result <- EnvCpt:::check_data(mat)
  expect_true(is.matrix(result))
})

# ═══════════════════════════════════════════════════════════════════════════════
# 4. ChangepointRegression
# ═══════════════════════════════════════════════════════════════════════════════
context("ChangepointRegression dispatch")

test_that("ChangepointRegression AMOC returns list", {
  mat <- make_ar1_mat(ar1_data)
  datai <- EnvCpt:::check_data(as.array(mat))
  result <- EnvCpt:::ChangepointRegression(
    datai, method = "AMOC", cpts.only = FALSE
  )
  expect_type(result, "list")
})

test_that("ChangepointRegression PELT returns list", {
  mat <- make_ar1_mat(ar1_data)
  datai <- EnvCpt:::check_data(as.array(mat))
  result <- EnvCpt:::ChangepointRegression(
    datai, method = "PELT", cpts.only = FALSE
  )
  expect_type(result, "list")
})

test_that("ChangepointRegression cpts.only=TRUE returns sorted vector", {
  mat <- make_ar1_mat(ar1_data)
  datai <- EnvCpt:::check_data(as.array(mat))
  result <- EnvCpt:::ChangepointRegression(
    datai, method = "PELT", cpts.only = TRUE
  )
  expect_true(is.numeric(result))
  expect_equal(result, sort(result))
})

test_that("ChangepointRegression errors on unrecognised method", {
  mat <- make_ar1_mat(ar1_data)
  datai <- EnvCpt:::check_data(as.array(mat))
  expect_error(
    EnvCpt:::ChangepointRegression(datai, method = "CROPS"),
    "method not recognised"
  )
})

# ═══════════════════════════════════════════════════════════════════════════════
# 5. envcpt — SUBSET MODELS
# ═══════════════════════════════════════════════════════════════════════════════
context("envcpt subset models")

if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  
  test_that("envcpt meanar1cpt detects AR1 changepoint", {
    out <- envcpt(ar1_data, models = "meanar1cpt")
    expect_s3_class(out, "envcpt")
    expect_s4_class(out$meanar1cpt, "cpt.reg")
    expect_true(any(abs(cpts(out$meanar1cpt) - 100) <= 10))
  })
  
  test_that("envcpt meanar2cpt runs AR2 changepoint model", {
    out <- envcpt(ar2_data, models = "meanar2cpt")
    expect_s3_class(out, "envcpt")
    expect_s4_class(out$meanar2cpt, "cpt.reg")
  })
  
  test_that("envcpt trendar1cpt works", {
    out <- envcpt(as.numeric(trend_data), 
                  models = "trendar1cpt", minseglen = 5)
    expect_s3_class(out, "envcpt")
  })
  
  test_that("envcpt trendar2cpt works", {
    out <- envcpt(as.numeric(trend_data), 
                  models = "trendar2cpt", minseglen = 5)
    expect_s3_class(out, "envcpt")
  })
  
  test_that("envcpt with numeric model indices 5 and 6", {
    out <- envcpt(ar1_data, models = c(5, 6))
    expect_s3_class(out, "envcpt")
  })
  
  test_that("envcpt trendcpt works", {
    out <- envcpt(as.numeric(trend_data), models = "trendcpt", 
                  minseglen = 5)
    expect_s4_class(out$trendcpt, "cpt.reg")
  })
  
  test_that("envcpt returns NA for skipped models", {
    out <- envcpt(ar1_data, models = "meanar1cpt")
    expect_true(is.na(out$mean))
  })
  
  test_that("envcpt minseglen respected in AR models", {
    out <- envcpt(ar1_data, models = "meanar1cpt", minseglen = 15)
    expect_s4_class(out$meanar1cpt, "cpt.reg")
  })
}

# ═══════════════════════════════════════════════════════════════════════════════
# 6. AICweights and BIC
# ═══════════════════════════════════════════════════════════════════════════════
context("AICweights and BIC")

test_that("AICweights returns numeric vector", {
  out <- envcpt(ar1_data)
  w   <- AICweights(out)
  expect_true(is.numeric(w))
})

test_that("AICweights sum to 1", {
  out <- envcpt(ar1_data)
  w   <- AICweights(out)
  expect_equal(sum(w, na.rm = TRUE), 1, tolerance = 1e-6)
})

test_that("AICweights all non-negative", {
  out <- envcpt(ar1_data)
  w   <- AICweights(out)
  expect_true(all(w >= 0, na.rm = TRUE))
})

test_that("AICweights errors on non-list", {
  tmp <- rnorm(100)
  class(tmp) <- "envcpt"
  expect_error(AICweights(tmp), "object argument must be a list")
})

test_that("AICweights errors on non-matrix summary", {
  tmp <- list(summary = rnorm(100))
  class(tmp) <- "envcpt"
  expect_error(AICweights(tmp), 
               "first element in the object list must be a matrix")
})

test_that("AICweights errors on non-numeric matrix", {
  tmp <- list(summary = matrix(LETTERS, nrow = 2))
  class(tmp) <- "envcpt"
  expect_error(AICweights(tmp), 
               "First two rows in matrix")
})

test_that("BIC returns numeric vector", {
  out <- envcpt(ar1_data)
  b   <- BIC(out)
  expect_true(is.numeric(b))
})

test_that("BIC length matches number of models", {
  out <- envcpt(ar1_data)
  b   <- BIC(out)
  expect_equal(length(b), 12)
})

test_that("BIC minimum is a valid model index", {
  out <- envcpt(ar1_data)
  b   <- BIC(out)
  expect_true(which.min(b) %in% 1:12)
})

test_that("BIC errors on non-matrix summary", {
  tmp <- list(summary = rnorm(10))
  class(tmp) <- "envcpt"
  expect_error(BIC(tmp), 
               "first element in the object list must be a matrix")
})

# ═══════════════════════════════════════════════════════════════════════════════
# 7. CptReg_AMOC_Normal and CptReg_PELT_Normal
# ═══════════════════════════════════════════════════════════════════════════════
context("CptReg_AMOC_Normal and CptReg_PELT_Normal")

test_that("CptReg_AMOC_Normal errors on invalid shape", {
  mat   <- make_ar1_mat(ar1_data)
  datai <- EnvCpt:::check_data(as.array(mat))
  expect_error(
    EnvCpt:::CptReg_AMOC_Normal(datai, shape = c(1, 2)),
    "Argument 'shape' is invalid"
  )
})

test_that("CptReg_AMOC_Normal errors on invalid dimensions", {
  bad <- as.matrix(data.frame(y = 1:5))
  expect_error(
    EnvCpt:::CptReg_AMOC_Normal(bad),
    "Invalid data dimensions"
  )
})

test_that("CptReg_AMOC_Normal returns list with cpts", {
  mat   <- make_ar1_mat(ar1_data)
  datai <- EnvCpt:::check_data(as.array(mat))
  result <- EnvCpt:::CptReg_AMOC_Normal(datai)
  expect_type(result, "list")
  expect_true("cpts" %in% names(result))
})

test_that("CptReg_PELT_Normal errors on invalid shape", {
  mat   <- make_ar1_mat(ar1_data)
  datai <- EnvCpt:::check_data(as.array(mat))
  expect_error(
    EnvCpt:::CptReg_PELT_Normal(datai, shape = c(1, 2)),
    "Argument 'shape' is invalid"
  )
})

test_that("CptReg_PELT_Normal returns list with cpts", {
  mat   <- make_ar1_mat(ar1_data)
  datai <- EnvCpt:::check_data(as.array(mat))
  result <- EnvCpt:::CptReg_PELT_Normal(datai)
  expect_type(result, "list")
  expect_true("cpts" %in% names(result))
})