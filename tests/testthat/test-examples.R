context("man file example tests")

# From envcpt.Rd & plot.fit.envcpt.Rd (same examples)
if(identical(Sys.getenv("NOT_CRAN"), "true")) {
  set.seed(1)
  x=c(rnorm(100,0,1),rnorm(100,5,1))
  out=envcpt(x) # run the 8 models with default values
  test_that('envcpt1',expect_equal_to_reference(out[[1]],file='envcpt1.rds'))
  
  test_that("envcpt2",expect_equal_to_reference(AIC(out),file="envcpt2.rds"))
  
  test_that("envcpt3",expect_equivalent(which.min(AIC(out)),2))
  
  test_that("envcpt4",expect_is(out[[3]],"cpt"))
  test_that("envcpt5",expect_equivalent(cpts(out[[3]]),100))
  
  test_that("envcpt5plot",expect_silent(plot(out,type='fit')))
  test_that("envcpt5aicplot",expect_silent(plot(out,type="aic")))
  
  test_that("messages1",expect_message(envcpt(x),"Fitting 12 models"))
}





if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  set.seed(10)
  x=c(0.01*(1:100),1.5-0.02*((101:250)-101))+rnorm(250,0,0.2)
  out=envcpt(x,minseglen=10) # run the 8 models with a minimum of 10 observations between changes
test_that("envcpt6",expect_equal_to_reference(AIC(out),file='envcpt6.rds'))

test_that("envcpt7",expect_equivalent(which.min(AIC(out)),8))

test_that("envcpt8",expect_is(out[[9]],"cpt.reg"))
test_that("envcpt9",expect_equivalent(cpts(out[[9]]),100))

test_that("envcpt9plot",expect_silent(plot(out,type='fit')))
test_that("envcpt9aicplot",expect_silent(plot(out,type="aic")))

test_that("messages2",expect_message(envcpt(x),"Fitting 12 models"))
}







if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  set.seed(65312)
  x=arima.sim(model=list(ar=c(0.7,0.2)),n=500)+0.01*(1:500)
  out=envcpt(x) # run the 8 models with 
test_that("envcpt10",expect_equal_to_reference(AIC(out),file='envcpt10.rds'))

test_that("envcpt11",expect_equivalent(which.min(AIC(out)),10))

test_that("envcpt12",expect_is(out[[11]],"lm"))
test_that("envcpt13",expect_equivalent(out[[11]]$coefficients, c(-0.016366180,     0.001629478,     0.639355022,     0.227675457)))

test_that("envcpt13plot",expect_silent(plot(out,type='fit')))
test_that("envcpt13aicplot",expect_silent(plot(out,type="aic")))

test_that("messages3",expect_message(envcpt(x),"Fitting 12 models"))
}




