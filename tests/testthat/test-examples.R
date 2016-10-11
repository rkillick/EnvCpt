context("man file example tests")

# From envcpt.Rd & plot.fit.envcpt.Rd (same examples)
if(identical(Sys.getenv("NOT_CRAN"), "true")) {
  set.seed(1)
  x=c(rnorm(100,0,1),rnorm(100,5,1))
  out=envcpt(x) # run the 8 models with default values
  truth=rbind(c(949.25169085, 535.4859516, 683.8641559, 532.9723591, 735.41930460, 535.4139881, 648.1919606,532.7957423),
              c(2,5,3,7,3,7,4,9))
  
  test_that('envcpt1',expect_equivalent(out[[1]],truth))
  
  truth=c(953.25169085,545.4859516,689.8641559,546.9723591,741.41930460,549.4139881,656.1919606,550.7957423)
  test_that("envcpt2",expect_equivalent(AIC(out),truth))
  
  test_that("envcpt3",expect_equivalent(which.min(AIC(out)),2))
  
  test_that("envcpt4",expect_is(out[[3]],"cpt"))
  test_that("envcpt5",expect_equivalent(cpts(out[[3]]),100))
  
  test_that("envcpt5plot",expect_silent(plot(out,type='fit')))
  test_that("envcpt5aicplot",expect_silent(plot(out,type="aic")))
  
  test_that("messages1",expect_message(envcpt(x),"Fitting 8 models"))
}





if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  set.seed(10)
  x=c(0.01*(1:100),1.5-0.02*((101:250)-101))+rnorm(250,0,0.2)
  out=envcpt(x,minseglen=10) # run the 8 models with a minimum of 10 observations between changes
  truth=c(578.94707,  -67.51376,   27.29077,   27.29077,  468.36581, -113.41437,   20.92805, -109.94452)
test_that("envcpt6",expect_equivalent(AIC(out),truth))

test_that("envcpt7",expect_equivalent(which.min(AIC(out)),6))

test_that("envcpt8",expect_is(out[[7]],"cpt.reg"))
test_that("envcpt9",expect_equivalent(cpts(out[[7]]),100))

test_that("envcpt9plot",expect_silent(plot(out,type='fit')))
test_that("envcpt9aicplot",expect_silent(plot(out,type="aic")))

test_that("messages2",expect_message(envcpt(x),"Fitting 8 models"))
}







if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  set.seed(100)
  x=arima.sim(model=list(ar=0.8),n=100)+5
  out=envcpt(x) # run the 8 models with 
  truth=c(330.121571,330.121571,291.371568,291.371568,324.688444,324.688444,287.383172,287.383172)
test_that("envcpt10",expect_equivalent(AIC(out),truth))

test_that("envcpt11",expect_equivalent(which.min(AIC(out)),7))

test_that("envcpt12",expect_is(out[[8]],"lm"))
test_that("envcpt13",expect_equivalent(out[[8]]$coefficients, c(2.538210310,-0.006158614,0.539862825)))

test_that("envcpt13plot",expect_silent(plot(out,type='fit')))
test_that("envcpt13aicplot",expect_silent(plot(out,type="aic")))

test_that("messages3",expect_message(envcpt(x),"Fitting 8 models"))
}




