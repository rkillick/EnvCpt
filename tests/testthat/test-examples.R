context("man file example tests")

# From envcpt.Rd
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,5,1))
out=envcpt(x) # run the 8 models with default values
truth=rbind(c(949.25169085, 535.4859516, 683.8641559, 532.9723591, 735.41930460, 535.4139881, 648.1919606,532.7957423),
            c(2,5,3,7,3,7,4,9))
test_that('envcpt1',expect_equivalent(out[[1]],truth))

truth=c(953.25169085,545.4859516,689.8641559,546.9723591,741.41930460,549.4139881,656.1919606,550.7957423)
test_that("envcpt2",expect_equivalent(aic.envcpt(out),truth))

test_that("envcpt3",expect_equivalent(which.min(aic.envcpt(out)),2))

test_that("envcpt4",expect_is(out[[3]],"cpt"))
test_that("envcpt5",expect_equivalent(cpts(out[[3]]),100))



set.seed(10)
x=c(0.01*(1:100),1.5-0.02*((101:250)-101))+rnorm(250,0,0.2)
out=envcpt(x,minseglen=10) # run the 8 models with a minimum of 10 observations between changes
truth=c(578.94707,  -67.51376,   27.29077,   25.89202,  468.36581, -113.41437,   20.92805, -109.94452)
test_that("envcpt6",expect_equivalent(aic.envcpt(out),truth))

test_that("envcpt7",expect_equivalent(which.min(aic.envcpt(out)),6))

test_that("envcpt8",expect_is(out[[7]],"cpt.reg"))
test_that("envcpt9",expect_equivalent(cpts(out[[7]]),100))



set.seed(1)
x=arima.sim(model=list(ar=0.6),n=100)+5
out=envcpt(x) # run the 8 models with 
truth=c(308.0158490,308.0158490,269.1805846,262.4384534,309.2885536,262.881711,264.4162925,264.4213773)
test_that("envcpt10",expect_equivalent(aic.envcpt(out),truth))

test_that("envcpt11",expect_equivalent(which.min(aic.envcpt(out)),4))

test_that("envcpt12",expect_is(out[[5]],"cpt.reg"))
test_that("envcpt13",expect_equivalent(cpts(out[[5]]),numeric()))






# From AIC.Rd
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,5,1))
out=envcpt(x) # run the 8 models with default values
truth=rbind(c(949.25169085, 535.4859516, 683.8641559, 532.9723591, 735.41930460, 535.4139881, 648.1919606,532.7957423),
            c(2,5,3,7,3,7,4,9))
test_that('envcpt14',expect_equivalent(out[[1]],truth))

truth=c(953.25169085,545.4859516,689.8641559,546.9723591,741.41930460,549.4139881,656.1919606,550.7957423)
test_that("envcpt15",expect_equivalent(aic.envcpt(out),truth))

test_that("envcpt16",expect_equivalent(which.min(aic.envcpt(out)),2))

test_that("envcpt17",expect_is(out[[3]],"cpt"))
test_that("envcpt18",expect_equivalent(cpts(out[[3]]),100))



set.seed(10)
x=c(0.01*(1:100),1.5-0.02*((101:250)-101))+rnorm(250,0,0.2)
out=envcpt(x,minseglen=10) # run the 8 models with a minimum of 10 observations between changes
truth=c(578.94707,  -67.51376,   27.29077,   25.89202,  468.36581, -113.41437,   20.92805, -109.94452)
test_that("envcpt19",expect_equivalent(aic.envcpt(out),truth))

test_that("envcpt20",expect_equivalent(which.min(aic.envcpt(out)),6))

test_that("envcpt21",expect_is(out[[7]],"cpt.reg"))
test_that("envcpt22",expect_equivalent(cpts(out[[7]]),100))



set.seed(1)
x=arima.sim(model=list(ar=0.6),n=100)+5
out=envcpt(x) # run the 8 models with 
truth=c(308.0158490,308.0158490,269.1805846,262.4384534,309.2885536,262.881711,264.4162925,264.4213773)
test_that("envcpt23",expect_equivalent(aic.envcpt(out),truth))

test_that("envcpt24",expect_equivalent(which.min(aic.envcpt(out)),4))

test_that("envcpt25",expect_is(out[[5]],"cpt.reg"))
test_that("envcpt26",expect_equivalent(cpts(out[[5]]),numeric()))

