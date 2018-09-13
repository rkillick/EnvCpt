\name{EnvCpt-package}
\alias{EnvCpt-package}
\alias{EnvCpt}
\docType{package}
\title{
Detection of Structural Changes in Climate and Environment Time Series
}
\description{
Tools for automatic model selection and diagnostics for Climate and Environmental data.  In particular the \code{envcpt()} function does automatic model selection between a variety of trend, changepoint and autocorrelation models.  The \code{envcpt()} function should be your first port of call.
}
\details{
\tabular{ll}{
Package: \tab EnvCpt\cr
Type: \tab Package\cr
Version: \tab 1.1.1 \cr
Date: \tab 2018-09-13\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Rebecca Killick <r.killick@lancs.ac.uk>, Claudie Beaulieu <c.beaulieu@soton.ac.uk>, Simon Taylor <s.taylor2@lancs.ac.uk>, Harjit Hullait <h.hullait@lancs.ac.uk>.

Maintainer: Rebecca Killick <r.killick@lancs.ac.uk>
}

\references{
EnvCpt Algorithm: Beaulieu, C, Killick, R (2018+) Distinguishing trends and shifts from memory in climate data. 

PELT Algorithm: Killick R, Fearnhead P, Eckley IA (2012) Optimal detection of changepoints with a linear computational cost, \emph{JASA} \bold{107(500)}, 1590--1598
}
\keyword{changepoint}
\keyword{segmentation}
\keyword{environment}
\keyword{climate}
\seealso{
\code{\link{envcpt}},\code{\link{AIC.envcpt}},\code{\link{plot.envcpt}}
}
\examples{
\dontrun{
set.seed(1)
x=c(rnorm(100,0,1),rnorm(100,5,1))
out=envcpt(x) # run all models with default values
out[[1]] # first row is twice the negative log-likelihood for each model
         # second row is the number of parameters
AIC(out) # returns AIC for each model.
which.min(AIC(out)) # gives meancpt (model 2) as the best model fit.
out$meancpt # gives the model fit for the meancpt model.
AICweights(out) # gives the AIC weights for each model
BIC(out) # returns the BIC for each model.
which.min(BIC(out)) # gives meancpt (model 2) as the best model fit too.
plot(out,type='fit') # plots the fits
plot(out,type="aic") # plots the aic values
plot(out,type="bic") # plots the bic values

set.seed(10)
x=c(0.01*(1:100),1.5-0.02*((101:250)-101))+rnorm(250,0,0.2)
out=envcpt(x,minseglen=10) # run all models with a minimum of 10 observations between changes
AIC(out) # returns the AIC for each model
which.min(AIC(out)) # gives trendcpt (model 8) as the best model fit.
out$trendcpt # gives the model fit for the trendcpt model.
AICweights(out) # gives the AIC weights for each model
BIC(out) # returns the BIC for each model.
which.min(BIC(out)) # gives trendcpt (model 8) as the best model fit too.
plot(out,type='fit') # plots the fits
plot(out,type="aic") # plots the aic values
plot(out,type="bic") # plots the bic values

set.seed(100)
x=arima.sim(model=list(ar=c(0.7,0.2)),n=500)+0.01*(1:500)
out=envcpt(x,models=c(3:6,9:12)) # runs a subset of models (those with AR components) 
AIC(out) # returns the AIC for each model
which.min(AIC(out)) # gives trendar2 (model 10) as the best model fit.
out$trendar2 # gives the model fit for the trendar2 model. Notice that the trend is tiny but does 
# produce a significantly better fit than the meanar2 model.
AICweights(out) # gives the AIC weights for each model
BIC(out) # returns the BIC for each model.
which.min(BIC(out)) # best fit is trendar2 (model 10) again.
plot(out,type='fit') # plots the fits
plot(out,type="aic") # plots the aic values
plot(out,type="bic") # plots the bic values
}
}
