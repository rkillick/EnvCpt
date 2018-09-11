\name{AIC.envcpt}
\alias{AIC.envcpt}
\alias{BIC.envcpt}
\title{
Computes the AIC/BIC measure for output from the envcpt function.
}
\description{
Uses the likelihood and number of parameters from the output of the \code{envcpt} function and calculates the AIC/BIC measure for each model.
}
\usage{
\method{AIC}{envcpt}(object,...,k=2)
\method{BIC}{envcpt}(object,...)
}
\arguments{
  \item{object}{
	A list produced as output from the \code{envcpt} function.  In essence any list where the first element contains a matrix whose first row is twice the negative log-likelihood and second row is the number of parameters.  Columns are different models to compare.
}
  \item{...}{Used to pass the length of the data to the BIC function as argument \code{n}.}
  \item{k}{numeric, the penalty per parameter to be used: the default \code{k=2} is the classical AIC. The value log(n) is hard coded into BIC.}
}
\details{
  Calculates the AIC defined as -2*<log-likelihood> + 2*<number of parameters> or the BIC as -2*<log-likelihood> + log(n)*<number of parameters>.  When comparing models the smaller the AIC/BIC the better the fit.
}
\value{
  Vector of AIC/BIC values the same length as the number of columns in the first entry to the input list (length 12 if output from envcpt is used).  The column names from the envcpt output are preserved to give clear indication on models.
}
\author{
Simon Taylor and Rebecca Killick
}


\seealso{
\code{\link{envcpt}}
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

\keyword{methods}
\keyword{univar}
\keyword{models}