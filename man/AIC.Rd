\name{AIC.envcpt}
\alias{AIC.envcpt}
\title{
Computes the AIC measure for output from the envcpt function.
}
\description{
Uses the likelihood and number of parameters from the output of the \code{envcpt} function and calculates the AIC measure for each model.
}
\usage{
\method{AIC}{envcpt}(object,...,k=2)
}
\arguments{
  \item{object}{
	A list produced as output from the \code{envcpt} function.  In essence any list where the first element contains a matrix whose first row is twice the negative log-likelihood and second row is the number of parameters.  Columns are different models to compare.
}
  \item{...}{Not used}
  \item{k}{numeric, the penalty per parameter to be used: the default \code{k=2} is the classical AIC.}
}
\details{
  Calculates the AIC defined as -2*<log-likelihood> + 2*<number of parameters>.  When comparing models the smaller the AIC the better the fit.
}
\value{
  Vector of AIC values the same length as the number of columns in the first entry to the input list (length 8 if output from envcpt is used).  The column names from the envcpt output are preserved to give clear indication on models.
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
out=envcpt(x) # run the 8 models with default values
out[[1]] # first row is twice the negative log-likelihood for each model
         # second row is the number of parameters
AIC(out) # returns AIC for each model.
which.min(AIC(out)) # gives meancpt (model 2) as the best model fit.
out[[3]] # gives the model fit for the meancpt model.
plot(out,type='fit') # plots the fits
plot(out,type="aic") # plots the aic values

set.seed(10)
x=c(0.01*(1:100),1.5-0.02*((101:250)-101))+rnorm(250,0,0.2)
out=envcpt(x,minseglen=10) # run the 8 models with a minimum of 10 observations between changes
AIC(out) # returns the AIC for each model
which.min(AIC(out)) # gives trendcpt (model 6) as the best model fit.
out[[7]] # gives the model fit for the trendcpt model.
plot(out,type='fit') # plots the fits
plot(out,type="aic") # plots the aic values

set.seed(100)
x=arima.sim(model=list(ar=0.8),n=100)+5
out=envcpt(x) # run the 8 models with 
AIC(out) # returns the AIC for each model
which.min(AIC(out)) # gives trendar (model 7) as the best model fit.
out[[7]] # gives the model fit for the trendar model. Notice that the trend is tiny but does 
# produce a significantly better fit than the meanar model.
plot(out,type='fit') # plots the fits
plot(out,type="aic") # plots the aic values
}
}

\keyword{methods}
\keyword{univar}
\keyword{models}