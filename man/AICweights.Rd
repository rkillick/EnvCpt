\name{AICweights}
\alias{AICweights.envcpt}
\alias{AICweights}
\title{
Computes the AIC weights for output from the envcpt function.
}
\description{
Uses the likelihood and number of parameters from the output of the \code{envcpt} function and calculates the AICweights for each model fitted relative to the model with the minimum AIC.
}
\usage{
\method{AICweights}{envcpt}(object)
}
\arguments{
  \item{object}{
	A list produced as output from the \code{envcpt} function.  In essence any list where the first element contains a matrix whose first row is twice the negative log-likelihood and second row is the number of parameters.  Columns are different models to compare.
}
}
\details{
  Calculates the AICweights defined as \deqn{w_i=\frac{\exp{-0.5 \Delta_i}}{\sum_r{\exp{-0.5\Delta_r}}}}{w_i=exp(-0.5 \Delta_i)  /  \sum exp(-0.5\Delta_r)} where the summation over \eqn{r} is across all models considered and \eqn{\Delta_i} is the difference between the AIC value for model \eqn{i} and the best model.
}
\value{
  Vector of AICweights the same length as the number of columns in the first entry to the input list (length 12 if output from envcpt where all models are considered).  The column names from the envcpt output are preserved to give clear indication on models.
}
\author{
Rebecca Killick
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