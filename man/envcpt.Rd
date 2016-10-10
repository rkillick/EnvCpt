\name{envcpt}
\alias{envcpt}
\title{
Assesses whether an environmental time series contains trend, autocorrelation and/or changes
}
\description{
Evaluates 8 different models (see details) and returns the model fits as well as a summary of the likelihood for each model.
}
\usage{
envcpt(data,minseglen=5,...,verbose=TRUE)
}
\arguments{
  \item{data}{
	A vector or ts object containing the data to fit the models to.
}
  \item{minseglen}{
  Positive integer giving the minimum segment length (no. of observations between changes) for the changepoint models, default is the minimum allowed by theory (for the largest model).
}
  \item{...}{
  Additional arguments to pass to the changepoint functions, if none are specified defaults with PELT multiple changepoint algorithm are used. See \code{\link[changepoint]{cpt.meanvar}} for options.
}
  \item{verbose}{
  If TRUE (default), prints to the console an progress bar indicating progression through the fit of the 8 models.
}
}
\details{
	This function is used to automatically fit 8 different models all with Normal distribution for errors:
	1. A constant mean and variance (using \code{\link[MASS]{fitdistr}})
	2. A piecewise constant mean and variance (using \code{\link[changepoint]{cpt.meanvar}})
	3. A constant mean with AR (1) errors (using \code{\link[forecast]{auto.arima}})
	4. A piecewise constant mean with AR(1) errors (using \code{cpt.reg} in this package - not exported)
	5. A linear trend over time (using \code{\link[stats]{lm}})
	6. A piecewise linear trend over time (using \code{cpt.reg} in this package - not exported)
	7. A linear trend over time with AR(1) errors (using \code{\link[stats]{lm}})
	8. A piecewise linear trend over time with AR(1) errors (using \code{cpt.reg} in this package - not exported)
	The default values for each function are used, except for the changepoint functions where multiple changes are identified and thus the PELT algorithm is used (see references or changepoint package for details).
}
\value{
  \code{envcpt} outputs a list of 9 elements.  The first element is a 2x8 matrix where the first row contains the likelihood for each model fit and the second row contains the number of parameters fit in each model.  The 8 columns are for the 8 different models in the order above (headings given).  If any element is NA then there was an error fitting this type of model (typically the AR models and this implies nonstationarity).
  
  Elements 2-9 of the list are the fits for each individual model, these are the direct output from the respective functions so see the individual functions for formats.  The first model fit is in element 2 and the eighth model fit is in element 9, but are named for convenience.
}
\references{
PELT Algorithm: Killick R, Fearnhead P, Eckley IA (2012) Optimal detection of changepoints with a linear computational cost, \emph{JASA} \bold{107(500)}, 1590--1598

MBIC: Zhang, N. R. and Siegmund, D. O. (2007) A Modified Bayes Information Criterion with Applications to the Analysis of Comparative Genomic Hybridization Data. \emph{Biometrics} \bold{63}, 22-32.
}
\author{
Rebecca Killick
}


\seealso{
\code{\link[changepoint]{cpt.meanvar}},\code{\link[changepoint]{plot-methods}},\code{\linkS4class{cpt}}, \code{\link{plot.envcpt}}
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
\keyword{ts}
