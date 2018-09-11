\name{plot.envcpt}
\alias{plot.envcpt}
\title{
Plots optionally either ("fit") the data and fits from each of the 8 models or ("aic") the aic for each of the 8 models.
}
\description{
Uses the output of the \code{envcpt} function and plots optionally ("fit") the original data and the fit from each of the 8 models or ("aic") the aic for each of the 8 models.
}
\usage{
\method{plot}{envcpt}(x,type=c('fit','bic','aic'),lwd=3,...,data=NA)
}
\arguments{
  \item{x}{
	A list produced as output from the \code{envcpt} function.  In essence a list with the named elements: "mean","meancpt","meanar1","meanar2","meanar1cpt","meanar2cpt","trend","trendcpt","trendar1","trendar2","trendar1cpt","trendar2cpt" where each element is the output from the appropriate function call.
}
  \item{type}{character vector. \code{fit} returns a plot of the data and fitted models. \code{aic} returns a bar chart of the aic values for each of the 8 models}
  \item{lwd}{Line width graphical parameter, see \code{\link[graphics]{par}} for further details.}
  \item{...}{Extra graphical parameters, passed to the original plot and the individual calls to \code{lines} or \code{barchar}.}
  \item{data}{This argument is only required when \code{x} doesn't contain any fits that provide the data which is required for the "fit" \code{type} option.  Default is NA as most of the time this is likely not required.}
}
\details{
  If \code{type="fit"}, the function plots the data at the bottom and stacks the different fits for the 8 models from the \code{\link{envcpt}} function on top.  No scale is given as all data and fits are scaled to be in (0,1).  This is designed as an initial visualization tool for the fits only.
  If \code{type="aic"} the function uses the \code{\link{AIC.envcpt}} function to calculate the AIC values for the \code{envcpt} output \code{x}.  Then barcharts the AIC values in the same order as the \code{type="fit"} option.  The minimum AIC is the preferred model and this is highlight by a solid block.  This is designed as an initial visualization tool for the AIC values only.
  If \code{type="bic"} the function uses the \code{\link{BIC.envcpt}} function to calculate the BIC values for the \code{envcpt} output \code{x}.  Then barcharts the BIC values in the same order as the \code{type="fit"} option.  The minimum BIC is the preferred model and this is highlight by a solid block.  This is designed as an initial visualization tool for the BIC values only.
}
\value{
  Returns the printed graphic to the active device.
}
\author{
Rebecca Killick & Claudie Beaulieu.
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
\keyword{ts}
