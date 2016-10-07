\name{plot.fit.envcpt}
\alias{plot.fit.envcpt}
\title{
Plots the data and fits from each of the 8 models.
}
\description{
Uses the output of the \code{envcpt} function and plots the original data and the fit from each of the 8 models.
}
\usage{
plot.fit.envcpt(x,mar=c(2,8,1,1),lwd=3,...)
}
\arguments{
  \item{x}{
	A list produced as output from the \code{envcpt} function.  In essence a list with the named elements: "mean","meancpt","meanar","meanarcpt","trend","trendcpt","trendar","trendarcpt" where each element is the output from the appropriate function call.
}
  \item{mar}{Margins for the plot, see \code{\link[graphics]{par}} for further details.  Defaults give small margins around default content.}
  \item{lwd}{Line width graphical parameter, see \code{\link[graphics]{par}} for further details.}
  \item{...}{Extra graphical parameters, passed to the original plot and the individual calls to \code{lines}}
}
\details{
  Plots the data at the bottom and stacks the different fits for the 8 models from the \code{\link{envcpt}} function on top.  No scale is given as all data and fits are scaled to be in (0,1).  This is designed as an initial visualization tool for the fits only.
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
out=envcpt(x) # run the 8 models with default values
out[[1]] # first row is twice the negative log-likelihood for each model
         # second row is the number of parameters
aic.envcpt(out) # returns AIC for each model.
which.min(aic.envcpt(out)) # gives meancpt (model 2) as the best model fit.
out[[3]] # gives the model fit for the meancpt model.
plot.fits.envcpt(out) # plots the fits
plot.aic.envcpt(out) # plots the aic values

set.seed(10)
x=c(0.01*(1:100),1.5-0.02*((101:250)-101))+rnorm(250,0,0.2)
out=envcpt(x,minseglen=10) # run the 8 models with a minimum of 10 observations between changes
aic.envcpt(out) # returns the AIC for each model
which.min(aic.envcpt(out)) # gives trendcpt (model 6) as the best model fit.
out[[7]] # gives the model fit for the trendcpt model.
plot.fits.envcpt(out) # plots the fits
plot.aic.envcpt(out) # plots the aic values

set.seed(100)
x=arima.sim(model=list(ar=0.8),n=100)+5
out=envcpt(x) # run the 8 models with 
aic.envcpt(out) # returns the AIC for each model
which.min(aic.envcpt(out)) # gives trendar (model 7) as the best model fit.
out[[7]] # gives the model fit for the trendar model. Notice that the trend is tiny but does 
# produce a significantly better fit than the meanar model.
plot.fits.envcpt(out) # plots the fits
plot.aic.envcpt(out) # plots the aic values
}
}

\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
