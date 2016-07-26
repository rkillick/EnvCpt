envcpt=function(data,minseglen=5,...,verbose=TRUE){
  # 8 models: mean, mean+cpt, mean+AR, mean+AR+cpt, trend, trend+AR, trend+cpt, trend+AR+cpt
  # assume normal errors throughout
  if(verbose==TRUE){
    print("Fitting 8 models")
    pb <- txtProgressBar(min = 0, max =8, style = 3)
  }
  # mean 
  mean.fit=fitdistr(data,densfun='normal')# function from MASS
  mean.loglik=-2*mean.fit$loglik 
  # literally the loglikelihood so need to do -2*
  if(verbose==TRUE){setTxtProgressBar(pb, 1)}

  # mean+cpt
  meancpt.fit=cpt.meanvar(data,method='PELT',minseglen=minseglen,...) # default MBIC penalty
  meancpt.loglik=logLik(meancpt.fit)[1] # function from changepoint
  # gives -2*loglik
  if(verbose==TRUE){setTxtProgressBar(pb, 2)}
  
  # mean+AR
  # use logLik(auto.arima(data,d=0,D=0,seasonal=FALSE))
  meanar.fit=try(arima(data,order=c(1,0,0),method="CSS-ML"))#function from stats
  if(class(meanar.fit)=='try-error'){
    meanar.loglik=NA
  }
  else{meanar.loglik=-2*logLik(meanar.fit)}
  # again this gives loglikelihood so need to do -2*
  if(verbose==TRUE){setTxtProgressBar(pb, 3)}
  
  # mean+AR+cpt
  meanarcpt.fit=cpt.reg(cbind(data[-1],rep(1,length(data)-1),data[-length(data)]),method="PELT",minseglen=minseglen,...) # default MBIC penalty
  meanarcpt.loglik=logLik(meanarcpt.fit)[1] # function from changepoint
  # gives -2*log likelihood
  # replace with faster version
  if(verbose==TRUE){setTxtProgressBar(pb, 4)}
  
  # trend
  trend.fit=lm(data~c(1:length(data))) # function from stats
  trend.loglik=-2*logLik(trend.fit)
  # again this gives loglikelihood so need to do -2*
  if(verbose==TRUE){setTxtProgressBar(pb, 5)}
  
  # trend+cpt
  trendcpt.fit=cpt.reg(cbind(data,rep(1,length(data)),1:length(data)),method='PELT',minseglen=minseglen,...) # default MBIC penalty
  trendcpt.loglik=logLik(trendcpt.fit)[1] # function from changepoint
  # gives -2*log likelihood
  if(verbose==TRUE){setTxtProgressBar(pb, 6)}
  
  # trend+AR
  trendar.fit=lm(data[-1]~c(1:(length(data)-1))+data[-length(data)]) # function from stats
  trendar.loglik=-2*logLik(trendar.fit)
  # again this gives loglikelihood so need to do -2*
  if(verbose==TRUE){setTxtProgressBar(pb, 7)}
  
  # trend+AR+cpt
  trendarcpt.fit=cpt.reg(cbind(data[-1],rep(1,length(data)-1),1:(length(data)-1),data[-length(data)]),method="PELT",minseglen=minseglen,...) # default MBIC penalty
  trendarcpt.loglik=logLik(trendarcpt.fit)[1] # function from changepoint
  # gives -2*log likelihood
  # replace with faster version
  if(verbose==TRUE){setTxtProgressBar(pb, 8)}
  
  out=list()
  logLik=c(mean.loglik,meancpt.loglik,meanar.loglik,meanarcpt.loglik,trend.loglik,trendcpt.loglik,trendar.loglik,trendarcpt.loglik)
  nparam=c(2,ncpts(meancpt.fit)+nseg(meancpt.fit)*2,3,ncpts(meanarcpt.fit)+nseg(meanarcpt.fit)*3,
           3,ncpts(trendcpt.fit)+nseg(trendcpt.fit)*3,4,ncpts(trendarcpt.fit)+nseg(trendarcpt.fit)*4
           )
  out[[1]]=rbind(logLik=logLik,nparam=nparam)
  colnames(out[[1]])=c("mean","meancpt","meanar","meanarcpt","trend","trendcpt","trendar","trendarcpt")
  out[[2]]=mean.fit
  out[[3]]=meancpt.fit
  out[[4]]=meanar.fit
  out[[5]]=meanarcpt.fit
  out[[6]]=trend.fit
  out[[7]]=trendcpt.fit
  out[[8]]=trendar.fit
  out[[9]]=trendarcpt.fit
  if(verbose==TRUE){close(pb)}
  return(out)
}