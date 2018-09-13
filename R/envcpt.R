envcpt=function(data,models=c("mean","meancpt","meanar1","meanar2","meanar1cpt","meanar2cpt","trend","trendcpt","trendar1","trendar2","trendar1cpt","trendar2cpt"),minseglen=5,...,verbose=TRUE){
  # 8 models: mean, mean+cpt, mean+AR, mean+AR+cpt, trend, trend+AR, trend+cpt, trend+AR+cpt
  # assume normal errors throughout
  if(any(!complete.cases(data))){stop("data has missing values, this function cannot handle missing values")}
  if(any(!is.numeric(data))){stop("data must be a numeric vector")}
  n=length(data)
  
  if(verbose==TRUE){
    message("Fitting ",length(models)," models")
    pb <- txtProgressBar(min = 0, max =length(models), style = 3)
  }
  # mean 
  if(any(models=="mean")|any(models==1)){
    mean.fit=fitdistr(data,densfun='normal')# function from MASS
    mean.loglik=-2*mean.fit$loglik 
    # literally the loglikelihood so need to do -2*
    if(verbose==TRUE){setTxtProgressBar(pb, 1)}
  }
  else(mean.fit=mean.loglik=NA)

  # mean+cpt
  if(any(models=="meancpt")|any(models==2)){
    meancpt.fit=cpt.meanvar(data,method='PELT',minseglen=minseglen,...) # default MBIC penalty
    meancpt.loglik=logLik(meancpt.fit)[1] # function from changepoint
    # gives -2*loglik
    if(verbose==TRUE){setTxtProgressBar(pb, 2)}
  }
  else{meancpt.fit=new("cpt");meancpt.loglik=NA}
  
  # mean+AR1
  if(any(models=="meanar1")|any(models==3)){
    # use logLik(auto.arima(data,d=0,D=0,seasonal=FALSE))
    meanar.fit=try(arima(data,order=c(1,0,0),method="CSS-ML"))#function from stats
    if(any(class(meanar.fit)=='try-error')){
      meanar.loglik=NA
    }
    else{meanar.loglik=-2*logLik(meanar.fit)}
    # again this gives loglikelihood so need to do -2*
    if(verbose==TRUE){setTxtProgressBar(pb, 3)}
  }
  else{meanar.fit=meanar.loglik=NA}
  
  # mean+AR2
  if(any(models=="meanar2")|any(models==4)){
    # use logLik(auto.arima(data,d=0,D=0,seasonal=FALSE))
    meanar2.fit=try(arima(data,order=c(2,0,0),method="CSS-ML"))#function from stats
    if(any(class(meanar2.fit)=='try-error')){
      meanar2.loglik=NA
    }
    else{meanar2.loglik=-2*logLik(meanar2.fit)}
    # again this gives loglikelihood so need to do -2*
    if(verbose==TRUE){setTxtProgressBar(pb, 4)}
  }
  else{meanar2.fit=meanar2.loglik=NA}

  # mean+AR1+cpt
  if(any(models=="meanar1cpt")|any(models==5)){
    meanarcpt.fit=cpt.reg(cbind(data[-1],rep(1,n-1),data[-n]),method="PELT",minseglen=minseglen,...) # default MBIC penalty
    if(ncpts(meanarcpt.fit)==0){
      if(any(models=="meanar1")|any(models==3)){meanarcpt.loglik=meanar.loglik}
      else{meanarcpt.loglik=logLik(meanarcpt.fit)[1]} # function from changepoint
    # used as cpt fit is not the same as no cpt fit due to regression not using the first data
    }
    else{
      meanarcpt.loglik=logLik(meanarcpt.fit)[1] # function from changepoint
      # gives -2*log likelihood
    }
    # replace with faster version (maybe trimmed version of auto.arima instead of lm)
    if(verbose==TRUE){setTxtProgressBar(pb, 5)}
  }
  else{meanarcpt.fit=new("cpt");meanarcpt.loglik=NA}
  
  # mean+AR2+cpt
  if(any(models=="meanar2cpt")|any(models==6)){
    meanar2cpt.fit=cpt.reg(cbind(data[-c(1:2)],rep(1,n-2),data[2:(n-1)],data[1:(n-2)]),method="PELT",minseglen=minseglen,...) # default MBIC penalty
    if(ncpts(meanar2cpt.fit)==0){
      if(any(models=="meanar2")|any(models==4)){meanar2cpt.loglik=meanar.loglik}
      else{meanar2cpt.loglik=logLik(meanar2cpt.fit)[1]} # function from changepoint
    # used as cpt fit is not the same as no cpt fit due to regression not using the first data
    }
    else{
      meanar2cpt.loglik=logLik(meanar2cpt.fit)[1] # function from changepoint
      # gives -2*log likelihood
    }
    # replace with faster version (maybe trimmed version of auto.arima instead of lm)
    if(verbose==TRUE){setTxtProgressBar(pb, 6)}
  }
  else{meanar2cpt.fit=new("cpt");meanar2cpt.loglik=NA}

  # trend
  if(any(models=="trend")|any(models==7)){
    trend.fit=lm(data~c(1:n)) # function from stats
    trend.loglik=-2*logLik(trend.fit)
    # again this gives loglikelihood so need to do -2*
    if(verbose==TRUE){setTxtProgressBar(pb, 7)}
  }
  else{trend.fit=trend.loglik=NA}
  
  # trend+cpt
  if(any(models=="trendcpt")|any(models==8)){
    trendcpt.fit=cpt.reg(cbind(data,rep(1,n),1:n),method='PELT',minseglen=minseglen,...) # default MBIC penalty
    if(ncpts(trendcpt.fit)==0){
      if(any(models=="trend")|any(models==7)){trendcpt.loglik=trend.loglik}
      else{trendcpt.loglik=logLik(trendcpt.fit)[1]} # function from changepoint
    }
    else{
      trendcpt.loglik=logLik(trendcpt.fit)[1] # function from changepoint
      # gives -2*log likelihood
    }
    if(verbose==TRUE){setTxtProgressBar(pb, 8)}
  }
  else{trendcpt.fit=new("cpt");trendcpt.loglik=NA}
  
  # trend+AR1
  if(any(models=="trendar1")|any(models==9)){
    trendar.fit=lm(data[-1]~c(1:(n-1))+data[-n]) # function from stats
    trendar.loglik=-2*logLik(trendar.fit)
    # again this gives loglikelihood so need to do -2*
    if(verbose==TRUE){setTxtProgressBar(pb, 9)}
  }
  else{trendar.fit=trendar.loglik=NA}

  # trend+AR2
  if(any(models=="trendar2")|any(models==10)){
    trendar2.fit=lm(data[-c(1:2)]~c(1:(n-2))+data[2:(n-1)]+data[1:(n-2)]) # function from stats
    trendar2.loglik=-2*logLik(trendar2.fit)
    # again this gives loglikelihood so need to do -2*
    if(verbose==TRUE){setTxtProgressBar(pb, 10)}
  }
  else{trendar2.fit=trendar2.loglik=NA}
  
  # trend+AR1+cpt
  if(any(models=="trendar1cpt")|any(models==11)){
    trendarcpt.fit=cpt.reg(cbind(data[-1],rep(1,n-1),1:(n-1),data[-n]),method="PELT",minseglen=minseglen,...) # default MBIC penalty
    if(ncpts(trendarcpt.fit)==0){
      if(any(models=="trendar1")|any(models==9)){trendarcpt.loglik=trendar.loglik}
      else{trendarcpt.loglik=logLik(trendarcpt.fit)[1]} # function from changepoint
    # used as cpt fit is not the same as no cpt fit due to regression not using the first data
    }
    else{
      trendarcpt.loglik=logLik(trendarcpt.fit)[1] # function from changepoint
      # gives -2*log likelihood
    }
    # replace with faster version, maybe trimmed auto.arima instead of lm
    if(verbose==TRUE){setTxtProgressBar(pb, 11)}
  }
  else{trendarcpt.fit=new("cpt");trendarcpt.loglik=NA}

  # trend+AR2+cpt
  if(any(models=="trendar2cpt")|any(models==12)){
    trendar2cpt.fit=cpt.reg(cbind(data[-c(1:2)],rep(1,n-2),1:(n-2),data[2:(n-1)],data[1:(n-2)]),method="PELT",minseglen=minseglen,...) # default MBIC penalty
    if(ncpts(trendar2cpt.fit)==0){
      if(any(models=="trendar2")|any(models==10)){trendar2cpt.loglik=trendar2.loglik}
      else{trendar2cpt.loglik=logLik(trendar2cpt.fit)[1]} # function from changepoint
    # used as cpt fit is not the same as no cpt fit due to regression not using the first data
    }
    else{
      trendar2cpt.loglik=logLik(trendar2cpt.fit)[1] # function from changepoint
      # gives -2*log likelihood
    }
    # replace with faster version, maybe trimmed auto.arima instead of lm
    if(verbose==TRUE){setTxtProgressBar(pb, 12)}
  }
  else{trendar2cpt.fit=new("cpt");trendar2cpt.loglik=NA}
  
  out=list()
  logLik=c(mean.loglik,meancpt.loglik,meanar.loglik,meanar2.loglik,meanarcpt.loglik,meanar2cpt.loglik,trend.loglik,trendcpt.loglik,trendar.loglik,trendar2.loglik,trendarcpt.loglik,trendar2cpt.loglik)
  nparam=c(2,ncpts(meancpt.fit)+nseg(meancpt.fit)*2,3,4,ncpts(meanarcpt.fit)+nseg(meanarcpt.fit)*3,ncpts(meanarcpt.fit)+nseg(meanarcpt.fit)*4,
           3,ncpts(trendcpt.fit)+nseg(trendcpt.fit)*3,4,5,ncpts(trendarcpt.fit)+nseg(trendarcpt.fit)*4,ncpts(trendarcpt.fit)+nseg(trendarcpt.fit)*5
           )
  out[[1]]=rbind(logLik=logLik,nparam=nparam)
  colnames(out[[1]])=c("mean","meancpt","meanar1","meanar2","meanar1cpt","meanar2cpt","trend","trendcpt","trendar1","trendar2","trendar1cpt","trendar2cpt")
  out[[2]]=mean.fit
  out[[3]]=meancpt.fit
  out[[4]]=meanar.fit
  out[[5]]=meanar2.fit
  out[[6]]=meanarcpt.fit
  out[[7]]=meanar2cpt.fit
  out[[8]]=trend.fit
  out[[9]]=trendcpt.fit
  out[[10]]=trendar.fit
  out[[11]]=trendar2.fit
  out[[12]]=trendarcpt.fit
  out[[13]]=trendar2cpt.fit
  names(out)=c("summary","mean","meancpt","meanar1","meanar2","meanar1cpt","meanar2cpt","trend","trendcpt","trendar1","trendar2","trendar1cpt","trendar2cpt")
  if(verbose==TRUE){close(pb)}
  class(out)="envcpt"
  return(out)
}