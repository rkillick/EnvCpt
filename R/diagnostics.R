AIC.envcpt = function(object,...,k=2){
  if(class(object)!="envcpt"){stop("object must be of class envcpt")}
  if(!is.list(object)){stop("object argument must be a list")}
  if(!is.matrix(object[[1]])){stop("first element in the object list must be a matrix.")}
  if(any(!is.numeric(object[[1]][c(1:2),]),na.rm=TRUE)){stop("First two rows in matrix in first element of object list must be numeric")}
  ##Calculate AIC for a secure list object
  return(object[[1]][1,] + k*object[[1]][2,])
}



plot.envcpt=function(x,type=c('fit','aic'),lwd=3,...){
  if(class(x)!="envcpt"){stop("x must be an object with class envcpt")}
  colors=rainbow(12)

  extra.args=list(...)
  
  if(any(type=="fit")){
    # data [0,1]
    if(!("mar"%in%names(extra.args))){
      par(mar=c(2,8,1,1))
    }
    else{par(mar=extra.args$mar)}
    
    data=data.set(x$meancpt)
    data=(data-min(data))/diff(range(data)) # scales to [0,1]
    
    plot(1:length(data),type='n',yaxt='n',ylim=c(0,13),ylab="",...)
    axis(2,at=seq(0.5,12.5,by=1),labels=c("Data","Mean","Mean + AR(1)","Mean + AR(2)","Trend","Trend + AR(1)","Trend + AR(2)","Mean cpt","Mean cpt + AR(1)","Mean cpt + AR(2)","Trend cpt","Trend cpt + AR(1)","Trend cpt + AR(2)"),col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
    
    # get the data and fitted values for each of the 8 models and scale to be within 0-1
    lines(data,col='black',lwd=lwd,...)
    
    # mean [1,2]
    mean=1.5 # always flat and we don't put a scale on it
    segments(1,mean,length(data),mean,col=colors[12],lwd=lwd,...)
    
    # meanar1 [2,3]
    meanar1= fitted.values(x$meanar1)
    meanar1=(meanar1-min(meanar1))/diff(range(meanar1))
    lines(meanar1+2,col=colors[11],lwd=lwd,...)
    
    # meanar2 [3,4]
    meanar2= fitted.values(x$meanar2)
    meanar2=(meanar2-min(meanar2))/diff(range(meanar2))
    lines(meanar2+3,col=colors[10],lwd=lwd,...)
    
    # trend [4,5]
    trend = x$trend$fitted.values
    trend=(trend-min(trend))/diff(range(trend))
    lines(trend+4,col=colors[9],lwd=lwd,...)
    
    # trendar1 [5,6]
    trendar1 = x$trendar1$fitted.values
    trendar1=(trendar1-min(trendar1))/diff(range(trendar1))
    lines(trendar1+5,col=colors[8],lwd=lwd,...)
    
    # trendar2 [6,7]
    trendar2 = x$trendar2$fitted.values
    trendar2=(trendar2-min(trendar2))/diff(range(trendar2))
    lines(trendar2+6,col=colors[7],lwd=lwd,...)

    # meancpt [7,8]
    means=param.est(x$meancpt)$mean
    cpts=c(0,x$meancpt@cpts)
    meancpt=rep(means,diff(cpts))
    meancpt=(meancpt-min(meancpt))/diff(range(meancpt))
    lines(meancpt+7,col=colors[6],lwd=lwd,...)
    
    # meanar1cpt [8,9]
    meanar1cpt=NULL
    cpts=c(0,x$meanar1cpt@cpts)
    betas=param.est(x$meanar1cpt)$beta
    for(i in 1:nseg(x$meanar1cpt)){
      meanar1cpt=c(meanar1cpt,betas[i,]%*%t(data.set(x$meanar1cpt)[(cpts[i]+1):cpts[i+1],-1]))
    }
    meanar1cpt=(meanar1cpt-min(meanar1cpt))/diff(range(meanar1cpt))
    lines(meanar1cpt+8,col=colors[5],lwd=lwd,...)
    
    # meanar2cpt [9,10]
    meanar2cpt=NULL
    cpts=c(0,x$meanar2cpt@cpts)
    betas=param.est(x$meanar2cpt)$beta
    for(i in 1:nseg(x$meanar2cpt)){
      meanar2cpt=c(meanar2cpt,betas[i,]%*%t(data.set(x$meanar2cpt)[(cpts[i]+1):cpts[i+1],-1]))
    }
    meanar2cpt=(meanar2cpt-min(meanar2cpt))/diff(range(meanar2cpt))
    lines(meanar2cpt+9,col=colors[4],lwd=lwd,...)

    # trendcpt [10,11]
    trendcpt=NULL
    cpts=c(0,x$trendcpt@cpts)
    betas=param.est(x$trendcpt)$beta
    for(i in 1:nseg(x$trendcpt)){
      trendcpt=c(trendcpt,betas[i,]%*%t(data.set(x$trendcpt)[(cpts[i]+1):cpts[i+1],-1]))
    }
    trendcpt=(trendcpt-min(trendcpt))/diff(range(trendcpt))
    lines(trendcpt+10,col=colors[3],lwd=lwd,...)
    
    # trendar1cpt [11,12]
    trendar1cpt=NULL
    cpts=c(0,x$trendar1cpt@cpts)
    betas=param.est(x$trendar1cpt)$beta
    for(i in 1:nseg(x$trendar1cpt)){
      trendar1cpt=c(trendar1cpt,betas[i,]%*%t(data.set(x$trendar1cpt)[(cpts[i]+1):cpts[i+1],-1]))
    }
    trendar1cpt=(trendar1cpt-min(trendar1cpt))/diff(range(trendar1cpt))
    lines(trendar1cpt+11,col=colors[2],lwd=lwd,...)

    # trendar2cpt [12,13]
    trendar2cpt=NULL
    cpts=c(0,x$trendar2cpt@cpts)
    betas=param.est(x$trendar2cpt)$beta
    for(i in 1:nseg(x$trendar2cpt)){
      trendar2cpt=c(trendar2cpt,betas[i,]%*%t(data.set(x$trendar2cpt)[(cpts[i]+1):cpts[i+1],-1]))
    }
    trendar2cpt=(trendar2cpt-min(trendar2cpt))/diff(range(trendar2cpt))
    lines(trendar2cpt+12,col=colors[1],lwd=lwd,...)
  }
  
  else if(any(type=="aic")){
    if(!("mar"%in%names(extra.args))){
      par(mar=c(4,8,1,1))
    }
    else{par(mar=extra.args$mar)}
    
    aic=AIC(x)
    col=rep("white",12)
    col[which.min(aic)]=rev(colors)[which.min(aic)]
    
    xlim=c(min(aic)-0.1*diff(range(aic)),max(aic)+0.1*diff(range(aic)))
    if(("xlim"%in%names(extra.args))&("main"%in%names(extra.args))){
      barplot(aic[c(1,3,4,7,9,10,2,5,6,8,11,12)],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,14),border=rev(colors)[c(1,3,4,7,9,10,2,5,6,8,11,12)],col=col[c(1,3,4,7,9,10,2,5,6,8,11,12)],xlab='<--  More likely      Less likely  -->',xpd=FALSE,...)
    }
    else if("xlim"%in%names(extra.args)){
      barplot(aic[c(1,3,4,7,9,10,2,5,6,8,11,12)],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,14),border=rev(colors)[c(1,3,4,7,9,10,2,5,6,8,11,12)],col=col[c(1,3,4,7,9,10,2,5,6,8,11,12)],xlab='<--  More likely      Less likely  -->',xpd=FALSE,main="AIC",...)
    }
    else if("main"%in%names(extra.args)){
      barplot(aic[c(1,3,4,7,9,10,2,5,6,8,11,12)],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,14),border=rev(colors)[c(1,3,4,7,9,10,2,5,6,8,11,12)],col=col[c(1,3,4,7,9,10,2,5,6,8,11,12)],xlab='<--  More likely      Less likely  -->',xpd=FALSE,xlim=xlim,...)
    }
    else{
      barplot(aic[c(1,3,4,7,9,10,2,5,6,8,11,12)],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,14),border=rev(colors)[c(1,3,4,7,9,10,2,5,6,8,11,12)],col=col[c(1,3,4,7,9,10,2,5,6,8,11,12)],xlab='<--  More likely      Less likely  -->',xpd=FALSE,xlim=xlim,main="AIC",...)
    }
    if("yaxt"%in%names(extra.args)){
      if(extra.args$yaxt!='n'){
        labels=c("Mean","Mean + AR(1)","Mean + AR(2)","Trend","Trend + AR(1)","Trend + AR(2)","Mean cpt","Mean cpt + AR(1)","Mean cpt + AR(2)","Trend cpt","Trend cpt + AR(1)","Trend cpt + AR(2)")
        if(ncpts(x$meancpt)==0){labels[7]="Mean cpt*"}
        if(ncpts(x$meanar1cpt)==0){labels[8]="Mean cpt + AR(1)*"}
        if(ncpts(x$meanar2cpt)==0){labels[9]="Mean cpt + AR(2)*"}
        if(ncpts(x$trendcpt)==0){labels[10]="Trend cpt*"}
        if(ncpts(x$trendar1cpt)==0){labels[11]="Trend cpt + AR(1)*"}
        if(ncpts(x$trendar2cpt)==0){labels[12]="Trend cpt + AR(2)*"}
        axis(2,at=seq(0.6,14,len=12),labels=labels,las=1)
      }
    }
    else{
      labels=c("Mean","Mean + AR(1)","Mean + AR(2)","Trend","Trend + AR(1)","Trend + AR(2)","Mean cpt","Mean cpt + AR(1)","Mean cpt + AR(2)","Trend cpt","Trend cpt + AR(1)","Trend cpt + AR(2)")
      if(ncpts(x$meancpt)==0){labels[7]="Mean cpt*"}
      if(ncpts(x$meanar1cpt)==0){labels[8]="Mean cpt + AR(1)*"}
      if(ncpts(x$meanar2cpt)==0){labels[9]="Mean cpt + AR(2)*"}
      if(ncpts(x$trendcpt)==0){labels[10]="Trend cpt*"}
      if(ncpts(x$trendar1cpt)==0){labels[11]="Trend cpt + AR(1)*"}
      if(ncpts(x$trendar2cpt)==0){labels[12]="Trend cpt + AR(2)*"}
      axis(2,at=seq(0.6,14,len=12),labels=labels,las=1)
    }
  }
  else{stop("type supplied can only be 'aic' or 'fit'.")}
}
