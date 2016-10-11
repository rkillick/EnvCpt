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
  colors=rainbow(8)

  extra.args=list(...)
  
  if(any(type=="fit")){
    # data [0,1]
    if(!("mar"%in%names(extra.args))){
      par(mar=c(2,8,1,1))
    }
    else{par(mar=extra.args$mar)}
    
    data=data.set(x$meancpt)
    data=(data-min(data))/diff(range(data)) # scales to [0,1]
    
    plot(1:length(data),type='n',yaxt='n',ylim=c(0,9),ylab="",...)
    axis(2,at=seq(0.5,8.5,by=1),labels=c("Data","Mean","Mean + AR(1)","Trend","Trend + AR(1)","Mean cpt","Mean cpt + AR(1)","Trend cpt","Trend cpt + AR(1)"),col=c("black",colors[c(8,7,5,6,3,4,2,1)]),las=1)
    
    # get the data and fitted values for each of the 8 models and scale to be within 0-1
    lines(data,col='black',lwd=lwd,...)
    
    # mean [1,2]
    mean=1.5 # always flat and we don't put a scale on it
    segments(1,mean,length(data),mean,col=colors[8],lwd=lwd,...)
    
    # meanar [2,3]
    meanar= fitted.values(x$meanar)
    meanar=(meanar-min(meanar))/diff(range(meanar))
    lines(meanar+2,col=colors[7],lwd=lwd,...)
    
    # trend [3,4]
    trend = x$trend$fitted.values
    trend=(trend-min(trend))/diff(range(trend))
    lines(trend+3,col=colors[5],lwd=lwd,...)
    
    # trendar [4,5]
    trendar = x$trendar$fitted.values
    trendar=(trendar-min(trendar))/diff(range(trendar))
    lines(trendar+4,col=colors[6],lwd=lwd,...)
    
    # meancpt [5,6]
    means=param.est(x$meancpt)$mean
    cpts=c(0,x$meancpt@cpts)
    meancpt=rep(means,diff(cpts))
    meancpt=(meancpt-min(meancpt))/diff(range(meancpt))
    lines(meancpt+5,col=colors[3],lwd=lwd,...)
    
    # meanarcpt [6,7]
    meanarcpt=NULL
    cpts=c(0,x$meanarcpt@cpts)
    betas=param.est(x$meanarcpt)$beta
    for(i in 1:nseg(x$meanarcpt)){
      meanarcpt=c(meanarcpt,betas[i,]%*%t(data.set(x$meanarcpt)[(cpts[i]+1):cpts[i+1],-1]))
    }
    meanarcpt=(meanarcpt-min(meanarcpt))/diff(range(meanarcpt))
    lines(meanarcpt+6,col=colors[4],lwd=lwd,...)
    
    # trendcpt [7,8]
    trendcpt=NULL
    cpts=c(0,x$trendcpt@cpts)
    betas=param.est(x$trendcpt)$beta
    for(i in 1:nseg(x$trendcpt)){
      trendcpt=c(trendcpt,betas[i,]%*%t(data.set(x$trendcpt)[(cpts[i]+1):cpts[i+1],-1]))
    }
    trendcpt=(trendcpt-min(trendcpt))/diff(range(trendcpt))
    lines(trendcpt+7,col=colors[2],lwd=lwd,...)
    
    # trendarcpt [8,9]
    trendarcpt=NULL
    cpts=c(0,x$trendarcpt@cpts)
    betas=param.est(x$trendarcpt)$beta
    for(i in 1:nseg(x$trendarcpt)){
      trendarcpt=c(trendarcpt,betas[i,]%*%t(data.set(x$trendarcpt)[(cpts[i]+1):cpts[i+1],-1]))
    }
    trendarcpt=(trendarcpt-min(trendarcpt))/diff(range(trendarcpt))
    lines(trendarcpt+8,col=colors[1],lwd=lwd,...)
  }
  
  else if(any(type=="aic")){
    if(!("mar"%in%names(extra.args))){
      par(mar=c(4,8,1,1))
    }
    else{par(mar=extra.args$mar)}
    
    aic=AIC(x)
    col=rep("white",8)
    col[which.min(aic)]=colors[c(8,3,7,4,5,2,6,1)][which.min(aic)]
    
    xlim=c(min(aic)-0.1*diff(range(aic)),max(aic)+0.1*diff(range(aic)))
    if(("xlim"%in%names(extra.args))&("main"%in%names(extra.args))){
      barplot(aic[c(1,3,5,7,2,4,6,8)],yaxt='n',ylab="",horiz=TRUE,ylim=c(-1,9.5),border=c(colors[c(8,7,5,6,3,4,2,1)]),col=col[c(1,3,5,7,2,4,6,8)],xlab='<--  More likely      Less likely  -->',xpd=FALSE,...)
    }
    else if("xlim"%in%names(extra.args)){
      barplot(aic[c(1,3,5,7,2,4,6,8)],yaxt='n',ylab="",horiz=TRUE,ylim=c(-1,9.5),border=c(colors[c(8,7,5,6,3,4,2,1)]),col=col[c(1,3,5,7,2,4,6,8)],xlab='<--  More likely      Less likely  -->',xpd=FALSE,main="AIC",...)
    }
    else if("main"%in%names(extra.args)){
      barplot(aic[c(1,3,5,7,2,4,6,8)],yaxt='n',ylab="",horiz=TRUE,ylim=c(-1,9.5),border=c(colors[c(8,7,5,6,3,4,2,1)]),col=col[c(1,3,5,7,2,4,6,8)],xlab='<--  More likely      Less likely  -->',xpd=FALSE,xlim=xlim,...)
    }
    else{
      barplot(aic[c(1,3,5,7,2,4,6,8)],yaxt='n',ylab="",horiz=TRUE,ylim=c(-1,9.5),border=c(colors[c(8,7,5,6,3,4,2,1)]),col=col[c(1,3,5,7,2,4,6,8)],xlab='<--  More likely      Less likely  -->',xpd=FALSE,xlim=xlim,main="AIC",...)
    }
    if("yaxt"%in%names(extra.args)){
      if(extra.args$yaxt!='n'){
        labels=c("Mean","Mean + AR(1)","Trend","Trend + AR(1)","Mean cpt","Mean cpt + AR(1)","Trend cpt","Trend cpt + AR(1)")
        if(ncpts(x$meancpt)==0){labels[5]="Mean cpt*"}
        if(ncpts(x$meanarcpt)==0){labels[6]="Mean cpt + AR(1)*"}
        if(ncpts(x$trendcpt)==0){labels[7]="Trend cpt*"}
        if(ncpts(x$trendarcpt)==0){labels[8]="Trend cpt + AR(1)*"}
        axis(2,at=seq(0.75,9,len=8),labels=labels,las=1)
      }
    }
    else{
      labels=c("Mean","Mean + AR(1)","Trend","Trend + AR(1)","Mean cpt","Mean cpt + AR(1)","Trend cpt","Trend cpt + AR(1)")
      if(ncpts(x$meancpt)==0){labels[5]="Mean cpt*"}
      if(ncpts(x$meanarcpt)==0){labels[6]="Mean cpt + AR(1)*"}
      if(ncpts(x$trendcpt)==0){labels[7]="Trend cpt*"}
      if(ncpts(x$trendarcpt)==0){labels[8]="Trend cpt + AR(1)*"}
      axis(2,at=seq(0.75,9,len=8),labels=labels,las=1)
    }
  }
  else{stop("type supplied can only be 'aic' or 'fit'.")}
}
