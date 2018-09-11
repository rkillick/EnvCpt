AIC.envcpt = function(object,...,k=2){
  if(class(object)!="envcpt"){stop("object must be of class envcpt")}
  if(!is.list(object)){stop("object argument must be a list")}
  if(!is.matrix(object[[1]])){stop("first element in the object list must be a matrix.")}
  if(any(!is.numeric(object[[1]][c(1:2),]),na.rm=TRUE)){stop("First two rows in matrix in first element of object list must be numeric")}
  ##Calculate AIC for a secure list object
  return(object[[1]][1,] + k*object[[1]][2,])
}

AICweights=function(object){
  UseMethod("AICweights")
}
AICweights.default=function(object){
  return("No default method created for S3 class AICweights.")
}
AICweights.envcpt=function(object){
  aic=AIC(object)
  deltaAIC=aic-min(aic,na.rm=T)
  expdeltaAIC=exp(-0.5*deltaAIC)
  weights=expdeltaAIC/(sum(expdeltaAIC,na.rm=T))
  return(weights)
}


BIC.envcpt = function(object,...){
  if(class(object)!="envcpt"){stop("object must be of class envcpt")}
  if(!is.list(object)){stop("object argument must be a list")}
  if(!is.matrix(object[[1]])){stop("first element in the object list must be a matrix.")}
  if(any(!is.numeric(object[[1]][c(1:2),]),na.rm=TRUE)){stop("First two rows in matrix in first element of object list must be numeric")}
  ##Calculate BIC for a secure list object
  z=list(...)
  if(!any(names(z)=="n")){
    models=which(!is.na(object[[1]][1,]))
    if(any(models==2)){n=length(data.set(object$meancpt))}
    else if(any(models==8)){n=length(data.set(object$trendcpt)[,1])}
    else if(any(models==5)){n=length(data.set(object$meanar1cpt)[,1])}
    else if(any(models==11)){n=length(data.set(object$trendar1cpt)[,1])}
    else if(any(models=6)){n=length(data.set(object$meanar2cpt)[,1])}
    else if(any(models==12)){n=length(data.set(object$trendar2cpt)[,1])}
    else{stop("Length of data is unknown, please pass an argument n to the function call.")}
  }
  return(object[[1]][1,] + log(n)*object[[1]][2,])
}



plot.envcpt=function(x,type=c('fit','bic','aic'),lwd=3,...,data=NA){
  if(class(x)!="envcpt"){stop("x must be an object with class envcpt")}
  colors=rainbow(12)

  extra.args=list(...)
  
  models=which(!is.na(x[[1]][1,]))
  reorder=c(1,3,4,7,9,10,2,5,6,8,11,12)
  reorder.index=c(1,7,2,3,8,9,4,10,5,6,11,12)[models]
  plot.order=sort(reorder.index)
  
  if(any(type=="fit")){
    # data [0,1]
    if(!("mar"%in%names(extra.args))){
      par(mar=c(2,8,1,1))
    }
    else{par(mar=extra.args$mar)}
    
    if(any(models==2)){data=data.set(x$meancpt)}
    else if(any(models==8)){data=data.set(x$trendcpt)[,1]}
    else if(any(models==5)){data=data.set(x$meanar1cpt)[,1];data=c(data[1],data)} # repeat first data point
    else if(any(models==11)){data=data.set(x$trendar1cpt)[,1];data=c(data[1],data)}
    else if(any(models=6)){data=data.set(x$meanar2cpt)[,1];data=c(data[1:2],data)}
    else if(any(models==12)){data=data.set(x$trendar2cpt)[,1];data=c(data[1:2],data)}
    else if(is.na(data)){stop("No changepoint models fit so no dataset is identifiable.  Please provide the data argument so the fits can be plotted.")}
    
    data=(data-min(data))/diff(range(data)) # scales to [0,1]
    
    plot(1:length(data),type='n',yaxt='n',ylim=c(0,length(models)+1),ylab="",...)
    # axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c("Data","Mean","Mean + AR(1)","Mean + AR(2)","Trend","Trend + AR(1)","Trend + AR(2)","Mean cpt","Mean cpt + AR(1)","Mean cpt + AR(2)","Trend cpt","Trend cpt + AR(1)","Trend cpt + AR(2)")[c(1,models+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
    
    # get the data and fitted values for each of the models and scale to be within 0-1
    lines(data,col='black',lwd=lwd,...)
    axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c("Data",rep("",12))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
    offset=1
    
    # mean [1,2]
    if(any(models==1)){
      axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c("","Mean",rep("",11))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
      mean=1.5 # always flat and we don't put a scale on it
      segments(1,mean,length(data),mean,col=colors[12],lwd=lwd,...)
      offset=offset+1
    }
    
    # meanar1 [2,3]
    if(any(models==3)){
      axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",2),"Mean + AR(1)",rep("",10))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
      meanar1= data-x$meanar1$residuals
      meanar1=(meanar1-min(meanar1))/diff(range(meanar1))
      lines(meanar1+offset,col=colors[11],lwd=lwd,...)
      offset=offset+1
    }
    
    # meanar2 [3,4]
    if(any(models==4)){
      axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",3),"Mean + AR(2)",rep("",9))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
      meanar2= data-x$meanar2$residuals
      meanar2=(meanar2-min(meanar2))/diff(range(meanar2))
      lines(meanar2+offset,col=colors[10],lwd=lwd,...)
      offset=offset+1
    }
    
    # trend [4,5]
    if(any(models==7)){
      axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",4),"Trend",rep("",8))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
      trend = x$trend$fitted.values
      trend=(trend-min(trend))/diff(range(trend))
      lines(trend+offset,col=colors[9],lwd=lwd,...)
      offset=offset+1
    }
    
    # trendar1 [5,6]
    if(any(models==9)){
      axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",5),"Trend + AR(1)",rep("",7))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
      trendar1 = x$trendar1$fitted.values
      trendar1=(trendar1-min(trendar1))/diff(range(trendar1))
      lines(trendar1+offset,col=colors[8],lwd=lwd,...)
      offset=offset+1
    }
    
    # trendar2 [6,7]
    if(any(models==10)){
      axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",6),"Trend + AR(2)",rep("",6))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
      trendar2 = x$trendar2$fitted.values
      trendar2=(trendar2-min(trendar2))/diff(range(trendar2))
      lines(trendar2+offset,col=colors[7],lwd=lwd,...)
      offset=offset+1
    }

    # meancpt [7,8]
    if(any(models==2)){
      means=param.est(x$meancpt)$mean
      cpts=c(0,x$meancpt@cpts)
      if(length(cpts)>2){ # there are cpts
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",7),"Mean cpt",rep("",5))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
        meancpt=rep(means,diff(cpts))
        meancpt=(meancpt-min(meancpt))/diff(range(meancpt))
        segments(x0=cpts[-c(1,length(cpts))],y0=offset,y1=offset+1,col=colors[6],...)
        lines(meancpt+offset,col=colors[6],lwd=lwd,...)
      }
      else{
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",7),"Mean cpt*",rep("",5))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
        meancpt=0.5 # flat mean, value doesn't matter
        segments(1,meancpt+offset,length(data),meancpt+offset,col=colors[6],lwd=lwd,...)
      }
      offset=offset+1
    }
    
    # meanar1cpt [8,9]
    if(any(models==5)){
      meanar1cpt=NULL
      cpts=c(0,x$meanar1cpt@cpts)
      betas=param.est(x$meanar1cpt)$beta
      for(i in 1:nseg(x$meanar1cpt)){
        meanar1cpt=c(meanar1cpt,betas[i,]%*%t(data.set(x$meanar1cpt)[(cpts[i]+1):cpts[i+1],-1]))
      }
      meanar1cpt=(meanar1cpt-min(meanar1cpt))/diff(range(meanar1cpt))
      if(length(cpts)>2){ # there are cpts
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",8),"Mean cpt + AR(1)",rep("",4))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
        segments(x0=cpts[-c(1,length(cpts))],y0=offset,y1=offset+1,col=colors[5],...)
      }
      else{
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",8),"Mean cpt + AR(1)*",rep("",4))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
      }
      lines(meanar1cpt+offset,col=colors[5],lwd=lwd,...)
      offset=offset+1
    }
    
    # meanar2cpt [9,10]
    if(any(models==6)){
      meanar2cpt=NULL
      cpts=c(0,x$meanar2cpt@cpts)
      betas=param.est(x$meanar2cpt)$beta
      for(i in 1:nseg(x$meanar2cpt)){
        meanar2cpt=c(meanar2cpt,betas[i,]%*%t(data.set(x$meanar2cpt)[(cpts[i]+1):cpts[i+1],-1]))
      }
      meanar2cpt=(meanar2cpt-min(meanar2cpt))/diff(range(meanar2cpt))
      if(length(cpts)>2){ # there are cpts
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",9),"Mean cpt + AR(2)",rep("",3))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
        segments(x0=cpts[-c(1,length(cpts))],y0=offset,y1=offset+1,col=colors[4],...)
      }
      else{
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",9),"Mean cpt + AR(2)*",rep("",3))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
      }
      lines(meanar2cpt+offset,col=colors[4],lwd=lwd,...)
      offset=offset+1
    }

    # trendcpt [10,11]
    if(any(models==8)){
      trendcpt=NULL
      cpts=c(0,x$trendcpt@cpts)
      betas=param.est(x$trendcpt)$beta
      for(i in 1:nseg(x$trendcpt)){
        trendcpt=c(trendcpt,betas[i,]%*%t(data.set(x$trendcpt)[(cpts[i]+1):cpts[i+1],-1]))
      }
      trendcpt=(trendcpt-min(trendcpt))/diff(range(trendcpt))
      if(length(cpts)>2){ # there are cpts
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",10),"Trend cpt",rep("",2))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
        segments(x0=cpts[-c(1,length(cpts))],y0=offset,y1=offset+1,col=colors[3],...)
      }
      else{
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",10),"Trend cpt*",rep("",2))[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
      }
      lines(trendcpt+offset,col=colors[3],lwd=lwd,...)
      offset=offset+1
    }
    
    # trendar1cpt [11,12]
    if(any(models==11)){
      trendar1cpt=NULL
      cpts=c(0,x$trendar1cpt@cpts)
      betas=param.est(x$trendar1cpt)$beta
      for(i in 1:nseg(x$trendar1cpt)){
        trendar1cpt=c(trendar1cpt,betas[i,]%*%t(data.set(x$trendar1cpt)[(cpts[i]+1):cpts[i+1],-1]))
      }
      trendar1cpt=(trendar1cpt-min(trendar1cpt))/diff(range(trendar1cpt))
      if(length(cpts)>2){ # there are cpts
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",11),"Trend cpt + AR(1)","")[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
        segments(x0=cpts[-c(1,length(cpts))],y0=offset,y1=offset+1,col=colors[2],...)
      }
      else{
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",11),"Trend cpt + AR(1)*","")[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
      }
      lines(trendar1cpt+offset,col=colors[2],lwd=lwd,...)
      offset=offset+1
    }

    # trendar2cpt [12,13]
    if(any(models==12)){
      trendar2cpt=NULL
      cpts=c(0,x$trendar2cpt@cpts)
      betas=param.est(x$trendar2cpt)$beta
      for(i in 1:nseg(x$trendar2cpt)){
        trendar2cpt=c(trendar2cpt,betas[i,]%*%t(data.set(x$trendar2cpt)[(cpts[i]+1):cpts[i+1],-1]))
      }
      trendar2cpt=(trendar2cpt-min(trendar2cpt))/diff(range(trendar2cpt))
      if(length(cpts)>2){ # there are cpts
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",12),"Trend cpt + AR(2)")[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
        segments(x0=cpts[-c(1,length(cpts))],y0=offset,y1=offset+1,col=colors[1],...)
      }
      else{
        axis(2,at=seq(0.5,length(models)+0.5,by=1),labels=c(rep("",12),"Trend cpt + AR(2)*")[c(1,plot.order+1)],col=c("black",colors[c(8,7,5,6,3,4,2,1,1,1,1,1)]),las=1)
      }
      lines(trendar2cpt+offset,col=colors[1],lwd=lwd,...)
      offset=offset+1
    }
  }
  
  else if(any(type=="bic")){
    if(!("mar"%in%names(extra.args))){
      par(mar=c(4,8,1,1))
    }
    else{par(mar=extra.args$mar)}
    
    bic=BIC(x)
    col=rep("white",12)
    col[which(reorder==which.min(bic))]=rev(colors)[which(reorder==which.min(bic))]
    col=col[plot.order]
    colors=rev(colors)[plot.order]
    
    xlim=c(min(bic,na.rm=T)-0.1*diff(range(bic,na.rm=T)),max(bic,na.rm=T)+0.1*diff(range(bic,na.rm=T)))
    if(("xlim"%in%names(extra.args))&("main"%in%names(extra.args))){
      barplot(bic[reorder][plot.order],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,length(models)*7/6),border=colors,col=col,xlab='<--  More likely      Less likely  -->',xpd=FALSE,...)
    }
    else if("xlim"%in%names(extra.args)){
      barplot(bic[reorder][plot.order],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,length(models)*7/6),border=colors,col=col,xlab='<--  More likely      Less likely  -->',xpd=FALSE,main="BIC",...)
    }
    else if("main"%in%names(extra.args)){
      barplot(bic[reorder][plot.order],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,length(models)*7/6),border=colors,col=col,xlab='<--  More likely      Less likely  -->',xpd=FALSE,xlim=xlim,...)
    }
    else{
      barplot(bic[reorder][plot.order],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,length(models)*7/6),border=colors,col=col,xlab='<--  More likely      Less likely  -->',xpd=FALSE,xlim=xlim,main="BIC",...)
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
        axis(2,at=seq(25/36,by=43/36,len=length(models)),labels=labels[plot.order],las=1)
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
      axis(2,at=seq(25/36,by=43/36,len=length(models)),labels=labels[plot.order],las=1)
    }
  }
  
  else if(any(type=="aic")){
    if(!("mar"%in%names(extra.args))){
      par(mar=c(4,8,1,1))
    }
    else{par(mar=extra.args$mar)}
    
    aic=AIC(x)
    col=rep("white",12)
    col[which(reorder==which.min(aic))]=rev(colors)[which(reorder==which.min(aic))]
    col=col[plot.order]
    colors=rev(colors)[plot.order]
    
    xlim=c(min(aic,na.rm=T)-0.1*diff(range(aic,na.rm=T)),max(aic,na.rm=T)+0.1*diff(range(aic,na.rm=T)))
    if(("xlim"%in%names(extra.args))&("main"%in%names(extra.args))){
      barplot(aic[reorder][plot.order],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,length(models)*7/6),border=colors,col=col,xlab='<--  More likely      Less likely  -->',xpd=FALSE,...)
    }
    else if("xlim"%in%names(extra.args)){
      barplot(aic[reorder][plot.order],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,length(models)*7/6),border=colors,col=col,xlab='<--  More likely      Less likely  -->',xpd=FALSE,main="AIC",...)
    }
    else if("main"%in%names(extra.args)){
      barplot(aic[reorder][plot.order],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,length(models)*7/6),border=colors,col=col,xlab='<--  More likely      Less likely  -->',xpd=FALSE,xlim=xlim,...)
    }
    else{
      barplot(aic[reorder][plot.order],yaxt='n',ylab="",horiz=TRUE,ylim=c(0,length(models)*7/6),border=colors,col=col,xlab='<--  More likely      Less likely  -->',xpd=FALSE,xlim=xlim,main="AIC",...)
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
        axis(2,at=seq(25/36,by=43/36,len=length(models)),labels=labels[plot.order],las=1)
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
      axis(2,at=seq(25/36,by=43/36,len=length(models)),labels=labels[plot.order],las=1)
    }
  }
  else{stop("type supplied can only be 'aic', 'bic' or 'fit'.")}
}
