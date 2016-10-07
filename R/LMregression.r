# functions to do the change in regression
single.reg.norm.calc <-
function(data,extrainf=TRUE){
  singledim=function(data,extrainf=TRUE){
    n=nrow(data)
    p=ncol(data)-1
    if(p==0){stop("Dimension of data is 1, no regressors found")}
    formula='-1+data[,2]'
  	if(p>1){
  		for(i in 2:p){
        formula=paste(formula,'+data[,',i+1,']',sep='')
      }
    }
    null=-2*logLik(eval(parse(text=paste('lm(data[,1]~',formula,')',sep=''))))
    tmp=NULL
    if((n-p)<p){stop("Too many regressors / not enough data to test for a change")}
    for(taustar in p:(n-p)){
      formula1='-1+data[1:taustar,2]'
      formulaT='-1+data[-c(1:taustar),2]'
      if(p>1){
    		for(i in 2:p){
          formula1=paste(formula1,'+data[1:taustar,',i+1,']',sep='')
          formulaT=paste(formulaT,'+data[-c(1:taustar),',i+1,']',sep='')
        }
      }
      tmp[taustar-1]=-2*logLik(eval(parse(text=paste('lm(data[1:taustar,1]~',formula1,')',sep=''))))-2*logLik(eval(parse(text=paste('lm(data[-c(1:taustar),1]~',formulaT,')',sep=''))))
    }
    tau=which(tmp==min(tmp))[1]+1
    taulike=tmp[tmp==min(tmp)]
    if(extrainf==TRUE){
      return(c(tau,null,taulike))
    }
    else{
      return(tau)
    }
  }


  if(length(dim(data))==2){
    # single data set
    cpt=singledim(data,extrainf)
    return(cpt)
  }
  else{
    rep=dim(data)[1]
    n=dim(data)[2]
    cpt=NULL
    if(extrainf==FALSE){
      for(i in 1:rep){
        cpt[i]=singledim(data[i,,],extrainf)
      }
    }
    else{
      cpt=matrix(0,ncol=3,nrow=rep)
      for(i in 1:rep){
        cpt[i,]=singledim(data[i,,],extrainf)
      }
    }
    return(cpt)
  }
}


single.reg.norm<-function(data,penalty="SIC",value=0,class=TRUE,param.estimates=TRUE){
	if(length(dim(data))==2){
		n=nrow(data)
		d=ncol(data)
		if(penalty=="Asymptotic"){
			alpha=value
      top=-(log(log((1 - alpha + exp(-2*exp(2*(log(log(n)))+(d/2)*(log(log(log(n))))- log(gamma(d/2)))))^(-1/2))))  +  2*(log(log(n)))+(d/2)*(log(log(log(n))))- log(gamma(d/2))
			bottom=(2*log(log(n)))^(1/2)
			value=(top/bottom)^2
		}
		tmp=single.reg.norm.calc(data,extrainf=TRUE)
		ans=decision(tmp[1],tmp[2],tmp[3],penalty,n,diffparam=d,value)
		if(class==TRUE){
			out=new("cpt.reg")
			data.set(out)=data;method(out)="AMOC";distribution(out)="Normal";pen.type(out)=penalty;pen.value(out)=ans$pen;ncpts.max(out)=1
			if(ans$cpt != n){cpts(out)=c(ans$cpt,n)}
			else{cpts(out)=ans$cpt}
			if(param.estimates==TRUE){
				out=param(out)
			}

			return(out)
		}
		else{ return(ans$cpt)}
	}
	else{
		n=dim(data)[2]
		d=dim(data)[3]
		if(penalty=="Asymptotic"){
			alpha=value
			top=-(log(log((1 - alpha + exp(-2*exp(2*(log(log(n)))+(d/2)*(log(log(log(n))))- log(gamma(d/2)))))^(-1/2))))  +  2*(log(log(n)))+(d/2)*(log(log(log(n))))- log(gamma(d/2))
			bottom=(2*log(log(n)))^(1/2)
			value=(top/bottom)^2
		}
		tmp=single.reg.norm.calc(data,extrainf=TRUE)
		ans=decision(tmp[,1],tmp[,2],tmp[,3],penalty,n,diffparam=d,value)
		if(class==TRUE){
			rep=nrow(data)
			out=list()
			for(i in 1:rep){
				out[[i]]=new("cpt.reg")
				data.set(out[[i]])=data[i,,];method(out[[i]])="AMOC";distribution(out[[i]])="Normal"; pen.type(out[[i]])=penalty;pen.value(out[[i]])=ans$pen;ncpts.max(out[[i]])=1
				if(ans$cpt[i] != n){cpts(out[[i]])=c(ans$cpt[i],n)}
				else{cpts(out[[i]])=ans$cpt[i]}
				if(param.estimates==TRUE){
					out[[i]]=param(out[[i]])
				}
			}
			return(out)
		}
		else{ return(ans$cpt)}
	}
}

decision<-function(tau,null,alt,penalty="SIC",n=0,diffparam=1,value=0){
	if((length(tau)!=length(null))||(length(tau)!=length(alt))){
		stop("Lengths of tau, null and alt do not match")
	}
  if((penalty=="SIC0") || (penalty=="BIC0")){
    value=diffparam*log(n)
  }
  else if((penalty=="SIC") || (penalty=="BIC")){
    value=(diffparam+1)*log(n)
  }
  else if(penalty=="MBIC"){
    value=(diffparam+2)*log(n)
  }
  else if(penalty=="AIC0"){
    value=2*diffparam
  }
  else if(penalty=="AIC"){
    value=2*(diffparam+1)
  }
  else if(penalty=="Hannan-Quinn0"){
    value=2*diffparam*log(log(n))
  }
  else if(penalty=="Hannan-Quinn"){
    value=2*(diffparam+1)*log(log(n))
  }
	else if(penalty=="None"){
		value=0
	}
	else if((penalty!="Manual")&&(penalty!="Asymptotic")){
		stop('Unknown Penalty')
	}
	if((penalty=="Manual")&&(is.numeric(value)==FALSE)){
		value=try(eval(parse(text=paste(value))),silent=TRUE)
		if(class(value)=='try-error'){
			stop('Your manual penalty cannot be evaluated')
		}
	}
	single.decision=function(tau,null,alt,n=0,diffparam=1,value=0){
		teststat=null-alt
		if(teststat>=value){return(tau)}
		else{return(n)}
	}
	if(length(tau)==1){
		out=single.decision(tau,null,alt,n,diffparam,value)
		names(out)="cpt"
		return(list(cpt=out,pen=value))
	}
	else{
		rep=length(tau)
		out=NULL
		for(i in 1:rep){
			out[i]=single.decision(tau[i],null[i],alt[i],n,diffparam,value)
		}
		names(out)=rep("cpt",rep)
		return(list(cpt=out,pen=value))
	}
}

PELT.reg.norm=function(data,pen=0,minseglen=3){

  mll.reg=function(x,p){
    x=matrix(x,ncol=p+1)
    formula='-1+x[,2]'
  	if(p>1){
  		for(i in 2:p){
        formula=paste(formula,'+x[,',i+1,']',sep='')
      }
    }
    like=-2*logLik(eval(parse(text=paste('lm(x[,1]~',formula,')',sep=''))))

		return(like)
  }
  n=nrow(data)
	p=ncol(data)-1
  if(p==0){stop("Dimension of data is 1, no regressors found")}
	if(minseglen<p){warning(paste("minseglen too small, set to:",p));minseglen=p}
  
  lastchangecpts=NULL
  lastchangelike=-pen
#	ncpts=NULL
  if(n<(2*(minseglen+1))){stop("Too many regressors / minseglen too large so not enough data to test for a change")}# minseglen+1 for the variance
	for(i in (minseglen+1):(2*(minseglen+1))){# p+1 for the variance
		lastchangelike[i+1]=mll.reg(data[1:i,],p)
		lastchangecpts[i]=0
	}
#	ncpts[1:p]=0
#	ncpts[(p+1):(2*p+2)]=1

  checklist=c(0,minseglen+1)
  for(tstar in (2*minseglen+3):n){
    tmplike=NULL
		for(i in 1:length(checklist)){
	    tmplike[i]=lastchangelike[checklist[i]+1]+mll.reg(data[(checklist[i]+1):tstar,],p)+pen
		}

    lastchangelike[tstar+1]=min(tmplike,na.rm=TRUE)
    lastchangecpts[tstar]=checklist[tmplike==lastchangelike[tstar+1]][1]
#		ncpts[tstar]=ncpts[lastchangecpts[tstar]]+1

    checklist=c(checklist[tmplike<=lastchangelike[tstar+1]+pen],tstar-minseglen)
  }
  last=n
  fcpt=NULL
  while(last!=0){
    fcpt=c(fcpt,last)
    last=lastchangecpts[last]
  }
  return(cpt=sort(fcpt))
}

multiple.reg.norm=function(data,mul.method="PELT",penalty="SIC",value=0,Q=5,class=TRUE,param.estimates=TRUE,minseglen=3){
	if(!((mul.method=="PELT"))){
		stop("Multiple Method is not recognised")
	}
	if(is.null(dim(data))==TRUE){
		stop("Data must be a matrix with first column being the response and the remaining columns covariates")
	}
	if(length(dim(data))==2){
		# single dataset
		n=nrow(data)
		diffparam=ncol(data)
	}
	else{
		n=dim(data)[2]
		diffparam=dim(data)[3]
	}
  if((penalty=="SIC0") || (penalty=="BIC0")){
    value=diffparam*log(n)
  }
  else if((penalty=="SIC") || (penalty=="BIC")){
    value=(diffparam+1)*log(n)
  }
  else if(penalty=="MBIC"){
    value=(diffparam+2)*log(n)
  }
  else if(penalty=="AIC0"){
    value=2*diffparam
  }
  else if(penalty=="AIC"){
    value=2*(diffparam+1)
  }
  else if(penalty=="Hannan-Quinn0"){
    value=2*diffparam*log(log(n))
  }
  else if(penalty=="Hannan-Quinn"){
    value=2*(diffparam+1)*log(log(n))
  }
  else if(penalty=="None"){
		value=0
	}
	else if((penalty!="Manual")&&(penalty!="Asymptotic")){
		stop('Unknown Penalty')
	}
	if((penalty=="Manual")&&(is.numeric(value)==FALSE)){
		value=try(eval(parse(text=paste(value))),silent=TRUE)
		if(class(value)=='try-error'){
			stop('Your manual penalty cannot be evaluated')
		}
	}

	if(penalty=="Asymptotic"){
		alpha=value
	  top= -(log(log((1-alpha+exp(-2*exp(2*(log(log(n)))+(diffparam/2)*(log(log(log(n))))-log(gamma(diffparam/2)))))^(-1/2))))+2*(log(log(n)))+(diffparam/2)*(log(log(log(n))))-log(gamma(diffparam/2))
	  bottom=(2*log(log(n)))^(1/2)
		value=(top/bottom)^2 + 2*(log(n))
	}
	if(length(dim(data))==2){
		# single dataset
		if(mul.method=="PELT"){
			out=PELT.reg.norm(data,value,minseglen)
			cpts=out
		}
		if(class==TRUE){
			ans=new("cpt.reg")
			data.set(ans)=data;cpttype(ans)="regression";method(ans)=mul.method; distribution(ans)="Normal";pen.type(ans)=penalty;pen.value(ans)=value;cpts(ans)=cpts
			if(mul.method=="PELT"){
				ncpts.max(ans)=Inf
			}
			if(param.estimates==TRUE){
				ans=param(ans)
			}
			return(ans)
		}
		else{ return(out)}
	}
	else{
		rep=dim(data)[1]
		out=list()
		if(class==TRUE){cpts=list()}
		if(mul.method=="PELT"){
			for(i in 1:rep){
				out=c(out,list(PELT.reg.norm(data[i,,],value)))
			}
			cpts=out
		}
		if(class==TRUE){
			ans=list()
			for(i in 1:rep){
				ans[[i]]=new("cpt.reg")
				data.set(ans[[i]])=data[i,,];cpttype(ans[[i]])="regression"; method(ans[[i]])=mul.method;distribution(ans[[i]])="Normal";pen.type(ans[[i]])=penalty;pen.value(ans[[i]])=value;cpts(ans[[i]])=cpts[[i]]
				if(mul.method=="PELT"){
					ncpts.max(ans[[i]])=Inf
				}
				if(param.estimates==TRUE){
					ans[[i]]=param(ans[[i]])
				}
			}
			return(ans)
		}
		else{return(out)}
	}
}

cpt.reg=function(data,penalty="MBIC",pen.value=0,method="AMOC",dist="Normal",class=TRUE,param.estimates=TRUE,minseglen=3){
	if(dist !="Normal"){ stop("Invalid distribution, must be Normal") }
	if(method=="AMOC"){
		return(single.reg.norm(data,penalty,pen.value,class,param.estimates))
	}
  else if(method=="PELT"){
    return(multiple.reg.norm(data,mul.method="PELT",penalty,pen.value,class=TRUE,param.estimates=param.estimates,minseglen=minseglen))
  }
	else{
		stop("Invalid Method, must be AMOC or PELT")
	}
}
