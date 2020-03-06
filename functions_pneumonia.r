read_data = function(file)
{
  x = read.csv(file)
  i=which(!is.na(x$ILI))
  output = list(t=x$Day.of.year[i],s1 = x$ILI[i], s2 = x$Pneumonia[i])
  return(output)
}

parsepar = function(parname,x)
{
  out = list(value=-1,fixed=TRUE)
  i=which(x$Parameter==parname)
  if(length(i)!=1)
  {
    if(parname=='n0')return(list(value=1,fixed=FALSE))
    if(parname=='b')return(list(value=0.11,fixed=FALSE))
    if(parname=='L1')return(list(value=5.5,fixed=TRUE))
    if(parname=='L2')return(list(value=7,fixed=TRUE))
    if(parname=='p1')return(list(value=0.017,fixed=TRUE))
    if(parname=='p2')return(list(value=0.19,fixed=TRUE))
  }
  if(length(i)==1)
  {
    temp = x$Value[i]
    if(temp=='unknown')
    {
      if(parname=='n0')return(list(value=1,fixed=FALSE))
      if(parname=='b')return(list(value=0.11,fixed=FALSE))
      if(parname=='L1')return(list(value=5.5,fixed=FALSE))
      if(parname=='L2')return(list(value=7,fixed=FALSE))
      if(parname=='p1')return(list(value=0.017,fixed=FALSE))
      if(parname=='p2')return(list(value=0.19,fixed=FALSE))
    }
    if(temp!='unknown')
    {
      return(list(value=as.numeric(temp),fixed=TRUE))
    }
  }
}

read_pars = function(file)
{
  x = read.csv(file,as.is=TRUE)
  output = list(n0=c(),b=c(),L1=c(),L2=c(),p1=c(),p2=c(),
                fixed_n0=FALSE,
                fixed_b=FALSE,
                fixed_L1=TRUE,
                fixed_L2=TRUE,
                fixed_p1=TRUE,
                fixed_p2=TRUE)
  temp = parsepar('n0',x);output$n0 = temp$value;output$fixed_n0 = temp$fixed
  temp = parsepar('b',x);output$b = temp$value;output$fixed_b = temp$fixed
  temp = parsepar('L1',x);output$L1 = temp$value;output$fixed_L1 = temp$fixed
  temp = parsepar('L2',x);output$L2 = temp$value;output$fixed_L2 = temp$fixed
  temp = parsepar('p1',x);output$p1 = temp$value;output$fixed_p1 = temp$fixed
  temp = parsepar('p2',x);output$p2 = temp$value;output$fixed_p2 = temp$fixed
  return(output)
}



modelmean = function(pars,TMAX)
{
  tv= 1:TMAX
  mu1 = pars$n0*pars$p1*exp(pars$b*(tv-pars$L1))
  mu2 = pars$n0*pars$p2*exp(pars$b*(tv-pars$L2))
  return(list(mu1=mu1,mu2=mu2))
}

loglikelihood = function(pars,surv=surveillance)
{
  mus = modelmean(pars,length(surv$s1))
  output = sum(dpois(surv$s2,mus$mu2,log=TRUE))
  return(output)
}
casestodate = function(pars,TMAX)
{
  tv= 1:TMAX
  nt = pars$n0*exp(pars$b*tv)
  return(sum(nt))
}
casestoday = function(pars,TMAX)
{
  tv= TMAX
  nt = pars$n0*exp(pars$b*tv)
  return(nt)
}

mh = function(oldpars,newpars,surv=surveillance)
{
  reject=FALSE
  if(newpars$n0<0)reject=TRUE
  if(newpars$b<0)reject=TRUE
  if(newpars$L1<0)reject=TRUE
  if(newpars$L2<0)reject=TRUE
  if(newpars$p1<0)reject=TRUE
  if(newpars$p1>1)reject=TRUE
  if(newpars$p2<0)reject=TRUE
  if(newpars$p2>1)reject=TRUE
  if(!reject)
  {
    newpars$LL = loglikelihood(newpars,surv)
    logaccprob = newpars$LL-oldpars$LL
    lu = -rexp(1)
    if(lu>logaccprob)reject=TRUE
  }
  if(reject)return(oldpars)
  return(newpars)
}

summariser = function(storepars,file)
{
  cat('Estimated cases to date: ',round(mean(storepars$size)),
      ' (95%CI ',round(quantile(storepars$size,0.025)),
      ', ',round(quantile(storepars$size,0.975)),')\n',file=file,sep='')
  cat('Estimated cases today: ',round(mean(storepars$sizetoday)),
      ' (95%CI ',round(quantile(storepars$sizetoday,0.025)),
      ', ',round(quantile(storepars$sizetoday,0.975)),')\n',file=file,append=TRUE,sep='')
  cat('\nWarning! This is only appropriate during the exponential growth phase!\n',
      file=file,append=TRUE,sep='')
}


mcmc = function(pars,surveillance,MCMCits=10000,SUBits=100)
{
  pars$LL = loglikelihood(pars,surveillance)
  storepars = list(n0=rep(0,MCMCits),b=rep(0,MCMCits),
                   L1=rep(0,MCMCits),L2=rep(0,MCMCits),
                   p1 =rep(0,MCMCits),p2 =rep(0,MCMCits),
                   LL=rep(0,MCMCits),size=rep(0,MCMCits),sizenow=rep(0,MCMCits))
  
  for(iteration in 1:MCMCits)
  {
    if(iteration%%100==0)cat('Iteration',iteration,'of',MCMCits,'log-L:',round(pars$LL),'\n')
    for(subit in 1:SUBits)
    {
      if(!pars$fixed_n0){oldpars = pars; pars$n0 = rnorm(1,pars$n0,0.1); pars = mh(oldpars,pars)}
      if(!pars$fixed_b){oldpars = pars; pars$b = rnorm(1,pars$b,0.01); pars = mh(oldpars,pars)}
      if(!pars$fixed_L1){oldpars = pars; pars$L1 = rnorm(1,pars$L1,0.1); pars = mh(oldpars,pars)}
      if(!pars$fixed_L2){oldpars = pars; pars$L2 = rnorm(1,pars$L2,0.1); pars = mh(oldpars,pars)}
      if(!pars$fixed_p1){oldpars = pars; pars$p1 = rnorm(1,pars$p1,0.01); pars = mh(oldpars,pars)}
      if(!pars$fixed_p2){oldpars = pars; pars$p2 = rnorm(1,pars$p2,0.01); pars = mh(oldpars,pars)}
    }
    
    
    storepars$n0[iteration] = pars$n0
    storepars$b[iteration] = pars$b
    storepars$L1[iteration] = pars$L1
    storepars$L2[iteration] = pars$L2
    storepars$p1[iteration] = pars$p1
    storepars$p2[iteration] = pars$p2
    storepars$LL[iteration] = pars$LL
    storepars$size[iteration]=casestodate(pars,length(surveillance$s1))
    storepars$sizetoday[iteration]=casestoday(pars,length(surveillance$s1))
  }
  
  return(storepars)
}