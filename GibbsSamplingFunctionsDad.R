library(mlbench)
library(MASS)
library(lars)
library(Matrix)
library(lattice)
library(truncnorm)
library(MCMCpack)

# Gibb Sampling Functions -------------------------------------------------

reportMCMC = function(obj)
{
  draws = obj
  if(!is.null(dim(draws)))
  {
    for (i in 1:dim(draws)[2])
    {
      ts.plot(draws[,i],xlab="iterations",ylab="",main=colnames(draws)[i]) 
      abline(h=mean(draws[,i]),col=2,lwd=2)
      acf(draws[,i],main="") 
      hist(draws[,i],prob=T,main="",xlab="") 
      abline(v=mean(draws[,i]),col=2,lwd=2)
    }
  }
  else
  {
    ts.plot(draws,xlab="iterations",ylab="",main=names(draws)[1]) 
    abline(h=mean(draws),col=2,lwd=2)
    acf(draws,main="") 
    hist(draws,prob=T,main="",xlab="") 
    abline(v=mean(draws),col=2,lwd=2)
    
  }
}

sampleB = function(y,X,sig2,tau2,mu)
{
  Tmat = diag(tau2)
  Tmatinv = diag(1/tau2)
  SIG1inv = (1/sig2)*((t(X))%*%X) + Tmatinv
  # SIG1 = ginv(as.matrix(SIG1inv))
  SIG1 = solve(SIG1inv)
  m = (1/sig2)*SIG1%*%t(X)%*%y + SIG1%*%Tmatinv%*%mu
  draw = mvrnorm(1,m,SIG1)
  return(draw)
}

samplesig2 = function(y,X,B)
{
  n = length(y)
  p1 = (n+4)/2
  p2 = as.numeric( ( (t(y - (X%*%B)))%*%(y - (X%*%B)) )/2 )
  draw = rgamma(1,p1,p2)
  return(1/draw)
}

sampletau2 = function(B)
{
  NS = length(B)
  p1 = (NS+3)/2
  p2 = ( (t(B))%*%B )/2
  draw = rgamma(1,p1,p2)
  return(1/draw)
}

samplemu = function(B,tau2,m,v2)
{
  NS = length(B)
  thetabar = mean(B)
  sigmu2 = 1/(NS/tau2 + 1/v2)
  mmu = sigmu2*((NS/tau2)*thetabar + (1/v2)*m)
  draw = rnorm(1,mmu,sqrt(sigmu2))
  return(draw)
}

samplealpha = function(ymxB,sigalpha2,sig2)
{  
  N = length(ymxB)
  ymxBbar = mean(ymxB)
  sa2 = 1/(N/sig2 + 1/sigalpha2)
  malpha = sa2*((N/sig2)*ymxBbar)
  draw = rnorm(1,malpha,sqrt(sa2))
}

samplesigmalpha2 = function(mu)
{
  N = length(mu)
  p1 = (N+3)/2
  p2 = ( (t(mu))%*%mu )/2
  draw = rgamma(1,p1,p2)
  return(1/draw)
}

# The wrapper function ---------------------------------------------------

Gibbswrapper = function(loops,y,X,numi,numj,alphaIDlist,BPrior)
{
  # prior on mui
  m = 0
  sig02 = 50
  
  # other stuff
  size = dim(X)[2]
  numpred = size / (numi*numj)
  
  BMCMC = matrix(0,size,loops)
  sig2MCMC = rep(0,loops)
  names(sig2MCMC) = rep('sig2',length(sig2MCMC))
  tau2MCMC = matrix(0,numpred,loops)
  muMCMC = matrix(0,numpred,loops)
  alphaMCMC = matrix(0,numi,loops)
  sigmalpha2MCMC = rep(0,loops)
  
  rownames(tau2MCMC) = paste('tau2',1:numpred)
  rownames(muMCMC) = paste('mu',1:numpred)
  
  # starting points for MCMC
  BMCMC[,1] = rep(1,size)
  sig2MCMC[1] = 1
  sigmalpha2MCMC[1] = 50
  tau2MCMC[,1] = rep(1,numpred)
  muMCMC[,1] = rep(1,numpred)
  alphaMCMC[,1] = rep(1,numi)
  loopind = seq(from = 1,to = size,by = numpred)
  
  #*************
  ## THE MCMC ##
  #*************
  for(i in 2:loops)
  {
    if(i%%10==0){print(noquote(paste('MCMC iter =',i)))}
    
    # constructing fixed effect vector
    alphas = rep(0,length(y))
    for(j in 1:numi)
    {
      ind = alphaIDlist[[j]]
      alphas[ind] = alphaMCMC[j,i-1]
    }
    
    # SAMPLE betas
    BMCMC[,i] = as.numeric(sampleB(y-alphas,X,sig2MCMC[i-1],rep(tau2MCMC[,i-1],numi*numj),rep(muMCMC[,i-1],numi*numj)))
    
    # SAMPLE sig2
    sig2MCMC[i] = samplesig2(y-alphas,X,BMCMC[,i])
    
    if(BPrior==TRUE)
    {
      # SAMPLE the prior parameters on the betas .. mu and tau2 (both vectors)
      for(j in 1:numpred)
      {
        predin = (loopind+j-1)
        tau2MCMC[j,i] = sampletau2((BMCMC[predin,i]-muMCMC[j,i-1]))
        muMCMC[j,i] = samplemu(BMCMC[predin,i],tau2MCMC[j,i],m,sig02)    
      }  
    }
    
    if(BPrior==FALSE)
    {
      tau2MCMC[,i] = 50
      muMCMC[,i] = 0
    }
    
    # sample prior parameter sigalpha2
    sigmalpha2MCMC[i] = samplesigmalpha2(alphaMCMC[,i-1])
    # sigmalpha2MCMC[i] = 50
    
    # SAMPLE the country i level fixed effects (alphas)
    for(j in 1:numi)
    {
      ind = alphaIDlist[[j]]
      ymxB = y[ind] - X[ind,]%*%BMCMC[,i]
      alphaMCMC[j,i] = samplealpha(ymxB,sigmalpha2MCMC[i],sig2MCMC[i])
    }
  }
  return(list(BMCMC,sig2MCMC,tau2MCMC,muMCMC,alphaMCMC,sigmalpha2MCMC))
}
