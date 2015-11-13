## playing around with data ##
library(foreign)
library(zoo)
library(plm)
library(Matrix)

setwd("~/Documents/Dad/CarmelRealEstate/")
source('GibbsSamplingFunctionsDad.R')
carmel = read.csv('Carmel real estate Nov 2013 to Nov 2015.csv',na.strings = '--')

# Data --------------------------------------------------------------------

N = dim(carmel)[1]

# unique values
uneigh = unique(carmel$Location.description)

# Building design matrix --------------------------------------------------

# removing missing data
nas = unique(c(which(is.na(carmel$Beds)),which(is.na(carmel$Baths.Total)),which(is.na(carmel$Baths.Full)),which(is.na(carmel$Lotsqft)),which(is.na(carmel$DOM))))
carmel2 = carmel[-nas,]
M = dim(carmel2)[1]

# construct X
numcovar = 5
size = numcovar*length(uneigh)
loopind = seq(from = 1,to = size,by = numcovar)
X = Matrix(matrix(0,M,size),sparse=T) #using sparse matrix class

# indicies for each i and j pair
indlist = list()
k=1
for(i in 1:length(uneigh))
{
  indlist[[k]] = which(carmel2$Location.description==uneigh[i])
  k=k+1
}

# constructing the covariates
neighnum = 1
for(i in loopind)
{
  ind = indlist[[neighnum]]
  X[ind,i] = carmel2$Beds[ind]
  X[ind,i+1] = carmel2$Baths.Full[ind]
  X[ind,i+2] = carmel2$Baths.Total[ind]
  X[ind,i+3] = carmel2$Lotsqft[ind]
  X[ind,i+4] = carmel2$DOM[ind]
  neighnum=neighnum+1
}

# saving the response as y
y = as.numeric(carmel2$Closeprice) / as.numeric(carmel2$Homesqft)

# Estimating Model via Gibbs Sampler --------------------------------------

alphaIDlist = list()
for(i in 1:length(uneigh))
{
  alphaIDlist[[i]] = which(carmel2$Location.description==uneigh[i])
}

source('GibbsSamplingFunctionsDad.R')
loops = 8000
numi = length(uneigh)
results = Gibbswrapper(loops,y,X,numi,alphaIDlist,BPrior=TRUE)

# Summarizing results -----------------------------------------------------

BMCMC = results[[1]]
sig2MCMC = results[[2]]
tau2MCMC = results[[3]]
muMCMC = results[[4]]
alphaMCMC = results[[5]]
sigmalpha2MCMC = results[[6]]

range=3000:loops

# naming betas
rownames(BMCMC) = rep("NA",dim(BMCMC)[1])
count=1
for(i in 1:length(ucounti))
{
  for(j in 1:length(ucountj))
  {
    ii = loopind[count]
    rownames(BMCMC)[ii] = paste(paste(ucounti[i],"/",ucountj[j],sep=''),paste('B',1,sep=''),sep=' ')
    rownames(BMCMC)[ii+1] = paste(paste(ucounti[i],"/",ucountj[j],sep=''),paste('B',2,sep=''),sep=' ')
    rownames(BMCMC)[ii+2] = paste(paste(ucounti[i],"/",ucountj[j],sep=''),paste('B',3,sep=''),sep=' ')
    count=count+1
  }
}

# naming alphas
rownames(alphaMCMC) = rep("NA",dim(alphaMCMC)[1])
for(i in 1:length(ucounti))
{
  rownames(alphaMCMC)[i] = paste(ucounti[i],'alpha',sep=' ')
}

pdf('prelimresults-BETAS.pdf',paper = "a4", width = 0, height = 0)
par(mfrow = c(6,3))
reportMCMC(t(BMCMC[,1000:2000]))
dev.off()

pdf('prelimresults-alphas.pdf',paper = "a4", width = 0, height = 0)
par(mfrow = c(3,3))
reportMCMC(t(alphaMCMC[,1000:2000]))
dev.off()


plot(sig2MCMC[range],type='l',col=2)
plot(sigmalpha2MCMC[range],type='l',col=5)

plot(tau2MCMC[1,range],type='l')
plot(tau2MCMC[2,range],type='l')
plot(tau2MCMC[3,range],type='l')

plot(muMCMC[1,range],type='l')
plot(muMCMC[2,range],type='l')
plot(muMCMC[3,range],type='l')

plot(BMCMC[1,range],type='l')
plot(BMCMC[2,range],type='l')
plot(BMCMC[3,range],type='l')
plot(BMCMC[4,range],type='l')

plot(BMCMC[10,range],type='l')
plot(BMCMC[11,range],type='l')
plot(BMCMC[12,range],type='l')
plot(BMCMC[13,range],type='l')

