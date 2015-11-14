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
y = log(as.numeric(carmel2$Closeprice) / as.numeric(carmel2$Homesqft))

# Estimating Model via Gibbs Sampler --------------------------------------

alphaIDlist = list()
for(i in 1:length(uneigh))
{
  alphaIDlist[[i]] = which(carmel2$Location.description==uneigh[i])
}

source('GibbsSamplingFunctionsDad.R')
loops = 10000
numi = length(uneigh)
results = Gibbswrapper(loops,y,X,numi,alphaIDlist,BPrior=FALSE)

# Summarizing results -----------------------------------------------------

BMCMC = results[[1]]
sig2MCMC = results[[2]]
tau2MCMC = results[[3]]
muMCMC = results[[4]]
alphaMCMC = results[[5]]
sigmalpha2MCMC = results[[6]]

range=(loops/2):loops

# naming betas

rownames(BMCMC) = rep("NA",dim(BMCMC)[1])
count=1
for(i in 1:length(uneigh))
{
  ii = loopind[count]
  rownames(BMCMC)[ii] = paste(paste(uneigh[i],'Beds',paste('B',1,sep=''),sep=' '))
  rownames(BMCMC)[ii+1] = paste(paste(uneigh[i],'Full Baths',paste('B',2,sep=''),sep=' '))
  rownames(BMCMC)[ii+2] = paste(paste(uneigh[i],'Total Baths',paste('B',3,sep=''),sep=' '))
  rownames(BMCMC)[ii+3] = paste(paste(uneigh[i],'Lot sqft',paste('B',4,sep=''),sep=' '))
  rownames(BMCMC)[ii+4] = paste(paste(uneigh[i],'DOM',paste('B',5,sep=''),sep=' '))
  count=count+1
}

# naming alphas
rownames(alphaMCMC) = rep("NA",dim(alphaMCMC)[1])
for(i in 1:length(uneigh))
{
  rownames(alphaMCMC)[i] = paste(uneigh[i],'alpha',sep=' ')
}

pdf('prelimresults-BETAS.pdf',paper = "a4", width = 0, height = 0)
par(mfrow = c(5,3))
reportMCMC(t(BMCMC[,range]))
dev.off()

pdf('prelimresults-alphas.pdf',paper = "a4", width = 0, height = 0)
par(mfrow = c(3,3))
reportMCMC(t(alphaMCMC[,range]))
dev.off()

pdf('prelimresults-mu.pdf',paper = "a4", width = 0, height = 0)
par(mfrow = c(3,3))
reportMCMC(t(muMCMC[,range]))
dev.off()

pdf('prelimresults-tau.pdf',paper = "a4", width = 0, height = 0)
par(mfrow = c(3,3))
reportMCMC(t(tau2MCMC[,range]))
dev.off()


