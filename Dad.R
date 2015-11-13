## playing around with data ##
library(foreign)
library(zoo)
library(plm)
library(Matrix)

setwd("~/Documents/Dad/CarmelRealEstate/")
source('GibbsSamplingFunctionsDad.R')
carmel = read.csv('Carmel real estate Nov 2013 to Nov 2015.csv',na.strings = '--')

# Functions ---------------------------------------------------------------

L <- function(x, k) 
{
  c(rep(NA, k), x)[1 : length(x)] 
}

deseasonalizeQ <- function (x)
{
  x <- ts(x)
  #Step1: Centered moving averages: create cma time series having the same length with the original time series x
  # cma has 2 NAs on both ends.
  cma <- filter(x, filter = c(1/8, 1/4, 1/4, 1/4, 1/8), sides=2)
  
  #Step2: Ratios = Original time series / centered moving averages
  ratio <- x/cma
  
  #Step3: Unadjusted 4 seasonal indexes
  unadj4si <- ts(1:4)
  # floor((length(x)-4)/4)  #"-4" is 4 NA at both ends; below "-1" is due to starting "0:" in multiplication
  
  unadj4si[1] <- mean(ratio[3+4*(0:(floor((length(x)-4)/4) - 1))])
  unadj4si[2] <- mean(ratio[4+4*(0:(floor((length(x)-4)/4) - 1))])
  unadj4si[3] <- mean(ratio[5+4*(0:(floor((length(x)-4)/4) - 1))])
  unadj4si[4] <- mean(ratio[6+4*(0:(floor((length(x)-4)/4) - 1))])
  
  #Step4: Adjusted 4 seasonal indexes
  adj4si <- ts(1:4)
  adj4si[1] <- unadj4si[1]/mean(c(unadj4si[1],unadj4si[2],unadj4si[3],unadj4si[4]))
  adj4si[2] <- unadj4si[2]/mean(c(unadj4si[1],unadj4si[2],unadj4si[3],unadj4si[4]))
  adj4si[3] <- unadj4si[3]/mean(c(unadj4si[1],unadj4si[2],unadj4si[3],unadj4si[4]))
  adj4si[4] <- unadj4si[4]/mean(c(unadj4si[1],unadj4si[2],unadj4si[3],unadj4si[4]))
  
  #Step5: Propogated adjusted seasonal indexes
  propadjsi <- ts(1:length(x))
  
  propadjsi[3+4*(0:(floor((length(x)-4)/4) - 1))] <- adj4si[1]
  propadjsi[4+4*(0:(floor((length(x)-4)/4) - 1))] <- adj4si[2]
  propadjsi[5+4*(0:(floor((length(x)-4)/4) - 1))] <- adj4si[3]
  propadjsi[6+4*(0:(floor((length(x)-4)/4) - 1))] <- adj4si[4]
  
  propadjsi[1] <- adj4si[3]
  propadjsi[2] <- adj4si[4]
  propadjsi[length(x)-1] <- adj4si[1]
  propadjsi[length(x)] <- adj4si[2]
  
  #Step6: Deseasonalized values
  out <- x/propadjsi  # deseasonalized = x/propadjsi
  return(out)
}

deseasonalizeDavidQ = function(x)
{
  Time = 1:length(x)
  N = rep(1:4,length(x)/4)
  
  Q1 = (N==1)
  Q2 = (N==2)
  Q3 = (N==3)
  
  fit = lm(x~Time+Q1+Q2+Q3)
  
  alpha = fit$coefficients[1]
  beta = fit$coefficients[2]
  deseason = alpha + beta*Time + fit$residuals
  return(deseason)
}

# Data --------------------------------------------------------------------

N = dim(carmel)[1]

# unique values
uneigh = unique(carmel$Location.description)

# c=1
# for(i in 1:length(ucounti))
# {
#   for(j in 1:length(ucountj))
#   {
#     ind = intersect(which(newmac$country_j==ucountj[j]),which(newmac$country_i==ucounti[i]))
#     if(length(ind)!=136) { print(c); print(length(ind)); print(ind)}
#     # export share change for all i and j
#     exp_share = newmac$exp_share[ind]
#     newmac$dexp_share[ind] = log(exp_share) - log(L(exp_share,4))
#     
#     # fx change for all i and j
#     fx = newmac$fx[ind]
#     newmac$dfx[ind] = log( (1 - ((fx)/(L(fx,4)))) + 1 )
#     
#     # comp change for all i and j
#     comp = newmac$comp[ind]
#     newmac$dcomp[ind] = log(comp) - log(L(comp,4))
#     
#     # gdp change for all i and j
#     gdp = newmac$gdp_j[ind]
#     newmac$dgdp[ind] = log( ((gdp)/(L(gdp,4)) - 1) + 1 )
#     c=c+1
#   }
# }

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
for(i in 1:length(ucounti))
{
  for(j in 1:length(ucountj))
  {
    indlist[[k]] = intersect(which(newmac2$country_j==ucountj[j]),which(newmac2$country_i==ucounti[i]))
    k=k+1
  }
}

# constructing the covariates
countrypairnum = 1
for(i in loopind)
{
  ind = indlist[[countrypairnum]]
  X[ind,i] = newmac2$dfx[ind]
  X[ind,i+1] = newmac2$dcomp[ind]
  X[ind,i+2] = newmac2$dgdp[ind]
  countrypairnum=countrypairnum+1
}

# saving the response as y
y = newmac2$dexp_share

# Estimating Model via Gibbs Sampler --------------------------------------

numi = length(ucounti)
numj = length(ucountj)
loops = 8000

alphaIDlist = list()
for(i in 1:length(ucounti))
{
  alphaIDlist[[i]] = which(newmac2$country_i==ucounti[i])
}

source('GibbsSamplingFunctionsKent.R')
results = Gibbswrapper(loops,y,X,numi,numj,alphaIDlist,BPrior=TRUE)

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

