###############################################################################
# This supplementary R script provides methods to
# 
#   simulate test positive/negative bovine PTBC data,
# 
#   fit the PTBC prevalence model described in the main text 
#   of Lang et al. (2024), using maximum likelihood methods,
# 
#   calculate the factor of proportionality using numerical differentiation,
#
#   diagnose the quality of fit of the model using model based simulation.
#
###############################################################################

# Packages needed
library(MASS)    # to draw sample from negative binomial distribution
library(mvtnorm) # to draw sample from multivariate normal distribution
library(gtools)  # we need the R function inv.logit
library(DPQ)     # we need the R function qbeta.R
library(maxLik)  # maximum likelihood estimation
library(rmutil)  # this package contains the function gauss.hermite
library(geepack) # GEE models are applied to obtain initial estimates
library(pracma)  # Numerical derivatives

##########################################################################
# 1. Simulation of test positive/negative bovine PTBC data 
# 
##########################################################################

# Set your working directory. The data set with the simulated MAP ELISA 
# test results will be created in it.
setwd("...")

##########################################################################
# Set the parameters of the simulated data set
#
##########################################################################

nHERD         = 40   # number of herds
prp_primi     = 0.4  # overall proportion of primiparous cows
M_herdsize    = 50   # minimum value of herd size
E_herdsize    = 250  # expected value of herd size
SD_herdsize   = 100  # standard deviation of herd size  
M_age_primi   = 1.8  # minimum value of age at first calving of primipara 
E_age_primi   = 2.3  # expected value of age at first calving of primipara
SD_age_primi  = 0.5  # standard deviation of age at first calving of primipara
D_age_primi   = 0.5  # expected value of age at ELISA test minus age at first calving of primipara
M_age_multi   = 2.8  # minimum value of age at ELISA test of multipara 
E_age_multi   = 4.8  # expected value of age at ELISA test of multipara
SD_age_multi  = 1.5  # standard deviation of age at ELISA test of multipara

##########################################################################
# Simulation of the data frame
# 
# HERD_ID      identifier of a herd (sequential number)
# COW_ID       identifier of a cow (sequential number)
# MULTIPAR     equals 1 for multiparous, 0 for primiparous cows, respectively
# AGE	         age in years
# CALVING_AGE  age in years at first calving for primiparous cows
# POS          indicator of test positive cows (1 positive, 0 negative)

herdsize = M_herdsize + rnegbin(nHERD, 
           mu=E_herdsize-M_herdsize, 
           theta=(E_herdsize-M_herdsize)^2/(SD_herdsize^2-E_herdsize+M_herdsize))

COW_ID   = 10^(ceiling(log10(sum(herdsize))))
HERD_ID  = 10^(ceiling(log10(nHERD)))

shrdsize = sum(herdsize)
PTBCdata = data.frame(HERD_ID=rep(NA_real_,shrdsize), COW_ID=NA_real_, MULTIPAR=NA_real_, 
                      AGE=NA_real_, CALVING_AGE=NA_real_, POS=NA_real_)

lastcow=0
for (n in 1:nHERD) {
  multi = sort(rbinom(herdsize[n],1,1-prp_primi))
  for (s in 1:herdsize[n]) {
    if (multi[s]==0) {
      PTBCdata[lastcow+s, "MULTIPAR"]    = 0
      PTBCdata[lastcow+s, "CALVING_AGE"] = M_age_primi + 
        rgamma(1,shape=(E_age_primi-M_age_primi)^2/SD_age_primi^2, 
               scale=SD_age_primi^2/(E_age_primi-M_age_primi))
      
      PTBCdata[lastcow+s, "AGE"] = PTBCdata[lastcow+s, "CALVING_AGE"] + 
        rgamma(1,shape=D_age_primi*(E_age_primi-M_age_primi)/SD_age_primi^2, 
               scale=SD_age_primi^2/(E_age_primi-M_age_primi))
      
    } else {
      PTBCdata[lastcow+s, "MULTIPAR"] = 1
      PTBCdata[lastcow+s, "AGE"]      = M_age_multi + 
        rgamma(1,shape=(E_age_multi-M_age_multi)^2/SD_age_multi^2, 
               scale=SD_age_multi^2/(E_age_multi-M_age_multi))
    }
    PTBCdata[lastcow+s, "COW_ID"] = COW_ID+lastcow+s
    PTBCdata[lastcow+s,"HERD_ID"] = HERD_ID+n
  }
  lastcow = lastcow + herdsize[n]
}

# View(PTBCdata)
# by(PTBCdata$AGE, PTBCdata$MULTIPAR, mean)

##########################################################################
# Simulation of ELISA test results (1=positive, 0=negative)  
# in the variable POS for each cow in each herd
#
##########################################################################

PTBCdata$one     = 1 # Ancillary variable
PTBCdata$HERD_ID = factor(PTBCdata$HERD_ID)
# str(PTBCdata)
# View(PTBCdata)

# Herd level aggregated data
PTBCHerdDat = aggregate(
  data.frame(
    nCOW = PTBCdata$one),
  list(
    MULTIPAR = PTBCdata$MULTIPAR,
    HERD_ID  = PTBCdata$HERD_ID),
  sum)

# str(PTBCHerdDat)
# View(PTBCHerdDat)

# Cumulated numbers of records (cows)
cnCOW = c(0,cumsum(PTBCHerdDat$nCOW))
# cnCOW

# Parameters of age-specific sensitivity of the MAP ELISA test (Meyer et al. 2018)
parA = 1.2 
parB = 3.0
parC = 0.30

# Specificity (Meyer et al. 2018)
Sp   = 0.995

# Standard deviation of the herd level random effect for primiparous cows
sigma1   = 0.20
# Standard deviation of the herd level random effect for multiparous cows
sigma2   = 0.40
# Pearson correlation coefficient between the herd level random effects 
# for primiparous and multiparous cows
rPearson = 0.50

# Dispersion model
psi1 = sigma1^(-2) 
psi2 = sigma2^(-2)

# Coefficients in the expressions for the marginal means (44) and (45)
# in the main text of Lang et al. (2024)
beta10 = -4.00 # primiparous cows, intercept
beta11 =  2.00 # primiparous cows, coefficiet of ln(age at first calving)
beta20 = -2.00 # multiparous cows, intercept

for (n in 1:nHERD) {
  # Herd level standardized random effects
  Eta = rmvnorm(n=1, mean=c(0,0), sigma=matrix(c(1,rPearson,rPearson,1), ncol=2))
  
  # Primiparous cows
  for (k in (cnCOW[2*n-1]+1):cnCOW[2*n]) {
    # Age dependent sensitivity of the kth cow (Meyer et al. 2018)
    Se  = inv.logit(parA-parB*exp(-parC*PTBCdata[k,"AGE"]))
    
    # Mean of conditional within herd animal level true prevalence for primiparous cows,
    # according to equation (44)
    mu1 = inv.logit(beta10+beta11*log(PTBCdata[k,"CALVING_AGE"]))
    
    # Conditional within herd animal level true prevalence for primiparous cows
    P1  = qbeta.R(pnorm(Eta[1],0,1), mu1*psi1, (1-mu1)*psi1)
    
    # Apparent prevalence
    AP1 = Se*P1 + (1-Sp)*(1-P1)
    
    # Simulated MAP ELISA test result for the kth cow 
    PTBCdata[k,"POS"] = rbinom(1, 1, AP1)
  }
  
  # Multiparous cows
  for (k in (cnCOW[2*n]+1):cnCOW[2*n+1]) {
    # Age dependent sensitivity of the kth cow (Meyer et al. 2018)
    Se  = inv.logit(parA-parB*exp(-parC*PTBCdata[k,"AGE"]))
    
    # Mean of conditional within herd animal level true prevalence for multiparous cows,
    # according to equation (45)
    mu2 = inv.logit(beta20)
    
    # Conditional within herd animal level true prevalence for multiparous cows
    P2  = qbeta.R(pnorm(Eta[2],0,1), mu2*psi2, (1-mu2)*psi2)
    
    # Apparent prevalence
    AP2 = Se*P2 + (1-Sp)*(1-P2)
    
    # Simulated MAP ELISA test result for the kth cow 
    PTBCdata[k,"POS"] = rbinom(1, 1, AP2)
  }
}

# Number of simulated test positive cases
# sum(PTBCdata$POS)

# Simulated apparent prevalence for primiparous and multiparous cows
# with(PTBCdata, by(POS, MULTIPAR, mean))

# Save the simulated data for further analysis in the current working directory
write.table(PTBCdata, "PTBC_simulated_data.csv", sep=";", dec=".",
            append=F, col.names=T, row.names=F)

##########################################################################
# 2. Fit the PTBC prevalence model using frequentist maximum likelihood
# method
# 
##########################################################################

# Packages needed
library(MASS)    # to draw sample from negative binomial distribution
library(mvtnorm) # to draw sample from multivariate normal distribution
library(gtools)  # we need the R function inv.logit
library(DPQ)     # we need the R function qbeta.R
library(maxLik)  # maximum likelihood estimation
library(rmutil)  # this package contains the function gauss.hermite
library(geepack) # GEE models are applied to obtain initial estimates
library(pracma)  # Numerical derivatives

# Set your working directory
setwd("...")

# Simulated data
PTBCdata = read.table("PTBC_simulated_data.csv", sep=";", dec=".", header=T)
PTBCdata$HERD_ID = factor(PTBCdata$HERD_ID)

# str(PTBCdata)
# View(PTBCdata)

# We group the ages of the animals into increasing categories: 
# 1.75, 2.00, ..., 2.75 years (increment=0.25), 
# 3.00, 3.33, ..., 4.33 years (increment=0.33), and 
# 4.50, 5.00, 5.50, ..., 16 years (increment=0.5).

PTBCdata[PTBCdata$AGE<=3,"AGEG"] =
  round(4*PTBCdata[PTBCdata$AGE<=3,"AGE"])/4

PTBCdata[PTBCdata$AGE>3 & PTBCdata$AGE<=4.5,"AGEG"] =
  round(3*PTBCdata[PTBCdata$AGE>3 & PTBCdata$AGE<=4.5,"AGE"])/3

PTBCdata[PTBCdata$AGE>4.5,"AGEG"] =
  round(2*PTBCdata[PTBCdata$AGE>4.5,"AGE"])/2

PTBCdata$AGEG = round(PTBCdata$AGEG,3)

PTBCdata$lCAGE = round(log(PTBCdata$CALVING_AGE),2)
PTBCdata[PTBCdata$MULTIPAR==1,"lCAGE"] = 0

# str(PTBCdata)
# View(PTBCdata)

# Herd level aggregated data
PTBCHerdDat = aggregate(
  data.frame(
    POS  = PTBCdata$POS,
    nCOW = PTBCdata$one),
  list(
    MULTIPAR = PTBCdata$MULTIPAR,
    HERD_ID  = PTBCdata$HERD_ID),
  sum)

# str(PTBCHerdDat)
# View(PTBCHerdDat)
# by(PTBCHerdDat$nCOW, INDICES=list(PTBCHerdDat$MULTIPAR), mean)

# Herd level aggregated data, 2nd version
PTBCHerd2 = aggregate(
  data.frame(
    POS  = PTBCdata$POS,
    nCOW = PTBCdata$one),
  list(
    AGEG     = PTBCdata$AGEG,
    lCAGE    = PTBCdata$lCAGE,
    MULTIPAR = PTBCdata$MULTIPAR,
    HERD_ID  = PTBCdata$HERD_ID),
  sum)

# str(PTBCHerd2)
# View(PTBCHerd2)

PTBCHerd2$one = 1 # ancillary variable

# Herd level aggregated data, 3rd version
PTBCHerd3 = aggregate(
  data.frame(
    nRec = PTBCHerd2$one),
  list(
    MULTIPAR = PTBCHerd2$MULTIPAR,
    HERD_ID  = PTBCHerd2$HERD_ID),
  sum)

# str(PTBCHerd3)
# View(PTBCHerd3)

# Cumulated numbers of records in PTBCHerd3
cnRec = c(0,cumsum(PTBCHerd3$nRec))
cnRec

##########################################################################
# Frequentist maximum likelihood analysis
#
# List of input data
LData =
  list(
    nHERD   = nlevels(PTBCdata$HERD_ID), # Number of all herds
    cnRec   = cnRec, # Cumulated numbers of records in PTBCHerd3
    nCOW    = PTBCHerd2$nCOW, # Number of cows in the records of PTBCHerd2
    POS     = PTBCHerd2$POS,  # Number of positive cows in the records of PTBCHerd2
    AgeG    = PTBCHerd2$AGEG, # Age groups of cows in PTBCHerd2, measured in years
    lCAGE   = PTBCHerd2$lCAGE ) # log-CALVING_Age groups of cows in PTBCHerd2, measured in years

# str(LData)
# LData

# Parameters of the model of age-dependent sensitivity (Meyer et al. 2018)
parA <<- 1.2 
parB <<- 3.0
parC <<- 0.30

# Specificity (Meyer et al. 2018)
Sp   <<- 0.995

##################################################################
# Multivariate Gauss-Hermite quadrature
#
##################################################################
## The function mgauss.hermite computes nodes (points) and weigths for
## multivariate Gauss-Hermite quadrature according to Jaeckel 2005. 
##
## It is available at
## http://biostatmatt.com/archives/2754 or at
## https://www.r-bloggers.com/2015/09/notes-on-multivariate-gaussian-quadrature-with-r-code/
##
## n     - number of points along each dimension before pruning
## mu    - mean vector
## sigma - covariance matrix
## prune - NULL - no pruning; [0-1] - fraction to prune
##################################################################
mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")
  
  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  
  # idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }
  
  ## rotate, scale, translate points
  eig <- eigen(sigma) 
  if (ncol(eig$vectors)>1.1) {
    rot <- eig$vectors %*% diag(sqrt(eig$values))
  } else {
    rot <- eig$vectors %*% diag(matrix(sqrt(eig$values),1,1))
  }
  pts <- t(rot %*% t(pts) + mu)
  
  return (list(points=pts, weights=wts))
}

###################################################################
# Marginal log-likelihood function
#
###################################################################
Marg_llik = function(
fixPar, # vector of fixed parameters, i.e. beta10, beta11, beta20, 
        # lsigm1, lsigm2, fisherz
nAGQ=1, # number of sample points
prune=NULL
) {

  beta10 <<- fixPar[1]  # intercept in (44)
  beta11 <<- fixPar[2]  # coefficient of log calving age
  beta20 <<- fixPar[3]  # intercept in (45)

  lsigm1  = fixPar[4]   # log standard deviation of random herd effects - primiparous cows
  lsigm2  = fixPar[5]   # log standard deviation of random herd effects - multiparous cows
  fisherz = fixPar[6]   # Fisher's z transformation of the Pearson correlation coefficient 
                        # between the random herd effects related to primiparous and 
                        # multiparous cows

  sigma1 = exp(lsigm1)
  sigma2 = exp(lsigm2)

  # Correlation is expressed from Fisher's z transformation
  rPearson = (exp(fisherz)-exp(-fisherz))/(exp(fisherz)+exp(-fisherz))

  # Dispersion model  
  psi1 <<- sigma1^(-2.0)
  psi2 <<- sigma2^(-2.0)

  # Note that the quadrature applied to standardized random effects yields
  # numerically more stable estimates.
  sig   = matrix(c(1,rPearson,rPearson,1),2,2)
  pts   = mgauss.hermite(nAGQ, mu=c(0,0), sigma=sig, prune=prune)

  # The components of the loglikelihood are calculated for each herd 
  llik = 0
  for (iHerd in 1:LData$nHERD) {

    # Marginal likelihood component
    MargLik = sum(apply(X=pts$points, MARGIN=1, FUN=Lik_Comp, iHerd) * pts$weights)

    llik = llik + log(MargLik)
  }

  iterac <<- iterac+1
  print(iterac)
  print(llik, digits=10)
  
  return (llik)
}

###################################################################
# The component of the likelihood is calculated for each herd 
#
###################################################################
Lik_Comp = function(
eta,  # Vector of random effects
iHerd # The index of the actual herd
) {
  t1 = 0 # Loglikelihood component 

  # Primiparous cows
  for (k in (LData$cnRec[2*iHerd-1]+1):LData$cnRec[2*iHerd]) {
    # Conditional within-herd animal-level true prevalence among primiparous cows
    mu1   = 1/(1+exp(-(beta10+beta11*LData$lCAGE[k])))
    CWHP1 = qbeta.R(pnorm(eta[1]), mu1*psi1, (1-mu1)*psi1) 
    
    # Age-dependent sensitivity of the kth cow
    Se = 1/(1+exp(-(parA-parB*exp(-parC*LData$AgeG[k]))))

    # Apparent prevalence
    pi1   = Se*CWHP1 + (1-Sp)*(1-CWHP1)

    t1 = t1 + dbinom(x=LData$POS[k], size=LData$nCOW[k], prob=pi1, log=T)
  }

  # Multiparous cows
  for (k in (LData$cnRec[2*iHerd]+1):LData$cnRec[2*iHerd+1]) {
    # Conditional within-herd animal-level true prevalence among multiparous cows
    mu2   = 1/(1+exp(-(beta20)))
    CWHP2 = qbeta.R(pnorm(eta[2]), mu2*psi2, (1-mu2)*psi2) 
    
    # Age-dependent sensitivity of the kth cow
    Se = 1/(1+exp(-(parA-parB*exp(-parC*LData$AgeG[k]))))

    # Apparent prevalence
    pi2  = Se*CWHP2 + (1-Sp)*(1-CWHP2)

    t1 = t1 + dbinom(x=LData$POS[k], size=LData$nCOW[k], prob=pi2, log=T)
  }

  t1 = exp(t1)
  
  return (t1)
}

###################################################################
# Here we run the maximum likelihood estimation
#
###################################################################

# GEE model to obtain initial estimates
geemod = geeglm(cbind(POS,nCOW) ~ -1+factor(MULTIPAR)+I((1-MULTIPAR)*lCAGE),
                 id=HERD_ID, family=binomial, data=PTBCHerd2)
summary(geemod)

# Constraints
A = rbind(diag(6),-diag(6)) 
A  
B = matrix(10,12,1)  
B  

# Iteration is started from the parameters of the GEE model. 
# The initial standard deviations and the 
# correlation between the two random effects is set to 0.  
coefstart = as.numeric(c(coef(geemod)[1],coef(geemod)[3],coef(geemod)[2],0,0,0))
iterac <<- 0
ml  = maxLik(Marg_llik, start=coefstart, 
      method="BFGS", constraints=list(ineqA=A, ineqB=B), nAGQ=3, prune=0.2)

sml = summary(ml)
sml

# Improvement
iterac <<- 0
ml  = maxLik(Marg_llik, start=coef(ml), 
             method="BFGS", constraints=list(ineqA=A, ineqB=B), nAGQ=5, prune=NULL)

sml = summary(ml)  
sml
# str(sml)

med0  = median(PTBCdata[PTBCdata$MULTIPAR==0,"CALVING_AGE"])
meda0 = median(PTBCdata[PTBCdata$MULTIPAR==0,"AGE"])
meda1 = median(PTBCdata[PTBCdata$MULTIPAR==1,"AGE"])

with(PTBCdata[PTBCdata$MULTIPAR==0,], hist(AGE-CALVING_AGE))
with(PTBCdata[PTBCdata$MULTIPAR==0,], median(AGE-CALVING_AGE))

min(PTBCdata[PTBCdata$MULTIPAR==0,"AGE"])
max(PTBCdata[PTBCdata$MULTIPAR==0,"AGE"])

min(PTBCdata[PTBCdata$MULTIPAR==1,"AGE"])
max(PTBCdata[PTBCdata$MULTIPAR==1,"AGE"])

1/(1+exp(-(sml$estimate[1,"Estimate"]))) # beta10
sml$estimate[2,"Estimate"] # beta11
mu1 = 1/(1+exp(-(sml$estimate[1,"Estimate"]+log(med0)*sml$estimate[2,"Estimate"]))) 
mu1

1/(1+exp(-(sml$estimate[3,"Estimate"]))) # beta20
mu2 = 1/(1+exp(-(sml$estimate[3,"Estimate"])))
mu2

exp(sml$estimate[4,"Estimate"]) # sigma1
exp(sml$estimate[5,"Estimate"]) # sigma2 

psi1 = exp(-2*sml$estimate[4,"Estimate"]) # psi1
psi2 = exp(-2*sml$estimate[5,"Estimate"]) # psi2 

# SD
sqrt(mu1*(1-mu1)/(1+psi1))
sqrt(mu2*(1-mu2)/(1+psi2))

# rPearson
f_z = sml$estimate[6,"Estimate"]
(exp(f_z)-exp(-f_z))/(exp(f_z)+exp(-f_z))

# 95% CI for rPearson
f_z_l = sml$estimate[6,"Estimate"]-qnorm(0.975)*sml$estimate[6,"Std. error"]
(exp(f_z_l)-exp(-f_z_l))/(exp(f_z_l)+exp(-f_z_l))

f_z_u = sml$estimate[6,"Estimate"]+qnorm(0.975)*sml$estimate[6,"Std. error"]
(exp(f_z_u)-exp(-f_z_u))/(exp(f_z_u)+exp(-f_z_u))

setwd("c:/_d/egyetem_ATE/cikkiras/2022/marg_multilev_models/_biostatistics/scripts/")
save(ml, file = "ml_v1_simu_calving.RData")
# load("ml_v1_simu_calving.RData")

#############################################################
# Calculate the factor of proportionality rho defined in (22)

# Reparameterized beta distribution function
mupbeta = function(
    mu,   # expected value
    psi,  # precision parameter
    q     # quantile
) {
  retv = pbeta(q, mu*psi,(1-mu)*psi)
  
  return (retv)
}

mupbeta(mu=mu1, psi=psi1, q=0.05)

# Ratio of differences of transformed conditional means to 
# differences of transformed marginal means (rho in (22))
CondMargRatio = function(
    mu, # expected value
    psi # precision parameter
) {
  
  numerator   = (-1)*grad(f=mupbeta, x0=mu, psi=psi, q=mu)
  denominator = dbeta(mu, mu*psi, (1-mu)*psi)
  
  retv = numerator/denominator
  
  return (retv)
}

rho = CondMargRatio(mu1, psi1)
rho

###########################################################################
# 3. Diagnose the quality of fit of the model using model based simulation
#
###########################################################################

# ML parameter estimates
coefs = coef(ml)

beta10  <<- coefs[1]  # intercept in (44)
beta11  <<- coefs[2]  # coefficient of log calving age in (44)
beta20  <<- coefs[3]  # intercept in (45)

lsigm1  <<- coefs[4]  # log standard deviation of random herd effects - primiparous cows
lsigm2  <<- coefs[5]  # log standard deviation of random herd effects - multiparous cows
fisherz <<- coefs[6]  # Fisher's z transformation of correlation coefficient between
                      # random herd effects related to primiparous and multiparous cows
# Dispersion model
psi1 <<- exp(-2*lsigm1)
psi2 <<- exp(-2*lsigm2)

# Correlation is expressed from Fisher's z transformation 
rPearson  <<- (exp(fisherz)-exp(-fisherz))/(exp(fisherz)+exp(-fisherz))

# Model based replications of the number of positive test results among 
# primiparous and/or multiparous cows for each herd
PTBC_Rep = function (
sel ) # 1: POS1, 2: POS2, 3: c(POS1,POS2) are returned, respectively
{
  POS1  = rep(NA_real_,LData$nHERD)
  POS2  = rep(NA_real_,LData$nHERD)
   
  for (n in 1:LData$nHERD) {
    # Number of simulated positive cases among primiparous cows
    POS1[n] = 0
    # Number of simulated positive cases among multiparous cows
    POS2[n] = 0

    # Herd specific random effect for primiparous and multiparous cows, respectively
    eta = as.numeric(rmvnorm(1, mean=c(0,0), sigma=matrix(c(1,rPearson,rPearson,1),2,2), 
                             checkSymmetry=F))
    # Primiparous cows
    for (k in (cnRec[2*n-1]+1):cnRec[2*n]) {
      # Conditional within-herd animal-level true prevalence among primiparous cows
      mu1   = 1/(1+exp(-(beta10+beta11*LData$lCAGE[k])))
      CWHP1 = qbeta.R(pnorm(eta[1]), mu1*psi1, (1-mu1)*psi1) 
      
      # Age-dependent sensitivity of the kth cow
      Se = 1/(1+exp(-(parA-parB*exp(-parC*LData$AgeG[k]))))
      
      # Apparent prevalence
      pi1     = Se*CWHP1 + (1-Sp)*(1-CWHP1)
      POS1[n] = POS1[n] + rbinom(1, LData$nCOW[k], pi1)
    }
    # Multiparous cows
    for (k in (cnRec[2*n]+1):cnRec[2*n+1]) {
      # Conditional within-herd animal-level true prevalence among multiparous cows
      mu2   = 1/(1+exp(-(beta20)))
      CWHP2 = qbeta.R(pnorm(eta[2]), mu2*psi2, (1-mu2)*psi2) 
      
      # Age-dependent sensitivity of the kth cow
      Se = 1/(1+exp(-(parA-parB*exp(-parC*LData$AgeG[k]))))
      
      # Apparent prevalence
      pi2     = Se*CWHP2 + (1-Sp)*(1-CWHP2)
      POS2[n] = POS2[n] + rbinom(1, LData$nCOW[k], pi2)
    }
  }
  
  # Return values: number of positive cases by parity and herd
  if (sel==1) {
    aRet = POS1
  } else if (sel==2) {
    aRet = POS2
  } else {
    aRet = cbind(POS1, POS2)
  }
  
  return (aRet)
}

# Simulated replications of test results
PTBCRep = PTBC_Rep( 3 )

PTBC_HD = PTBCHerdDat[order(PTBCHerdDat$MULTIPAR),]
rownames(PTBC_HD) = NULL

PTBC_HD[PTBC_HD$MULTIPAR==0,"sPOS"] = PTBCRep[,"POS1"]
PTBC_HD[PTBC_HD$MULTIPAR==1,"sPOS"] = PTBCRep[,"POS2"]

# Quantile-quantile plots are used to compare distributions 
# of the number of original and replicated positive cases for both
# primiparous and multiparous cows
qqplot(PTBC_HD[PTBC_HD$MULTIPAR==0,"POS"]+rnorm(LData$nHERD,0,0.1),
       PTBC_HD[PTBC_HD$MULTIPAR==0,"sPOS"]+rnorm(LData$nHERD,0,0.1))
abline(0,1)

qqplot(PTBC_HD[PTBC_HD$MULTIPAR==1,"POS"]+rnorm(LData$nHERD,0,0.1),
       PTBC_HD[PTBC_HD$MULTIPAR==1,"sPOS"]+rnorm(LData$nHERD,0,0.1))
abline(0,1)

###############################################################################


