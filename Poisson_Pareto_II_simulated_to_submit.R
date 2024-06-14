###############################################################################
# This supplementary R script provides methods to
# 
#   simulate repeated counts of events from Poisson-Pareto II mixture,
# 
#   fit a Bayesian model to the simulated data, as described in the main text 
#   of Lang et al. (2024),
# 
#   investigate the effect of prior distributions on the convergence and
#   accuracy of the posterior estimates of the parameters.
#
###############################################################################

# Set your working directory
setwd("...")

library(rmutil)  # This package contains [d,p,q,r]pareto

##############################################
# simulated data
# set.seed(1000)

# Number of patients
nSubj = 100 # 100, 300, 500, 1000

# Square root of the dispersion parameter
Sigma = 0.2 # 0.1, 0.2, 0.5, 1.0 

beta  = c(2.0,1.0)
tau   = 1.0

kRep  = length(beta)
kRep

# Factor of proportionality
rho = log(exp(beta[2])/(exp(beta[2])+Sigma^(-2)))-log(exp(beta[1])/(exp(beta[1])+Sigma^(-2)))
rho = rho/(log(log(1+exp(beta[2])*Sigma^2))-log(log(1+exp(beta[1])*Sigma^2))) 
rho

# Simulated Poisson-Pareto data set
ppdata = data.frame(ID = rep(1:nSubj, each=kRep),
                    IRep = rep(1:kRep, nSubj))
ppdata$events = NA_real_
# View(ppdata)

for (n in 1:nSubj) {
  ranEf = runif(1)
  for (k in 1:kRep) {
    mu  = exp(beta[k])

    # m=mu=lambda/(alpha-1) the expected value, s=alpha=(lambda+mu)/mu=1+Sigma^(-2)*mu^(-1) 
    muc = tau + qpareto(p=ranEf, m=mu, s=1+Sigma^(-2)*mu^(-1))

    ppdata[(n-1)*kRep+k,"events"] = rpois(1, muc)
  }
}

str(ppdata)
# View(ppdata)

write.table(ppdata, "ppdata.csv", sep=";", dec=".", 
            append=F, col.names=T, row.names=F)

ppdata = read.table("ppdata.csv", sep=";", dec=".", header=T, stringsAsFactors=T)
str(ppdata)

##########################################################################
# Bayesian analysis
#
##########################################################################

# List of input data
LData =
  list(
    nSubj   = nSubj,   # Number of all patients
    kRep    = kRep,    # Number of repeated observations of patients
    Events  = ppdata$events) # Number of events recorded at each observation

str(LData)
LData

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')

gc()
bayes_test = stan('Poisson_Pareto_II_sim_to_submit.stan', 
                  data=LData, chains=4, iter=20000, verbose=F,
                  control=list(adapt_delta=0.99, max_treedepth=15))

save(bayes_test, file="Poisson_Pareto_II_sim_to_submit_n100_beta3_S1_4x20000.RData")
load("Poisson_Pareto_II_sim_to_submit_n100_beta3_S1_4x20000.RData")

check_hmc_diagnostics(bayes_test)

stan_diag(bayes_test)
stan_diag(bayes_test, info="treedepth")
stan_diag(bayes_test, info="divergence")
stan_diag(bayes_test, info="stepsize")

stan_trace(bayes_test,pars="beta",inc_warmup=F) 
stan_trace(bayes_test,pars="dbeta",inc_warmup=F) 
stan_trace(bayes_test,pars="tau",inc_warmup=F) 
stan_trace(bayes_test,pars="Sigma",inc_warmup=F) 
stan_trace(bayes_test,pars="rho",inc_warmup=F) 
stan_trace(bayes_test,pars="lp__",inc_warmup=F) 

print(summary(bayes_test, pars=c("beta"))$summary,   digits=3)
print(summary(bayes_test, pars=c("dbeta"))$summary,   digits=3)
print(summary(bayes_test, pars=c("tau"))$summary,   digits=3)
print(summary(bayes_test, pars=c("Sigma"))$summary,   digits=3)
print(summary(bayes_test, pars=c("rho"))$summary,   digits=3)

pairs(bayes_test, pars="beta")
pairs(bayes_test, pars="dbeta")
pairs(bayes_test, pars=c("beta","tau","Sigma","rho"))

library(ggplot2) # This package is required to run stan_ac

stan_ac(bayes_test,pars="beta",inc_warmup=F) 
stan_ac(bayes_test,pars="dbeta",inc_warmup=F) 
stan_ac(bayes_test,pars="tau",inc_warmup=F) 
stan_ac(bayes_test,pars="Sigma",inc_warmup=F) 
stan_ac(bayes_test,pars="rho",inc_warmup=F) 
