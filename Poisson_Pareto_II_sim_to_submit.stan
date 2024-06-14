/* STAN code of the Poisson-Pareto II Bayesian model */
data {
  int<lower=1> nSubj; // Number of all patients
  int<lower=1> kRep;  // Number of repeated observations of patients
  int<lower=0> Events[nSubj*kRep]; // Number of events recorded at each observation
}
parameters {
  real beta[kRep];           // Effects of the repeated observations
  real<lower=0.0001> tau;    // Threshold parameter of the Pareto II distribution
  real eta[nSubj];           // Patient level random effects
  real<lower=0.0001> Sigma2; // Dispersion parameter of the patient level random effects
}
transformed parameters {
}
model {
  real mu;
  real muc;

  real t1;

  beta   ~ normal(0, 3); // Vectorized, weakly informative prior
//beta   ~ normal(0, 5); // Vectorized, weakly informative prior
  Sigma2 ~ inv_gamma(2.04, 1.04); // Weakly informative prior; mean=1, SD=5
  tau    ~ inv_gamma(2.04, 1.04); // Weakly informative prior; mean=1, SD=5
  eta    ~ normal(0, 1); // Vectorized

  for (iSubj in 1:nSubj) {
    t1  = 0; // Loglikelihood component of the iSubj-th patient

    for (iRep in 1:kRep) {

      // Marginal expected value of the number of events, 
      // minus the threshold parameter tau
      mu = exp(beta[iRep]);

      // Conditional expected value of the number of events
      // muc is the Phi(eta[iSubj])-quantile of the Pareto II distribution
      muc = tau + 1/Sigma2*((1-Phi(eta[iSubj]))^(-mu/(mu+1/Sigma2))-1);

      // Increment of the loglikelihood component      
      t1 += poisson_lpmf(Events[(iSubj-1)*kRep+iRep] | muc);    
    }
    target += t1;
  }
}
generated quantities {
  real Sigma; // Square root of the dispersion parameter
  real dbeta[kRep]; // beta[1], beta[2]-beta[1], beta[3]-beta[1], ...
  real rho; // Factor of proportionality 
  
  Sigma  = sqrt(Sigma2);
  rho    = (log(exp(beta[2])/(exp(beta[2])+Sigma^(-2)))-
           log(exp(beta[1])/(exp(beta[1])+Sigma^(-2)))) /
           (log(log(1+exp(beta[2])*Sigma^2))-log(log(1+exp(beta[1])*Sigma^2))); 
  
  for (iRep in 1:kRep) {
    if (iRep==1) {
      dbeta[1] = beta[1];
    } else {
      dbeta[iRep] = beta[iRep]-beta[1];
    }
  }
}
