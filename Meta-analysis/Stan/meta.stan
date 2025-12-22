data{
  int n;
  int n_group;
  vector[n] Proportion;
  vector[n] Day;
  array[n] int group;
}

parameters{
  // Group parameters
  vector<lower=0>[n_group] k;
  vector<lower=0>[n_group] mu;

  // Likelihood uncertainty parameter
  real<lower=0> P_sigma;
}

model{
  // Group priors
  k ~ normal( 5 , 4 ) T[0,];
  mu ~ normal( 365 , 300 ) T[0,];
  
  // Constraint on k and mu (joint prior)
  target += normal_lpdf( log(k) + log(mu) | log(5) , 0.01 );
  
  // Likelihood uncertainty prior
  P_sigma ~ exponential( 1 );

  // Model
  vector[n] P_mu = 1 / ( 1 + exp( k[group] .* ( Day - mu[group] ) ) );

  // Likelihood
  Proportion ~ normal( P_mu , P_sigma );
}