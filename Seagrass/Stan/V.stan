data{
  int n;
  vector[n] Vs;
}

parameters{
  real Vsmu;
  real<lower=0> sigma;
}

model{
  // Priors
  Vsmu ~ normal( 0, 1 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Vsmu;
  }

  // Likelihood
  Vs ~ normal( mu , sigma );
}