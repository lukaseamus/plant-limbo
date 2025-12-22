data{
  int n;
  vector[n] Ts;
}

parameters{
  real Tsmu;
  real<lower=0> sigma;
}

model{
  // Priors
  Tsmu ~ normal( 0, 1 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Tsmu;
  }

  // Likelihood
  Ts ~ normal( mu , sigma );
}