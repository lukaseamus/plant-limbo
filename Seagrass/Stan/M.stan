data{
  int n;
  int n_Species;
  vector<lower=0>[n] Mass;
  vector[n] Day_c;
  array[n] int Species;
}

parameters{
  vector<lower=0>[n_Species] alpha;
  vector[n_Species] beta;
  real<lower=0> sigma;
}

model{
  // Priors
  alpha ~ gamma( 0.52^2 / 0.25^2, 0.52 / 0.25^2 ); // reparameterised with mean and sd
  beta ~ normal( 0, 0.01 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = alpha[Species[i]] + beta[Species[i]] * Day_c[i];
  }

  // Likelihood
  Mass ~ normal( mu , sigma );
}