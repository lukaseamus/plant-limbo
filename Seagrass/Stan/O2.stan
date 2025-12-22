data{
  int n;
  vector<lower=0>[n] Value;
  vector[n] delta_t_c;
}

parameters{
  real<lower=0> alpha;
  real beta;
  real<lower=0> sigma;
}

model{
  // Priors (from Rose et al. 2012, doi 10.5194/os-8-545-2012 and Woo & Pattiaratchi 2008, doi 10.1016/j.dsr.2008.05.005)
  alpha ~ gamma( 227.9194^2 / 20^2, 227.9194 / 20^2 ); // reparameterised with mean and sd
  beta ~ normal( 0, 5 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = alpha + beta * delta_t_c[i];
  }

  // Likelihood
  Value ~ normal( mu , sigma );
}