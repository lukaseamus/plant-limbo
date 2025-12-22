data{
  int n;
  int n_Species;
  vector[n] Pl_mean;
  vector[n] Pl_sd;
  vector[n] Day;
  array[n] int Species;
}

parameters{
  // Estimate of Pl for measurment error
  vector[n] Pl;

  // Species parameters
  vector<lower=0>[n_Species] alpha; // no pooling on alpha
  vector<lower=0>[n_Species] tau; // no pooling on tau
  vector<lower=0>[n_Species] mu;

  // Pooled parameters
  real<lower=0> mu_mu;
  real<lower=0> kmu; // complete pooling on kmu

  real<lower=0> mu_theta;

  // Likelihood uncertainty parameter
  real<lower=0> Pl_sigma;
}

transformed parameters{
  vector<lower=0>[n_Species] k = kmu ./ mu;
}

model{
  // Pooled priors
  mu_mu ~ gamma( 25^2 / 10^2 , 25 / 10^2 ); // estimated as 5 / 0.2
  kmu ~ gamma( 5^2 / 1^2 , 5 / 1^2 );

  mu_theta ~ exponential( 1 ); 

  // Species priors
  alpha ~ gamma( 16^2 / 4^2 , 16 / 4^2 );
  tau ~ gamma( 1^2 / 0.8^2 , 1 / 0.8^2 );
  mu ~ gamma( mu_mu / mu_theta , 1 / mu_theta );

  // Likelihood uncertainty prior
  Pl_sigma ~ exponential( 1 );

  // Model
  vector[n] Pl_mu = ( alpha[Species] + tau[Species] ) ./
    ( 1 + exp( k[Species] .* ( Day - mu[Species] ) ) ) 
    - tau[Species];

  // Likelihood incorporating measurement error
  Pl ~ normal( Pl_mu , Pl_sigma );
  Pl_mean ~ normal( Pl , Pl_sd );
}

generated quantities {
  real mu_new = gamma_rng( mu_mu / mu_theta , 1 / mu_theta );
  real k_new = kmu / mu_new;
}