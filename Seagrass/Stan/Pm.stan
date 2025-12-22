data{
  int n;
  int n_Species;
  vector[n] Pm_mean;
  vector[n] Pm_sd;
  vector[n] Day;
  array[n] int Species;
}

parameters{
  // Estimate of Pm for measurment error
  vector[n] Pm;

  // Species parameters
  vector<lower=0>[n_Species] alpha;
  vector<lower=0>[n_Species] tau;
  vector<lower=0>[n_Species] mu;

  // Pooled parameters
  real<lower=0> alpha_mu;
  real<lower=0> tau_mu;
  real<lower=0> mu_mu;
  real<lower=0> kmu; // complete pooling on kmu

  real<lower=0> alpha_theta;
  real<lower=0> tau_theta;
  real<lower=0> mu_theta;

  // Likelihood uncertainty parameter
  real<lower=0> Pm_sigma;
}

transformed parameters{
  vector<lower=0>[n_Species] k = kmu ./ mu;
}

model{
  // Pooled priors
  alpha_mu ~ gamma( 30^2 / 5^2 , 30 / 5^2 );
  tau_mu ~ gamma( 2^2 / 1^2 , 2 / 1^2 ); // reduced sd because variation is added by theta
  mu_mu ~ gamma( 25^2 / 10^2 , 25 / 10^2 ); // estimated as 5/0.2
  kmu ~ gamma( 5^2 / 1^2 , 5 / 1^2 );

  alpha_theta ~ exponential( 1 );
  tau_theta ~ exponential( 1 );
  mu_theta ~ exponential( 1 ); 

  // Species priors
  alpha ~ gamma( alpha_mu / alpha_theta , 1 / alpha_theta );
  tau ~ gamma( tau_mu / tau_theta , 1 / tau_theta );
  mu ~ gamma( mu_mu / mu_theta , 1 / mu_theta );

  // Likelihood uncertainty prior
  Pm_sigma ~ exponential( 1 );

  // Model
  vector[n] Pm_mu = ( alpha[Species] + tau[Species] ) ./
    ( 1 + exp( k[Species] .* ( Day - mu[Species] ) ) ) 
    - tau[Species];

  // Likelihood incorporating measurement error
  Pm ~ normal( Pm_mu , Pm_sigma );
  Pm_mean ~ normal( Pm , Pm_sd );
}

generated quantities {
  real alpha_new = gamma_rng( alpha_mu / alpha_theta , 1 / alpha_theta );
  real tau_new = gamma_rng( tau_mu / tau_theta , 1 / tau_theta );
  real mu_new = gamma_rng( mu_mu / mu_theta , 1 / mu_theta );
  real k_new = kmu / mu_new;
}