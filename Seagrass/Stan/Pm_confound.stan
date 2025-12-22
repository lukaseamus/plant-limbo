data{
  int n;
  vector[n] Pm_mean;
  vector[n] Pm_sd;
  vector[n] O2_mean_std;
  vector[n] O2_sd_std;
  vector[n] T_mean_std;
  vector[n] T_sd_std;
  vector[n] P_mean_std;
  vector[n] S_std;
  vector[n] M_std;
}

parameters{
  // True estimates for measurement error
  vector[n] Pm;
  vector[n] O2;
  vector[n] T;

  // Intercept parameter
  real<lower=0> alpha;

  // Slope parameters
  real beta_O2;
  real beta_T;
  real beta_P;
  real beta_S;
  real beta_M;

  // Likelihood uncertainty parameter
  real<lower=0> Pm_sigma;
}

model{
  // Predictor measurement error
  O2 ~ normal( 0 , 1 );
  O2_mean_std ~ normal( O2 , O2_sd_std );
  T ~ normal( 0 , 1 );
  T_mean_std ~ normal( T , T_sd_std );

  // Intercept and slope priors
  alpha ~ gamma( 20^2 / 10^2 , 20 / 10^2 );
  beta_O2 ~ normal( 0 , 1 );
  beta_T ~ normal( 0 , 1 );
  beta_P ~ normal( 0 , 1 );
  beta_S ~ normal( 0 , 1 );
  beta_M ~ normal( 0 , 1 );

  // Likelihood uncertainty prior
  Pm_sigma ~ exponential( 1 );

  // Model
  vector[n] Pm_mu;
  for ( i in 1:n ) {
    Pm_mu[i] = alpha + beta_O2 * O2[i] + beta_T * T[i] +
    beta_P * P_mean_std[i] + beta_S * S_std[i] + beta_M * M_std[i];
  }

  // Likelihood incorporating measurement error
  Pm ~ normal( Pm_mu , Pm_sigma );
  Pm_mean ~ normal( Pm , Pm_sd );
}