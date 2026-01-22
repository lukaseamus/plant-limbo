data{
  int n;
  vector[n] Pm_mean;
  vector[n] Pm_sd;
  array[n] int Species;
  int n_Species;
  vector[n] O2_mean_std;
  vector[n] O2_sd_std;
  vector[n] T_mean_std;
  vector[n] T_sd_std;
  vector[n] P_mean_std;
  vector[n] S_std;
  vector[n] M_std;
}

parameters{
  // Unobserved/true/latent variables
  vector[n] Pm;
  vector[n] O2;
  vector[n] T;

  // Intercept parameter
  vector[n_Species] alpha;

  // Slope parameters
  vector[n_Species] beta_O2;
  vector[n_Species] beta_T;
  vector[n_Species] beta_P;
  vector[n_Species] beta_S;
  vector[n_Species] beta_M;

  // Likelihood standard deviation
  vector<lower=0>[n_Species] Pm_sigma;
}

model{
  // Predictor measurement error
  O2 ~ normal( 0 , 1 ); // All confounders are standardised, hence the standard normal
  O2_mean_std ~ normal( O2 , O2_sd_std );
  T ~ normal( 0 , 1 );
  T_mean_std ~ normal( T , T_sd_std );

  // Intercept and slope priors
  alpha ~ normal( 21 , 10 );
  beta_O2 ~ normal( 0 , 1 ); // these slopes are the change per standard deviation
  beta_T ~ normal( 0 , 1 );
  beta_P ~ normal( 0 , 1 );
  beta_S ~ normal( 0 , 1 );
  beta_M ~ normal( 0 , 1 );

  // Likelihood standard deviation prior
  Pm_sigma ~ exponential( 1 );

  // Model
  vector[n] Pm_mu = alpha[Species] + beta_O2[Species] .* O2 + 
  beta_T[Species] .* T + beta_P[Species] .* P_mean_std + 
  beta_S[Species] .* S_std + beta_M[Species] .* M_std;

  // Normal likelihood with normal measurement error
  Pm ~ normal( Pm_mu , Pm_sigma[Species] );
  Pm_mean ~ normal( Pm , Pm_sd );
}