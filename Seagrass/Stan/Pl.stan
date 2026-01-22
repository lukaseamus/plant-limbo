data{
  int n;
  vector[n] Day; // Main continuous predictor
  vector[n] Pl_mean; // Response
  vector[n] Pl_sd; // Response is measured with error
  array[n] int Species; // Main categorical predictor
  int n_Species;
  vector[n] O2_mean_std; // Confounders
  vector[n] O2_sd_std; // O2 and temperature are
  vector[n] T_mean_std; // measured with error
  vector[n] T_sd_std;
  vector[n] P_mean_std;
  vector[n] S_std;
  vector[n] M_std;
}

parameters{
  // Unobserved/true/latent variables
  vector[n] Pl;
  vector[n] O2;
  vector[n] T;

  // Species parameters
  vector[n_Species] log_alpha_mu;
  vector[n_Species] log_mu_mu;
  vector[n_Species] log_tau_mu;

  // Confounder effects
  vector[n_Species] log_alpha_O2;
  vector[n_Species] log_alpha_T;
  vector[n_Species] log_alpha_P;
  vector[n_Species] log_alpha_S;
  vector[n_Species] log_alpha_M;
  vector[n_Species] log_mu_O2;
  vector[n_Species] log_mu_T;
  vector[n_Species] log_mu_P;
  vector[n_Species] log_mu_S;
  vector[n_Species] log_mu_M;
  vector[n_Species] log_tau_O2;
  vector[n_Species] log_tau_T;
  vector[n_Species] log_tau_P;
  vector[n_Species] log_tau_S;
  vector[n_Species] log_tau_M;

  // Likelihood standard deviation
  vector<lower=0>[n_Species] Pl_sigma;
}

model{
  // Species priors
  log_alpha_mu ~ normal( log(7) , 1 );
  log_mu_mu ~ normal( log(14) , 0.6 );
  log_tau_mu ~ normal( log(1) , 0.4 );
  
  // Confounder priors
  O2 ~ normal( 0 , 1 ); // All confounders are standardised, hence the standard normal
  O2_mean_std ~ normal( O2 , O2_sd_std ); // Normal measurement error
  T ~ normal( 0 , 1 );
  T_mean_std ~ normal( T , T_sd_std );
  
  log_alpha_O2 ~ normal( 0 , 1 ); // these slopes are the change per standard deviation
  log_alpha_T ~ normal( 0 , 1 ); 
  log_alpha_P ~ normal( 0 , 1 );
  log_alpha_S ~ normal( 0 , 1 );
  log_alpha_M ~ normal( 0 , 1 );
  log_mu_O2 ~ normal( 0 , 1 );
  log_mu_T ~ normal( 0 , 1 );
  log_mu_P ~ normal( 0 , 1 );
  log_mu_S ~ normal( 0 , 1 );
  log_mu_M ~ normal( 0 , 1 );
  log_tau_O2 ~ normal( 0 , 1 );
  log_tau_T ~ normal( 0 , 1 );
  log_tau_P ~ normal( 0 , 1 );
  log_tau_S ~ normal( 0 , 1 );
  log_tau_M ~ normal( 0 , 1 );
  
  // Likelihood standard deviation prior
  Pl_sigma ~ exponential( 1 );

  // Model
  /// Parameters
  vector[n] alpha = exp(
      log_alpha_mu[Species] + 
      log_alpha_O2[Species] .* O2 + 
      log_alpha_T[Species] .* T + 
      log_alpha_P[Species] .* P_mean_std + 
      log_alpha_S[Species] .* S_std + 
      log_alpha_M[Species] .* M_std
  );
  
  vector[n] mu = exp(
      log_mu_mu[Species] + 
      log_mu_O2[Species] .* O2 + 
      log_mu_T[Species] .* T + 
      log_mu_P[Species] .* P_mean_std + 
      log_mu_S[Species] .* S_std + 
      log_mu_M[Species] .* M_std
  );
  
  vector[n] tau = exp(
      log_tau_mu[Species] + 
      log_tau_O2[Species] .* O2 + 
      log_tau_T[Species] .* T + 
      log_tau_P[Species] .* P_mean_std + 
      log_tau_S[Species] .* S_std + 
      log_tau_M[Species] .* M_std
  );
  
  /// Function
  vector[n] Pl_mu = ( alpha + tau ) .* inv_logit( -5 / mu .* ( Day - mu ) ) - tau;
  
  // Normal likelihood with normal measurement error
  Pl ~ normal( Pl_mu , Pl_sigma[Species] );
  Pl_mean ~ normal( Pl , Pl_sd );
}