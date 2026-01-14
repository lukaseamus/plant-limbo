data{
  int n;
  vector[n] Day;
  vector[n] Ratio;
  vector[n] Light_effect; // dark (-0.5) vs. light (+0.5)
  vector[n] Chl_effect; // photosynthesis (-0.5) vs. chlorophyll (+0.5)
  array[n] int Group; // non-taxonomic plant group
  int n_Group;
  array[n] int Species;
  int n_Species;
  array[n] int Experiment;
  int n_Experiment;
}

parameters{
  // Global mu parameters
  /// Intercepts
  real alpha_mu;
  real<lower=0> alpha_sigma_g;
  real<lower=0> alpha_sigma_s;
  real<lower=0> alpha_sigma_e;
  
  /// Effects
  real beta_mu_l; // mean effect of light
  real beta_mu_c; // mean effect of chlorophyll
  real<lower=0> beta_sigma_l;
  real<lower=0> beta_sigma_c;
  
  // Partially pooled mu parameters
  /// Intercepts
  vector[n_Group] alpha_z_g; // z-scores
  vector[n_Species] alpha_z_s;
  vector[n_Experiment] alpha_z_e;
  
  /// Effects
  vector[n_Group] beta_z_l; // z-scores
  vector[n_Group] beta_z_c;
  
  // Likelihood standard deviation
  real<lower=0> sigma;
}

transformed parameters{
  vector[n_Group] alpha_g = alpha_z_g * alpha_sigma_g + alpha_mu;
  vector[n_Species] alpha_s = alpha_z_s * alpha_sigma_s + 0;
  vector[n_Experiment] alpha_e = alpha_z_e * alpha_sigma_e + 0;
  
  vector[n_Group] beta_l = beta_z_l * beta_sigma_l + beta_mu_l;
  vector[n_Group] beta_c = beta_z_c * beta_sigma_c + beta_mu_c;
}

model{
  // Priors
  /// Global parameters
  alpha_mu ~ normal( log(14) , 0.3 );
  beta_mu_l ~ normal( 0 , 0.4 );
  beta_mu_c ~ normal( 0 , 0.4 );
  
  alpha_sigma_g ~ normal( 0 , 0.3 ) T[0,]; // half-normal priors
  alpha_sigma_s ~ normal( 0 , 0.3 ) T[0,];
  alpha_sigma_e ~ normal( 0 , 0.3 ) T[0,];
  beta_sigma_l ~ normal( 0 , 0.4 ) T[0,];
  beta_sigma_c ~ normal( 0 , 0.4 ) T[0,];
  
  sigma ~ exponential( 1 );
  
  /// Partially pooled parameters
  alpha_z_g ~ normal( 0 , 1 ); // standard normal
  alpha_z_s ~ normal( 0 , 1 );
  alpha_z_e ~ normal( 0 , 1 );
  beta_z_l ~ normal( 0 , 1 );
  beta_z_c ~ normal( 0 , 1 );

  // Model
  /// Parameters
  vector[n] mu = exp(
      alpha_g[Group] + 
      beta_l[Group] .* Light_effect + 
      beta_c[Group] .* Chl_effect +
      alpha_s[Species] + alpha_e[Experiment]
  );
  
  /// Function
  vector[n] Ratio_mu = log_inv_logit(
    -5 / mu .* ( Day - mu )
  );
  
  // Lognormal likelihood
  Ratio ~ lognormal( Ratio_mu , sigma );
}