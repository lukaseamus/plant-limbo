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
  vector[n_Group] alpha_g;
  vector[n_Species] alpha_s;
  vector[n_Experiment] alpha_e;
  
  /// Effects
  vector[n_Group] beta_l; // effect of light
  vector[n_Group] beta_c; // effect of chlorophyll
  
  // Likelihood standard deviation
  real<lower=0> sigma;
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
  alpha_g ~ normal( alpha_mu , alpha_sigma_g );
  alpha_s ~ normal( 0 , alpha_sigma_s );
  alpha_e ~ normal( 0 , alpha_sigma_e );
  beta_l ~ normal( beta_mu_l , beta_sigma_l );
  beta_c ~ normal( beta_mu_c , beta_sigma_c );
  
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