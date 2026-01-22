data{
  int n;
  vector[n] Halflife;
  array[n] int Group;
  int n_Group;
  array[n] int Reference;
  int n_Reference;
}

parameters{
  // Global parameters
  real alpha_mu; // intercept on log scale
  real<lower=0> alpha_sigma_g;
  real<lower=0> alpha_sigma_r;
  real<lower=0> sigma; // likelihood standard deviation (log scale)
  
  // Partially pooled parameters
  vector[n_Group] alpha_z_g; // z-scores
  vector[n_Reference] alpha_z_r;
}

transformed parameters{
  // Convert z-scores
  vector[n_Group] alpha_g = alpha_z_g * alpha_sigma_g + alpha_mu;
  vector[n_Reference] alpha_r = alpha_z_r * alpha_sigma_r + 0;
}

model{
  // Priors
  /// Global parameters
  alpha_mu ~ normal( log(14) , 1 );
  alpha_sigma_g ~ normal( 0 , 1 ) T[0,]; // half-normal prior
  alpha_sigma_r ~ normal( 0 , 1 ) T[0,];
  sigma ~ exponential( 1 );
  
  /// Partially pooled parameters
  alpha_z_g ~ normal( 0 , 1 ); // standard normal for z-scores
  alpha_z_r ~ normal( 0 , 1 );

  // Model
  vector[n] mu = alpha_g[Group] + alpha_r[Reference];
  Halflife ~ lognormal( mu , sigma );
}