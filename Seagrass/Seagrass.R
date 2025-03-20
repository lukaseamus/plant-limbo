# Load oxygen files
directory <- "~/Desktop/Projects/Seagrass/Data/Oxygen"
files <- list.files(path = directory, pattern = "\\.csv$", full.names = TRUE)

require(tidyverse)
O2 <- files %>%
  map(~ read.csv(., skip = 1, header = TRUE) %>%
        drop_na(Value) %>%
        mutate(delta_t = as.numeric(delta_t),
               delta_t_c = delta_t - mean(delta_t))) %>%
  set_names(str_remove(basename(files), "\\.csv$") %>% make.names)
str(O2$X220915_B1)

require(tidybayes)
O2_list <- O2 %>%
  map(~ select(., Value, delta_t_c) %>%
        compose_data())

O2_stan <- "
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
"

require(cmdstanr)
O2_mod <- cmdstan_model(stan_file = write_stan_file(code = O2_stan))

O2_samples <- O2_list %>%
  map(~ O2_mod$sample(data = .,
                      seed = 100,
                      chains = 8,
                      parallel_chains = parallel::detectCores(),
                      iter_warmup = 1e4,
                      iter_sampling = 1e4))

# check Rhat, effective sample size and chains
O2_summary <- O2_samples %>%
  map(~ .$summary())

O2_summary %>%
  bind_rows() %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

O2_summary %>%
  bind_rows() %>%
  summarise(rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# decent Rhat and effective sample size

O2_draws <- O2_samples %>%
  map(~ .$draws(format = "df"))

require(bayesplot)
require(patchwork)
O2_draws %>% 
  map(~ mcmc_rank_overlay(.)) %>%
  wrap_plots() %>%
  ggsave(filename = "rank.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")

O2_draws %>% 
  sample(6) %>%
  map(~ mcmc_pairs(., pars = c("alpha", "beta"))) %>%
  wrap_plots() %>%
  ggsave(filename = "pairs.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 50, height = 50, units = "cm")


# Prior-posterior comparison
O2_prior_posterior <- O2_draws %>%
  map(~ spread_draws(., alpha, beta, sigma) %>%
        ungroup() %>%
        mutate(
          sigma_prior = rexp(length(.draw), 1), 
          beta_prior = rnorm(length(.draw), 0, 5),
          alpha_prior = rgamma(length(.draw), 227.9194^2 / 20^2, 227.9194 / 20^2)
        )
      )

O2_prior_posterior %>%
  map(~ pivot_longer(., cols = c("beta", "beta_prior", "alpha", "alpha_prior"),
                        names_to = c("parameter", "distribution"),
                        names_sep = "_", # this will produce NAs and throw a warning message 
                        values_to = "samples") %>%
        mutate(parameter = fct(parameter),
               distribution = fct(ifelse(is.na(distribution), # here the NAs are dealt with
                                         "posterior", distribution))) %>%
        ggplot(aes(samples, fill = distribution)) +
          facet_wrap(~ parameter, scales = "free", nrow = 1) +
          geom_density(colour = NA, alpha = 0.5) +
          theme_minimal() +
          theme(legend.position = "top", legend.justification = 0)) %>%
  wrap_plots() %>%
  ggsave(filename = "ppcomp.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")

O2_mu <- O2_prior_posterior %>% 
  map2(O2,
       ~ select(.x, alpha, beta) %>%
         expand_grid(delta_t_c = seq(.y %>% pull(delta_t_c) %>% min(),
                                     .y %>% pull(delta_t_c) %>% max(),
                                     length.out = 50)) %>%
         mutate(mu = alpha + beta * delta_t_c) %>%
         select(mu, delta_t_c)
       )

O2_mu_summary <- O2_mu %>%
  map(~ group_by(., delta_t_c) %>%
        mean_qi(mu, .width = c(.5, .8, .9))
      )
         
O2 %>%
  map2(O2_mu_summary,
       ~ ggplot() +
          geom_point(data = .x, aes(delta_t_c, Value)) +
          geom_line(data = .y, aes(delta_t_c, mu)) +
          geom_ribbon(data = .y, aes(delta_t_c, ymin = .lower, ymax = .upper,
                                     alpha = factor(.width))) +
          scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
          theme_minimal()) %>%
  imap(~ .x + ggtitle(.y)) %>%
  wrap_plots() %>%
  ggsave(filename = "fit.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")
# fit looks good

O2_initial <- O2_mu %>%
  map(~ slice_min(., delta_t_c))

O2_initial %>%
  map(~ ggplot(., aes(mu)) +
          geom_density() +
          theme_minimal()) %>%
  wrap_plots() %>%
  ggsave(filename = "initial.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")


V <- read.csv("~/Desktop/Projects/Seagrass/Data/Volume.csv") %>%
  mutate(Vs = (Volume - mean(Volume)) / sd(Volume)) # standardise Volume

V_stan <- "
data{
  int n;
  vector[n] Vs;
}

parameters{
  real Vsmu;
  real<lower=0> sigma;
}

model{
  // Priors
  Vsmu ~ normal( 0, 1 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Vsmu;
  }

  // Likelihood
  Vs ~ normal( mu , sigma );
}
"

V_mod <- cmdstan_model(stan_file = write_stan_file(code = V_stan))

V_samples <- V_mod$sample(data = compose_data(V),
                          seed = 100,
                          chains = 8,
                          parallel_chains = parallel::detectCores(),
                          iter_warmup = 1e4,
                          iter_sampling = 1e4)

# check Rhat, effective sample size and chains
V_summary <- V_samples$summary()

V_summary %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

V_draws <- V_samples$draws(format = "df")

V_draws %>% mcmc_rank_overlay()

V_draws %>% mcmc_pairs()

# Prior-posterior comparison
V_prior_posterior <- V_draws %>%
  spread_draws(Vsmu, sigma) %>%
  ungroup() %>%
  mutate(
    sigma_prior = rexp(length(.draw), 1), 
    Vsmu_prior = rnorm(length(.draw), 0, 1)
    )

V_prior_posterior %>%
  pivot_longer(., cols = c("Vsmu", "Vsmu_prior", "sigma", "sigma_prior"),
                     names_to = c("parameter", "distribution"),
                     names_sep = "_", # this will produce NAs and throw a warning message 
                     values_to = "samples") %>%
  mutate(parameter = fct(parameter),
         distribution = fct(ifelse(is.na(distribution), # here the NAs are dealt with
                                   "posterior", distribution))) %>%
  ggplot(aes(samples, fill = distribution)) +
    facet_wrap(~ parameter, scales = "free", nrow = 1) +
    geom_density(colour = NA, alpha = 0.5) +
    theme_minimal() +
    theme(legend.position = "top", legend.justification = 0)

require(magrittr)
V_prior_posterior %<>%
  mutate(Vmu = V %$% (Vsmu * sd(Volume) + mean(Volume))) %>%
  select(Vmu)


O2 %>%
  map(~ summarise(., Tmean = mean(Temp),
                  Tsd = sd(Temp),
                  Pmean = mean(Pressure),
                  Psd = sd(Pressure))) %>%
  bind_rows() %>%
  filter(Psd == 0 | Tsd == 0)
# in several cases pressure did not vary at all across the incubation,
# so measurement error in pressure can not be included in the model

O2 %<>%
  map(~ mutate(., Ts = ( Temp - mean(Temp) ) / sd(Temp)))

T_list <- O2 %>%
  map(~ select(., Ts) %>%
        compose_data())

T_stan <- "
data{
  int n;
  vector[n] Ts;
}

parameters{
  real Tsmu;
  real<lower=0> sigma;
}

model{
  // Priors
  Tsmu ~ normal( 0, 1 );
  sigma ~ exponential( 1 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = Tsmu;
  }

  // Likelihood
  Ts ~ normal( mu , sigma );
}
"

T_mod <- cmdstan_model(stan_file = write_stan_file(code = T_stan))

T_samples <- T_list %>%
  map(~ T_mod$sample(data = .,
                     seed = 100,
                     chains = 8,
                     parallel_chains = parallel::detectCores(),
                     iter_warmup = 1e4,
                     iter_sampling = 1e4))

# check Rhat, effective sample size and chains
T_summary <- T_samples %>%
  map(~ .$summary())

T_summary %>%
  bind_rows() %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

T_summary %>%
  bind_rows() %>%
  summarise(rhat_mean = mean(rhat),
            rhat_sd = sd(rhat),
            ess_mean = mean(ess_bulk),
            ess_sd = sd(ess_bulk))
# decent Rhat and effective sample size

T_draws <- T_samples %>%
  map(~ .$draws(format = "df"))

T_draws %>% 
  map(~ mcmc_rank_overlay(.)) %>%
  wrap_plots() %>%
  ggsave(filename = "Trank.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")

T_draws %>% 
  sample(6) %>%
  map(~ mcmc_pairs(., pars = c("Tsmu", "sigma"))) %>%
  wrap_plots() %>%
  ggsave(filename = "Tpairs.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 50, height = 50, units = "cm")


# Prior-posterior comparison
T_prior_posterior <- T_draws %>%
  map(~ spread_draws(., Tsmu, sigma) %>%
        ungroup() %>%
        mutate(
          sigma_prior = rexp(length(.draw), 1), 
          Tsmu_prior = rnorm(length(.draw), 0, 1)
        )
  )

T_prior_posterior %>%
  map(~ pivot_longer(., cols = c("Tsmu", "Tsmu_prior", "sigma", "sigma_prior"),
                     names_to = c("parameter", "distribution"),
                     names_sep = "_", # this will produce NAs and throw a warning message 
                     values_to = "samples") %>%
        mutate(parameter = fct(parameter),
               distribution = fct(ifelse(is.na(distribution), # here the NAs are dealt with
                                         "posterior", distribution))) %>%
        ggplot(aes(samples, fill = distribution)) +
        facet_wrap(~ parameter, scales = "free", nrow = 1) +
        geom_density(colour = NA, alpha = 0.5) +
        theme_minimal() +
        theme(legend.position = "top", legend.justification = 0)) %>%
  wrap_plots() %>%
  ggsave(filename = "Tppcomp.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 100, height = 100, units = "cm")

T_prior_posterior %<>%
  map2(O2,
       ~ mutate(.x, Tmu = .y %$% (Tsmu * sd(Temp) + mean(Temp))) %>%
         select(Tmu)
       )

# combine posterior distributions into one tibble

P <- O2_prior_posterior %>%
  map2(O2_initial,
       ~ bind_cols(.x, .y)) %>%
  map2(T_prior_posterior,
       ~ bind_cols(.x, .y, V_prior_posterior)) %>%
  map2(O2 %>% 
         map(~ summarise(., Salinity = mean(Salinity),
                         Pressure = mean(Pressure))),
       ~ bind_cols(.x, .y)) %>%
  map2(O2 %>% 
         map(~ select(., Date, Time) %>%
               slice(1)),
       ~ bind_cols(.x, .y)) %>%
  imap(~ mutate(.x, ID = .y)) %>%
  bind_rows()

# clean up, build explanatory variables from ID

P %<>%
  select(.draw, ID, Date, Time, beta, mu, Tmu, Vmu, Salinity, Pressure) %>%
  rename(O2 = mu, Temperature = Tmu, Volume = Vmu) %>%
  mutate(Date = Date %>% mdy(),
         Time = Time %>% hms(),
         Species = case_when(
           str_extract(ID, "(?<=_)\\p{L}(?=\\d*)") == "B" ~ "Blank",
           str_extract(ID, "(?<=_)\\p{L}(?=\\d*)") == "A" ~ "Amphibolis antarctica",
           str_extract(ID, "(?<=_)\\p{L}(?=\\d*)") == "H" ~ "Halophila ovalis"
           ),
         Leaf = str_extract(ID, "(?<=_[A-Za-z])\\d") %>% as.numeric())

# load and explore mass data
M <- read.csv("~/Desktop/Projects/Seagrass/Data/Mass.csv") %>%
  mutate(Species = fct_relevel(Species, "Halophila ovalis"),
         Date = Date %>% dmy()) %>% # convert Date to a usable date
  group_by(Species) %>%
  mutate(Day = min(Date) %--% Date %>%
           time_length("day")) %>%
  ungroup() %>% # very important
  mutate(Day_c = Day - mean(Day)) # centre Day

M %>%
  ggplot(aes(Day, Mass, colour = Species)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_minimal()

mean(M$Mass) # 0.5195062 g

M_prior <- 
  tibble(n = 1:1e3,
         alpha = rgamma(n = 1e3, shape = 0.52^2 / 0.25^2, rate = 0.52 / 0.25^2),
         beta = rnorm(n = 1e3, mean = 0, sd = 0.01)) %>% 
  expand_grid(Day = M %$% seq(min(Day_c), max(Day_c), length.out = 50)) %>%
  mutate(mu = alpha + beta * Day)

M_prior %>%
  ggplot(aes(x = Day, y = mu, group = n)) +
  geom_hline(yintercept = M %$% c(min(Mass), max(Mass))) +
  geom_line(alpha = 0.05) +
  coord_cartesian(expand = F, clip = "off") +
  theme_minimal()
# prior simulation looks good


M_stan <- "
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
"

M_mod <- cmdstan_model(stan_file = write_stan_file(code = M_stan))

M_samples <- M_mod$sample(data = M %>%
                            select(Mass, Day_c, Species) %>%
                            compose_data(),
                          seed = 100,
                          chains = 8,
                          parallel_chains = parallel::detectCores(),
                          iter_warmup = 1e4,
                          iter_sampling = 1e4)


M_summary <- M_samples$summary()
M_summary %>%
  filter(rhat > 1.001)

M_draws <- M_samples$draws(format = "df")
M_draws %>% mcmc_rank_overlay() # chains look good
M_draws %>% mcmc_pairs(pars = c("alpha[1]", "beta[1]"))
M_draws %>% mcmc_pairs(pars = c("alpha[2]", "beta[2]"))

# intercept and slope estimates for the text
M_samples %>%
  recover_types(M %>% select(Species)) %>%
  spread_draws(alpha[Species], beta[Species]) %>%
  mutate(alpha = if_else(Species == "Halophila ovalis",
                         alpha / 10, alpha),
         beta = if_else(Species == "Halophila ovalis",
                        beta / 10, beta)) %>%
  mutate(alpha_nc = alpha + beta * -mean(M$Day), # calculate alpha for non-centred predictor
         beta_p = beta / alpha_nc * 100, # calculate beta in % d^-1
         beta_mg = beta * 1e3) %>% # calculate beta in mg leaf^-1 d^-1
  group_by(Species) %>%
  summarise(alpha_mu = mean(alpha_nc),
            alpha_sd = sd(alpha_nc),
            beta_mu = mean(beta_mg),
            beta_sd = sd(beta_mg),
            beta_p_mu = mean(beta_p),
            beta_p_sd = sd(beta_p),
            P_beta = mean(beta_mg < 0))
  

# correct photosynthesis for blank, volume and mass (µM min^-1 to µmol g^-1 h^-1)
P %<>%
  filter(Species != "Blank") %>%
  left_join(P %>% filter(Species == "Blank") %>%
              select(.draw, Date, Time, beta) %>%
              rename(blank = beta),
            by = c(".draw", "Date", "Time"),
            relationship = "many-to-one") %>%
  left_join(M %>% mutate(Date = Date %>% dmy()), 
            by = c("Date", "Species", "Leaf"),
            relationship = "many-to-one") %>%
  mutate(Pm = ( beta - blank ) * Volume * 60 / Mass,
         Pl = if_else(Species == "Halophila ovalis", 
                      ( beta - blank ) * Volume * 60 / 10, # there were 10 H. ovalis leaves per sample
                      ( beta - blank ) * Volume * 60)) %>%
  select(-c(beta, Volume, blank))

# calculate detrital age (days) from dates
P %<>%
  group_by(Species) %>%
  mutate(Day = min(Date) %--% Date %>%
           time_length("day"))

# summarise P for modelling purposes
P_summary <- P %>%
  group_by(ID, Day, Species, Leaf) %>%
  summarise(Pm_mean = mean(Pm),
            Pm_sd = sd(Pm),
            Pl_mean = mean(Pl),
            Pl_sd = sd(Pl),
            O2_mean = mean(O2),
            O2_sd = sd(O2),
            T_mean = mean(Temperature),
            T_sd = sd(Temperature),
            P_mean = mean(Pressure),
            S = mean(Salinity),
            M = mean(Mass),
            n = length(.draw)) %>%
  ungroup() %>%
  mutate(Species = fct_relevel(Species, "Halophila ovalis"))

P_summary %>%
  mutate(Incubation = if_else(nchar(ID) == 12, "Second", "First")) %>%
  ggplot() +
    geom_pointrange(aes(Day, Pm_mean, 
                        ymin = Pm_mean - Pm_sd, 
                        ymax = Pm_mean + Pm_sd,
                        colour = Species,
                        shape = Incubation)) +
    geom_smooth(aes(Day, Pm_mean, colour = Species),
                se = F) +
    theme_minimal()

P_summary %>%
  filter(nchar(ID) == 10) %>%
  ggplot() +
  geom_pointrange(aes(Day, Pm_mean, 
                      ymin = Pm_mean - Pm_sd, 
                      ymax = Pm_mean + Pm_sd,
                      colour = Species)) +
  geom_smooth(aes(Day, Pm_mean, colour = Species),
              se = F) +
  theme_minimal()

P_summary %>%
  mutate(Incubation = if_else(nchar(ID) == 12, "Second", "First")) %>%
  ggplot() +
    geom_pointrange(aes(Day, Pl_mean, 
                        ymin = Pl_mean - Pl_sd, 
                        ymax = Pl_mean + Pl_sd,
                        colour = Species,
                        shape = Incubation)) +
    geom_smooth(aes(Day, Pl_mean, colour = Species),
                se = F) +
    facet_wrap(~Species, scales = "free") +
    theme_minimal()

P_summary %>%
  filter(nchar(ID) == 10) %>%
  ggplot() +
    geom_pointrange(aes(Day, Pl_mean, 
                        ymin = Pl_mean - Pl_sd, 
                        ymax = Pl_mean + Pl_sd,
                        colour = Species)) +
    geom_smooth(aes(Day, Pl_mean, colour = Species),
                se = F) +
    facet_wrap(~Species, scales = "free") +
    theme_minimal()

# clearly the second incubation skews the data at 14 d, perhaps because of the higher initial O2
# -> proceed with only the first incubation
P %<>%
  filter(nchar(ID) == 10)
P_summary %<>%
  filter(nchar(ID) == 10)

# clean up
rm(O2, O2_draws, O2_initial, O2_list, O2_mod, O2_mu, O2_mu_summary,
   O2_prior_posterior, O2_samples, O2_summary, T_draws, T_list, T_mod,
   T_prior_posterior, T_samples, T_summary, V, V_draws, V_mod, V_samples,
   V_prior_posterior, V_summary, directory, files, O2_stan, T_stan, V_stan)

# save data
write.csv(M, file = "~/Desktop/M.csv", row.names = F)
write.csv(P, file = "~/Desktop/P.csv", row.names = F)
write.csv(P_summary, file = "~/Desktop/P_summary.csv", row.names = F)

# Confounders
# create standardised explanatory variables
P_summary %<>%
  mutate(O2_mean_std = ( O2_mean - mean(O2_mean) ) / sd(O2_mean),
         O2_sd_std = O2_sd / sd(O2_mean),
         T_mean_std = ( T_mean - mean(T_mean) ) / sd(T_mean),
         T_sd_std = T_sd / sd(T_mean),
         P_mean_std = ( P_mean - mean(P_mean) ) / sd(P_mean),
         S_std = ( S - mean(S) ) / sd(S),
         M_std = ( M - mean(M) ) / sd(M))

# Prior simulation
Prior <- read.csv("~/Desktop/Projects/Seagrass/Data/Prior.csv")

# calculate as species-specific means and sds
Prior %>%
  left_join(M %>% select(Species, Mass) %>%
              mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                        Mass / 10, Mass)) %>%
              group_by(Species) %>%
              summarise(Leafmass = mean(Leafmass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Fl = Flux * Leafmass) %>%
  group_by(Species, Variable) %>%
  summarise(Fm_mean = mean(Flux),
            Fm_sd = sd(Flux),
            Fm_median = median(Flux),
            Fl_mean = mean(Fl),
            Fl_sd = sd(Fl),
            Fl_median = median(Fl),
            n = length(Flux))

# calculate overall mean and sd 
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  left_join(M %>% select(Species, Mass) %>%
              mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                        Mass / 10, Mass)) %>%
              group_by(Species) %>%
              summarise(Leafmass = mean(Leafmass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Fl = Flux * Leafmass) %>%
  summarise(Pm_mean = mean(Flux),
            Pm_sd = sd(Flux),
            Pm_median = median(Flux),
            Pl_mean = mean(Fl),
            Pl_sd = sd(Fl),
            Pl_median = median(Fl),
            n = length(Flux))


# based on the mean-median comparison photosynthesis and 
# detrital respiration are right-skewed

# prior simulation for mass-based estimates
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  ggplot() +
    geom_density(aes(Flux), fill = "red", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 35.2173, colour = "red") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 35.2173^2 / 10^2, rate = 35.2173 / 10^2)), aes(x)) + # arbitrary sd
    theme_minimal()
# shift mean to capture the most probable peak
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  ggplot() +
    geom_density(aes(Flux), fill = "red", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 35.2173, colour = "red") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 20^2 / 10^2, rate = 20 / 10^2)), aes(x)) +
    theme_minimal()
# looks better and 20 is closer to median (20.67893)


# prior for mass-based confounder model
Pm_prior <- 
  tibble(n = 1:1e3,
         alpha = rgamma(n = 1e3, shape = 20^2 / 10^2, rate = 20 / 10^2), # Pm at mean of predictor is positive
         beta = rnorm(n = 1e3, mean = 0, sd = 1)) %>% 
  expand_grid(Predictor = seq(-3, 3, length.out = 50)) %>%
  mutate(P_mu = alpha + beta * Predictor)

Pm_prior %>%
  ggplot(aes(x = Predictor, y = P_mu, group = n)) +
  geom_hline(yintercept = P_summary %$% c(min(Pm_mean), 0, max(Pm_mean))) +
  geom_line(alpha = 0.05) +
  coord_cartesian(expand = F, clip = "off") +
  theme_minimal()
# covers all probable scenarios

Pm_confound_stan <- "
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
"

Pm_confound_mod <- cmdstan_model(stan_file = write_stan_file(code = Pm_confound_stan))

Pm_confound_samples <- Pm_confound_mod$sample(data = P_summary %>%
                                                select(Pm_mean, Pm_sd, O2_mean_std, O2_sd_std,
                                                       T_mean_std, T_sd_std, P_mean_std, S_std, M_std) %>%
                                                compose_data(),
                                              seed = 100,
                                              chains = 8,
                                              parallel_chains = parallel::detectCores(),
                                              iter_warmup = 1e4,
                                              iter_sampling = 1e4)


Pm_confound_summary <- Pm_confound_samples$summary()
Pm_confound_summary %>%
  filter(rhat > 1.001)

Pm_confound_draws <- Pm_confound_samples$draws(format = "df")

Pm_confound_draws %>% mcmc_rank_overlay() # chains look good

Pm_confound_draws %>% mcmc_pairs(pars = c("alpha", "beta_O2"))
Pm_confound_draws %>% mcmc_pairs(pars = c("alpha", "beta_T"))
Pm_confound_draws %>% mcmc_pairs(pars = c("alpha", "beta_P"))
Pm_confound_draws %>% mcmc_pairs(pars = c("alpha", "beta_S"))
Pm_confound_draws %>% mcmc_pairs(pars = c("alpha", "beta_M"))

# Prior-posterior comparison
Pm_confound_posterior <- Pm_confound_samples %>%
  gather_draws(alpha, beta_O2, beta_T, beta_P, beta_S, beta_M) %>%
  ungroup() %>%
  mutate(Distribution = "Posterior")

Pm_confound_prior <- tibble(.chain = 1:8 %>% rep(each = 1e4),
                            .iteration = 1:1e4 %>% rep(8),
                            .draw = 1:8e4,
                            alpha = rgamma(8e4, 20^2 / 10^2 , 20 / 10^2),
                            beta_O2 = rnorm(8e4, 0 , 1),
                            beta_T = rnorm(8e4, 0 , 1),
                            beta_P = rnorm(8e4, 0 , 1),
                            beta_S = rnorm(8e4, 0 , 1),
                            beta_M = rnorm(8e4, 0 , 1),
                            Distribution = "Prior") %>%
  pivot_longer(cols = starts_with(c("alpha", "beta")), values_to = ".value", names_to = ".variable")


Pm_confound_prior_posterior <- Pm_confound_posterior %>%
  bind_rows(Pm_confound_prior) %>%
  mutate(.value = P_summary %$% case_when( # reverse standardisation
    .variable == "beta_O2" ~ .value / sd(O2_mean), # this typically involves multiplying by sd
    .variable == "beta_T" ~ .value / sd(T_mean), # but because this is a slope, it's inverse
    .variable == "beta_P" ~ .value / sd(P_mean),
    .variable == "beta_S" ~ .value / sd(S),
    .variable == "beta_M" ~ .value / sd(M),
    .variable == "alpha" ~ .value
    ))

Pm_confound_prior_posterior %>%
  ggplot(aes(.value, fill = Distribution)) +
    facet_wrap(~ .variable, scales = "free") +
    geom_density(colour = NA, alpha = 0.5) +
    theme_minimal()
# looks like priors were reasonable

# calculate P_mu from parameters
Pm_confound_mu <- Pm_confound_samples %>%
  spread_draws(alpha, beta_O2, beta_T, beta_P, beta_S, beta_M) %>%
  ungroup() %>%
  pivot_longer(cols = c(beta_O2, beta_T, beta_P, beta_S, beta_M),
               names_to = ".variable", values_to = "beta", names_prefix = "beta_") %>%
  mutate(Predictor = P_summary %$% case_when(
    .variable == "O2" ~ list( seq(min(O2_mean_std), max(O2_mean_std), length.out = 100) ),
    .variable == "T" ~ list( seq(min(T_mean_std), max(T_mean_std), length.out = 100) ),
    .variable == "P" ~ list( seq(min(P_mean_std), max(P_mean_std), length.out = 100) ),
    .variable == "S" ~ list( seq(min(S_std), max(S_std), length.out = 100) ),
    .variable == "M" ~ list( seq(min(M_std), max(M_std), length.out = 100) )
  )) %>%
  unnest(Predictor) %>% # expand the list column Predictor
  mutate(P_mu = alpha + beta * Predictor) # since predictor variables were standardised,
# excluding the other slopes from the equation is equivalent to setting the other predictors
# to their mean and multiplying by their slope since mean = 0 for standardised variables

Pm_confound_mu_summary <- Pm_confound_mu %>%
  group_by(.variable, Predictor) %>%
  mean_qi(P_mu, .width = c(.5, .8, .9)) %>%
  mutate(Predictor_original = P_summary %$% case_when(
    .variable == "O2" ~ Predictor * sd(O2_mean) + mean(O2_mean), # reverse standardisation
    .variable == "T" ~ Predictor * sd(T_mean) + mean(T_mean),
    .variable == "P" ~ Predictor * sd(P_mean) + mean(P_mean),
    .variable == "S" ~ Predictor * sd(S) + mean(S),
    .variable == "M" ~ Predictor * sd(M) + mean(M)
  ))

rm(Pm_confound_draws, Pm_confound_mu, Pm_confound_posterior, Pm_confound_prior, 
   Pm_confound_summary, Pm_confound_mod, Pm_confound_samples)

mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(0, 0.5, 0.2, 0, unit = "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 12, hjust = 0),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black", lineend = "square"),
                 legend.key = element_blank(),
                 legend.key.width = unit(.25, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.spacing.x = unit(.5, "cm"),
                 legend.key.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.position = "top",
                 legend.justification = 0,
                 legend.text = element_text(size = 12, hjust = 0),
                 legend.title = element_blank(),
                 legend.margin = margin(0, 0, 0, 0, unit = "cm"),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 12, hjust = 0),
                 panel.spacing = unit(0.6, "cm"),
                 text = element_text(family = "Futura"))

require(ggh4x)
Fig_S1_beta <- Pm_confound_prior_posterior %>%
               filter(.variable != "alpha") %>%
               mutate(.variable = fct_relevel(.variable, "beta_O2", "beta_T", "beta_P", "beta_S"),
                      Distribution = fct_relevel(Distribution, "Prior")) %>%
                  ggplot(aes(.value, alpha = Distribution)) +
                    geom_vline(xintercept = 0) +
                    geom_density(fill = "black", colour = NA, adjust = 1.8) +
                    scale_alpha_manual(values = c(0.2, 0.5)) +
                    facet_wrap(~ .variable, scales = "free", nrow = 1,
                               labeller = labeller(.variable = as_labeller(
                                                      c("beta_O2" = "italic('b')['O'[2]]*' (µM'^-1*')'",
                                                        "beta_T" = "italic('b')['T']*' (°C'^-1*')'",
                                                        "beta_P" = "italic('b')['Pressure']*' (hPa'^-1*')'",
                                                        "beta_S" = "italic('b')['Salinity']*' (‰'^-1*')'",
                                                        "beta_M" = "italic('b')['Mass']*' (g'^-1*')'"),
                                                      label_parsed))) +
                    facetted_pos_scales(x = list(
                      .variable == "beta_O2" ~ scale_x_continuous(limits = c(-0.5, 1),
                                                                  oob = scales::oob_keep,
                                                                  breaks = seq(-0.5, 1, by = 0.5),
                                                                  labels = scales::label_number(accuracy = c(0.1, 1, 0.1, 1),
                                                                                                style_negative = "minus")),
                      .variable == "beta_T" ~ scale_x_continuous(limits = c(-4, 4),
                                                                 oob = scales::oob_keep,
                                                                 breaks = seq(-4, 4, by = 4),
                                                                 labels = scales::label_number(style_negative = "minus")),
                      .variable == "beta_P" ~ scale_x_continuous(limits = c(-0.5, 1),
                                                                 oob = scales::oob_keep,
                                                                 breaks = seq(-0.5, 1, by = 0.5),
                                                                 labels = scales::label_number(accuracy = c(0.1, 1, 0.1, 1),
                                                                                               style_negative = "minus")),
                      .variable == "beta_S" ~ scale_x_continuous(limits = c(-6, 6),
                                                                 oob = scales::oob_keep,
                                                                 breaks = seq(-6, 6, by = 6),
                                                                 labels = scales::label_number(style_negative = "minus")),
                      .variable == "beta_M" ~ scale_x_continuous(limits = c(-12, 12),
                                                                 oob = scales::oob_keep,
                                                                 breaks = seq(-12, 12, by = 12),
                                                                 labels = scales::label_number(style_negative = "minus"))
                    )) +
                    coord_cartesian(expand = F, clip = "off") +
                    mytheme +
                    theme(axis.line.y = element_blank(),
                          axis.title = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.text.y = element_blank(),
                          legend.position = c(-0.1, 0.72), # ignore warning since legend.position.inside will not deliver
                          plot.margin = margin(0, 0.5, 0.5, 0, unit = "cm"))

require(ggdensity)
Fig_S1_O2 <- ggplot() +
                geom_hdr(data = P %>% group_by(ID) %>% mutate(O2 = sample(O2)), # randomly reshuffle within ID
                         aes(O2, Pm, group = ID), n = 500, method = "mvnorm", probs = 0.999,
                         position = "identity") +
                geom_line(data = Pm_confound_mu_summary %>%
                            filter(.variable == "O2"),
                          aes(Predictor_original, P_mu)) +
                geom_ribbon(data = Pm_confound_mu_summary %>%
                              filter(.variable == "O2"),
                            aes(Predictor_original, ymin = .lower,
                                ymax = .upper, alpha = factor(.width)),
                            colour = NA) +
                scale_alpha_manual(values = c(0.2, 0.4, 0.3, 0.2), guide = "none") + # first alpha is for geom_hdr
                scale_y_continuous(breaks = seq(-10, 50, 10), 
                                   labels = scales::label_number(style_negative = "minus")) +
                scale_x_continuous(breaks = seq(210, 250, 20)) +
                labs(x = expression("Initial O"[2]*" (µM)"),
                     y = expression(italic(P)["max"]*" (µmol O"[2]*" g"^-1*" h"^-1*")")) +
                coord_cartesian(xlim = c(210, 250), ylim = c(-10, 50), 
                                expand = F, clip = "off") +
                mytheme +
                theme(plot.margin = margin(0, 0.6, 0.2, 0, unit = "cm"))

Fig_S1_T <- ggplot() +
                geom_hdr(data = P %>% group_by(ID) %>% mutate(Temperature = sample(Temperature)),
                         aes(Temperature, Pm, group = ID), n = 500, method = "mvnorm", probs = 0.999,
                         position = "identity") +
                geom_line(data = Pm_confound_mu_summary %>%
                            filter(.variable == "T"),
                          aes(Predictor_original, P_mu)) +
                geom_ribbon(data = Pm_confound_mu_summary %>%
                              filter(.variable == "T"),
                            aes(Predictor_original, ymin = .lower,
                                ymax = .upper, alpha = factor(.width)),
                            colour = NA) +
                scale_alpha_manual(values = c(0.2, 0.4, 0.3, 0.2), guide = "none") +
                scale_y_continuous(breaks = seq(-10, 50, 10), 
                                   labels = scales::label_number(style_negative = "minus")) +
                scale_x_continuous(breaks = seq(16, 20, 2)) +
                labs(x = "T (°C)") +
                coord_cartesian(xlim = c(16, 20), ylim = c(-10, 50), 
                                expand = F, clip = "off") +
                mytheme +
                theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      plot.margin = margin(0, 0.6, 0.2, 0, unit = "cm"))

Fig_S1_P <- ggplot() +
                geom_violin(data = P,
                            aes(Pressure, Pm, group = ID), position = "identity",
                            alpha = 0.2, width = 6, colour = NA, fill = "black") +
                geom_line(data = Pm_confound_mu_summary %>%
                            filter(.variable == "P"),
                          aes(Predictor_original, P_mu)) +
                geom_ribbon(data = Pm_confound_mu_summary %>%
                              filter(.variable == "P"),
                            aes(Predictor_original, ymin = .lower,
                                ymax = .upper, alpha = factor(.width)),
                            colour = NA) +
                scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                scale_y_continuous(breaks = seq(-10, 50, 10), 
                                   labels = scales::label_number(style_negative = "minus")) +
                scale_x_continuous(breaks = seq(1000, 1030, 15)) +
                labs(x = "Pressure (hPa)") +
                coord_cartesian(xlim = c(1000, 1030), ylim = c(-10, 50),
                                expand = F, clip = "off") +
                mytheme +
                theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      plot.margin = margin(0, 0.6, 0.2, 0, unit = "cm"))

Fig_S1_S <- ggplot() +
                geom_violin(data = P,
                            aes(Salinity, Pm, group = ID), position = "identity",
                            alpha = 0.2, width = 0.4, colour = NA, fill = "black") +
                geom_line(data = Pm_confound_mu_summary %>%
                            filter(.variable == "S"),
                          aes(Predictor_original, P_mu)) +
                geom_ribbon(data = Pm_confound_mu_summary %>%
                              filter(.variable == "S"),
                            aes(Predictor_original, ymin = .lower,
                                ymax = .upper, alpha = factor(.width)),
                            colour = NA) +
                scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                scale_y_continuous(breaks = seq(-10, 50, 10), 
                                   labels = scales::label_number(style_negative = "minus")) +
                scale_x_continuous(breaks = seq(34, 36, 1)) +
                labs(x = "Salinity (‰)") +
                coord_cartesian(xlim = c(34, 36), ylim = c(-10, 50),
                                expand = F, clip = "off") +
                mytheme +
                theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      plot.margin = margin(0, 0.6, 0.2, 0, unit = "cm"))

Fig_S1_M <- ggplot() +
                geom_violin(data = P,
                            aes(Mass, Pm, group = ID), position = "identity",
                            alpha = 0.2, width = 0.3, colour = NA, fill = "black") +
                geom_line(data = Pm_confound_mu_summary %>%
                            filter(.variable == "M"),
                          aes(Predictor_original, P_mu)) +
                geom_ribbon(data = Pm_confound_mu_summary %>%
                              filter(.variable == "M"),
                            aes(Predictor_original, ymin = .lower,
                                ymax = .upper, alpha = factor(.width)),
                            colour = NA) +
                scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                scale_y_continuous(breaks = seq(-10, 50, 10), 
                                   labels = scales::label_number(style_negative = "minus")) +
                scale_x_continuous(breaks = seq(0, 1.5, 0.5),
                                   labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1))) +
                labs(x = "Mass (g)") +
                coord_cartesian(xlim = c(0, 1.5), ylim = c(-10, 50),
                                expand = F, clip = "off") +
                mytheme +
                theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      plot.margin = margin(0, 0.6, 0.2, 0, unit = "cm"))

require(patchwork)
Fig_S1 <- Fig_S1_beta / ( Fig_S1_O2 | Fig_S1_T | Fig_S1_P | Fig_S1_S | Fig_S1_M ) +
  plot_layout(heights = c(0.4, 1))

ggsave(plot = Fig_S1, filename = "Fig_S1.pdf", device = cairo_pdf, path = "~/Desktop",
       width = 20, height = 10, units = "cm")

rm(Fig_S1, Fig_S1_beta, Fig_S1_O2, Fig_S1_T, Fig_S1_P, Fig_S1_S, Fig_S1_M,
   Pm_confound_mu_summary, Pm_confound_prior_posterior, Pm_confound_stan)

# Pm model
# Prior simulation
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  left_join(M %>% select(Species, Mass) %>%
              mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                        Mass / 10, Mass)) %>%
              group_by(Species) %>%
              summarise(Leafmass = mean(Leafmass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Fl = Flux * Leafmass) %>%
  summarise(Pm_mean = mean(Flux),
            Pm_sd = sd(Flux),
            Pm_median = median(Flux),
            Pl_mean = mean(Fl),
            Pl_sd = sd(Fl),
            Pl_median = median(Fl),
            n = length(Flux))

Prior %>%
  filter(Variable == "Detrital respiration") %>%
  mutate(Fl = Flux * M %>% 
           mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                     Mass / 10, Mass)) %>%
           pull(Leafmass) %>% mean()) %>%
  summarise(Rm_mean = mean(Flux),
            Rm_sd = sd(Flux),
            Rm_median = median(Flux),
            Rl_mean = mean(Fl),
            Rl_sd = sd(Fl),
            Rl_median = median(Fl),
            n = length(Flux))

# I focussed on the low mode for average Pmax across detrital age
# alpha is Pmax at day 0 so needs to be higher, so something in between mean and median should be right
# I also need to reduce sd for this more complicated model to run 
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  ggplot() +
    geom_density(aes(Flux), fill = "red", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 35.2173, colour = "red") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 30^2 / 5^2, rate = 30 / 5^2)), aes(x)) +
    theme_minimal()

Prior %>%
  filter(Variable == "Detrital respiration") %>%
  ggplot() +
    geom_density(aes(Flux), fill = "red", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 5.355247, colour = "red") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 5.355247^2 / 5^2, rate = 5.355247 / 5^2)), aes(x)) + # arbitrary sd
    theme_minimal()
# reduce mean and sd to capture peak
Prior %>%
  filter(Variable == "Detrital respiration") %>%
  ggplot() +
    geom_density(aes(Flux), fill = "red", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 5.355247, colour = "red") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 2^2 / 1.5^2, rate = 2 / 1.5^2)), aes(x)) +
    theme_minimal()

Pm_prior <-
  tibble(n = 1:1e3,
         alpha = rgamma(n = 1e3, shape = 30^2 / 5^2, rate = 30 / 5^2),
         tau = rgamma(n = 1e3, shape = 2^2 / 1.5^2, rate = 2 / 1.5^2),
         k = rgamma(n = 1e3, shape = 0.2^2 / 0.1^2, rate = 0.2 / 0.1^2), # see Fig. 2b in doi 10.1093/aob/mcad167
         mu = rgamma(n = 1e3, shape = 17^2 / 10^2, rate = 17 / 10^2)) %>% # half the experimental duration
  filter(k * mu > 4) %>% # simulates constraining k and mu with a joint prior
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 50)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau)

Pm_prior %>%
  ggplot(aes(x = Day, y = P_mu, group = n)) +
    geom_hline(yintercept = P_summary %$% c(min(Pm_mean), 0, max(Pm_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()

# prior simulation will also need to be done with the joint prior in Stan
Pm_prior_stan <- "
parameters{
  real<lower=0> alpha;
  real<lower=0> tau;
  real<lower=0> k;
  real<lower=0> mu;
}

model{
  alpha ~ gamma( 30^2 / 5^2 , 30 / 5^2 );
  tau ~ gamma( 2^2 / 1.5^2 , 2 / 1.5^2 );
  k ~ gamma( 0.2^2 / 0.1^2 , 0.2 / 0.1^2 );
  mu ~ gamma( 17^2 / 10^2 , 17 / 10^2 );
  target += gamma_lpdf( k * mu | 4^2 / 1^2 , 4 / 1^2 );
}
"

Pm_prior_mod <- cmdstan_model(stan_file = write_stan_file(code = Pm_prior_stan))
Pm_prior_samples <- Pm_prior_mod$sample(data = list(), # no data to condition on
                                        seed = 100,
                                        chains = 8,
                                        parallel_chains = parallel::detectCores(),
                                        iter_warmup = 1e4,
                                        iter_sampling = 1e4)

Pm_prior_draws <- Pm_prior_samples$draws(format = "df")

Pm_prior_draws %>%
  ggplot() +
  geom_density(aes(alpha), fill = "black", alpha = 0.5) +
  geom_density(aes(rgamma( 8e4, 30^2 / 5^2 , 30 / 5^2 )),
               fill = "red", alpha = 0.5) +
  theme_minimal()
# unchanged

Pm_prior_draws %>%
  ggplot() +
  geom_density(aes(tau), fill = "black", alpha = 0.5) +
  geom_density(aes(rgamma( 8e4, 2^2 / 1.5^2 , 2 / 1.5^2 )),
               fill = "red", alpha = 0.5) +
  theme_minimal()
# unchanged

Pm_prior_draws %>%
  ggplot() +
    geom_density(aes(k), fill = "black", alpha = 0.5) +
    geom_density(aes(rgamma( 8e4, 0.2^2 / 0.1^2 , 0.2 / 0.1^2 )),
                 fill = "red", alpha = 0.5) +
    scale_x_continuous(limits = c(0, 0.5), oob = scales::oob_keep) +
    theme_minimal()
# shifted by joint prior

Pm_prior_draws %>%
  ggplot() +
    geom_density(aes(mu), fill = "black", alpha = 0.5) +
    geom_density(aes(rgamma( 8e4, 17^2 / 10^2 , 17 / 10^2 )),
                 fill = "red", alpha = 0.5) +
    scale_x_continuous(limits = c(0, 50), oob = scales::oob_keep) +
    theme_minimal()
# shifted by joint prior

Pm_prior_draws %>%
  slice_sample(n = 1e3) %>% # subsample
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 50)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau) %>%
  ggplot(aes(x = Day, y = P_mu, group = .draw)) +
    geom_hline(yintercept = P_summary %$% c(min(Pm_mean), 0, max(Pm_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()

Pm_stan <- "
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
  vector<lower=0>[n_Species] k;
  vector<lower=0>[n_Species] mu;

  // Pooled parameters
  real<lower=0> alpha_mu;
  real<lower=0> tau_mu;
  real<lower=0> k_mu;
  real<lower=0> mu_mu;

  real<lower=0> alpha_sigma;
  real<lower=0> tau_sigma;
  real<lower=0> k_sigma;
  real<lower=0> mu_sigma;

  // Likelihood uncertainty parameter
  real<lower=0> Pm_sigma;
}

model{
  // Pooled priors
  alpha_mu ~ gamma( 30^2 / 5^2 , 30 / 5^2 );
  tau_mu ~ gamma( 2^2 / 1.5^2 , 2 / 1.5^2 );
  k_mu ~ gamma( 0.2^2 / 0.1^2 , 0.2 / 0.1^2 );
  mu_mu ~ gamma( 17^2 / 10^2 , 17 / 10^2 );

  alpha_sigma ~ exponential( 1 );
  tau_sigma ~ exponential( 1 );
  k_sigma ~ exponential( 1 );
  mu_sigma ~ exponential( 1 );

  // Constraint on k_mu and mu_mu (joint prior)
  target += gamma_lpdf( k_mu * mu_mu | 4^2 / 1^2 , 4 / 1^2 );

  // Species priors
  alpha ~ gamma( alpha_mu^2 / alpha_sigma^2 , alpha_mu / alpha_sigma^2 );
  tau ~ gamma( tau_mu^2 / tau_sigma^2 , tau_mu / tau_sigma^2 );
  k ~ gamma( k_mu^2 / k_sigma^2 , k_mu / k_sigma^2 );
  mu ~ gamma( mu_mu^2 / mu_sigma^2 , mu_mu / mu_sigma^2 );

  // Constraint on k and mu (joint prior)
  target += gamma_lpdf( k .* mu | 4^2 / 1^2 , 4 / 1^2 );

  // Likelihood uncertainty prior
  Pm_sigma ~ exponential( 1 );

  // Model
  vector[n] Pm_mu;
  for ( i in 1:n ) {
    Pm_mu[i] =
    ( alpha[Species[i]] + tau[Species[i]] ) /
    ( 1 + exp( k[Species[i]] * ( Day[i] - mu[Species[i]] ) ) )
    - tau[Species[i]];
  }

  // Likelihood incorporating measurement error
  Pm ~ normal( Pm_mu , Pm_sigma );
  Pm_mean ~ normal( Pm , Pm_sd );
}
"
Pm_mod <- cmdstan_model(stan_file = write_stan_file(code = Pm_stan))

Pm_samples <- Pm_mod$sample(data = P_summary %>%
                              select(Pm_mean, Pm_sd, Day, Species) %>%
                              compose_data(),
                            seed = 100,
                            chains = 8,
                            parallel_chains = parallel::detectCores(),
                            iter_warmup = 1e4,
                            iter_sampling = 1e4)
# the lowest percentage of divergent transitions attained is 14%

Pm_summary <- Pm_samples$summary()
Pm_summary %>%
  filter(rhat > 1.001) # 16 rhat above 1.001
Pm_summary %>%
  filter(rhat > 1.01) # none above 1.01

Pm_draws <- Pm_samples$draws(format = "df")

Pm_draws %>% mcmc_rank_overlay() # chains look fine

Pm_draws %>% mcmc_pairs(pars = c("alpha[1]", "tau[1]", "k[1]", "mu[1]"))
Pm_draws %>% mcmc_pairs(pars = c("alpha[2]", "tau[2]", "k[2]", "mu[2]"))
# this is chosen as the optimal model


# Prior-posterior comparison
Pm_posterior <- Pm_samples %>%
  recover_types(P_summary %>% select(Species)) %>%
  gather_draws(alpha[Species], tau[Species], k[Species], mu[Species],
               alpha_mu, tau_mu, k_mu, mu_mu) %>%
  ungroup() %>%
  mutate(Species = if_else(is.na(Species), "Seagrasses", Species) %>%
                    fct_relevel("Seagrasses", "Halophila ovalis"),
         .variable = case_when(
                      .variable %in% c("alpha", "alpha_mu") ~ "alpha",
                      .variable %in% c("tau", "tau_mu") ~ "tau",
                      .variable %in% c("k", "k_mu") ~ "k",
                      .variable %in% c("mu", "mu_mu") ~ "mu"
                      ) %>% fct_relevel("alpha", "tau", "k"))

# Combine with prior so that Distribution and Species are merged into one grouping variable
Pm_prior_posterior <- Pm_posterior %>%
  bind_rows(Pm_prior_samples %>%
              gather_draws(alpha, tau, k, mu) %>%
              ungroup() %>%
              mutate(Species = NA)) %>%
  mutate(Group = if_else(is.na(Species), "Prior", Species) %>%
           fct_relevel("Seagrasses", "Amphibolis antarctica", "Halophila ovalis"),
         .variable = fct_relevel(.variable, "alpha", "tau", "k"))

require(ggridges)
Fig_1a_top <- Pm_prior_posterior %>%
  ggplot() +
    geom_density_ridges(aes(.value, Group, colour = Group, fill = Group),
                        quantile_lines = TRUE, quantiles = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95),
                        alpha = 0.5, scale = 2, rel_min_height = 0.001, bandwidth = c(0.2, 0.2, 0.005, 0.5),
                        from = c(20, rep(0, 3)), to = c(50, 9, 0.6, 42)) +
    scale_colour_manual(values = c("#363538", "#4a7518", "#bdd268", "#b5b8ba"),
                        labels = c("Seagrasses", expression(italic("Amphibolis antarctica")),
                                   expression(italic("Halophila ovalis")), "Prior"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_fill_manual(values = c("#363538", "#4a7518", "#bdd268", "#b5b8ba"),
                      labels = c("Seagrasses", expression(italic("Amphibolis antarctica")),
                                 expression(italic("Halophila ovalis")), "Prior"),
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(~ .variable, scales = "free_x",
               labeller = labeller(.variable = as_labeller(
                 c("alpha" = "italic('a')*' (µmol g'^-1*' h'^-1*')'",
                   "tau" = "italic('t')*' (µmol g'^-1*' h'^-1*')'",
                   "k" = "italic('k')*' (d'^-1*')'",
                   "mu" = "italic('µ')*' (d)'"),
                 label_parsed))
               ) +
    facetted_pos_scales(x = list(
      .variable == "alpha" ~ scale_x_continuous(limits = c(20, 50),
                                                breaks = seq(20, 50, by = 10)),
      .variable == "tau" ~ scale_x_continuous(limits = c(0, 9),
                                              breaks = seq(0, 9, by = 3)),
      .variable == "k" ~ scale_x_continuous(limits = c(0, 0.6),
                                            breaks = seq(0, 0.6, by = 0.2),
                                            labels = scales::label_number(accuracy = c(1, 0.1, 0.1, 0.1))),
      .variable == "mu" ~ scale_x_continuous(limits = c(0, 42),
                                             breaks = seq(0, 42, by = 14))
    )) +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.spacing = unit(1, "cm"))

# Parameters and differences
# Parameters for Table 1
Pm_posterior %>%
  group_by(Species, .variable) %>%
  summarise(mean = mean(.value),
            sd = sd(.value),
            n = length(.value)) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(mean, digits = 2), signif(mean, digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        )
        ) %>%
  pivot_wider(names_from = Species, values_from = rounded) %>%
  write.csv("~/Desktop/Pm_para.csv", row.names = FALSE)

# Calculate differences
Pm_diff <- Pm_posterior %>%
  filter(Species != "Seagrasses") %>%
  pivot_wider(names_from = c(Species, .variable), values_from = .value) %>%
  mutate(delta_alpha = `Halophila ovalis_alpha` - `Amphibolis antarctica_alpha`,
         delta_tau = `Halophila ovalis_tau` - `Amphibolis antarctica_tau`,
         delta_k = `Halophila ovalis_k` - `Amphibolis antarctica_k`,
         delta_mu = `Halophila ovalis_mu` - `Amphibolis antarctica_mu`) %>%
  select(.chain, .iteration, .draw, delta_alpha, delta_tau, delta_k, delta_mu) %>%
  pivot_longer(cols = c(delta_alpha, delta_tau, delta_k, delta_mu),
               names_to = "Parameter", values_to = "Difference", names_prefix = "delta_") %>%
  mutate(Parameter = fct_relevel(Parameter, "alpha", "tau", "k"))

# Differences and probabilities for Table 1
Pm_diff %>%
  group_by(Parameter) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P_more = mean(Difference > 0),
            P_less = mean(Difference < 0),
            n = length(Difference)) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(abs(mean), digits = 2), signif(abs(mean), digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        ),
         P = signif(pmax(P_less, P_more), digits = 2)) %>%
  write.csv("~/Desktop/Pm_diff.csv", row.names = FALSE)


# calculate P_mu from parameters
Pm_mu <- Pm_prior_posterior %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau)

Pm_mu_summary <- Pm_mu %>%
  group_by(Day, Group) %>%
  median_qi(P_mu, .width = c(.5, .8, .9))

Fig_1a_bottom <- ggplot() +
                  geom_hline(yintercept = 0) +
                  geom_violin(data = P,
                              aes(Day, Pm, fill = Species,
                                  colour = Species, group = ID),
                              alpha = 0.2, position = "identity", width = 2) +
                  # too much overplotting, so I reduced the prior to the 0.05 and 0.95 quantiles
                  geom_ribbon(data = Pm_mu_summary %>% filter(Group == "Prior", .width == 0.9),
                              aes(Day, ymin = .lower, ymax = .upper), fill = NA, colour = "#b5b8ba") +
                  geom_line(data = Pm_mu_summary %>% filter(Group != "Prior"),
                            aes(Day, P_mu, colour = fct_relevel(Group, "Halophila ovalis",
                                                                "Amphibolis antarctica"))) +
                  geom_ribbon(data = Pm_mu_summary %>% filter(Group != "Prior"),
                              aes(Day, ymin = .lower, ymax = .upper,
                                  fill = fct_relevel(Group, "Halophila ovalis",
                                                     "Amphibolis antarctica"),
                              alpha = factor(.width)), colour = NA) +
                  labs(y = expression(italic(P)["max"]*" (µmol O"[2]*" g"^-1*" h"^-1*")"),
                       x = "Detrital age (d)") +
                  scale_colour_manual(values = c("#bdd268", "#4a7518", "#363538"),
                                      guide = "none") +
                  scale_fill_manual(values = c("#bdd268", "#4a7518", "#363538"),
                                    guide = "none") +
                  scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
                  scale_y_continuous(breaks = seq(-10, 50, by = 10),
                                     labels = scales::label_number(style_negative = "minus")) +
                  coord_cartesian(xlim = c(-0.8, 35), ylim = c(-10, 50), expand = FALSE, clip = "off") +
                  mytheme

# leaf-based estimates
# based on Pm I know alpha needs to be constrained for posteriors to converge on one mode
# for this to be possible for leaf-based estimates, alpha needs to be made as similar as
# possible, by adjusting the number of H. ovalis leaves used in the model
# calculate overall mean and sd 

M %>% group_by(Species) %>%
  summarise(Mass = mean(Mass)) %>%
  pivot_wider(names_from = Species, values_from = Mass) %>%
  mutate(Ratio = `Halophila ovalis` / `Amphibolis antarctica`)
# the sample mass of 10 H. ovalis leaves was on average only 
# 71.5% of that of the single A. antarctica leaf
# this means one A. antarctica leaf is equivalent in mass to
10 / 0.715
# about 14 H. ovalis leaves
# however, this is similar enough to analyse on a sample basis

Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  left_join(M %>% select(Species, Mass) %>%
              mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                        Mass / 10, Mass)) %>%
              group_by(Species) %>%
              summarise(Leafmass = mean(Leafmass),
                        Samplemass = mean(Mass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Fl = Flux * Leafmass,
         Fs = Flux * Samplemass) %>%
  summarise(Pl_mean = mean(Fl),
            Pl_sd = sd(Fl),
            Pl_median = median(Fl),
            Ps_mean = mean(Fs),
            Ps_sd = sd(Fs),
            Ps_median = median(Fs),
            n = length(Flux))

Prior %>%
  filter(Variable == "Detrital respiration") %>%
  mutate(Fl = Flux * M %>% 
           mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                     Mass / 10, Mass)) %>%
           pull(Leafmass) %>% mean(),
         Fs = Flux * M %>%  pull(Mass) %>% mean()) %>%
  summarise(Rl_mean = mean(Fl),
            Rl_sd = sd(Fl),
            Rl_median = median(Fl),
            Rs_mean = mean(Fs),
            Rs_sd = sd(Fs),
            Rs_median = median(Fs),
            n = length(Flux))

# prior simulation for sample-based estimates
Prior %>%
  filter(Variable == "Light-saturated net photosynthesis") %>%
  left_join(M %>% select(Species, Mass) %>%
              group_by(Species) %>%
              summarise(Samplemass = mean(Mass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Flux = Flux * Samplemass) %>%
  ggplot() +
    geom_density(aes(Flux), fill = "red", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 15.85307, colour = "red") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 16^2 / 4^2, rate = 16 / 4^2)), aes(x)) +
    theme_minimal()

# focus on peak
Prior %>%
  filter(Variable == "Detrital respiration") %>%
  mutate(Flux = Flux * M %>% pull(Mass) %>% mean()) %>%
  ggplot() +
    geom_density(aes(Flux), fill = "red", alpha = 0.5, colour = NA) +
    geom_vline(xintercept = 2.782084, colour = "red") +
    geom_density(data = tibble(x = rgamma(n = 1e3, shape = 1^2 / 0.8^2, rate = 1 / 0.8^2)), aes(x)) + 
    theme_minimal()


Pl_prior <-
  tibble(n = 1:1e3,
         alpha = rgamma(n = 1e3, shape = 16^2 / 4^2, rate = 16 / 4^2),
         tau = rgamma(n = 1e3, shape = 1^2 / 0.8^2, rate = 1 / 0.8^2),
         k = rgamma(n = 1e3, shape = 0.2^2 / 0.1^2, rate = 0.2 / 0.1^2), # same as before
         mu = rgamma(n = 1e3, shape = 17^2 / 10^2, rate = 17 / 10^2)) %>% # same as before
  filter(k * mu > 4) %>% # same as before
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 50)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau)

Pl_prior %>%
  ggplot(aes(x = Day, y = P_mu, group = n)) +
    geom_hline(yintercept = P_summary %$% c(min(Pl_mean), 0, max(Pl_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()

# prior simulation will also need to be done with the joint prior in Stan
Pl_prior_stan <- "
parameters{
  real<lower=0> alpha;
  real<lower=0> tau;
  real<lower=0> k;
  real<lower=0> mu;
}

model{
  alpha ~ gamma( 16^2 / 4^2 , 16 / 4^2 );
  tau ~ gamma( 1^2 / 0.8^2 , 1 / 0.8^2 );
  k ~ gamma( 0.2^2 / 0.1^2 , 0.2 / 0.1^2 );
  mu ~ gamma( 17^2 / 10^2 , 17 / 10^2 );
  target += gamma_lpdf( k * mu | 4^2 / 1^2 , 4 / 1^2 );
}
"

Pl_prior_mod <- cmdstan_model(stan_file = write_stan_file(code = Pl_prior_stan))
Pl_prior_samples <- Pl_prior_mod$sample(data = list(), # no data to condition on
                                        seed = 100,
                                        chains = 8,
                                        parallel_chains = parallel::detectCores(),
                                        iter_warmup = 1e4,
                                        iter_sampling = 1e4)

Pl_prior_draws <- Pl_prior_samples$draws(format = "df")

Pl_prior_draws %>%
  ggplot() +
  geom_density(aes(alpha), fill = "black", alpha = 0.5) +
  geom_density(aes(rgamma( 8e4, 16^2 / 4^2 , 16 / 4^2 )),
               fill = "red", alpha = 0.5) +
  theme_minimal()
# unchanged

Pl_prior_draws %>%
  ggplot() +
  geom_density(aes(tau), fill = "black", alpha = 0.5) +
  geom_density(aes(rgamma( 8e4, 1^2 / 0.8^2 , 1 / 0.8^2 )),
               fill = "red", alpha = 0.5) +
  theme_minimal()
# unchanged

Pl_prior_draws %>%
  ggplot() +
    geom_density(aes(k), fill = "black", alpha = 0.5) +
    geom_density(aes(rgamma( 8e4, 0.2^2 / 0.1^2 , 0.2 / 0.1^2 )),
                 fill = "red", alpha = 0.5) +
    scale_x_continuous(limits = c(0, 0.5), oob = scales::oob_keep) +
    theme_minimal()
# shifted by joint prior

Pl_prior_draws %>%
  ggplot() +
    geom_density(aes(mu), fill = "black", alpha = 0.5) +
    geom_density(aes(rgamma( 8e4, 17^2 / 10^2 , 17 / 10^2 )),
                 fill = "red", alpha = 0.5) +
    scale_x_continuous(limits = c(0, 50), oob = scales::oob_keep) +
    theme_minimal()
# shifted by joint prior

Pl_prior_draws %>%
  slice_sample(n = 1e3) %>% # subsample
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 50)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau) %>%
  ggplot(aes(x = Day, y = P_mu, group = .draw)) +
    geom_hline(yintercept = P_summary %$% c(min(Pl_mean), 0, max(Pl_mean))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()

Pl_stan <- "
data{
  int n;
  int n_Species;
  vector[n] Pl_mean;
  vector[n] Pl_sd;
  vector[n] Day;
  array[n] int Species;
}

parameters{
  // Estimate of Pm for measurment error
  vector[n] Pl;

  // Species parameters
  vector<lower=0>[n_Species] alpha;
  vector<lower=0>[n_Species] tau;
  vector<lower=0>[n_Species] k;
  vector<lower=0>[n_Species] mu;

  // Pooled parameters
  real<lower=0> k_mu;
  real<lower=0> mu_mu;

  real<lower=0> k_sigma;
  real<lower=0> mu_sigma;

  // Likelihood uncertainty parameter
  real<lower=0> Pl_sigma;
}

model{
  // Pooled priors
  k_mu ~ gamma( 0.2^2 / 0.1^2 , 0.2 / 0.1^2 );
  mu_mu ~ gamma( 17^2 / 10^2 , 17 / 10^2 );

  k_sigma ~ exponential( 1 );
  mu_sigma ~ exponential( 1 );

  // Constraint on k_mu and mu_mu (joint prior)
  target += gamma_lpdf( k_mu * mu_mu | 4^2 / 1^2 , 4 / 1^2 );

  // Species priors
  alpha ~ gamma( 16^2 / 4^2 , 16 / 4^2 );
  tau ~ gamma( 1^2 / 0.8^2 , 1 / 0.8^2 );
  k ~ gamma( k_mu^2 / k_sigma^2 , k_mu / k_sigma^2 );
  mu ~ gamma( mu_mu^2 / mu_sigma^2 , mu_mu / mu_sigma^2 );

  // Constraint on k and mu (joint prior)
  target += gamma_lpdf( k .* mu | 4^2 / 1^2 , 4 / 1^2 );

  // Likelihood uncertainty prior
  Pl_sigma ~ exponential( 1 );

  // Model
  vector[n] Pl_mu;
  for ( i in 1:n ) {
    Pl_mu[i] =
    ( alpha[Species[i]] + tau[Species[i]] ) /
    ( 1 + exp( k[Species[i]] * ( Day[i] - mu[Species[i]] ) ) )
    - tau[Species[i]];
  }

  // Likelihood incorporating measurement error
  Pl ~ normal( Pl_mu , Pl_sigma );
  Pl_mean ~ normal( Pl , Pl_sd );
}
"
Pl_mod <- cmdstan_model(stan_file = write_stan_file(code = Pl_stan))

Pl_samples <- Pl_mod$sample(data = P_summary %>%
                              select(Pl_mean, Pl_sd, Day, Species) %>%
                              # revert H. ovalis leaf-based estimates to sample mass
                              # to make them comparable to A. antarctica estimates
                              mutate(Pl_mean = if_else(Species == "Halophila ovalis",
                                                       Pl_mean * 10, Pl_mean),
                                     Pl_sd = if_else(Species == "Halophila ovalis",
                                                     Pl_sd * 10, Pl_sd)) %>%
                              compose_data(),
                            seed = 100,
                            chains = 8,
                            parallel_chains = parallel::detectCores(),
                            iter_warmup = 1e4,
                            iter_sampling = 1e4)
# the lowest percentage of divergent transitions attained is 5%

Pl_summary <- Pl_samples$summary()
Pl_summary %>%
  filter(rhat > 1.001) # one rhat above 1.001
Pl_summary %>%
  filter(rhat > 1.01) # none above 1.001

Pl_draws <- Pl_samples$draws(format = "df")

Pl_draws %>% mcmc_rank_overlay() # chains look fine

Pl_draws %>% mcmc_pairs(pars = c("alpha[1]", "tau[1]", "k[1]", "mu[1]"))
Pl_draws %>% mcmc_pairs(pars = c("alpha[2]", "tau[2]", "k[2]", "mu[2]"))
# this is chosen as the optimal model

# Prior-posterior comparison
Pl_posterior <- Pl_samples %>%
  recover_types(P_summary %>% select(Species)) %>%
  gather_draws(alpha[Species], tau[Species], k[Species], mu[Species],
               k_mu, mu_mu) %>%
  ungroup() %>%
  mutate(Species = if_else(is.na(Species), "Seagrasses", Species) %>%
                    fct_relevel("Seagrasses", "Halophila ovalis"),
         .variable = case_when(
                      .variable %in% c("k", "k_mu") ~ "k",
                      .variable %in% c("mu", "mu_mu") ~ "mu",
                      TRUE ~ .variable
                      ) %>% fct_relevel("alpha", "tau", "k"))

# Combine with prior so that Distribution and Species are merged into one grouping variable
Pl_prior_posterior <- Pl_posterior %>%
  bind_rows(Pl_prior_samples %>%
              gather_draws(alpha, tau, k, mu) %>%
              ungroup() %>%
              mutate(Species = NA)) %>%
  mutate(Group = if_else(is.na(Species), "Prior", Species) %>%
           fct_relevel("Seagrasses", "Amphibolis antarctica", "Halophila ovalis"),
         .variable = fct_relevel(.variable, "alpha", "tau", "k"))

ggsave(
Pl_prior_posterior %>%
  ggplot() +
    geom_density_ridges(aes(.value, Group, colour = Group, fill = Group),
                        quantile_lines = TRUE, quantiles = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95),
                        alpha = 0.5, scale = 2, rel_min_height = 0.001, bandwidth = c(0.2, 0.2, 0.005, 0.5),
                        from = c(10, rep(0, 3)), to = c(30, 6, 0.6, 42)) +
    scale_colour_manual(values = c("#363538", "#4a7518", "#bdd268", "#b5b8ba"),
                        labels = c("Seagrasses", expression(italic("Amphibolis antarctica")),
                                   expression(italic("Halophila ovalis")), "Prior"),
                        guide = guide_legend(reverse = TRUE)) +
    scale_fill_manual(values = c("#363538", "#4a7518", "#bdd268", "#b5b8ba"),
                      labels = c("Seagrasses", expression(italic("Amphibolis antarctica")),
                                 expression(italic("Halophila ovalis")), "Prior"),
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(~ .variable, scales = "free_x",
               labeller = labeller(.variable = as_labeller(
                 c("alpha" = "italic('a')*' (µmol sample'^-1*' h'^-1*')'",
                   "tau" = "italic('t')*' (µmol sample'^-1*' h'^-1*')'",
                   "k" = "italic('k')*' (d'^-1*')'",
                   "mu" = "italic('µ')*' (d)'"),
                 label_parsed))
               ) +
    facetted_pos_scales(x = list(
      .variable == "alpha" ~ scale_x_continuous(limits = c(10, 30),
                                                breaks = seq(10, 30, by = 10)),
      .variable == "tau" ~ scale_x_continuous(limits = c(0, 6),
                                              breaks = seq(0, 6, by = 2)),
      .variable == "k" ~ scale_x_continuous(limits = c(0, 0.6),
                                            breaks = seq(0, 0.6, by = 0.2),
                                            labels = scales::label_number(accuracy = c(1, 0.1, 0.1, 0.1))),
      .variable == "mu" ~ scale_x_continuous(limits = c(0, 42),
                                             breaks = seq(0, 42, by = 14))
    )) +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()),
filename = "Pl_prior_posterior.pdf", device = cairo_pdf, path = "~/Desktop",
         width = 20, height = 5, units = "cm")

# plot without alpha and tau since they are too divergent on the leaf scale
Fig_1b_top_right <- Pl_prior_posterior %>%
  filter(.variable %in% c("k", "mu")) %>%
  ggplot() +
    geom_density_ridges(aes(.value, Group, colour = Group, fill = Group),
                        quantile_lines = TRUE, quantiles = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95),
                        alpha = 0.5, scale = 2, rel_min_height = 0.001, bandwidth = c(0.005, 0.5),
                        from = c(0, 0), to = c(0.6, 42)) +
    scale_colour_manual(values = c("#363538", "#4a7518", "#bdd268", "#b5b8ba"),
                        labels = c("Seagrasses", expression(italic("Amphibolis antarctica")),
                                   expression(italic("Halophila ovalis")), "Prior"),
                        guide = "none") +
    scale_fill_manual(values = c("#363538", "#4a7518", "#bdd268", "#b5b8ba"),
                      labels = c("Seagrasses", expression(italic("Amphibolis antarctica")),
                                 expression(italic("Halophila ovalis")), "Prior"),
                      guide = "none") +
    facet_grid(~ .variable, scales = "free_x",
               labeller = labeller(.variable = as_labeller(
                 c("k" = "italic('k')*' (d'^-1*')'",
                   "mu" = "italic('µ')*' (d)'"),
                 label_parsed))
               ) +
    facetted_pos_scales(x = list(
      .variable == "k" ~ scale_x_continuous(limits = c(0, 0.6),
                                            breaks = seq(0, 0.6, by = 0.2),
                                            labels = scales::label_number(accuracy = c(1, 0.1, 0.1, 0.1))),
      .variable == "mu" ~ scale_x_continuous(limits = c(0, 42),
                                             breaks = seq(0, 42, by = 14))
    )) +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.spacing = unit(1, "cm"))

# Parameters and differences
# Parameters for Table 1
# correct scale for alpha and tau
Pl_posterior %<>%
  mutate(.value = if_else(Species == "Halophila ovalis" & .variable %in% c("alpha", "tau"),
                          .value / 10, .value))

Pl_posterior %>%
  group_by(Species, .variable) %>%
  summarise(mean = mean(.value),
            sd = sd(.value),
            n = length(.value)) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(mean, digits = 2), signif(mean, digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        )
        ) %>%
  pivot_wider(names_from = Species, values_from = rounded) %>%
  write.csv("~/Desktop/Pl_para.csv", row.names = FALSE)

# Calculate differences
Pl_diff <- Pl_posterior %>%
  filter(Species != "Seagrasses") %>%
  pivot_wider(names_from = c(Species, .variable), values_from = .value) %>%
  mutate(delta_alpha = `Halophila ovalis_alpha` - `Amphibolis antarctica_alpha`,
         delta_tau = `Halophila ovalis_tau` - `Amphibolis antarctica_tau`,
         delta_k = `Halophila ovalis_k` - `Amphibolis antarctica_k`,
         delta_mu = `Halophila ovalis_mu` - `Amphibolis antarctica_mu`) %>%
  select(.chain, .iteration, .draw, delta_alpha, delta_tau, delta_k, delta_mu) %>%
  pivot_longer(cols = c(delta_alpha, delta_tau, delta_k, delta_mu),
               names_to = "Parameter", values_to = "Difference", names_prefix = "delta_") %>%
  mutate(Parameter = fct_relevel(Parameter, "alpha", "tau", "k"))

# Differences and probabilities for Table 1
Pl_diff %>%
  group_by(Parameter) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P_more = mean(Difference > 0),
            P_less = mean(Difference < 0),
            n = length(Difference)) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(abs(mean), digits = 2), signif(abs(mean), digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        ),
         P = signif(pmax(P_less, P_more), digits = 2)) %>%
  write.csv("~/Desktop/Pl_diff.csv", row.names = FALSE)


# Differences and probabilities for text: comparing k and mu between mass- and leaf-based
Pm_Pl_diff <- Pm_posterior %>%
  mutate(Scale = "Mass") %>%
  bind_rows(Pl_posterior %>% mutate(Scale = "Leaf")) %>%
  pivot_wider(names_from = c(Scale, .variable), values_from = .value) %>%
  mutate(delta_k = Mass_k - Leaf_k,
         delta_mu = Mass_mu - Leaf_mu) %>%
  select(.chain, .iteration, .draw, Species, delta_k, delta_mu) %>%
  pivot_longer(cols = c(delta_k, delta_mu),
               names_to = "Parameter", values_to = "Difference", names_prefix = "delta_")

Pm_Pl_diff %>%
  group_by(Species, Parameter) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P_more = mean(Difference > 0),
            P_less = mean(Difference < 0),
            n = length(Difference)) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(abs(mean), digits = 2), signif(abs(mean), digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        ),
  P = signif(pmax(P_less, P_more), digits = 2)) %>%
  write.csv("~/Desktop/Pm_Pl_diff.csv", row.names = FALSE)

# calculate P_mu from parameters
Pl_mu <- Pl_posterior %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(P_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau)

Pl_mu_summary <- Pl_mu %>%
  group_by(Day, Species) %>%
  median_qi(P_mu, .width = c(.5, .8, .9))

Fig_1b_bottom_left <- ggplot() +
                  geom_violin(data = P %>%
                                filter(Species == "Halophila ovalis"),
                              aes(Day, Pl, group = ID), colour = "#bdd268",
                              fill = "#bdd268", alpha = 0.2,
                              position = "identity", width = 4) +
                  geom_line(data = Pl_mu_summary %>%
                              filter(Species == "Halophila ovalis"),
                            aes(Day, P_mu), colour = "#bdd268") +
                  geom_ribbon(data = Pl_mu_summary %>%
                                filter(Species == "Halophila ovalis"),
                              aes(Day, ymin = .lower, ymax = .upper,
                                  alpha = factor(.width)), fill = "#bdd268", colour = NA) +
                  labs(y = expression(italic(P)["max"]*" (µmol O"[2]*" leaf"^-1*" h"^-1*")"),
                       x = "Detrital age (d)") +
                  scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
                  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5),
                                     labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1, 0.1))) +
                  coord_cartesian(xlim = c(-1.8, 35), ylim = c(0, 2.5), expand = FALSE, clip = "off") +
                  mytheme

Fig_1b_bottom_right <- ggplot() +
                  geom_hline(yintercept = 0) +
                  geom_violin(data = P %>%
                                filter(Species == "Amphibolis antarctica"),
                              aes(Day, Pl, group = ID), colour = "#4a7518",
                              fill = "#4a7518", alpha = 0.2,
                              position = "identity", width = 4) +
                  geom_line(data = Pl_mu_summary %>%
                              filter(Species == "Amphibolis antarctica"),
                            aes(Day, P_mu), colour = "#4a7518") +
                  geom_ribbon(data = Pl_mu_summary %>%
                                filter(Species == "Amphibolis antarctica"),
                              aes(Day, ymin = .lower, ymax = .upper,
                              alpha = factor(.width)), fill = "#4a7518", colour = NA) +
                  labs(y = expression(italic(P)["max"]*" (µmol O"[2]*" leaf"^-1*" h"^-1*")"),
                       x = "Detrital age (d)") +
                  scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
                  scale_y_continuous(breaks = seq(-10, 40, by = 10),
                                     labels = scales::label_number(style_negative = "minus")) +
                  coord_cartesian(xlim = c(-1.8, 35), ylim = c(-10, 40), expand = FALSE, clip = "off") +
                  mytheme +
                  theme(axis.title = element_blank())

Fig_1 <- Fig_1a_top / Fig_1a_bottom / 
         ( ( wrap_elements() | Fig_1b_top_right ) + plot_layout(widths = c(1, 0.9345)) ) /
         ( Fig_1b_bottom_left | Fig_1b_bottom_right ) +
         plot_layout(heights = c(0.2, 1, 0.2, 1)) +
         plot_annotation(tag_levels = list(c("A", "", "B"))) &
         theme(plot.tag = element_text(family = "Futura", size = 15, face = "bold"))

ggsave(plot = Fig_1, filename = "Fig_1.pdf", device = cairo_pdf, path = "~/Desktop",
       width = 20, height = 20, units = "cm")

