#### Plant limbo: seagrass detrital photosynthesis ####
#### Luka Seamus Wright                            ####

# 1. Oxygen ####
# 1.1 Load data ####
require(here)
files <- list.files(path = here("Seagrass", "Oxygen"), 
                    pattern = "\\.csv$", full.names = TRUE)

require(tidyverse)
O2 <- files %>%
  map(~ read.csv(., skip = 1, header = TRUE) %>%
        drop_na(Value) %>%
        mutate(delta_t = as.numeric(delta_t),
               delta_t_c = delta_t - mean(delta_t))) %>%
  set_names(str_remove(basename(files), "\\.csv$") %>% make.names)
str(O2$X220915_B1)

# 1.2 Stan model ####
require(tidybayes)
O2_list <- O2 %>%
  map(~ select(., Value, delta_t_c) %>%
        compose_data())

require(cmdstanr)
O2_mod <- here("Seagrass", "Stan", "O2.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

O2_samples <- O2_list %>%
  map(~ O2_mod$sample(data = .,
                      seed = 100,
                      chains = 8,
                      parallel_chains = parallel::detectCores(),
                      iter_warmup = 1e4,
                      iter_sampling = 1e4))

# 1.3 Model checks ####
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
# good Rhat and effective sample size: 
# rhat = 1.00 ± 0.0000895
# ess = 56507 ± 27269

O2_draws <- O2_samples %>%
  map(~ .$draws(format = "df"))

require(bayesplot)
require(patchwork)
O2_draws %>% 
  map(~ mcmc_rank_overlay(.)) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_rank.pdf", device = cairo_pdf, 
         path = here("Seagrass", "Plots"),
         width = 100, height = 100, units = "cm")

O2_draws %>% 
  sample(6) %>%
  map(~ mcmc_pairs(., pars = c("alpha", "beta"))) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_pairs.pdf", device = cairo_pdf, 
         path = here("Seagrass", "Plots"),
         width = 50, height = 50, units = "cm")


# 1.4 Prior-posterior comparison #### 
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
  ggsave(filename = "O2_prior_posterior.pdf", device = cairo_pdf, 
         path = here("Seagrass", "Plots"),
         width = 100, height = 100, units = "cm")

# 1.5 Prediction #### 
# 1.5.1 Mu #### 
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
        mean_qi(mu, .width = c(.5, .8, .9)))
         
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
  ggsave(filename = "O2_fit.pdf", device = cairo_pdf, 
         path = here("Seagrass", "Plots"),
         width = 100, height = 100, units = "cm")
# fit looks good

# 1.5.2 Initial O2 #### 
O2_initial <- O2_mu %>%
  map(~ slice_min(., delta_t_c))

O2_initial %>%
  map(~ ggplot(., aes(mu)) +
          geom_density() +
          theme_minimal()) %>%
  wrap_plots() %>%
  ggsave(filename = "O2_initial.pdf", device = cairo_pdf, 
         path = here("Seagrass", "Plots"),
         width = 100, height = 100, units = "cm")

# 2. Volume #### 
# 2.1 Load data #### 
V <- here("Seagrass", "Volume.csv") %>% read.csv() %>%
  mutate(Vs = (Volume - mean(Volume)) / sd(Volume)) # standardise Volume

# 2.2 Stan model #### 
V_mod <- here("Seagrass", "Stan", "V.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

V_samples <- V_mod$sample(data = compose_data(V),
                          seed = 100,
                          chains = 8,
                          parallel_chains = parallel::detectCores(),
                          iter_warmup = 1e4,
                          iter_sampling = 1e4)

# 2.3 Model checks #### 
# check Rhat, effective sample size and chains
V_summary <- V_samples$summary()

V_summary %>%
  filter(rhat > 1.001)
# no Rhat above 1.001

V_draws <- V_samples$draws(format = "df")

V_draws %>% mcmc_rank_overlay()

V_draws %>% mcmc_pairs()

# 2.4 Prior-posterior comparison ####
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

# 2.5 Prediction ####
require(magrittr)
V_prior_posterior %<>%
  mutate(Vmu = V %$% (Vsmu * sd(Volume) + mean(Volume))) %>%
  select(Vmu)

# 3. Pressure ####
O2 %>%
  map(~ summarise(., Tmean = mean(Temp),
                  Tsd = sd(Temp),
                  Pmean = mean(Pressure),
                  Psd = sd(Pressure))) %>%
  bind_rows() %>%
  filter(Psd == 0 | Tsd == 0)
# in several cases pressure did not vary at all across the incubation,
# so measurement error in pressure can not be included in the model

# 4. Temperature ####
# 4.1 Prepare data ####
O2 %<>%
  map(~ mutate(., Ts = ( Temp - mean(Temp) ) / sd(Temp)))

T_list <- O2 %>%
  map(~ select(., Ts) %>%
        compose_data())

# 4.2 Stan model ####
T_mod <- here("Seagrass", "Stan", "T.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

T_samples <- T_list %>%
  map(~ T_mod$sample(data = .,
                     seed = 100,
                     chains = 8,
                     parallel_chains = parallel::detectCores(),
                     iter_warmup = 1e4,
                     iter_sampling = 1e4))

# 4.3 Model checks ####
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
# good Rhat and effective sample size:
# rhat = 1.00 ± 0.0000628
# ess = 59260 ± 15876

T_draws <- T_samples %>%
  map(~ .$draws(format = "df"))

T_draws %>% 
  map(~ mcmc_rank_overlay(.)) %>%
  wrap_plots() %>%
  ggsave(filename = "T_rank.pdf", device = cairo_pdf, 
         path = here("Seagrass", "Plots"),
         width = 100, height = 100, units = "cm")

T_draws %>% 
  sample(6) %>%
  map(~ mcmc_pairs(., pars = c("Tsmu", "sigma"))) %>%
  wrap_plots() %>%
  ggsave(filename = "T_pairs.pdf", device = cairo_pdf, 
         path = here("Seagrass", "Plots"),
         width = 50, height = 50, units = "cm")


# 4.4 Prior-posterior comparison ####
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
  ggsave(filename = "T_prior_posterior.pdf", device = cairo_pdf, 
         path = here("Seagrass", "Plots"),
         width = 100, height = 100, units = "cm")

T_prior_posterior %<>%
  map2(O2,
       ~ mutate(.x, Tmu = .y %$% (Tsmu * sd(Temp) + mean(Temp))) %>%
         select(Tmu)
       )

# 5. Synthesis of predictors ####
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

# 6. Mass ####
# 6.1 Load data ####
M <- here("Seagrass", "Mass.csv") %>% read.csv() %>%
  mutate(Species = fct_relevel(Species, "Halophila ovalis"),
         Date = Date %>% dmy()) %>% # convert Date to a usable date
  group_by(Species) %>%
  mutate(Day = min(Date) %--% Date %>%
           time_length("day")) %>%
  ungroup() %>% # very important
  mutate(Day_c = Day - mean(Day)) # centre Day

# 6.2 Visualise data ####
M %>%
  ggplot(aes(Day, Mass, colour = Species)) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_minimal()

# 6.3 Prior simulation ####
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

# 6.4 Stan model ####
M_mod <- here("Seagrass", "Stan", "M.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

M_samples <- M_mod$sample(data = M %>%
                            select(Mass, Day_c, Species) %>%
                            compose_data(),
                          seed = 100,
                          chains = 8,
                          parallel_chains = parallel::detectCores(),
                          iter_warmup = 1e4,
                          iter_sampling = 1e4)

# 6.5 Model checks ####
M_summary <- M_samples$summary()
M_summary %>%
  filter(rhat > 1.001)

M_draws <- M_samples$draws(format = "df")
M_draws %>% mcmc_rank_overlay() # chains look good
M_draws %>% mcmc_pairs(pars = c("alpha[1]", "beta[1]"))
M_draws %>% mcmc_pairs(pars = c("alpha[2]", "beta[2]"))

# 6.6 Prediction ####
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
  
# 7. Photosynthesis ####
# 7.1 Calculate response ####
# correct photosynthesis for blank, volume and mass (µM min^-1 to µmol g^-1 h^-1)
P %<>%
  filter(Species != "Blank") %>%
  left_join(P %>% filter(Species == "Blank") %>%
              select(.draw, Date, Time, beta) %>%
              rename(blank = beta),
            by = c(".draw", "Date", "Time"),
            relationship = "many-to-one") %>%
  left_join(M, by = c("Date", "Species", "Leaf"),
            relationship = "many-to-one") %>%
  mutate(Pm = ( beta - blank ) * Volume * 60 / Mass,
         Pl = if_else(Species == "Halophila ovalis", 
                      ( beta - blank ) * Volume * 60 / 10, # there were 10 H. ovalis leaves per sample
                      ( beta - blank ) * Volume * 60)) %>%
  select(-c(beta, Volume, blank))

# 7.2 Calculate predictor ####
# calculate detrital age (days) from dates
P %<>%
  group_by(Species) %>%
  mutate(Day = min(Date) %--% Date %>%
           time_length("day")) %>%
  ungroup() %T>%
  print()

# 7.3 Summarise response and predictors ####
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
  mutate(Species = fct_relevel(Species, "Halophila ovalis")) %T>%
  print()

# 7.4 Check duplicates ####
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

# 7.5 Clean up and back up ####
# clean up
rm(O2, O2_draws, O2_initial, O2_list, O2_mod, O2_mu, O2_mu_summary,
   O2_prior_posterior, O2_samples, O2_summary, T_draws, T_list, T_mod,
   T_prior_posterior, T_samples, T_summary, V, V_draws, V_mod, V_samples,
   V_prior_posterior, V_summary, files, O2_stan, T_stan, V_stan)

# save data
M %>% write_rds(file = here("Seagrass", "RDS", "M.rds"))
P %>% write_rds(file = here("Seagrass", "RDS", "P.rds"))
P_summary %>% write_rds(file = here("Seagrass", "RDS", "P_summary.rds"))

# 7.6 Mass-based photosynthesis ####
# 7.6.1 Standardise confounders ####
# Confounders are best standardised (z-scores) because they have
# varying scales that are hard to anticipate in the prior. This
# also helps when predicting because setting them to zero is
# equivalent to holding them at their mean.
P_summary %<>%
  mutate(O2_mean_std = ( O2_mean - mean(O2_mean) ) / sd(O2_mean),
         O2_sd_std = O2_sd / sd(O2_mean),
         T_mean_std = ( T_mean - mean(T_mean) ) / sd(T_mean),
         T_sd_std = T_sd / sd(T_mean),
         P_mean_std = ( P_mean - mean(P_mean) ) / sd(P_mean),
         S_std = ( S - mean(S) ) / sd(S),
         M_std = ( M - mean(M) ) / sd(M)) %T>%
  print()

# 7.6.2 Prior simulation ####
# Load data from literature
Prior <- here("Seagrass", "Prior", "Prior.csv") %>% 
  read.csv() %T>%
  print()

# Calculate leaf-based flux
Prior %<>%
  left_join(M %>% select(Species, Mass) %>%
              mutate(Leafmass = if_else(Species == "Halophila ovalis",
                                        Mass / 10, Mass)) %>% # There were 10 leaves for H. ovalis
              group_by(Species) %>%
              summarise(Leafmass = mean(Leafmass)),
            by = "Species", relationship = "many-to-one") %>%
  mutate(Flux_leaf = Flux * Leafmass) %T>%
  print()

# Calculate species-specific summary stats
Prior %>%
  group_by(Species, Variable) %>%
  summarise(Fm_mean = mean(Flux),
            Fm_sd = sd(Flux),
            Fm_median = median(Flux),
            Fl_mean = mean(Flux_leaf),
            Fl_sd = sd(Flux_leaf),
            Fl_median = median(Flux_leaf),
            n = n())

# Calculate overall summary stats
Prior %>%
  group_by(Variable) %>%
  summarise(Pm_mean = mean(Flux),
            Pm_sd = sd(Flux),
            Pm_median = median(Flux),
            Pl_mean = mean(Flux_leaf),
            Pl_sd = sd(Flux_leaf),
            Pl_median = median(Flux_leaf),
            n = n())
# As expected, based on the mean-median comparison,
# photosynthesis and respiration are right-skewed.
# Better proceed with the median.

# Simulate
tibble(n = 1:1e3, 
       # I am exponentiating here because that's how I'm including confounders
       alpha = exp( rnorm( 1e3 , log(21) , 0.4 ) ), # based on prior median of 20.7
       # Two weeks seems reasonable based on E. radiata (doi: 10.1093/aob/mcad167)
       mu = exp( rnorm( 1e3 , log(14) , 0.6 ) ), # more uncertainty around mu
       tau = exp( rnorm( 1e3 , log(2) , 0.4 ) ), # based on prior median of 1.93
       Pm_sigma = rexp( 1e3 , 1 )) %>%
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(
    Pm_mu = (alpha + tau) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    Pm = rnorm( n() , Pm_mu , Pm_sigma )
  ) %>%
  pivot_longer(cols = c(Pm_mu, Pm),
               names_to = "parameter") %>%
  ggplot(aes(Day, value, group = n)) +
    geom_line(alpha = 0.05) +
    geom_hline(yintercept = P_summary %$% range(Pm_mean)) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks very reasonable!

# 7.6.3 Stan model ####
Pm_mod <- here("Seagrass", "Stan", "Pm.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

Pm_samples <- Pm_mod$sample(
  data = P_summary %>%
    select(Day, Pm_mean, Pm_sd, Species,
           O2_mean_std, O2_sd_std, T_mean_std, T_sd_std,
           P_mean_std, S_std, M_std) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4
)

# Save draws
Pm_samples$draws() %>%
  write_rds(here("Seagrass", "RDS", "Pm_samples.rds"))
Pm_samples$draws(format = "df") %>%
  write_rds(here("Seagrass", "RDS", "Pm_samples_df.rds"))

# 7.6.4 Model checks ####
# R-hat
Pm_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000994.

# Chains
Pm_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(units = "cm", height = 50, width = 50, 
         device = cairo_pdf, filename = "Pm_chains.pdf",
         path = here("Seagrass", "Plots"))
# Looks good

# Pairs
Pm_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c(
      "log_alpha_mu[1]", "log_alpha_mu[2]",
      "log_mu_mu[1]", "log_mu_mu[2]",
      "log_tau_mu[1]", "log_tau_mu[2]",
      "Pm_sigma[1]", "Pm_sigma[2]"
    )
  ) %>%
  ggsave(units = "cm", height = 50, width = 50, 
         filename = "Pm_pairs.png",
         path = here("Seagrass", "Plots"))

Pm_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c(
      "log_alpha_mu[1]", "log_alpha_O2[1]",
      "log_alpha_T[1]", "log_alpha_M[1]",
      "log_mu_mu[1]", "log_mu_O2[1]",
      "log_mu_T[1]", "log_mu_M[1]",
      "log_tau_mu[1]", "log_tau_O2[1]",
      "log_tau_T[1]", "log_tau_M[1]"
    )
  ) %>%
  ggsave(units = "cm", height = 50, width = 50, 
         filename = "Pm_pairs_confound.png",
         path = here("Seagrass", "Plots"))
# No correlations. All looks well.

# 7.6.5 Prior-posterior comparison ####
# Sample prior
source("functions.R")
Pm_prior <- prior_samples(
  model = Pm_mod,
  data = P_summary %>%
    select(Day, Pm_mean, Pm_sd, Species,
           O2_mean_std, O2_sd_std, T_mean_std, T_sd_std,
           P_mean_std, S_std, M_std) %>%
    compose_data()
  )

Pm_prior %>% 
  prior_posterior_draws(
    posterior_samples = Pm_samples,
    group = P_summary %>% select(Species),
    parameters = c("log_alpha_mu[Species]", "log_alpha_O2[Species]",
                   "log_alpha_T[Species]", "log_alpha_P[Species]",
                   "log_alpha_S[Species]", "log_alpha_M[Species]",
                   "log_mu_mu[Species]", "log_mu_O2[Species]",
                   "log_mu_T[Species]", "log_mu_P[Species]",
                   "log_mu_S[Species]", "log_mu_M[Species]",
                   "log_tau_mu[Species]", "log_tau_O2[Species]",
                   "log_tau_T[Species]", "log_tau_P[Species]",
                   "log_tau_S[Species]", "log_tau_M[Species]",
                   "Pm_sigma[Species]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Species") %>%
  ggsave(units = "cm", height = 30, width = 50, 
         device = cairo_pdf, filename = "Pm_prior_posterior.pdf",
         path = here("Seagrass", "Plots"))

# 7.6.6 Parameters ####
Pm_prior_posterior <- Pm_prior %>% 
  prior_posterior_draws(
    posterior_samples = Pm_samples,
    group = P_summary %>% select(Species),
    parameters = c("log_alpha_mu[Species]", "log_alpha_O2[Species]",
                   "log_alpha_T[Species]", "log_alpha_P[Species]",
                   "log_alpha_S[Species]", "log_alpha_M[Species]",
                   "log_mu_mu[Species]", "log_mu_O2[Species]",
                   "log_mu_T[Species]", "log_mu_P[Species]",
                   "log_mu_S[Species]", "log_mu_M[Species]",
                   "log_tau_mu[Species]", "log_tau_O2[Species]",
                   "log_tau_T[Species]", "log_tau_P[Species]",
                   "log_tau_S[Species]", "log_tau_M[Species]",
                   "Pm_sigma[Species]"),
    format = "short"
  ) %>%
  filter(Species == "Halophila ovalis" & distribution == "prior" |
           distribution == "posterior") %>% # Remove redundant priors
  mutate( # Embed prior in Species
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct()
  ) %>%
  select(-distribution) %T>%
  print()

Pm_prior_posterior %>% # Save
  write_rds(here("Seagrass", "RDS", "Pm_prior_posterior.rds"))

Pm_contrast <- Pm_prior_posterior %>%
  filter(Species != "Prior") %>%
  mutate( # Parameters conditional on confounders held at their means (zero)
    alpha = exp(log_alpha_mu),
    mu = exp(log_mu_mu),
    tau = exp(log_tau_mu),
  ) %>%
  select(starts_with("."), Species, 
         alpha, mu, tau) %>%
  pivot_longer(cols = c(alpha, mu, tau),
               names_to = "Parameter") %>%
  pivot_wider(names_from = Species) %>%
  mutate(
    # Calculate contrast as difference and ratio
    diff = `Halophila ovalis` - `Amphibolis antarctica`,
    ratio = `Halophila ovalis` / `Amphibolis antarctica`,
  ) %T>%
  print()

# 7.6.7 Prediction ####
Pm_prediction <- Pm_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "Day") %>%
  mutate( # Predictions conditional on confounders held at their means (zero)
    alpha = exp(log_alpha_mu),
    mu = exp(log_mu_mu),
    tau = exp(log_tau_mu),
    Pm_mu = (alpha + tau) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    Pm = rnorm( n() , Pm_mu , Pm_sigma )
  ) %>% # Summarise predictions
  group_by(Day, Species) %>%
  median_qi(Pm_mu, Pm, .width = c(.5, .8, .9)) %T>%
  print()

Pm_prediction %>% # Save
  write_rds(here("Seagrass", "RDS", "Pm_prediction.rds"))

# 7.6.8 Confounder prediction ####
Pm_prediction_O2 <- Pm_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "O2_mean_std") %>%
  mutate( # The other confounders are held at their mean by exclusion
    log_alpha = log_alpha_mu + log_alpha_O2 * O2_mean_std,
    log_mu = log_mu_mu + log_mu_O2 * O2_mean_std,
    log_tau = log_tau_mu + log_tau_O2 * O2_mean_std
  ) %>%
  group_by(O2_mean_std, Species) %>%
  median_qi(log_alpha, log_mu, log_tau, .width = c(.5, .8, .9)) %T>%
  print()

Pm_prediction_T <- Pm_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "T_mean_std") %>%
  mutate(
    log_alpha = log_alpha_mu + log_alpha_T * T_mean_std,
    log_mu = log_mu_mu + log_mu_T * T_mean_std,
    log_tau = log_tau_mu + log_tau_T * T_mean_std
  ) %>%
  group_by(T_mean_std, Species) %>%
  median_qi(log_alpha, log_mu, log_tau, .width = c(.5, .8, .9)) %T>%
  print()

Pm_prediction_P <- Pm_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "P_mean_std") %>%
  mutate(
    log_alpha = log_alpha_mu + log_alpha_P * P_mean_std,
    log_mu = log_mu_mu + log_mu_P * P_mean_std,
    log_tau = log_tau_mu + log_tau_P * P_mean_std
  ) %>%
  group_by(P_mean_std, Species) %>%
  median_qi(log_alpha, log_mu, log_tau, .width = c(.5, .8, .9)) %T>%
  print()

Pm_prediction_S <- Pm_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "S_std") %>%
  mutate(
    log_alpha = log_alpha_mu + log_alpha_S * S_std,
    log_mu = log_mu_mu + log_mu_S * S_std,
    log_tau = log_tau_mu + log_tau_S * S_std
  ) %>%
  group_by(S_std, Species) %>%
  median_qi(log_alpha, log_mu, log_tau, .width = c(.5, .8, .9)) %T>%
  print()

Pm_prediction_M <- Pm_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "M_std") %>%
  mutate(
    log_alpha = log_alpha_mu + log_alpha_M * M_std,
    log_mu = log_mu_mu + log_mu_M * M_std,
    log_tau = log_tau_mu + log_tau_M * M_std
  ) %>%
  group_by(M_std, Species) %>%
  median_qi(log_alpha, log_mu, log_tau, .width = c(.5, .8, .9)) %T>%
  print()

# Combine confounder predictions
Pm_prediction_confound <- bind_rows(
  Pm_prediction_O2 %>%
    mutate(Confounder = "O2" %>% fct()) %>%
    rename(Predictor = O2_mean_std),
  Pm_prediction_T %>%
    mutate(Confounder = "T" %>% fct()) %>%
    rename(Predictor = T_mean_std),
  Pm_prediction_P %>%
    mutate(Confounder = "P" %>% fct()) %>%
    rename(Predictor = P_mean_std),
  Pm_prediction_S %>%
    mutate(Confounder = "S" %>% fct()) %>%
    rename(Predictor = S_std),
  Pm_prediction_M %>%
    mutate(Confounder = "M" %>% fct()) %>%
    rename(Predictor = M_std)
) %T>%
  print()

# Clean up
rm(Pm_prediction_O2, Pm_prediction_T, Pm_prediction_P,
   Pm_prediction_S, Pm_prediction_M)

# Reverse standardisation
Pm_prediction_confound %<>%
  mutate(
    Predictor = P_summary %$% case_when(
      Confounder == "O2" ~ Predictor * sd(O2_mean) + mean(O2_mean),
      Confounder == "T" ~ Predictor * sd(T_mean) + mean(T_mean),
      Confounder == "P" ~ Predictor * sd(P_mean) + mean(P_mean),
      Confounder == "S" ~ Predictor * sd(S) + mean(S),
      Confounder == "M" ~ Predictor * sd(M) + mean(M)
    )
  ) %T>%
  print()

Pm_prediction_confound %>% # Save
  write_rds(here("Seagrass", "RDS", "Pm_prediction_confound.rds"))

# 7.7 Leaf-based photosynthesis ####
# 7.7.1 Prior simulation ####
tibble(n = 1:1e3, 
       # I am picking a value from the prior medians 1 and 7 (see above) 
       # for alpha and increasing uncertainty because I know that leaves 
       # vary dramatically in size.
       alpha = exp( rnorm( 1e3 , log(7) , 1 ) ), # the upper end is more sensible
       mu = exp( rnorm( 1e3 , log(14) , 0.6 ) ), # same as before
       tau = exp( rnorm( 1e3 , log(1) , 0.4 ) ), # no prior data, but likely smaller
       Pl_sigma = rexp( 1e3 , 1 )) %>%
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(
    Pl_mu = (alpha + tau) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    Pl = rnorm( n() , Pl_mu , Pl_sigma )
  ) %>%
  pivot_longer(cols = c(Pl_mu, Pl),
               names_to = "parameter") %>%
  ggplot(aes(Day, value, group = n)) +
    geom_line(alpha = 0.05) +
    geom_hline(yintercept = P_summary %$% range(Pl_mean)) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks fine

# 7.7.2 Stan model ####
Pl_mod <- here("Seagrass", "Stan", "Pl.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

Pl_samples <- Pl_mod$sample(
  data = P_summary %>%
    select(Day, Pl_mean, Pl_sd, Species,
           O2_mean_std, O2_sd_std, T_mean_std, T_sd_std,
           P_mean_std, S_std, M_std) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4
)

# Save draws
Pl_samples$draws() %>%
  write_rds(here("Seagrass", "RDS", "Pl_samples.rds"))
Pl_samples$draws(format = "df") %>%
  write_rds(here("Seagrass", "RDS", "Pl_samples_df.rds"))

# 7.7.3 Model checks ####
# R-hat
Pl_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000101.

# Chains
Pl_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(units = "cm", height = 50, width = 50, 
         device = cairo_pdf, filename = "Pl_chains.pdf",
         path = here("Seagrass", "Plots"))
# Looks good

# Pairs
Pl_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c(
      "log_alpha_mu[1]", "log_alpha_mu[2]",
      "log_mu_mu[1]", "log_mu_mu[2]",
      "log_tau_mu[1]", "log_tau_mu[2]",
      "Pl_sigma[1]", "Pl_sigma[2]"
    )
  ) %>%
  ggsave(units = "cm", height = 50, width = 50, 
         filename = "Pl_pairs.png",
         path = here("Seagrass", "Plots"))

Pl_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c(
      "log_alpha_mu[1]", "log_alpha_O2[1]",
      "log_alpha_T[1]", "log_alpha_M[1]",
      "log_mu_mu[1]", "log_mu_O2[1]",
      "log_mu_T[1]", "log_mu_M[1]",
      "log_tau_mu[1]", "log_tau_O2[1]",
      "log_tau_T[1]", "log_tau_M[1]"
    )
  ) %>%
  ggsave(units = "cm", height = 50, width = 50, 
         filename = "Pl_pairs_confound.png",
         path = here("Seagrass", "Plots"))
# No or very little correlation. All looks well.

# 7.7.4 Prior-posterior comparison ####
# Sample prior
Pl_prior <- prior_samples(
  model = Pl_mod,
  data = P_summary %>%
    select(Day, Pl_mean, Pl_sd, Species,
           O2_mean_std, O2_sd_std, T_mean_std, T_sd_std,
           P_mean_std, S_std, M_std) %>%
    compose_data()
  )

Pl_prior %>% 
  prior_posterior_draws(
    posterior_samples = Pl_samples,
    group = P_summary %>% select(Species),
    parameters = c("log_alpha_mu[Species]", "log_alpha_O2[Species]",
                   "log_alpha_T[Species]", "log_alpha_P[Species]",
                   "log_alpha_S[Species]", "log_alpha_M[Species]",
                   "log_mu_mu[Species]", "log_mu_O2[Species]",
                   "log_mu_T[Species]", "log_mu_P[Species]",
                   "log_mu_S[Species]", "log_mu_M[Species]",
                   "log_tau_mu[Species]", "log_tau_O2[Species]",
                   "log_tau_T[Species]", "log_tau_P[Species]",
                   "log_tau_S[Species]", "log_tau_M[Species]",
                   "Pl_sigma[Species]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Species") %>%
  ggsave(units = "cm", height = 30, width = 50, 
         device = cairo_pdf, filename = "Pl_prior_posterior.pdf",
         path = here("Seagrass", "Plots"))

# 7.7.5 Parameters ####
Pl_prior_posterior <- Pl_prior %>% 
  prior_posterior_draws(
    posterior_samples = Pl_samples,
    group = P_summary %>% select(Species),
    parameters = c("log_alpha_mu[Species]", "log_alpha_O2[Species]",
                   "log_alpha_T[Species]", "log_alpha_P[Species]",
                   "log_alpha_S[Species]", "log_alpha_M[Species]",
                   "log_mu_mu[Species]", "log_mu_O2[Species]",
                   "log_mu_T[Species]", "log_mu_P[Species]",
                   "log_mu_S[Species]", "log_mu_M[Species]",
                   "log_tau_mu[Species]", "log_tau_O2[Species]",
                   "log_tau_T[Species]", "log_tau_P[Species]",
                   "log_tau_S[Species]", "log_tau_M[Species]",
                   "Pl_sigma[Species]"),
    format = "short"
  ) %>%
  filter(Species == "Halophila ovalis" & distribution == "prior" |
           distribution == "posterior") %>% # Remove redundant priors
  mutate( # Embed prior in Species
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct()
  ) %>%
  select(-distribution) %T>%
  print()

Pl_prior_posterior %>% # Save
  write_rds(here("Seagrass", "RDS", "Pl_prior_posterior.rds"))

Pl_contrast <- Pl_prior_posterior %>%
  filter(Species != "Prior") %>%
  mutate(
    alpha = exp(log_alpha_mu),
    mu = exp(log_mu_mu),
    tau = exp(log_tau_mu),
  ) %>%
  select(starts_with("."), Species, 
         alpha, mu, tau) %>%
  pivot_longer(cols = c(alpha, mu, tau),
               names_to = "Parameter") %>%
  pivot_wider(names_from = Species) %>%
  mutate(
    diff = `Halophila ovalis` - `Amphibolis antarctica`,
    ratio = `Halophila ovalis` / `Amphibolis antarctica`,
  ) %T>%
  print()

# 7.7.6 Prediction ####
Pl_prediction <- Pl_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "Day") %>%
  mutate( # Predictions conditional on confounders held at their means (zero)
    alpha = exp(log_alpha_mu),
    mu = exp(log_mu_mu),
    tau = exp(log_tau_mu),
    Pl_mu = (alpha + tau) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    Pl = rnorm( n() , Pl_mu , Pl_sigma )
  ) %>% # Summarise predictions
  group_by(Day, Species) %>%
  median_qi(Pl_mu, Pl, .width = c(.5, .8, .9)) %T>%
  print()

Pl_prediction %>% # Save
  write_rds(here("Seagrass", "RDS", "Pl_prediction.rds"))

# 7.8 Linear confounder model ####
# 7.8.1 Prior simulation ####
# I am only going to run this on mass-based photosynthesis, so the intercept
# prior will be centred on the prior median of 21.
tibble(n = 1:1e3, 
       alpha = rnorm( 1e3 , 21 , 10 ),
       beta_O2 = rnorm( 1e3 , 0 , 1 ), # These slopes represent the change per
       beta_T = rnorm( 1e3 , 0 , 1 ), # standard deviation in the predictor.
       beta_P = rnorm( 1e3 , 0 , 1 ),
       beta_S = rnorm( 1e3 , 0 , 1 ),
       beta_M = rnorm( 1e3 , 0 , 1 ),
       Pm_sigma = rexp( 1e3 , 1 )) %>%
  expand_grid(
    Predictor = P_summary %$% 
      seq(
        min( c(O2_mean_std, T_mean_std, P_mean_std, S_std, M_std) ), 
        max( c(O2_mean_std, T_mean_std, P_mean_std, S_std, M_std) ), 
        length.out = 100
      )
  ) %>%
  mutate(
    Pm_mu = alpha + beta_O2 * Predictor + beta_T * Predictor +
      beta_P * Predictor + beta_S * Predictor + beta_M * Predictor,
    Pm = rnorm( n() , Pm_mu , Pm_sigma )
  ) %>%
  pivot_longer(cols = c(Pm_mu, Pm),
               names_to = "parameter") %>%
  ggplot(aes(Predictor, value, group = n)) +
    geom_line(alpha = 0.05) +
    geom_hline(yintercept = P_summary %$% range(Pm_mean)) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Looks reasonable

# 7.8.2 Stan model ####
Pm_confound_mod <- here("Seagrass", "Stan", "Pm_confound.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

Pm_confound_samples <- Pm_confound_mod$sample(
  data = P_summary %>%
    select(Pm_mean, Pm_sd, Species,
           O2_mean_std, O2_sd_std, T_mean_std, T_sd_std,
           P_mean_std, S_std, M_std) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4
)

# Save draws
Pm_confound_samples$draws() %>%
  write_rds(here("Seagrass", "RDS", "Pm_confound_samples.rds"))
Pm_confound_samples$draws(format = "df") %>%
  write_rds(here("Seagrass", "RDS", "Pm_confound_samples_df.rds"))

# 7.8.3 Model checks ####
# R-hat
Pm_confound_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000100.

# Chains
Pm_confound_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(units = "cm", height = 50, width = 50, 
         device = cairo_pdf, filename = "Pm_confound_chains.pdf",
         path = here("Seagrass", "Plots"))
# Looks good

# Pairs
Pm_confound_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c(
      "alpha[1]", "alpha[2]",
      "beta_O2[1]", "beta_O2[2]",
      "beta_T[1]", "beta_T[2]",
      "beta_P[1]", "beta_P[2]",
      "beta_S[1]", "beta_S[2]",
      "beta_M[1]", "beta_M[2]"
    )
  ) %>%
  ggsave(units = "cm", height = 50, width = 50, 
         filename = "Pm_confound_pairs.png",
         path = here("Seagrass", "Plots"))
# No or very little correlation. All looks well.

# 7.8.4 Prior-posterior comparison ####
# Sample prior
Pm_confound_prior <- prior_samples(
  model = Pm_confound_mod,
  data = P_summary %>%
    select(Pm_mean, Pm_sd, Species,
           O2_mean_std, O2_sd_std, T_mean_std, T_sd_std,
           P_mean_std, S_std, M_std) %>%
    compose_data()
  )

Pm_confound_prior %>% 
  prior_posterior_draws(
    posterior_samples = Pm_confound_samples,
    group = P_summary %>% select(Species),
    parameters = c("alpha[Species]", "beta_O2[Species]",
                   "beta_T[Species]", "beta_P[Species]",
                   "beta_S[Species]", "beta_M[Species]",
                   "Pm_sigma[Species]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Species") %>%
  ggsave(units = "cm", height = 20, width = 40, 
         device = cairo_pdf, filename = "Pm_confound_prior_posterior.pdf",
         path = here("Seagrass", "Plots"))

# 7.8.5 Parameters ####
Pm_confound_prior_posterior <- Pm_confound_prior %>% 
  prior_posterior_draws(
    posterior_samples = Pm_confound_samples,
    group = P_summary %>% select(Species),
    parameters = c("alpha[Species]", "beta_O2[Species]",
                   "beta_T[Species]", "beta_P[Species]",
                   "beta_S[Species]", "beta_M[Species]",
                   "Pm_sigma[Species]"),
    format = "short"
  ) %>%
  filter(Species == "Halophila ovalis" & distribution == "prior" |
           distribution == "posterior") %>% # Remove redundant priors
  mutate( # Embed prior in Species
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct()
  ) %>%
  select(-distribution) %T>%
  print()

Pm_confound_prior_posterior %>% # Save
  write_rds(here("Seagrass", "RDS", "Pm_confound_prior_posterior.rds"))

# 7.8.6 Prediction ####
Pm_confound_prediction_O2 <- Pm_confound_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "O2_mean_std") %>%
  mutate( # The other confounders are held at their mean by exclusion
    Pm_mu = alpha + beta_O2 * O2_mean_std,
    Pm = rnorm( n() , Pm_mu , Pm_sigma )
  ) %>%
  group_by(O2_mean_std, Species) %>%
  median_qi(Pm_mu, Pm, .width = c(.5, .8, .9)) %T>%
  print()

Pm_confound_prediction_T <- Pm_confound_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "T_mean_std") %>%
  mutate(
    Pm_mu = alpha + beta_T * T_mean_std,
    Pm = rnorm( n() , Pm_mu , Pm_sigma )
  ) %>%
  group_by(T_mean_std, Species) %>%
  median_qi(Pm_mu, Pm, .width = c(.5, .8, .9)) %T>%
  print()

Pm_confound_prediction_P <- Pm_confound_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "P_mean_std") %>%
  mutate(
    Pm_mu = alpha + beta_P * P_mean_std,
    Pm = rnorm( n() , Pm_mu , Pm_sigma )
  ) %>%
  group_by(P_mean_std, Species) %>%
  median_qi(Pm_mu, Pm, .width = c(.5, .8, .9)) %T>%
  print()

Pm_confound_prediction_S <- Pm_confound_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "S_std") %>%
  mutate(
    Pm_mu = alpha + beta_S * S_std,
    Pm = rnorm( n() , Pm_mu , Pm_sigma )
  ) %>%
  group_by(S_std, Species) %>%
  median_qi(Pm_mu, Pm, .width = c(.5, .8, .9)) %T>%
  print()

Pm_confound_prediction_M <- Pm_confound_prior_posterior %>% 
  spread_continuous(data = P_summary,
                    group_name = "Species",
                    predictor_name = "M_std") %>%
  mutate(
    Pm_mu = alpha + beta_M * M_std,
    Pm = rnorm( n() , Pm_mu , Pm_sigma )
  ) %>%
  group_by(M_std, Species) %>%
  median_qi(Pm_mu, Pm, .width = c(.5, .8, .9)) %T>%
  print()

# Combine confounder predictions
Pm_confound_prediction <- bind_rows(
  Pm_confound_prediction_O2 %>%
    mutate(Confounder = "O2" %>% fct()) %>%
    rename(Predictor = O2_mean_std),
  Pm_confound_prediction_T %>%
    mutate(Confounder = "T" %>% fct()) %>%
    rename(Predictor = T_mean_std),
  Pm_confound_prediction_P %>%
    mutate(Confounder = "P" %>% fct()) %>%
    rename(Predictor = P_mean_std),
  Pm_confound_prediction_S %>%
    mutate(Confounder = "S" %>% fct()) %>%
    rename(Predictor = S_std),
  Pm_confound_prediction_M %>%
    mutate(Confounder = "M" %>% fct()) %>%
    rename(Predictor = M_std)
) %T>%
  print()

# Clean up
rm(Pm_confound_prediction_O2, Pm_confound_prediction_T, Pm_confound_prediction_P,
   Pm_confound_prediction_S, Pm_confound_prediction_M)

# Reverse standardisation
Pm_confound_prediction %<>%
  mutate(
    Predictor = P_summary %$% case_when(
      Confounder == "O2" ~ Predictor * sd(O2_mean) + mean(O2_mean),
      Confounder == "T" ~ Predictor * sd(T_mean) + mean(T_mean),
      Confounder == "P" ~ Predictor * sd(P_mean) + mean(P_mean),
      Confounder == "S" ~ Predictor * sd(S) + mean(S),
      Confounder == "M" ~ Predictor * sd(M) + mean(M)
    )
  ) %T>%
  print()

Pm_confound_prediction %>% # Save
  write_rds(here("Seagrass", "RDS", "Pm_confound_prediction.rds"))

# 8. Figures ####
# 8.1 Figure 1 ####
# Define custom theme
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = "cm"),
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
                 panel.spacing.x = unit(1, "cm"),
                 panel.spacing.y = unit(0.6, "cm"),
                 text = element_text(family = "Futura"))

Fig_1_Pm <- Pm_prediction %>%
  filter(Species != "Prior") %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_violin(data = P %>%
                  mutate(Species = Species %>% fct_relevel("Halophila ovalis")),
                aes(Day, Pm, fill = Species, colour = Species, group = ID),
                alpha = 0.2, position = "identity", width = 2) +
    geom_ribbon(aes(Day, ymin = Pm.lower, ymax = Pm.upper,
                    alpha = factor(.width), fill = Species)) +
    geom_line(aes(Day, Pm, colour = Species)) +
    stat_slab(data = Pm_prior_posterior %>%
                filter(Species != "Prior"),
              aes(-1, exp(log_alpha_mu), fill = Species),
              colour = NA, scale = -4, n = 2^10) +
    stat_slab(data = Pm_prior_posterior %>%
                filter(Species != "Prior"),
              aes(exp(log_mu_mu), -10, fill = Species),
              colour = NA, scale = 4, n = 2^10) +
    stat_slab(data = Pm_prior_posterior %>%
                filter(Species != "Prior"),
              aes(42, -exp(log_tau_mu), fill = Species),
              colour = NA, scale = 4, n = 2^10) +
    scale_colour_manual(values = c("#bdd268", "#4a7518"), guide = "none") +
    scale_fill_manual(values = c("#bdd268", "#4a7518"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_continuous(breaks = seq(0, 42, by = 7)) +
    scale_y_continuous(breaks = seq(-10, 50, by = 10),
                       labels = scales::label_number(style_negative = "minus")) +
    facet_grid(NA ~ Species, switch = "y") +
    labs(y = expression(italic(P)["max"]*" (µmol O"[2]*" g"^-1*" h"^-1*")"),
         x = "Detrital age (days)") +
    coord_cartesian(xlim = c(0, 42), ylim = c(-10, 50), 
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(strip.text = element_text(face = "italic"),
          strip.text.y = element_text(size = 20, colour = "transparent"),
          panel.spacing.x = unit(1.5, "cm"),
          plot.margin = margin(0.2, 1, 0.2, 0.2, unit = "cm"))

Fig_1_Pm # densities are not ideal because unbounded

require(ggh4x)
Fig_1_Pl <- Pl_prediction %>%
  filter(Species != "Prior") %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_violin(data = P %>%
                  mutate(Species = Species %>% fct_relevel("Halophila ovalis")),
                aes(Day, Pl, fill = Species, colour = Species, group = ID),
                alpha = 0.2, position = "identity", width = 2) +
    geom_ribbon(aes(Day, ymin = Pl.lower, ymax = Pl.upper,
                    alpha = factor(.width), fill = Species)) +
    geom_line(aes(Day, Pl, colour = Species)) +
    scale_colour_manual(values = c("#bdd268", "#4a7518"), guide = "none") +
    scale_fill_manual(values = c("#bdd268", "#4a7518"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    facet_wrap(~ Species, strip.position = "left", scales = "free") +
    facetted_pos_scales(
      y = list(
        Species == "Halophila ovalis" ~
          scale_y_continuous(limits = c(0, 2.5),
                             breaks = seq(0, 2.5, by = 0.5),
                             labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1, 0.1))),
        Species == "Amphibolis antarctica" ~
          scale_y_continuous(limits = c(-10, 50),
                             breaks = seq(-10, 50, by = 10),
                             labels = scales::label_number(style_negative = "minus"))
      )
    ) +
    scale_x_continuous(breaks = seq(0, 35, by = 7)) +
    labs(y = expression(italic(P)["max"]*" (µmol O"[2]*" leaf"^-1*" h"^-1*")"),
         x = "Detrital age (days)") +
    coord_cartesian(xlim = c(0, 35), expand = FALSE, clip = "off") +
    mytheme +
    theme(strip.text = element_text(face = "italic"),
          strip.text.y = element_text(size = 20, colour = "transparent"),
          plot.margin = margin(0.2, 1, 0.2, 0.2, unit = "cm"))

Fig_1_Pl
# I can't even remotely achieve the same effect with densities 
# as in Fig_1_Pm. The only solution would be separate plots, but 
# I do not think the additional insight gained from Pl is worthwhile. So
# I'll go with a reworked Fig_1_Pm.

# Build densities manually
alpha_dens <- Pm_prior_posterior %>%
  group_by(Species) %>%
  reframe(y = c(0, density(exp(log_alpha_mu), n = 2^10, from = 0, to = 50, bw = 60 * 0.02)$x, 50),
          x = c(0, density(exp(log_alpha_mu), n = 2^10, from = 0, to = 50, bw = 60 * 0.02)$y, 0)) %>%
  group_by(Species) %>% # Standardise area with Riemann sum (avoid manually added y[1]).
  mutate(x = x * 25 / ( sum(x) * ( y[3] - y[2] ) )) %>%
  ungroup() %T>%
  print()

mu_dens <- Pm_prior_posterior %>%
  group_by(Species) %>%
  reframe(x = c(0, density(exp(log_mu_mu), n = 2^10, from = 0, to = 49, bw = 42 * 0.02)$x, 49),
          y = c(0, density(exp(log_mu_mu), n = 2^10, from = 0, to = 49, bw = 42 * 0.02)$y, 0)) %>%
  group_by(Species) %>%
  mutate(y = y * 25 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

tau_dens <- Pm_prior_posterior %>%
  group_by(Species) %>%
  reframe(y = c(0, density(exp(log_tau_mu), n = 2^10, from = 0, to = 10, bw = 60 * 0.02)$x, 10),
          x = c(0, density(exp(log_tau_mu), n = 2^10, from = 0, to = 10, bw = 60 * 0.02)$y, 0)) %>%
  group_by(Species) %>%
  mutate(x = x * 25 / ( sum(x) * ( y[3] - y[2] ) )) %>%
  ungroup() %T>%
  print()

Fig_1 <- Pm_prediction %>%
  filter(Species != "Prior") %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_violin(data = P %>%
                  mutate(Species = Species %>% fct_relevel("Halophila ovalis")),
                aes(Day, Pm, fill = Species, group = ID),
                alpha = 0.7, position = "identity", 
                width = 3, colour = NA) +
    geom_ribbon(aes(Day, ymin = Pm.lower, ymax = Pm.upper,
                    alpha = factor(.width), fill = Species)) +
    geom_line(aes(Day, Pm, colour = Species)) +
    geom_polygon(data = alpha_dens %>% filter(Species != "Prior"),
                 aes(-(x+1.2), y, fill = Species), colour = NA) +
    geom_polygon(data = mu_dens %>% filter(Species != "Prior"),
                 aes(x, y-10, fill = Species), colour = NA) +
    geom_polygon(data = tau_dens %>% filter(Species != "Prior"),
                 aes(x+42, -y, fill = Species), colour = NA) +
    scale_colour_manual(values = c("#bdd268", "#4a7518"), guide = "none") +
    scale_fill_manual(values = c("#bdd268", "#4a7518"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_continuous(breaks = seq(0, 42, by = 7)) +
    scale_y_continuous(breaks = seq(-10, 50, by = 10),
                       labels = scales::label_number(style_negative = "minus")) +
    facet_grid(NA ~ Species, switch = "y") +
    labs(y = expression(italic(P)["max"]*" (µmol O"[2]*" g"^-1*" h"^-1*")"),
         x = "Detrital age (days)") +
    coord_cartesian(xlim = c(0, 42), ylim = c(-10, 50), 
                    expand = FALSE, clip = "off") +
    mytheme +
    theme(strip.text = element_text(face = "italic"),
          strip.text.y = element_text(colour = "transparent"),
          panel.spacing.x = unit(1.5, "cm"),
          plot.margin = margin(0.2, 1.5, 0.2, 0.2, unit = "cm"))

Fig_1

Fig_1 %>%
  ggsave(filename = "Fig_1.pdf", path = "Figures",
         device = cairo_pdf, height = 8, width = 20, units = "cm")

# 8.2 Figure S1 ####
# 8.2.1 Figure S1a ####
# Here I need to plot separately because of P
require(ggdensity) # bivariate density
Fig_S1a_O2 <- Pm_confound_prediction %>%
  filter(Species != "Prior" & Confounder == "O2") %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_hdr(data = P %>% group_by(ID) %>% mutate(O2 = sample(O2)), # randomly reshuffle within ID
             aes(O2, Pm, group = ID, fill = Species), colour = NA,
             n = 500, method = "mvnorm", probs = 0.999,
             position = "identity") +
    geom_ribbon(aes(Predictor, ymin = Pm_mu.lower, ymax = Pm_mu.upper,
                    alpha = factor(.width), fill = Species)) +
    geom_line(aes(Predictor, Pm_mu, colour = Species)) +
    scale_colour_manual(values = c("#bdd268", "#4a7518")) +
    scale_fill_manual(values = c("#bdd268", "#4a7518")) +
    scale_alpha_manual(values = c(0.7, 0.5, 0.4, 0.3), guide = "none") + # first alpha is for geom_hdr
    scale_x_continuous(breaks = seq(210, 250, 20)) +
    scale_y_continuous(breaks = seq(-10, 50, by = 10),
                       labels = scales::label_number(style_negative = "minus")) +
    # facet_grid(rows = vars(Species)) +
    labs(y = expression(italic(P)["max"]*" (µmol O"[2]*" g"^-1*" h"^-1*")"),
         x = expression("Initial O"[2]*" (µM)")) +
    coord_cartesian(ylim = c(-10, 50), xlim = c(210, 250), expand = FALSE, clip = "off") +
    mytheme +
    theme(plot.margin = margin(0.2, 1, 0.2, 0.2, unit = "cm"),
          legend.text = element_text(face = "italic"))

Fig_S1a_O2

Fig_S1a_T <- Pm_confound_prediction %>%
  filter(Species != "Prior" & Confounder == "T") %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_hdr(data = P %>% group_by(ID) %>% mutate(Temperature = sample(Temperature)), # randomly reshuffle within ID
             aes(Temperature, Pm, group = ID, fill = Species), colour = NA,
             n = 500, method = "mvnorm", probs = 0.999,
             position = "identity") +
    geom_ribbon(aes(Predictor, ymin = Pm_mu.lower, ymax = Pm_mu.upper,
                    alpha = factor(.width), fill = Species)) +
    geom_line(aes(Predictor, Pm_mu, colour = Species)) +
    scale_colour_manual(values = c("#bdd268", "#4a7518")) +
    scale_fill_manual(values = c("#bdd268", "#4a7518")) +
    scale_alpha_manual(values = c(0.7, 0.5, 0.4, 0.3), guide = "none") + # first alpha is for geom_hdr
    scale_x_continuous(breaks = seq(16, 20, 2)) +
    labs(x = "Temp. (°C)") +
    coord_cartesian(ylim = c(-10, 50), xlim = c(16, 20), expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.text.y = element_blank(), # Emulate faceting
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(0.2, 1, 0.2, 0, unit = "cm"),
          legend.text = element_text(face = "italic"))

Fig_S1a_T

Fig_S1a_P <- Pm_confound_prediction %>%
  filter(Species != "Prior" & Confounder == "P") %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_violin(data = P,
                aes(Pressure, Pm, fill = Species, group = ID),
                alpha = 0.7, position = "identity", 
                width = 7, colour = NA) +
    geom_ribbon(aes(Predictor, ymin = Pm_mu.lower, ymax = Pm_mu.upper,
                    alpha = factor(.width), fill = Species)) +
    geom_line(aes(Predictor, Pm_mu, colour = Species)) +
    scale_colour_manual(values = c("#bdd268", "#4a7518")) +
    scale_fill_manual(values = c("#bdd268", "#4a7518")) +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_continuous(breaks = seq(1000, 1030, 15)) +
    labs(x = "Pressure (hPa)") +
    coord_cartesian(ylim = c(-10, 50), xlim = c(1000, 1030), expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.text.y = element_blank(), # Emulate faceting
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(0.2, 1, 0.2, 0, unit = "cm"),
          legend.text = element_text(face = "italic"))

Fig_S1a_P

Fig_S1a_S <- Pm_confound_prediction %>%
  filter(Species != "Prior" & Confounder == "S") %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_violin(data = P,
                aes(Salinity, Pm, fill = Species, group = ID),
                alpha = 0.7, position = "identity", 
                width = 0.5, colour = NA) +
    geom_ribbon(aes(Predictor, ymin = Pm_mu.lower, ymax = Pm_mu.upper,
                    alpha = factor(.width), fill = Species)) +
    geom_line(aes(Predictor, Pm_mu, colour = Species)) +
    scale_colour_manual(values = c("#bdd268", "#4a7518")) +
    scale_fill_manual(values = c("#bdd268", "#4a7518")) +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_continuous(breaks = 34:36) +
    labs(x = "Salinity (‰)") +
    coord_cartesian(ylim = c(-10, 50), xlim = c(34, 36), expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.text.y = element_blank(), # Emulate faceting
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(0.2, 1, 0.2, 0, unit = "cm"),
          legend.text = element_text(face = "italic"))

Fig_S1a_S

Fig_S1a_M <- Pm_confound_prediction %>%
  filter(Species != "Prior" & Confounder == "M") %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_violin(data = P,
                aes(Mass, Pm, fill = Species, group = ID),
                alpha = 0.7, position = "identity", 
                width = 0.4, colour = NA) +
    geom_ribbon(aes(Predictor, ymin = Pm_mu.lower, ymax = Pm_mu.upper,
                    alpha = factor(.width), fill = Species)) +
    geom_line(aes(Predictor, Pm_mu, colour = Species)) +
    scale_colour_manual(values = c("#bdd268", "#4a7518")) +
    scale_fill_manual(values = c("#bdd268", "#4a7518")) +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_continuous(breaks = seq(0, 1.5, 0.5),
                       labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1))) +
    labs(x = "Mass (g)") +
    coord_cartesian(ylim = c(-10, 50), xlim = c(0, 1.5), expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.text.y = element_blank(), # Emulate faceting
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(0.2, 0.5, 0.2, 0, unit = "cm"),
          legend.text = element_text(face = "italic"))

Fig_S1a_M

# 8.2.2 Figure S1b ####
Pm_confound_beta <- Pm_confound_prior_posterior %>%
  select(-c(alpha, Pm_sigma)) %>%
  pivot_longer(cols = -c(starts_with("."), Species),
               names_to = "Parameter") %>%
  filter(Species != "Prior") %>%
  mutate(
    Species = Species %>% fct_drop() %>%
      fct_relevel("Amphibolis antarctica"),
    Parameter = Parameter %>% 
      fct_relevel("beta_O2", "beta_T", "beta_P", "beta_S", "beta_M"),
    value = P_summary %$% case_when( # Reverse standardisation
      Parameter == "beta_O2" ~ value / sd(O2_mean), # This typically involves multiplying by sd
      Parameter == "beta_T" ~ value / sd(T_mean), # but because this is a slope, it's inverse
      Parameter == "beta_P" ~ value / sd(P_mean),
      Parameter == "beta_S" ~ value / sd(S),
      Parameter == "beta_M" ~ value / sd(M)
    )
  ) %T>%
  print()

require(ggridges)
Fig_S1b <- Pm_confound_beta %>%
  ggplot() +
    stat_density_ridges(aes(value, Species, fill = Species),
                        colour = NA, n = 2^10, scale = 2, 
                        from = c(-0.8, -8, -0.8, -8, -15),
                        to = c(0.8, 8, 0.8, 8, 15),
                        bandwidth = c(1.6, 16, 1.6, 16, 30)*0.02) +
    geom_vline(xintercept = 0) +
    geom_text(
      data = tibble(
        Parameter = Pm_confound_beta %$% levels(Parameter) %>% fct(),
        label = c("italic(P)[max]", rep(NA, 4))
      ),
      aes(x = -1.84, y = 4, label = label),
      family = "Futura", size = 12, size.unit = "pt",
      hjust = 0, vjust = 1, parse = TRUE
    ) +
    scale_fill_manual(values = c("#4a7518", "#bdd268"),
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(
      ~ Parameter, scales = "free",
        labeller = labeller(
          Parameter = as_labeller(
            c("beta_O2" = "italic('β')['O'[2]]*' (µM'^-1*')'",
              "beta_T" = "italic('β')['T']*' (°C'^-1*')'",
              "beta_P" = "italic('β')['P']*' (hPa'^-1*')'",
              "beta_S" = "italic('β')['Sal']*' (‰'^-1*')'",
              "beta_M" = "italic('β')['Mass']*' (g'^-1*')'"),
            label_parsed
          )
        )
    ) +
    # facet_grid2(
    #   "Pmax" ~ Parameter, scales = "free", switch = "y",
    #   strip = strip_nested(
    #     text_y = element_text(angle = 0, hjust = 0, vjust = 1)
    #   ),
    #   labeller = labeller(
    #     .rows = as_labeller(
    #       c(Pmax = "italic(P)[max]"),
    #       label_parsed
    #     ),
    #     Parameter = as_labeller(
    #       c("beta_O2" = "italic('β')['O'[2]]*' (µM'^-1*')'",
    #         "beta_T" = "italic('β')['T']*' (°C'^-1*')'",
    #         "beta_P" = "italic('β')['P']*' (hPa'^-1*')'",
    #         "beta_S" = "italic('β')['Sal']*' (‰'^-1*')'",
    #         "beta_M" = "italic('β')['Mass']*' (g'^-1*')'"),
    #       label_parsed
    #     )
    #   )
    # ) +
    facetted_pos_scales(
      x = list(
        Parameter == "beta_O2" ~ 
          scale_x_continuous(limits = c(-0.8, 0.8),
                             oob = scales::oob_keep,
                             breaks = seq(-0.8, 0.8, by = 0.8),
                             labels = scales::label_number(style_negative = "minus",
                                                           accuracy = c(0.1, 1, 0.1))),
        Parameter == "beta_T" ~ 
          scale_x_continuous(limits = c(-8, 8),
                             oob = scales::oob_keep,
                             breaks = seq(-8, 8, by = 8),
                             labels = scales::label_number(style_negative = "minus")),
        Parameter == "beta_P" ~ 
          scale_x_continuous(limits = c(-0.8, 0.8),
                             oob = scales::oob_keep,
                             breaks = seq(-0.8, 0.8, by = 0.8),
                             labels = scales::label_number(style_negative = "minus",
                                                           accuracy = c(0.1, 1, 0.1))),
        Parameter == "beta_S" ~ 
          scale_x_continuous(limits = c(-8, 8),
                             oob = scales::oob_keep,
                             breaks = seq(-8, 8, by = 8),
                             labels = scales::label_number(style_negative = "minus")),
        Parameter == "beta_M" ~ 
          scale_x_continuous(limits = c(-15, 15),
                             oob = scales::oob_keep,
                             breaks = seq(-15, 15, by = 15),
                             labels = scales::label_number(style_negative = "minus"))
      )
    ) +
    coord_cartesian(ylim = c(1, 4), expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.text = element_text(face = "italic"))

Fig_S1b
# Warnings are due to deliberate NAs in geom_text and can be safely ignored

# 8.2.3 Figure S1c ####
Pm_confound <- Pm_prior_posterior %>%
  select(-c(ends_with("mu"), Pm_sigma)) %>%
  pivot_longer(cols = -c(starts_with("."), Species),
               names_to = "Parameter") %>%
  separate(Parameter, into = c("Parameter", "Confounder"),
           sep = "_(?=[^_]+$)") %>%
  filter(Species != "Prior") %>%
  mutate(
    Species = Species %>% fct_drop() %>%
      fct_relevel("Amphibolis antarctica"),
    Parameter = Parameter %>% 
      fct_relevel("log_alpha", "log_mu"),
    Confounder = Confounder %>%
      fct_relevel("O2", "T", "P", "S"),
    value = P_summary %$% case_when( # Reverse standardisation as before
      Confounder == "O2" ~ value / sd(O2_mean),
      Confounder == "T" ~ value / sd(T_mean),
      Confounder == "P" ~ value / sd(P_mean),
      Confounder == "S" ~ value / sd(S),
      Confounder == "M" ~ value / sd(M)
    )
  ) %T>%
  print()

# I cannot have completely free scales in facet_grid, so need
# to plot separately.
Fig_S1c_alpha <- Pm_confound %>%
  filter(Parameter == "log_alpha") %>%
  ggplot() +
    stat_density_ridges(aes(value, Species, fill = Species),
                        colour = NA, n = 2^10, scale = 2,
                        from = c(-0.04, -2.5, -0.3, -5, -1.5),
                        to = c(0.04, 2.5, 0.3, 5, 1.5),
                        bandwidth = c(0.08, 5, 0.6, 10, 3)*0.02) +
    geom_vline(xintercept = 0) +
    geom_text(
      data = tibble(
        Confounder = Pm_confound %$% levels(Confounder) %>% fct(),
        label = c("'ln '*italic('α')", rep(NA, 4))
      ),
      aes(x = -0.092, y = 4, label = label),
      family = "Futura", size = 12, size.unit = "pt",
      hjust = 0, vjust = 1, parse = TRUE
    ) +
    scale_fill_manual(values = c("#4a7518", "#bdd268"),
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(
      ~ Confounder, scales = "free",
        labeller = labeller(
          Confounder = as_labeller(
            c("O2" = "italic('β')['O'[2]]*' (µM'^-1*')'",
              "T" = "italic('β')['T']*' (°C'^-1*')'",
              "P" = "italic('β')['P']*' (hPa'^-1*')'",
              "S" = "italic('β')['Sal']*' (‰'^-1*')'",
              "M" = "italic('β')['Mass']*' (g'^-1*')'"),
            label_parsed
          )
        )
    ) +
    # facet_grid2(
    #   Parameter ~ Confounder, scales = "free", switch = "y",
    #   strip = strip_nested(
    #     text_y = element_text(angle = 0, hjust = 0, vjust = 1)
    #   ),
    #   labeller = labeller(
    #     Confounder = as_labeller(
    #       c("O2" = "italic('β')['O'[2]]*' (µM'^-1*')'",
    #         "T" = "italic('β')['T']*' (°C'^-1*')'",
    #         "P" = "italic('β')['P']*' (hPa'^-1*')'",
    #         "S" = "italic('β')['Sal]*' (‰'^-1*')'",
    #         "M" = "italic('β')['Mass']*' (g'^-1*')'"),
    #       label_parsed
    #     ),
    #     Parameter = as_labeller(
    #       c("log_alpha" = "'ln '*italic('α')",
    #         "log_mu" = "'ln '*italic('μ')",
    #         "log_tau" = "'ln '*italic('τ')"),
    #       label_parsed
    #     )
    #   )
    # ) +
    facetted_pos_scales(
      x = list(
        Confounder == "O2" ~ 
          scale_x_continuous(limits = c(-0.04, 0.04),
                             oob = scales::oob_keep,
                             breaks = seq(-0.04, 0.04, by = 0.04),
                             labels = scales::label_number(style_negative = "minus",
                                                           accuracy = c(0.01, 1, 0.01))),
        Confounder == "T" ~ 
          scale_x_continuous(limits = c(-2.5, 2.5),
                             oob = scales::oob_keep,
                             breaks = seq(-2.5, 2.5, by = 2.5),
                             labels = scales::label_number(style_negative = "minus",
                                                           accuracy = c(0.1, 1, 0.1))),
        Confounder == "P" ~ 
          scale_x_continuous(limits = c(-0.3, 0.3),
                             oob = scales::oob_keep,
                             breaks = seq(-0.3, 0.3, by = 0.3),
                             labels = scales::label_number(style_negative = "minus",
                                                           accuracy = c(0.1, 1, 0.1))),
        Confounder == "S" ~ 
          scale_x_continuous(limits = c(-5, 5),
                             oob = scales::oob_keep,
                             breaks = seq(-5, 5, by = 5),
                             labels = scales::label_number(style_negative = "minus")),
        Confounder == "M" ~ 
          scale_x_continuous(limits = c(-1.5, 1.5),
                             oob = scales::oob_keep,
                             breaks = seq(-1.5, 1.5, by = 1.5),
                             labels = scales::label_number(style_negative = "minus",
                                                           accuracy = c(0.1, 1, 0.1)))
      )
    ) +
    coord_cartesian(ylim = c(1, 4), expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.text = element_text(face = "italic"))

Fig_S1c_alpha

Fig_S1c_mu <- Pm_confound %>%
  filter(Parameter == "log_mu") %>%
  ggplot() +
    stat_density_ridges(aes(value, Species, fill = Species),
                        colour = NA, n = 2^10, scale = 2,
                        from = c(-0.15, -1.5, -0.2, -6, -3),
                        to = c(0.15, 1.5, 0.2, 6, 3),
                        bandwidth = c(0.3, 3, 0.4, 12, 6)*0.02) +
    geom_vline(xintercept = 0) +
    geom_text(
      data = tibble(
        Confounder = Pm_confound %$% levels(Confounder) %>% fct(),
        label = c("'ln '*italic('μ')", rep(NA, 4))
      ),
      aes(x = -0.345, y = 4, label = label),
      family = "Futura", size = 12, size.unit = "pt",
      hjust = 0, vjust = 1, parse = TRUE
    ) +
    scale_fill_manual(values = c("#4a7518", "#bdd268"),
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(
      ~ Confounder, scales = "free",
        labeller = labeller(
          Confounder = as_labeller(
            c("O2" = "italic('β')['O'[2]]*' (µM'^-1*')'",
              "T" = "italic('β')['T']*' (°C'^-1*')'",
              "P" = "italic('β')['P']*' (hPa'^-1*')'",
              "S" = "italic('β')['Sal']*' (‰'^-1*')'",
              "M" = "italic('β')['Mass']*' (g'^-1*')'"),
            label_parsed
          )
        )
    ) +
    facetted_pos_scales(
      x = list(
        Confounder == "O2" ~ 
          scale_x_continuous(limits = c(-0.15, 0.15),
                             oob = scales::oob_keep,
                             breaks = seq(-0.15, 0.15, by = 0.15),
                             labels = scales::label_number(style_negative = "minus",
                                                           accuracy = c(0.01, 1, 0.01))),
        Confounder == "T" ~ 
          scale_x_continuous(limits = c(-1.5, 1.5),
                             oob = scales::oob_keep,
                             breaks = seq(-1.5, 1.5, by = 1.5),
                             labels = scales::label_number(style_negative = "minus",
                                                           accuracy = c(0.1, 1, 0.1))),
        Confounder == "P" ~ 
          scale_x_continuous(limits = c(-0.2, 0.2),
                             oob = scales::oob_keep,
                             breaks = seq(-0.2, 0.2, by = 0.2),
                             labels = scales::label_number(style_negative = "minus",
                                                           accuracy = c(0.1, 1, 0.1))),
        Confounder == "S" ~ 
          scale_x_continuous(limits = c(-6, 6),
                             oob = scales::oob_keep,
                             breaks = seq(-6, 6, by = 6),
                             labels = scales::label_number(style_negative = "minus")),
        Confounder == "M" ~ 
          scale_x_continuous(limits = c(-3, 3),
                             oob = scales::oob_keep,
                             breaks = seq(-3, 3, by = 3),
                             labels = scales::label_number(style_negative = "minus"))
      )
    ) +
    coord_cartesian(ylim = c(1, 4), expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.text = element_text(face = "italic"))

Fig_S1c_mu

Fig_S1c_tau <- Pm_confound %>%
  filter(Parameter == "log_tau") %>%
  ggplot() +
    stat_density_ridges(aes(value, Species, fill = Species),
                        colour = NA, n = 2^10, scale = 2,
                        from = c(-0.5, -4, -0.5, -7, -12),
                        to = c(0.5, 4, 0.5, 7, 12),
                        bandwidth = c(1, 8, 1, 14, 24)*0.02) +
    geom_vline(xintercept = 0) +
    geom_text(
      data = tibble(
        Confounder = Pm_confound %$% levels(Confounder) %>% fct(),
        label = c("'ln '*italic('τ')", rep(NA, 4))
      ),
      aes(x = -1.15, y = 4, label = label),
      family = "Futura", size = 12, size.unit = "pt",
      hjust = 0, vjust = 1, parse = TRUE
    ) +
    scale_fill_manual(values = c("#4a7518", "#bdd268"),
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(
      ~ Confounder, scales = "free",
        labeller = labeller(
          Confounder = as_labeller(
            c("O2" = "italic('β')['O'[2]]*' (µM'^-1*')'",
              "T" = "italic('β')['T']*' (°C'^-1*')'",
              "P" = "italic('β')['P']*' (hPa'^-1*')'",
              "S" = "italic('β')['Sal']*' (‰'^-1*')'",
              "M" = "italic('β')['Mass']*' (g'^-1*')'"),
            label_parsed
          )
        )
    ) +
    facetted_pos_scales(
      x = list(
        Confounder == "O2" ~ 
          scale_x_continuous(limits = c(-0.5, 0.5),
                             oob = scales::oob_keep,
                             breaks = seq(-0.5, 0.5, by = 0.5),
                             labels = scales::label_number(style_negative = "minus",
                                                           accuracy = c(0.1, 1, 0.1))),
        Confounder == "T" ~ 
          scale_x_continuous(limits = c(-4, 4),
                             oob = scales::oob_keep,
                             breaks = seq(-4, 4, by = 4),
                             labels = scales::label_number(style_negative = "minus")),
        Confounder == "P" ~ 
          scale_x_continuous(limits = c(-0.5, 0.5),
                             oob = scales::oob_keep,
                             breaks = seq(-0.5, 0.5, by = 0.5),
                             labels = scales::label_number(style_negative = "minus",
                                                           accuracy = c(0.1, 1, 0.1))),
        Confounder == "S" ~ 
          scale_x_continuous(limits = c(-7, 7),
                             oob = scales::oob_keep,
                             breaks = seq(-7, 7, by = 7),
                             labels = scales::label_number(style_negative = "minus")),
        Confounder == "M" ~ 
          scale_x_continuous(limits = c(-12, 12),
                             oob = scales::oob_keep,
                             breaks = seq(-12, 12, by = 12),
                             labels = scales::label_number(style_negative = "minus"))
      )
    ) +
    coord_cartesian(ylim = c(1, 4), expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.text = element_text(face = "italic"))

Fig_S1c_tau

# 8.2.4 Combine panels ####
require(patchwork)
Fig_S1 <- (
  ( Fig_S1a_O2 | 
      (Fig_S1a_T + theme(legend.position = "none")) | 
      (Fig_S1a_P + theme(legend.position = "none")) | 
      (Fig_S1a_S + theme(legend.position = "none")) | 
      (Fig_S1a_M + theme(legend.position = "none")) ) /
    (Fig_S1b + theme(legend.position = "none")) / 
    (Fig_S1c_alpha + theme(legend.position = "none",
                           strip.text = element_blank())) / 
    (Fig_S1c_mu + theme(legend.position = "none",
                        strip.text = element_blank())) / 
    (Fig_S1c_tau + theme(legend.position = "none",
                         strip.text = element_blank()))
) +
  plot_layout(heights = c(1, rep(0.25, 4)))

Fig_S1 %>%
  ggsave(filename = "Fig_S1.pdf", path = "Figures",
         device = cairo_pdf, height = 18, width = 20, units = "cm")

# 9. Table 1 ####
require(glue)
Table_1 <- Pm_contrast %>%
  mutate(Response = "Mass-based" %>% fct()) %>%
  bind_rows(
    Pl_contrast %>%
      mutate(Response = "Leaf-based" %>% fct())
  ) %>%
  select(-starts_with(".")) %>%
  group_by(Response, Parameter) %>%
  summarise(
    across(
      everything(), 
      list(
        mean = mean, 
        sd = sd, 
        median = median
      )
    ),
    n = n(),
    P = pmax( mean( diff > 0 ) , mean( diff < 0 ) ) %>% signif(2)
  ) %>%
  ungroup() %>%
  mutate(
    `Halophila ovalis` = glue(
      "{signif(`Halophila ovalis_mean`, 2)} ± {signif(`Halophila ovalis_sd`, 2)} ({signif(`Halophila ovalis_median`, 2)})"
    ),
    `Amphibolis antarctica` = glue(
      "{signif(`Amphibolis antarctica_mean`, 2)} ± {signif(`Amphibolis antarctica_sd`, 2)} ({signif(`Amphibolis antarctica_median`, 2)})"
    ),
    diff = glue(
      "{signif(diff_mean, 2)} ± {signif(diff_sd, 2)} ({signif(diff_median, 2)})"
    ),
    ratio = glue(
      "{signif(ratio_mean, 2)} ± {signif(ratio_sd, 2)} ({signif(ratio_median, 2)})"
    )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median") | n)) %T>%
  print()

Table_1 %>%
  write_csv(here("Tables", "Table_1.csv"))

require(officer)
read_docx() %>%
  body_add_table(value = Table_1) %>%
  print(target = here("Tables", "Table_1.docx"))











# 7.6.7 Visualisation ####
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
                          legend.position = "inside",
                          legend.position.inside = c(-0.1, 0.72),
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

ggsave(plot = Fig_S1, filename = "Fig_S1.pdf", device = cairo_pdf, 
       path = "Figures", width = 20, height = 10, units = "cm")

rm(Fig_S1, Fig_S1_beta, Fig_S1_O2, Fig_S1_T, Fig_S1_P, Fig_S1_S, Fig_S1_M,
   Pm_confound_mu_summary, Pm_confound_prior_posterior, Pm_confound_stan)



# Reorder factors for visualisation
Pm_prior_posterior %<>%
  mutate(Group = Group %>%
           fct_relevel("Seagrasses", "Amphibolis antarctica", "Halophila ovalis"),
         .variable = .variable %>% fct_relevel("alpha", "tau", "k"))

Fig_1a_top <- Pm_prior_posterior %>%
  filter(.variable != "Pm_sigma") %>%
  ggplot() +
    geom_density_ridges(aes(.value, Group, colour = Group, fill = Group),
                        # quantile_lines = TRUE, quantiles = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95),
                        alpha = 0.5, scale = 2, rel_min_height = 0.002, 
                        bandwidth = c(40*0.02, 8*0.02, 0.8*0.02, 56*0.02),
                        from = c(10, rep(0, 3)), to = c(50, 8, 0.8, 56)) +
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
      .variable == "alpha" ~ scale_x_continuous(limits = c(10, 50),
                                                breaks = seq(10, 50, by = 20),
                                                oob = scales::oob_keep),
      .variable == "tau" ~ scale_x_continuous(limits = c(0, 8),
                                              breaks = seq(0, 8, by = 4),
                                              oob = scales::oob_keep),
      .variable == "k" ~ scale_x_continuous(limits = c(0, 0.8),
                                            breaks = seq(0, 0.8, by = 0.4),
                                            oob = scales::oob_keep,
                                            labels = scales::label_number(accuracy = c(1, rep(0.1, 2)))),
      .variable == "mu" ~ scale_x_continuous(limits = c(0, 56),
                                             breaks = seq(0, 56, by = 28),
                                             oob = scales::oob_keep)
    )) +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.spacing = unit(1, "cm"))
Fig_1a_top

# 7.7.8 Parameters and differences ####
# Parameters for Table 1
Pm_prior_posterior %>%
  filter(Group != "Prior" & .variable != "Pm_sigma") %>%
  group_by(Group, .variable) %>%
  summarise(mean = mean(.value),
            sd = sd(.value),
            n = n()) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(mean, digits = 2), signif(mean, digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        )
        ) %>%
  pivot_wider(names_from = Group, values_from = rounded) %>%
  write.csv(here("Tables", "Pm_para.csv"), row.names = FALSE)

# Calculate differences
Pm_diff <- Pm_prior_posterior %>%
  filter(!Group %in% c("Prior", "Seagrasses") & .variable != "Pm_sigma") %>%
  pivot_wider(names_from = c(Group, .variable), values_from = .value) %>%
  mutate(delta_alpha = `Halophila ovalis_alpha` - `Amphibolis antarctica_alpha`,
         delta_tau = `Halophila ovalis_tau` - `Amphibolis antarctica_tau`,
         delta_k = `Halophila ovalis_k` - `Amphibolis antarctica_k`,
         delta_mu = `Halophila ovalis_mu` - `Amphibolis antarctica_mu`) %>%
  select(starts_with("."), starts_with("delta")) %>%
  pivot_longer(cols = c(delta_alpha, delta_tau, delta_k, delta_mu),
               names_to = "Parameter", values_to = "Difference", names_prefix = "delta_") %>%
  mutate(Parameter = fct_relevel(Parameter, "alpha", "tau", "k")) %T>%
  print()

# Differences and probabilities for Table 1
Pm_diff %>%
  group_by(Parameter) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P_more = mean(Difference > 0),
            P_less = mean(Difference < 0),
            n = n()) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(abs(mean), digits = 2), signif(abs(mean), digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        ),
         P = signif(pmax(P_less, P_more), digits = 2)) %>%
  write.csv(here("Tables", "Pm_diff.csv"), row.names = FALSE)

# 7.7.9 Prediction ####
# calculate Pm_mu and Pm_obs from parameters
Pm_prediction <- Pm_prior_posterior %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(Pm_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau,
         Pm_obs = rnorm( n() , Pm_mu , Pm_sigma ))

Pm_prediction_summary <- Pm_prediction %>%
  group_by(Day, Group) %>%
  median_qi(Pm_mu, Pm_obs, .width = c(.5, .8, .9))

Fig_1a_bottom <- ggplot() +
                  geom_hline(yintercept = 0) +
                  geom_violin(data = P,
                              aes(Day, Pm, fill = Species,
                                  colour = Species, group = ID),
                              alpha = 0.2, position = "identity", width = 2) +
                  # too much overplotting, so I reduced the prior to the 0.05 and 0.95 quantiles
                  geom_ribbon(data = Pm_prediction_summary %>% filter(Group == "Prior", .width == 0.9),
                              aes(Day, ymin = Pm_mu.lower, ymax = Pm_mu.upper), fill = NA, colour = "#b5b8ba") +
                  geom_line(data = Pm_prediction_summary %>% filter(Group != "Prior"),
                            aes(Day, Pm_mu, colour = fct_relevel(Group, "Halophila ovalis",
                                                                "Amphibolis antarctica"))) +
                  geom_ribbon(data = Pm_prediction_summary %>% filter(Group != "Prior"),
                              aes(Day, ymin = Pm_mu.lower, ymax = Pm_mu.upper,
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
Fig_1a_bottom


# Reorder factors for visualisation
Pl_prior_posterior %<>%
  mutate(Group = Group %>%
           fct_relevel("Seagrasses", "Amphibolis antarctica", "Halophila ovalis"),
         .variable = .variable %>% fct_relevel("alpha", "tau", "k"))

# plot without alpha and tau since they are too divergent on the leaf scale
Fig_1b_top_right <- Pl_prior_posterior %>%
  filter(.variable %in% c("k", "mu")) %>%
  ggplot() +
    geom_density_ridges(aes(.value, Group, colour = Group, fill = Group),
                        # quantile_lines = TRUE, quantiles = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95),
                        alpha = 0.5, scale = 2, rel_min_height = 0.002, 
                        bandwidth = c(0.8*0.02, 56*0.02),
                        from = c(0, 0), to = c(0.8, 56)) +
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
      .variable == "k" ~ scale_x_continuous(limits = c(0, 0.8),
                                            breaks = seq(0, 0.8, by = 0.4),
                                            oob = scales::oob_keep,
                                            labels = scales::label_number(accuracy = c(1, rep(0.1, 2)))),
      .variable == "mu" ~ scale_x_continuous(limits = c(0, 56),
                                             breaks = seq(0, 56, by = 28),
                                             oob = scales::oob_keep)
    )) +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.spacing = unit(1, "cm"))
Fig_1b_top_right

# 7.8.5 Parameters and differences ####
# Parameters for Table 1
# correct scale for alpha, tau and Pl_sigma
Pl_prior_posterior %<>%
  mutate(.value = if_else(Group == "Halophila ovalis" & 
                            .variable %in% c("alpha", "tau", "Pl_sigma"),
                          .value / 10, .value))

Pl_prior_posterior %>%
  filter(Group != "Prior" & .variable != "Pl_sigma") %>%
  group_by(Group, .variable) %>%
  summarise(mean = mean(.value),
            sd = sd(.value),
            n = n()) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(mean, digits = 2), signif(mean, digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        )
        ) %>%
  pivot_wider(names_from = Group, values_from = rounded) %>%
  write.csv(here("Tables", "Pl_para.csv"), row.names = FALSE)

# Calculate differences
Pl_diff <- Pl_prior_posterior %>%
  filter(!Group %in% c("Prior", "Seagrasses") & .variable != "Pl_sigma") %>%
  pivot_wider(names_from = c(Group, .variable), values_from = .value) %>%
  mutate(delta_alpha = `Halophila ovalis_alpha` - `Amphibolis antarctica_alpha`,
         delta_tau = `Halophila ovalis_tau` - `Amphibolis antarctica_tau`,
         delta_k = `Halophila ovalis_k` - `Amphibolis antarctica_k`,
         delta_mu = `Halophila ovalis_mu` - `Amphibolis antarctica_mu`) %>%
  select(starts_with("."), starts_with("delta")) %>%
  pivot_longer(cols = c(delta_alpha, delta_tau, delta_k, delta_mu),
               names_to = "Parameter", values_to = "Difference", names_prefix = "delta_") %>%
  mutate(Parameter = fct_relevel(Parameter, "alpha", "tau", "k")) %T>%
  print()

# Differences and probabilities for Table 1
Pl_diff %>%
  group_by(Parameter) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P_more = mean(Difference > 0),
            P_less = mean(Difference < 0),
            n = n()) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(abs(mean), digits = 2), signif(abs(mean), digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        ),
         P = signif(pmax(P_less, P_more), digits = 2)) %>%
  write.csv(here("Tables", "Pl_diff.csv"), row.names = FALSE)


# Differences and probabilities for text: comparing k and mu between mass- and leaf-based
Pm_Pl_diff <- Pm_prior_posterior %>%
  filter(Group != "Prior" & .variable %in% c("k", "mu")) %>%
  mutate(Scale = "Mass") %>%
  bind_rows(Pl_prior_posterior %>%
              filter(Group != "Prior" & .variable %in% c("k", "mu")) %>%
              mutate(Scale = "Leaf")) %>%
  pivot_wider(names_from = c(Scale, .variable), values_from = .value) %>%
  mutate(delta_k = Mass_k - Leaf_k,
         delta_mu = Mass_mu - Leaf_mu) %>%
  select(starts_with("."), Group, starts_with("delta")) %>%
  pivot_longer(cols = c(delta_k, delta_mu),
               names_to = "Parameter", values_to = "Difference", 
               names_prefix = "delta_") %T>%
  print()

Pm_Pl_diff %>%
  group_by(Group, Parameter) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P_more = mean(Difference > 0),
            P_less = mean(Difference < 0),
            n = n()) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(abs(mean), digits = 2), signif(abs(mean), digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        ),
  P = signif(pmax(P_less, P_more), digits = 2)) %>%
  write.csv(here("Tables", "Pm_Pl_diff.csv"), row.names = FALSE)

# 7.8.6 Prediction ####
# calculate P_mu from parameters
Pl_prediction <- Pl_prior_posterior %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  expand_grid(Day = P_summary %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(Pl_mu = ( alpha + tau ) / ( 1 + exp( k * ( Day - mu ) ) ) - tau,
         Pl_obs = rnorm( n() , Pl_mu , Pl_sigma ))
# Ignore warning because NAs are known to be produced for new seagrasses since
# alpha and tau were not estimated for new seagrasses.

Pl_prediction_summary <- Pl_prediction %>%
  group_by(Day, Group) %>%
  median_qi(Pl_mu, Pl_obs, .width = c(.5, .8, .9))

Fig_1b_bottom_left <- ggplot() +
                  geom_violin(data = P %>%
                                filter(Species == "Halophila ovalis"),
                              aes(Day, Pl, group = ID), colour = "#bdd268",
                              fill = "#bdd268", alpha = 0.2,
                              position = "identity", width = 4) +
                  geom_line(data = Pl_prediction_summary %>%
                              filter(Group == "Halophila ovalis"),
                            aes(Day, Pl_mu), colour = "#bdd268") +
                  geom_ribbon(data = Pl_prediction_summary %>%
                                filter(Group == "Halophila ovalis"),
                              aes(Day, ymin = Pl_mu.lower, ymax = Pl_mu.upper,
                                  alpha = factor(.width)), fill = "#bdd268", colour = NA) +
                  labs(y = expression(italic(P)["max"]*" (µmol O"[2]*" leaf"^-1*" h"^-1*")"),
                       x = "Detrital age (d)") +
                  scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                  scale_x_continuous(breaks = seq(0, 35, by = 5)) +
                  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5),
                                     labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1, 0.1))) +
                  coord_cartesian(xlim = c(-1.8, 35), ylim = c(0, 2.5), expand = FALSE, clip = "off") +
                  mytheme
Fig_1b_bottom_left

Fig_1b_bottom_right <- ggplot() +
                  geom_hline(yintercept = 0) +
                  geom_violin(data = P %>%
                                filter(Species == "Amphibolis antarctica"),
                              aes(Day, Pl, group = ID), colour = "#4a7518",
                              fill = "#4a7518", alpha = 0.2,
                              position = "identity", width = 4) +
                  geom_line(data = Pl_prediction_summary %>%
                              filter(Group == "Amphibolis antarctica"),
                            aes(Day, Pl_mu), colour = "#4a7518") +
                  geom_ribbon(data = Pl_prediction_summary %>%
                                filter(Group == "Amphibolis antarctica"),
                              aes(Day, ymin = Pl_mu.lower, ymax = Pl_mu.upper,
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
Fig_1b_bottom_right

Fig_1 <- Fig_1a_top / Fig_1a_bottom / 
         ( ( wrap_elements() | Fig_1b_top_right ) + plot_layout(widths = c(1, 0.9345)) ) /
         ( Fig_1b_bottom_left | Fig_1b_bottom_right ) +
         plot_layout(heights = c(0.2, 1, 0.2, 1)) +
         plot_annotation(tag_levels = list(c("A", "", "B"))) &
         theme(plot.tag = element_text(family = "Futura", size = 15, face = "bold"))

ggsave(plot = Fig_1, filename = "Fig_1.pdf", 
       device = cairo_pdf, path = "Figures",
       width = 20, height = 20, units = "cm")