# 1. Data ####
# 1.1 Load data ####
# Load and expand data using draws from the Gaussian distribution
require(tidyverse)
require(here)
set.seed(100)
meta <- here("Meta-analysis", "Meta.csv") %>% read.csv() %>%
  rowwise() %>%
  mutate(Observation = if_else(
    !is.na(SEM),
    list( rnorm( N , Mean , SEM * sqrt(N) ) ),
    list( Mean )
    )) %>%
  unnest(Observation)
# warnings can be safely ignored (I tested this with a longer method using filter)

# Mean and Observation should be identical for rows with only one observation
meta %>%
  filter(is.na(SEM)) %>%
  mutate(diff = Mean - Observation) %>%
  pull(diff) %>%
  range()
# diff = 0 so Mean and Observation are identical

# 1.2 Calculate metadata ####
meta %>% nrow()
# 10566 observations

meta %>% group_by(Reference, Series) %>% n_groups()
# 535 measurement series

meta %>% group_by(Reference) %>% n_groups()
# 127 studies

meta %>% group_by(Species) %>% n_groups()
# 92 species

meta %>% group_by(Family) %>% n_groups()
# 37 families

meta %>% group_by(Order) %>% n_groups()
# 23 orders

meta %>% group_by(Phylum) %>% n_groups()
# 3 phyla

meta %>% 
  group_by(Group) %>% 
  summarise(Species = n_distinct(Species))
# 15 seaweeds
# 5 seagrasses
# 6 freshwater plants
# 66 terrestrial plants

# 1.3 Export table of species ####
meta %>% 
  group_by(Phylum, Order, Family, Species, Group) %>%
  summarise(Days = if_else(max(Day) < 100, 
                           signif(max(Day), digits = 2), 
                           signif(max(Day), digits = 3)),
            Studies = n_distinct(Reference),
            Observations = n()) %>%
  arrange(Phylum, Order, Family, Species) %>%
  write.csv(here("Tables", "Species.csv"), row.names = FALSE)

# Check that there are no highly improbable day ranges
meta %>% 
  group_by(Species) %>%
  summarise(Day_range = str_c(min(Day), max(Day), sep = "–"),
            Hour_range = str_c(min(Day*24), max(Day*24), sep = "–"),
            Minute_range = str_c(min(Day*24*60), max(Day*24*60), sep = "–"))
# all looks good

# 1.4 Transform data ####
# Re-express data as a proportion of the initial value in the measurement series
# or if there are multiple initial values as a proportion of their mean
require(magrittr)
meta %<>%
  group_by(Reference, Series) %>%
  mutate(Proportion = if_else(
          is.na(SEM),
          Observation / mean( Observation[ Day == min(Day) ] ),
          Observation / mean( Mean[ Day == min(Day) ] )
          )) %>%
  ungroup() %>%
  mutate(Group = fct_relevel(Group, "Terrestrial", "Freshwater", "Seagrass"))

# 1.5 Visually explore data ####
meta %>%
  ggplot(aes(Day, Proportion)) +
    geom_point(alpha = 0.2, shape = 16) +
    facet_grid(Response + Light ~ Group + Water, scales = "free_x") +
    coord_cartesian(ylim = c(0, 2)) +
    theme_minimal()
# too many predictors

meta %>%
  ggplot(aes(Day, Proportion, colour = Light)) +
  geom_point(alpha = 0.2, shape = 16) +
  facet_grid(Response ~ Group, scales = "free_x") +
  coord_cartesian(ylim = c(0, 2)) +
  theme_minimal()
# much clearer

# Create composite categorical predictor
meta %<>%
  mutate(group = fct((str_c(Group, Response, Light, sep = "_"))))

# 2. Prior simulation ####
# There is no information on k and mu for plants at large. Consideirng that
# the product of k and mu has to be aorund 5 for the logitsic curve to start
# near its maximum at t0, ranges for k and mu can be reciprocally estimated.
# For instance, I'd assume a reasonable range for mu is half a day to 3 years.
# This would mean the maximum for k has to be around 5 / 0.5 = 10. For simplicity
# the lower bound can be set to zero in both cases. To allow the posterior to
# go beyond these estimated bounds, a truncated normal prior is better than
# a uniform one.

# prior simulation will also need to be done with the joint prior in Stan
meta_prior_stan <- "
parameters{
  real<lower=0> k;
  real<lower=0> mu;
}

model{
  // k ~ uniform( 0 , 10 );
  // mu ~ uniform( 0 , 1095 );
  k ~ normal( 5 , 4 ) T[0,];
  mu ~ normal( 365 , 300 ) T[0,];
  target += normal_lpdf( log(k) + log(mu) | log(5) , 0.05 );
}
"
require(cmdstanr)
meta_prior_mod <- cmdstan_model(stan_file = write_stan_file(code = meta_prior_stan))
meta_prior_samples <- meta_prior_mod$sample(data = list(), # no data to condition on
                                            seed = 100,
                                            chains = 8,
                                            parallel_chains = parallel::detectCores(),
                                            iter_warmup = 1e4,
                                            iter_sampling = 1e4)

meta_prior_draws <- meta_prior_samples$draws(format = "df")

meta_prior_draws %>%
  ggplot() +
    geom_density(aes(k), fill = "black", alpha = 0.5, bw = 0.01) +
    scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep) +
    theme_minimal()

meta_prior_draws %>%
  ggplot() +
    geom_density(aes(mu), fill = "black", alpha = 0.5, bw = 1) +
    scale_x_continuous(limits = c(0, 100), oob = scales::oob_keep) +
    theme_minimal()

meta_prior_draws %>%
  slice_sample(n = 1e3) %>% # subsample
  expand_grid(Day = meta %$% seq(min(Day), max(Day), length.out = 1e4)) %>%
  mutate(P_mu = 1 / ( 1 + exp( k * ( Day - mu ) ) )) %>%
  ggplot(aes(x = Day, y = P_mu, group = .draw)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off", 
                    # xlim = c(0, 5)
                    ) +
    theme_minimal()

# 3. Stan model ####
meta_stan <- "
data{
  int n;
  int n_group;
  vector[n] Proportion;
  vector[n] Day;
  array[n] int group;
}

parameters{
  // Group parameters
  vector<lower=0>[n_group] k;
  vector<lower=0>[n_group] mu;

  // Likelihood uncertainty parameter
  real<lower=0> P_sigma;
}

model{
  // Group priors
  k ~ normal( 5 , 4 ) T[0,];
  mu ~ normal( 365 , 300 ) T[0,];
  
  // Constraint on k and mu (joint prior)
  target += normal_lpdf( log(k) + log(mu) | log(5) , 0.01 );
  
  // Likelihood uncertainty prior
  P_sigma ~ exponential( 1 );

  // Model
  vector[n] P_mu = 1 / ( 1 + exp( k[group] .* ( Day - mu[group] ) ) );

  // Likelihood
  Proportion ~ normal( P_mu , P_sigma );
}
"
meta_mod <- cmdstan_model(stan_file = write_stan_file(code = meta_stan))

require(tidybayes)
meta_samples <- meta_mod$sample(data = meta %>%
                                  select(Proportion, Day, group) %>%
                                  compose_data(),
                                seed = 100,
                                chains = 8,
                                parallel_chains = parallel::detectCores(),
                                iter_warmup = 1e4,
                                iter_sampling = 1e4)

# 4. Model checks ####
meta_summary <- meta_samples$summary()
meta_summary %>%
  filter(rhat > 1.001) # no rhat above 1.001

meta_draws <- meta_samples$draws(format = "df")

require(bayesplot)
meta_draws %>% mcmc_rank_overlay() # chains look good

meta_draws %>% mcmc_pairs(pars = c("k[1]", "mu[1]")) # strong correlation, but this is expected
meta_draws %>% mcmc_pairs(pars = c("k[13]", "mu[13]")) # due to joint prior

# 5. Prior-posterior comparison ####
meta_prior_posterior <- meta_samples %>%
  recover_types(meta %>% select(group)) %>%
  gather_draws(k[group], mu[group]) %>%
  ungroup() %>%
  separate(group, into = c("Group", "Response", "Treatment"), sep = "_") %>%
  bind_rows(
    meta_prior_samples %>%
      gather_draws(k, mu) %>%
      ungroup() %>%
      slice(rep( 1:n(), 4 )) %>%
      mutate(Group = c("Terrestrial", "Freshwater", "Seagrass", "Seaweed") %>%
               rep(each = 8e4 * 2),
             Response = "Prior", Treatment = "Prior")
  ) %>%
  mutate(Group = case_when(
                  Group == "Terrestrial" ~ "Terrestrial plants",
                  Group == "Freshwater" ~ "Freshwater plants",
                  Group == "Seagrass" ~ "Seagrasses",
                  Group == "Seaweed" ~ "Seaweeds"
                 ) %>% fct_relevel("Terrestrial plants", "Freshwater plants", "Seagrasses"),
         Response = fct_relevel(Response, "Chlorophyll", "Photosynthesis"),
         Treatment = case_when(
                      Treatment == "Yes" ~ "Light", 
                      Treatment == "No" ~ "Dark",
                      Treatment == "Prior" ~ "Prior"
                     ) %>% fct_relevel("Prior", "Light")
         )

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

require(ggridges)
require(ggh4x)
Fig_2_top <- meta_prior_posterior %>%
  ggplot() +
    geom_density_ridges(aes(.value, Response, colour = Treatment, fill = Treatment),
                        # quantile_lines = TRUE, quantiles = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95),
                        alpha = 0.5, scale = 2, rel_min_height = 0.002,
                        bandwidth = c(2*0.02, 10*0.02, 2*0.02, 10*0.02, 
                                      1.5*0.02, 40*0.02, 1.5*0.02, 1000*0.02),
                        from = rep(0, 8), to = c(2, 10, 2, 10, 1.5, 40, 1.5, 1000)) +
    scale_colour_manual(values = c("#b5b8ba", "#f5a54a", "#2e4a5b")) +
    scale_fill_manual(values = c("#b5b8ba", "#f5a54a", "#2e4a5b")) +
    scale_y_discrete(labels = c("Chlorophyll" = "Chl", "Photosynthesis" = "P", "Prior" = "")) +
    facet_nested(~ Group + .variable, scales = "free_x",
                 labeller = labeller(.variable = as_labeller(
                                     c("k" = "italic('k')*' (d'^-1*')'",
                                       "mu" = "italic('µ')*' (d)'"),
                                     label_parsed))
                 ) +
    facetted_pos_scales(x = list(
      Group == "Terrestrial plants" & .variable == "k" ~ scale_x_continuous(limits = c(0, 2),
                                                                            breaks = seq(0, 2, by = 1),
                                                                            oob = scales::oob_keep),
      Group == "Terrestrial plants" & .variable == "mu" ~ scale_x_continuous(limits = c(0, 10),
                                                                             breaks = seq(0, 10, by = 5)),
      Group == "Freshwater plants" & .variable == "k" ~ scale_x_continuous(limits = c(0, 2),
                                                                           breaks = seq(0, 2, by = 1),
                                                                           oob = scales::oob_keep),
      Group == "Freshwater plants" & .variable == "mu" ~ scale_x_continuous(limits = c(0, 10),
                                                                            breaks = seq(0, 10, by = 5),
                                                                            oob = scales::oob_keep),
      Group == "Seagrasses" & .variable == "k" ~ scale_x_continuous(limits = c(0, 1.5),
                                                                    breaks = seq(0, 1.5, by = 0.5),
                                                                    oob = scales::oob_keep,
                                                                    labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1))),
      Group == "Seagrasses" & .variable == "mu" ~ scale_x_continuous(limits = c(0, 40),
                                                                     breaks = seq(0, 40, by = 20),
                                                                     oob = scales::oob_keep),
      Group == "Seaweeds" & .variable == "k" ~ scale_x_continuous(limits = c(0, 1.5),
                                                                  breaks = seq(0, 1.5, by = 0.5),
                                                                  oob = scales::oob_keep,
                                                                  labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1))),
      Group == "Seaweeds" & .variable == "mu" ~ scale_x_continuous(limits = c(0, 1000),
                                                                   oob = scales::oob_keep,
                                                                   breaks = seq(0, 1000, by = 500))
    )) +
    coord_cartesian(expand = FALSE) +
    mytheme +
    theme(axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(face = "italic"),
          plot.margin = margin(0, 0.5, 0.5, 0, unit = "cm"))
Fig_2_top

# 6. Parameters and differences ####
# Export parameters for Table 2
meta_prior_posterior %>%
  filter(Response != "Prior") %>%
  mutate(Response = fct_relevel(Response, "Photosynthesis")) %>%
  group_by(Group, Response, Treatment, .variable) %>%
  summarise(mean = mean(.value),
            sd = sd(.value),
            n = n()) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(mean, digits = 2), signif(mean, digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        )
         ) %>%
  pivot_wider(names_from = c(Response, Treatment), values_from = c(mean, sd, rounded)) %>%
  write.csv(here("Tables", "meta_para.csv"), row.names = FALSE)

# Calculate differences for Table 2
meta_diff <- meta_samples %>%
  recover_types(meta %>% select(group)) %>%
  gather_draws(k[group], mu[group]) %>%
  ungroup() %>%
  pivot_wider(names_from = c(.variable, group), values_from = .value) %>%
  mutate(
         delta_k_Terrestrial_Photosynthesis = k_Terrestrial_Photosynthesis_Yes - k_Terrestrial_Photosynthesis_No,
         delta_mu_Terrestrial_Photosynthesis = mu_Terrestrial_Photosynthesis_Yes - mu_Terrestrial_Photosynthesis_No,
         delta_k_Seaweed_Photosynthesis = k_Seaweed_Photosynthesis_Yes - k_Seaweed_Photosynthesis_No,
         delta_mu_Seaweed_Photosynthesis = mu_Seaweed_Photosynthesis_Yes - mu_Seaweed_Photosynthesis_No,
         delta_k_Terrestrial_Chlorophyll = k_Terrestrial_Chlorophyll_Yes - k_Terrestrial_Chlorophyll_No,
         delta_mu_Terrestrial_Chlorophyll = mu_Terrestrial_Chlorophyll_Yes - mu_Terrestrial_Chlorophyll_No,
         delta_k_Freshwater_Chlorophyll = k_Freshwater_Chlorophyll_Yes - k_Freshwater_Chlorophyll_No,
         delta_mu_Freshwater_Chlorophyll = mu_Freshwater_Chlorophyll_Yes - mu_Freshwater_Chlorophyll_No,
         delta_k_Seagrass_Chlorophyll = k_Seagrass_Chlorophyll_Yes - k_Seagrass_Chlorophyll_No,
         delta_mu_Seagrass_Chlorophyll = mu_Seagrass_Chlorophyll_Yes - mu_Seagrass_Chlorophyll_No,
         delta_k_Terrestrial_Light = k_Terrestrial_Photosynthesis_Yes - k_Terrestrial_Chlorophyll_Yes,
         delta_mu_Terrestrial_Light = mu_Terrestrial_Photosynthesis_Yes - mu_Terrestrial_Chlorophyll_Yes,
         delta_k_Terrestrial_Dark = k_Terrestrial_Photosynthesis_No - k_Terrestrial_Chlorophyll_No,
         delta_mu_Terrestrial_Dark = mu_Terrestrial_Photosynthesis_No - mu_Terrestrial_Chlorophyll_No,
         delta_k_Freshwater_Dark = k_Freshwater_Photosynthesis_No - k_Freshwater_Chlorophyll_No,
         delta_mu_Freshwater_Dark = mu_Freshwater_Photosynthesis_No - mu_Freshwater_Chlorophyll_No,
         delta_k_Seagrass_Light = k_Seagrass_Photosynthesis_Yes - k_Seagrass_Chlorophyll_Yes,
         delta_mu_Seagrass_Light = mu_Seagrass_Photosynthesis_Yes - mu_Seagrass_Chlorophyll_Yes,
         delta_k_Seaweed_Light = k_Seaweed_Photosynthesis_Yes - k_Seaweed_Chlorophyll_Yes,
         delta_mu_Seaweed_Light = mu_Seaweed_Photosynthesis_Yes - mu_Seaweed_Chlorophyll_Yes
         ) %>%
  select(.chain, .iteration, .draw, starts_with("delta")) %>%
  pivot_longer(cols = -c(.chain, .iteration, .draw),
               names_to = "Contrast", values_to = "Difference", names_prefix = "delta_") %T>%
  print()

# Export differences and probabilities for Table 2
meta_diff %>%
  group_by(Contrast) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P_less = mean(Difference < 0),
            P_more = mean(Difference > 0),
            n = n()) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(abs(mean), digits = 2), signif(abs(mean), digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        ),
         P = signif(pmax(P_less, P_more), digits = 2)
  ) %>%
  select(-c(P_less, P_more)) %>%
  separate(Contrast, into = c("Parameter", "Group", "Contrast"), sep = "_") %>%
  mutate(Parameter = fct_relevel(Parameter, "k"),
         Group = fct_relevel(Group, "Terrestrial", "Freshwater", "Seagrass"),
         Contrast = fct_relevel(Contrast, "Photosynthesis", "Chlorophyll", "Light")) %>%
  arrange(Group, Contrast, Parameter) %>%
  write.csv(here("Tables", "meta_diff.csv"), row.names = FALSE)

# Calculate differences between plant groups for text
meta_diff_group <- meta_samples %>%
  recover_types(meta %>% select(group)) %>%
  gather_draws(k[group], mu[group]) %>%
  ungroup() %>%
  pivot_wider(names_from = c(.variable, group), values_from = .value) %>%
  mutate(
        delta_k_Photosynthesis_Light_Terrestrial_Seagrass = k_Terrestrial_Photosynthesis_Yes - k_Seagrass_Photosynthesis_Yes,
        delta_mu_Photosynthesis_Light_Terrestrial_Seagrass = mu_Terrestrial_Photosynthesis_Yes - mu_Seagrass_Photosynthesis_Yes,
        delta_k_Photosynthesis_Light_Terrestrial_Seaweed = k_Terrestrial_Photosynthesis_Yes - k_Seaweed_Photosynthesis_Yes,
        delta_mu_Photosynthesis_Light_Terrestrial_Seaweed = mu_Terrestrial_Photosynthesis_Yes - mu_Seaweed_Photosynthesis_Yes,
        delta_k_Photosynthesis_Dark_Terrestrial_Freshwater = k_Terrestrial_Photosynthesis_No - k_Freshwater_Photosynthesis_No,
        delta_mu_Photosynthesis_Dark_Terrestrial_Freshwater = mu_Terrestrial_Photosynthesis_No - mu_Freshwater_Photosynthesis_No,
        delta_k_Photosynthesis_Dark_Terrestrial_Seaweed = k_Terrestrial_Photosynthesis_No - k_Seaweed_Photosynthesis_No,
        delta_mu_Photosynthesis_Dark_Terrestrial_Seaweed = mu_Terrestrial_Photosynthesis_No - mu_Seaweed_Photosynthesis_No,
        delta_k_Chlorophyll_Light_Terrestrial_Freshwater = k_Terrestrial_Chlorophyll_Yes - k_Freshwater_Chlorophyll_Yes,
        delta_mu_Chlorophyll_Light_Terrestrial_Freshwater = mu_Terrestrial_Chlorophyll_Yes - mu_Freshwater_Chlorophyll_Yes,
        delta_k_Chlorophyll_Light_Terrestrial_Seagrass = k_Terrestrial_Chlorophyll_Yes - k_Seagrass_Chlorophyll_Yes,
        delta_mu_Chlorophyll_Light_Terrestrial_Seagrass = mu_Terrestrial_Chlorophyll_Yes - mu_Seagrass_Chlorophyll_Yes,
        delta_k_Chlorophyll_Light_Terrestrial_Seaweed = k_Terrestrial_Chlorophyll_Yes - k_Seaweed_Chlorophyll_Yes,
        delta_mu_Chlorophyll_Light_Terrestrial_Seaweed = mu_Terrestrial_Chlorophyll_Yes - mu_Seaweed_Chlorophyll_Yes,
        delta_k_Chlorophyll_Dark_Terrestrial_Freshwater = k_Terrestrial_Chlorophyll_No - k_Freshwater_Chlorophyll_No,
        delta_mu_Chlorophyll_Dark_Terrestrial_Freshwater = mu_Terrestrial_Chlorophyll_No - mu_Freshwater_Chlorophyll_No,
        delta_k_Chlorophyll_Dark_Terrestrial_Seagrass = k_Terrestrial_Chlorophyll_No - k_Seagrass_Chlorophyll_No,
        delta_mu_Chlorophyll_Dark_Terrestrial_Seagrass = mu_Terrestrial_Chlorophyll_No - mu_Seagrass_Chlorophyll_No,
        delta_k_Photosynthesis_Light_Seagrass_Seaweed = k_Seagrass_Photosynthesis_Yes - k_Seaweed_Photosynthesis_Yes,
        delta_mu_Photosynthesis_Light_Seagrass_Seaweed = mu_Seagrass_Photosynthesis_Yes - mu_Seaweed_Photosynthesis_Yes,
        delta_k_Photosynthesis_Dark_Freshwater_Seaweed = k_Freshwater_Photosynthesis_No - k_Seaweed_Photosynthesis_No,
        delta_mu_Photosynthesis_Dark_Freshwater_Seaweed = mu_Freshwater_Photosynthesis_No - mu_Seaweed_Photosynthesis_No,
        delta_k_Chlorophyll_Light_Freshwater_Seagrass = k_Freshwater_Chlorophyll_Yes - k_Seagrass_Chlorophyll_Yes,
        delta_mu_Chlorophyll_Light_Freshwater_Seagrass = mu_Freshwater_Chlorophyll_Yes - mu_Seagrass_Chlorophyll_Yes,
        delta_k_Chlorophyll_Light_Freshwater_Seaweed = k_Freshwater_Chlorophyll_Yes - k_Seaweed_Chlorophyll_Yes,
        delta_mu_Chlorophyll_Light_Freshwater_Seaweed = mu_Freshwater_Chlorophyll_Yes - mu_Seaweed_Chlorophyll_Yes,
        delta_k_Chlorophyll_Light_Seagrass_Seaweed = k_Seagrass_Chlorophyll_Yes - k_Seaweed_Chlorophyll_Yes,
        delta_mu_Chlorophyll_Light_Seagrass_Seaweed = mu_Seagrass_Chlorophyll_Yes - mu_Seaweed_Chlorophyll_Yes,
        delta_k_Chlorophyll_Dark_Freshwater_Seagrass = k_Freshwater_Chlorophyll_No - k_Seagrass_Chlorophyll_No,
        delta_mu_Chlorophyll_Dark_Freshwater_Seagrass = mu_Freshwater_Chlorophyll_No - mu_Seagrass_Chlorophyll_No
  ) %>%
  select(.chain, .iteration, .draw, starts_with("delta")) %>%
  pivot_longer(cols = -c(.chain, .iteration, .draw),
               names_to = "Contrast", values_to = "Difference", names_prefix = "delta_") %T>%
  print()

# Export differences and probabilities for text
meta_diff_group %>%
  group_by(Contrast) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P_less = mean(Difference < 0),
            P_more = mean(Difference > 0),
            n = n()) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(abs(mean), digits = 2), signif(abs(mean), digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        ),
         P = signif(pmax(P_less, P_more), digits = 2)
        ) %>%
  select(-c(P_less, P_more)) %>%
  separate(Contrast, into = c("Parameter", "Response", "Treatment", "First", "Second"), sep = "_") %>%
  mutate(Parameter = fct_relevel(Parameter, "k"),
         Response = fct_relevel(Response, "Photosynthesis"),
         Treatment = fct_relevel(Treatment, "Light"),
         First = fct_relevel(First, "Terrestrial", "Freshwater", "Seagrass"),
         Second = fct_relevel(Second, "Freshwater", "Seagrass", "Seaweed")) %>%
  arrange(Response, Treatment, First, Second, Parameter) %>%
  write.csv(here("Tables", "meta_diff_group.csv"), row.names = FALSE)

# 7. Prediction ####
# calculate P_mu from parameters
meta_renamed <- meta %>%
  rename("Treatment" = Light) %>%
  mutate(Group = case_when(
                  Group == "Terrestrial" ~ "Terrestrial plants",
                  Group == "Freshwater" ~ "Freshwater plants",
                  Group == "Seagrass" ~ "Seagrasses",
                  Group == "Seaweed" ~ "Seaweeds"
                 ) %>% fct_relevel("Terrestrial plants", "Freshwater plants", "Seagrasses"),
         Treatment = if_else(Treatment == "Yes", "Light", "Dark") %>%
                      fct_relevel("Light"))

meta_mu <- meta_prior_posterior %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  left_join(meta_renamed %>%
              group_by(Group, Response, Treatment) %>%
              summarise(Days = max(Day)),
            by = c("Group", "Response", "Treatment")) %>%
  mutate(Days = if_else(is.na(Days), 320, Days)) %>% # 320 days is the maximum required range
  rowwise() %>% # this is necessary because Days needs to be called separately for each row
  mutate(Day = case_when(
    Response == "Prior" ~ list( c(seq(0, Days * 0.01, length.out = 50), # Prior will be plotted on multiple scales
                                  seq(Days * 0.01 + 1, Days, length.out = 50)) ), # so needs higher resolution at lower values
    Group %in% c("Terrestrial plants", "Freshwater plants") &
      Response == "Chlorophyll" & Treatment == "Dark" ~ list( seq(0, Days * 0.2, length.out = 100) ),
    TRUE ~ list( seq(0, Days, length.out = 100) ) 
                       )) %>%
  unnest(Day) %>% # expand the list column Day
  mutate(P_mu = 1 / ( 1 + exp( k * ( Day - mu ) ) ) )

meta_mu_summary <- meta_mu %>%
  mutate(Group = fct_relevel(Group, "Terrestrial plants", "Freshwater plants", "Seagrasses"),
         Treatment = fct_relevel(Treatment, "Light")) %>%
  group_by(Day, Group, Response, Treatment) %>%
  median_qi(P_mu, .width = c(.5, .8, .9)) # mean is sometimes outside 50% interval, so use median

Fig_2_middle <- ggplot() +
                  geom_point(data = meta_renamed %>%
                               filter(Response == "Photosynthesis"),
                             aes(Day, Proportion, colour = Treatment),
                             alpha = 0.2, shape = 16) +
                  geom_ribbon(data = meta_mu_summary %>% filter(Response == "Prior", .width == 0.9),
                              aes(Day, ymin = .lower, ymax = .upper), fill = NA, colour = "#b5b8ba") +
                  geom_line(data = meta_mu_summary %>% filter(Response == "Photosynthesis"),
                            aes(Day, P_mu, colour = Treatment)) +
                  geom_ribbon(data = meta_mu_summary %>% filter(Response == "Photosynthesis"),
                              aes(Day, ymin = .lower, ymax = .upper,
                                  fill = Treatment, alpha = factor(.width)),
                              colour = NA) +
                  labs(y = expression(italic("P")*" (normalised)")) +
                  scale_colour_manual(values = c("#f5a54a", "#2e4a5b"),
                                      guide = "none") +
                  scale_fill_manual(values = c("#f5a54a", "#2e4a5b"),
                                    guide = "none") +
                  scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                  scale_y_continuous(labels = scales::label_number(accuracy = c(1, 0.01, 0.1, 0.01, 1))) +
                  facet_grid(~ Group, scales = "free_x") +
                  facetted_pos_scales(x = list(
                    Group == "Terrestrial plants" ~ scale_x_continuous(limits = c(0, 30), 
                                                                       breaks = seq(0, 30, by = 10),
                                                                       labels = NULL),
                    Group == "Freshwater plants" ~ scale_x_continuous(limits = c(0, 18), 
                                                                      breaks = seq(0, 18, by = 6),
                                                                      labels = NULL),
                    Group == "Seagrasses" ~ scale_x_continuous(limits = c(0, 60), 
                                                               breaks = seq(0, 60, by = 20),
                                                               labels = NULL),
                    Group == "Seaweeds" ~ scale_x_continuous(limits = c(0, 180), 
                                                             breaks = seq(0, 180, by = 60))
                  )) +
                  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
                  mytheme +
                  theme(axis.title.x = element_blank(),
                        strip.text = element_blank())
Fig_2_middle

Fig_2_bottom <- ggplot() +
                  geom_point(data = meta_renamed %>%
                               filter(Response == "Chlorophyll"),
                             aes(Day, Proportion, colour = Treatment),
                             alpha = 0.2, shape = 16) +
                  geom_ribbon(data = meta_mu_summary %>% filter(Response == "Prior", .width == 0.9),
                              aes(Day, ymin = .lower, ymax = .upper), fill = NA, colour = "#b5b8ba") +
                  geom_line(data = meta_mu_summary %>% filter(Response == "Chlorophyll"),
                            aes(Day, P_mu, colour = Treatment)) +
                  geom_ribbon(data = meta_mu_summary %>% filter(Response == "Chlorophyll"),
                              aes(Day, ymin = .lower, ymax = .upper, fill = Treatment, 
                                  alpha = factor(.width)),
                              colour = NA) +
                  labs(y = expression(italic("Chl")*" (normalised)"),
                       x = "Detrital age (d)") +
                  scale_colour_manual(values = c("#f5a54a", "#2e4a5b"),
                                      guide = "none") +
                  scale_fill_manual(values = c("#f5a54a", "#2e4a5b"),
                                    guide = "none") +
                  scale_alpha_manual(values = c(0.4, 0.3, 0.2), guide = "none") +
                  scale_y_continuous(labels = scales::label_number(accuracy = c(1, 0.01, 0.1, 0.01, 1))) +
                  facet_grid(~ Group, scales = "free_x") +
                  facetted_pos_scales(x = list(
                    Group == "Terrestrial plants" ~ scale_x_continuous(limits = c(0, 30), 
                                                                       breaks = seq(0, 30, by = 10)),
                    Group == "Freshwater plants" ~ scale_x_continuous(limits = c(0, 18), 
                                                                      breaks = seq(0, 18, by = 6)),
                    Group == "Seagrasses" ~ scale_x_continuous(limits = c(0, 60), 
                                                               breaks = seq(0, 60, by = 20)),
                    Group == "Seaweeds" ~ scale_x_continuous(limits = c(0, 320), 
                                                             breaks = seq(0, 320, by = 80))
                  )) +
                  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
                  mytheme +
                  theme(strip.text = element_blank())
Fig_2_bottom

require(patchwork)
Fig_2 <- ( (Fig_2_top + theme(legend.position = c(0.14, -1.35))) / Fig_2_middle / Fig_2_bottom ) +
  plot_layout(heights = c(0.4, 1, 1))
# ignore legend.position warning: legend.position.inside will not produce the same result

ggsave(plot = Fig_2, filename = "Fig_2.pdf", 
       device = cairo_pdf, path = "Figures",
       width = 20, height = 15, units = "cm")