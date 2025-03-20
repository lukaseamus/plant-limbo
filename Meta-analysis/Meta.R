# Load and expand data using draws from the Gaussian distribution
require(tidyverse)
set.seed(100)
meta <- read.csv("~/Desktop/Projects/Seagrass/Data/Meta.csv") %>%
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

# Calculate metadata
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

# Export table of species
meta %>% 
  group_by(Phylum, Order, Family, Species, Group) %>%
  summarise(Days = if_else(max(Day) < 100, 
                           signif(max(Day), digits = 2), 
                           signif(max(Day), digits = 3)),
            Studies = n_distinct(Reference),
            Observations = n()) %>%
  arrange(Phylum, Order, Family, Species) %>%
  write.csv("~/Desktop/Species.csv", row.names = FALSE)

# Check that there are no highly improbable day ranges
meta %>% 
  group_by(Species) %>%
  summarise(Day_range = str_c(min(Day), max(Day), sep = "–"),
            Hour_range = str_c(min(Day*24), max(Day*24), sep = "–"),
            Minute_range = str_c(min(Day*24*60), max(Day*24*60), sep = "–"))
# all looks good

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

# Visually explore data
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

# Prior simulation
# k is informed by mean for seagrasses (0.22 d^-1) and kelp (0.12 d^-1)
# if seagrasses are assumed to be intermediate between seaweeds and terrestrial plants,
# k for terrestrial plants can be assumed to be about 0.3 to 0.4 and it is sensible to
# pick the estimate for seagrasses as the central tendency for the prior
# mu is informed by half the mean experimental duration
meta %>%
  group_by(Reference, Series) %>%
  summarise(Days = max(Day)) %>%
  ungroup() %$%
  mean(Days)/2

meta_prior <- 
  tibble(n = 1:5e3,
         k = rgamma(n = 5e3, shape = 0.22^2 / 0.3^2, rate = 0.22 / 0.3^2), 
         mu = rgamma(n = 5e3, shape = 8^2 / 100^2, rate = 8 / 100^2)) %>%
  filter(k * mu > 4) %>% # simulates constraining k and mu with a joint prior
  expand_grid(Day = meta %$% seq(min(Day), max(Day), length.out = 1e3)) %>%
  mutate(P_mu = 1 / ( 1 + exp( k * ( Day - mu ) ) ))

meta_prior %>%
  ggplot(aes(x = Day, y = P_mu, group = n)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal()
# filtering method doesn't work well


# prior simulation will also need to be done with the joint prior in Stan
meta_prior_stan <- "
parameters{
  real<lower=0> k;
  real<lower=0> mu;
}

model{
  k ~ gamma( 0.22^2 / 0.3^2 , 0.22 / 0.3^2 );
  mu ~ gamma( 8^2 / 100^2 , 8 / 100^2 );
  target += gamma_lpdf( k * mu | 4^2 / 0.1^2 , 4 / 0.1^2 );
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
    geom_density(aes(k), fill = "black", alpha = 0.5) +
    geom_density(aes(rgamma( 8e4, 0.22^2 / 0.3^2 , 0.22 / 0.3^2 )),
                 fill = "red", alpha = 0.5) +
    scale_x_continuous(limits = c(0, 1), oob = scales::oob_keep) +
    theme_minimal()

meta_prior_draws %>%
  ggplot() +
    geom_density(aes(mu), fill = "black", alpha = 0.5) +
    # geom_density(aes(rgamma( 8e4, 8^2 / 100^2 , 8 / 100^2 )),
    #              fill = "red", alpha = 0.5) +
    scale_x_continuous(limits = c(0, 100), oob = scales::oob_keep) +
    theme_minimal()

meta_prior_draws %>%
  slice_sample(n = 1e3) %>% # subsample
  expand_grid(Day = meta %$% seq(min(Day), max(Day), length.out = 1e3)) %>%
  mutate(P_mu = 1 / ( 1 + exp( k * ( Day - mu ) ) )) %>%
  ggplot(aes(x = Day, y = P_mu, group = .draw)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off", 
                    # xlim = c(0, 10)
                    ) +
    theme_minimal()


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
  k ~ gamma( 0.22^2 / 0.3^2 , 0.22 / 0.3^2 );
  mu ~ gamma( 8^2 / 100^2 , 8 / 100^2 );
  
  // Constraint on k and mu (joint prior)
  target += gamma_lpdf( k .* mu | 4^2 / 0.1^2 , 4 / 0.1^2 );
  
  // Likelihood uncertainty prior
  P_sigma ~ exponential( 1 );

  // Model
  vector[n] P_mu;
  for ( i in 1:n ) {
    P_mu[i] = 1 / ( 1 + exp( k[group[i]] * ( Day[i] - mu[group[i]] ) ) );
  }

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

meta_summary <- meta_samples$summary()
meta_summary %>%
  filter(rhat > 1.001) # no rhat above 1.001

meta_draws <- meta_samples$draws(format = "df")

require(bayesplot)
meta_draws %>% mcmc_rank_overlay() # chains look good

meta_draws %>% mcmc_pairs(pars = c("k[1]", "mu[1]")) # high correlation, but this is expected
meta_draws %>% mcmc_pairs(pars = c("k[2]", "mu[2]")) # due to joint prior
meta_draws %>% mcmc_pairs(pars = c("k[3]", "mu[3]"))
meta_draws %>% mcmc_pairs(pars = c("k[4]", "mu[4]"))
meta_draws %>% mcmc_pairs(pars = c("k[5]", "mu[5]"))
meta_draws %>% mcmc_pairs(pars = c("k[6]", "mu[6]"))
meta_draws %>% mcmc_pairs(pars = c("k[7]", "mu[7]"))
meta_draws %>% mcmc_pairs(pars = c("k[8]", "mu[8]"))
meta_draws %>% mcmc_pairs(pars = c("k[9]", "mu[9]"))
meta_draws %>% mcmc_pairs(pars = c("k[10]", "mu[10]"))
meta_draws %>% mcmc_pairs(pars = c("k[11]", "mu[11]"))
meta_draws %>% mcmc_pairs(pars = c("k[12]", "mu[12]"))
meta_draws %>% mcmc_pairs(pars = c("k[13]", "mu[13]"))
# this is chosen as the optimal model


# Prior-posterior comparison
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
                        quantile_lines = TRUE, quantiles = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95),
                        alpha = 0.5, scale = 2, rel_min_height = 0.001,
                        bandwidth = c(0.01, 0.1, 0.01, 0.1, 0.01, 0.6, 0.01, 10),
                        from = rep(0, 8), to = c(1.5, 10, 1.5, 10, 1, 60, 1, 1000)) +
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
      Group == "Terrestrial plants" & .variable == "k" ~ scale_x_continuous(limits = c(0, 1.5),
                                                                            breaks = seq(0, 1.5, by = 0.5),
                                                                            labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1))),
      Group == "Terrestrial plants" & .variable == "mu" ~ scale_x_continuous(limits = c(0, 10),
                                                                             breaks = seq(0, 10, by = 5)),
      Group == "Freshwater plants" & .variable == "k" ~ scale_x_continuous(limits = c(0, 1.5),
                                                                           breaks = seq(0, 1.5, by = 0.5),
                                                                           labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1))),
      Group == "Freshwater plants" & .variable == "mu" ~ scale_x_continuous(limits = c(0, 10),
                                                                            breaks = seq(0, 10, by = 5)),
      Group == "Seagrasses" & .variable == "k" ~ scale_x_continuous(limits = c(0, 1),
                                                                    breaks = seq(0, 1, by = 0.5),
                                                                    labels = scales::label_number(accuracy = c(1, 0.1, 1))),
      Group == "Seagrasses" & .variable == "mu" ~ scale_x_continuous(limits = c(0, 60),
                                                                     breaks = seq(0, 60, by = 30)),
      Group == "Seaweeds" & .variable == "k" ~ scale_x_continuous(limits = c(0, 1),
                                                                  breaks = seq(0, 1, by = 0.5),
                                                                  labels = scales::label_number(accuracy = c(1, 0.1, 1))),
      Group == "Seaweeds" & .variable == "mu" ~ scale_x_continuous(limits = c(0, 1000),
                                                                   breaks = seq(0, 1000, by = 500))
    )) +
    coord_cartesian(expand = FALSE) +
    mytheme +
    theme(axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(face = "italic"),
          plot.margin = margin(0, 0.5, 0.5, 0, unit = "cm"))

ggsave(plot = Fig_2_top, filename = "meta_prior_posterior.pdf", device = cairo_pdf, path = "~/Desktop",
       width = 20, height = 6, units = "cm")


# Parameters and differences
# Export parameters for Table 2
meta_prior_posterior %>%
  filter(Response != "Prior") %>%
  mutate(Response = fct_relevel(Response, "Photosynthesis")) %>%
  group_by(Group, Response, Treatment, .variable) %>%
  summarise(mean = mean(.value),
            sd = sd(.value),
            n = length(.value)) %>%
  mutate(rounded = paste( # 2 significant figures except for very large numbers
    if_else(mean < 100, signif(mean, digits = 2), signif(mean, digits = 3)), "±",
    if_else(sd < 100, signif(sd, digits = 2), signif(sd, digits = 3))
                        )
         ) %>%
  pivot_wider(names_from = c(Response, Treatment), values_from = c(mean, sd, rounded)) %>%
  write.csv("~/Desktop/meta_para.csv", row.names = FALSE)

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
               names_to = "Contrast", values_to = "Difference", names_prefix = "delta_")

# Export differences and probabilities for Table 2
meta_diff %>%
  group_by(Contrast) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P_less = mean(Difference < 0),
            P_more = mean(Difference > 0),
            n = length(Difference)) %>%
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
  write.csv("~/Desktop/meta_diff.csv", row.names = FALSE)

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
               names_to = "Contrast", values_to = "Difference", names_prefix = "delta_")

# Export differences and probabilities for text
meta_diff_group %>%
  group_by(Contrast) %>%
  summarise(mean = mean(Difference),
            sd = sd(Difference),
            P_less = mean(Difference < 0),
            P_more = mean(Difference > 0),
            n = length(Difference)) %>%
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
  write.csv("~/Desktop/meta_diff_group.csv", row.names = FALSE)

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
  mutate(Day = if_else(Response == "Prior",
                       list( c(seq(0, Days * 0.035, length.out = 50), # Prior will be plotted on multiple scales
                               seq(Days * 0.035 + 1, Days, length.out = 50)) ), # so needs higher resolution at lower values
                       list( seq(0, Days, length.out = 100) ) 
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

require(patchwork)
Fig_2 <- ( (Fig_2_top + theme(legend.position = c(0.14, -1.35))) / Fig_2_middle / Fig_2_bottom ) +
  plot_layout(heights = c(0.4, 1, 1))
# ignore legend.position warning: legend.position.inside will not produce the same result

ggsave(plot = Fig_2, filename = "Fig_2.pdf", device = cairo_pdf, path = "~/Desktop",
       width = 20, height = 15, units = "cm")



