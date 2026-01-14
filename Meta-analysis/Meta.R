# 1. Data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)
meta <- here("Meta-analysis", "Meta.csv") %>% 
  read_csv(
    col_types = list(
      "f", "c", "f", "f", "f", "f",
      "f", "f", "f", "f", "d", "d",
      "d", "d", "f", "f", "c", "c"
    )
  ) %T>%
  print()

meta %>% distinct(Reference, Series)
# There are 551 unique reference-series combinations,
# that is 551 unique experiments across studies.

# Check that series are not repeated in reference
meta %>%
  group_by(Reference) %>%
  mutate(
    block = ( Series != Series %>% 
                lag(default = first(Series)) ) %>%
      cumsum()
  ) %>%
  distinct(Reference, Series, block) %>%
  count(Reference, Series) %>%
  filter(n > 1)
# Every series in reference matches one block, 
# so e.g. series 1 only exists once in any
# given reference.

# Check that there are no highly improbable day ranges
meta %>% 
  group_by(Species) %>%
  summarise(Day_range = str_c(min(Day), max(Day), sep = "–"),
            Hour_range = str_c(min(Day*24), max(Day*24), sep = "–"),
            Minute_range = str_c(min(Day*24*60), max(Day*24*60), sep = "–")) %>%
  print(n = 100)
# All looks reasonable.

meta %$% any(is.na(N)) # Sample size contains no NAs

# 1.2 Simulate observations ####
# Expand data using draws from probability distributions
meta %>% distinct(Response, Method)
# Most of these response variables are restricted to positive
# numbers, so are best simulated using a gamma distribution.
meta %>% distinct(Method, Unit) %>% arrange(Method) %>% print(n = 70)
# Some are also bounded above, such as Fv/Fm, so are best simulated using 
# a beta distribution. However, % of initial, % of total and proportion 
# of initial are not true proportions because they may exceed 100% and 1
# respectively when photosynthesis increases after excision:
meta %>% filter(Unit %>% str_detect("%")) %$% max(Mean)
meta %>% filter(Unit %>% str_detect("Prop")) %$% max(Mean)
# so they are best described by a gamma distribution. Gas exchange should
# be simulated with a normal distribution because cases of net gas exchange 
# can be negative.
# Now the problem is that Mean = 0 is reported for some cases of Fv/Fm:
meta %>% filter(!is.na(SEM) & N > 1 & Method != "Gas exchange" & Mean <= 0)
# This precludes using the beta distribution. In addition, beta only allows 
# SEM < sqrt( Mean * (1 - Mean ) / N ) which is violated in more cases:
meta %>%
  filter(!is.na(SEM) & N > 1 & Unit == "Fv/Fm") %>% 
  mutate(SEM_max = sqrt( Mean * (1 - Mean) / N )) %>%
  filter(SEM >= SEM_max)
# Hence for simulations from beta, I would need to enforce Mean > 0 
# and SEM < sqrt( Mean * (1 - Mean ) / N ). The alternative would be a
# truncated normal bounded by 0 and 1 but this does not correctly
# represent the generative distribution that could have resulted
# in the observed Mean and SEM. If I enforce a maximal SEM < SEM_max,
# but close to SEM_max, the beta distribution only predicts values 
# near 0 and 1, so the truncated normal distribution is likely still
# more representative of Fv/Fm observations.

# Draw from probability distributions
require(extraDistr) # R doesn't have a built-in truncated normal
set.seed(100) # Ensure reproducibility of simulation
meta %<>%
  rowwise() %>%
  mutate(
    Observation = if( !is.na(SEM) & N > 1 ){ # Cases with summaries of observations
      SD <- SEM * sqrt(N) # Calculate standard deviation
      if( Method == "Gas exchange" ){ # Cases of gas exchange that may be negative
        list( rnorm( N , Mean , SD ) )
      } else if( Unit == "Fv/Fm" ){ # Fv/Fm is bounded by 0 and 1
        list( rtnorm( N , Mean , SD , 0 , 1 ) )
      } else { # All other cases are positive but otherwise unbounded
        list( rgamma( N , Mean^2 / SD^2 , Mean / SD^2 ) )
      }
    } else { # Cases with raw observations
      list( Mean )
    }
  ) %>%
  unnest(Observation) %T>%
  print()

meta %$% all(is.finite(Observation)) # No NAs or NaNs or Infs

# Mean and Observation should be identical 
# for rows with only one observation
meta %>%
  filter(is.na(SEM) & N == 1) %>%
  filter(Mean != Observation)
# Mean and Observation are identical

# 1.3 Calculate metadata ####
meta %>% nrow()
# 10646 observations

meta %>% group_by(Reference) %>% 
  pull(Series) %>% as.numeric() %>% median()
# median of 3 timeseries/experiments per study

meta %>% group_by(Reference, Series) %>% n_groups()
# 551 timeseries/experiments in total

meta %>% group_by(Reference) %>% n_groups()
# 129 studies

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
  summarise(Species = n_distinct(Species)) %>%
  arrange(Species)
# 5 seagrasses
# 6 freshwater plants
# 15 seaweeds
# 66 terrestrial plants

meta %>% 
  group_by(Group) %>% 
  summarise(References = n_distinct(Reference)) %>%
  arrange(References)
# 4 studies on seagrasses
# 12 on freshwater plants
# 17 on seaweeds
# 98 on terrestrial plants


# 1.4 Export table of species ####
meta %>% 
  group_by(Phylum, Order, Family, Species, Group) %>%
  summarise(
    Days = if_else(
      max(Day) < 100, 
      signif(max(Day), digits = 2), 
      signif(max(Day), digits = 3)
    ),
    Studies = n_distinct(Reference),
    Experiments = n_distinct(Reference, Series),
    Observations = n()
  ) %>%
  arrange(Phylum, Order, Family, Species) %>%
  write_csv(here("Tables", "Species.csv"))

# 1.5 Transform data ####
# Data need to be re-expressed as the ratio of the final
# to the initial value in each measurement. This assumes
# that the initial value in the timeseries is representative
# of baseline photosynthesis.

meta %>%
  group_by(Reference, Series) %>%
  filter(Day == min(Day)) %>%
  distinct(Reference, Series, Day) %>%
  ungroup() %>%
  count(Day == 0)
# Most timeseries indeed start with t = 0.

meta %>%
  group_by(Reference, Series) %>%
  slice(1) %>%
  ungroup() %>%
  count(Day == 0)
# This yields the same result so all
# series start with the minimum value.

meta %>%
  group_by(Reference, Series) %>%
  slice(1) %>%
  ungroup() %>%
  filter(Day > 0) %>%
  distinct(Reference, Species, Day) %>%
  arrange(desc(Day))
# The greatest delay is 9 days for Frontier et al. 2021,
# but since this study is on Laminaria which was kelp under
# natural conditions with flowing seawater, I expect 9 days 
# after excision to be representative of the baseline.
# The only other studies that delayed measurement for over 1 day
# are Asim et al. 2022 and Karatas et al. 2010. These studies
# are perhaps more concerning because they are on terrestrial 
# plants which I expect lose photosynthesis more quickly. But
# given that the experimental duration was 56 days in the first
# case and 6-7 days in the second case, I am willing to accept 
# these data.

# There are many cases where multiple initial values were recorded
# either as the initial mean ± standard error or as observations.
# Hence it is important to divide each observation by the mean of
# observations at the initial time. Importantly the mean of simulated
# observations will not exactly match Mean in rnorm(N, Mean, SEM*sqrt(N)).
# So the ratio needs to be calculated separately for cases where
# the initial mean is provided and cases where several initial
# observations were provided. The mean of initial observations
# also covers cases where there is just a single initial observation.
meta %<>%
  group_by(Reference, Series) %>%
  mutate(
    Ratio = if_else( 
      is.na(SEM) & N == 1,
      Observation / mean( Observation[ Day == min(Day) ] ),
      Observation / first( Mean[ Day == min(Day) ] )
    )
  ) %>%
  ungroup() %T>%
  print()

meta %$% all(is.finite(Ratio))

# 1.6 Organise predictors ####
# Create a variable that indexes unique experiments
# across studies, since this will be easier to code
# the model with.
meta %<>%
  group_by(Reference, Series) %>%
  mutate(Experiment = cur_group_id() %>% 
           as_factor()) %>%
  ungroup() %T>%
  print()

# Relevel some factors
meta %<>%
  mutate(
    Group = Group %>%
      fct_relevel("Terrestrial", "Freshwater", "Seagrass"),
    Light = Light %>% fct_relevel("Yes"),
    Water = Water %>% fct_relevel("Yes"),
    Response = Response %>% fct_relevel("Photosynthesis")
  ) %T>%
  print()

# 1.7 Visually explore data ####
# Explore range
meta %>%
  ggplot(aes(Ratio, 0)) +
    geom_jitter(alpha = 0.2, shape = 16) +
    facet_grid(Response ~ Group) +
    theme_minimal()
# One clear negative outlier in seaweed photosynthesis.

meta %>%
  ggplot(aes(Ratio, 0)) +
    geom_jitter(alpha = 0.2, shape = 16) +
    facet_grid(Response ~ Group) +
    coord_cartesian(xlim = c(-1, 3)) +
    theme_minimal()
# Some negative values, as expected for net photosynthesis:
meta %>%
  filter(Ratio < 0) %>%
  ggplot(aes(Ratio, 0)) +
    geom_jitter(alpha = 0.2, shape = 16) +
    facet_grid(Response ~ Group) +
    coord_cartesian(xlim = c(-3, 1)) +
    theme_minimal()
# For the purposes of this meta-analysis it is more sensible 
# to treat net photosynthesis as gross photosynthesis because 
# all other ratios measure this; net gas exchange is just the 
# odd one out. Hence I am replacing all negative values and 
# zeros with the minimal observed positive ratio, assumed to be 
# the minimum observable ratio within measurement error of 0:
meta %<>%
  mutate(
    Ratio = if_else(
      Ratio <= 0, 
      min( Ratio[Ratio > 0] ), 
      Ratio 
    )
  ) %T>%
  print()

meta %$% range(Ratio)

meta %>%
  ggplot(aes(Ratio, 0)) +
    geom_jitter(alpha = 0.2, shape = 16) +
    facet_grid(Response ~ Group) +
    theme_minimal()
# Looks better.

# Explore relationship
meta %>%
  ggplot(aes(Day, Ratio)) +
    geom_point(alpha = 0.2, shape = 16) +
    facet_grid(Response + Light ~ Group + Water, scales = "free_x") +
    theme_minimal()
# Water availability does not make sense to explore because
# only studies on terrestrial plants include cases with no water.

meta %>%
  ggplot(aes(Day, Ratio, colour = Light)) +
    geom_point(alpha = 0.2, shape = 16) +
    facet_grid(Response ~ Group, scales = "free_x") +
    theme_minimal()
# Much clearer. I will explore the predictors Light and Response 
# across groups, species, studies and experiments. 

# Add effect coding
meta %<>%
  mutate( # -0.5/+0.5 effect coding as opposed to 0/1 indicator coding
    Light_effect = if_else(
      Light == "Yes", 0.5, -0.5 
    ),
    Chl_effect = if_else(
      Response == "Chlorophyll", 0.5, -0.5
    )
  ) %T>%
  print()

# 2. Model ####
# 2.1 Prior simulation ####
# I am using the same logistic model as in Seagrass.R but simplified
# because alpha and tau do not need to be estimated:
# 1 / ( 1 + exp( 5 / mu * ( t - mu ) ) ) 
# This causal model does not allow increases in photosynthesis after
# excision because it solely measures decay of photosynthesis, so
# since observed ratios exceed 1
meta %$% max(Ratio)
# these must be accounted for by the likelihood, meaning the likelihood
# cannot be beta (upper bound at 1). The most straightforward likelihood
# that fits is gamma because it is strictly positive and has additive
# variance which I assume for most response variables.

# There is no information on mu for plants at large, but based on the
# literature review and my kelp and seagrass experiments, 2 weeks
# seems reasonable. Nonetheless, prior simulation is necessary:
tibble(n = 1:1e3,
       alpha_mu = rnorm( 1e3 , log(14) , 0.3 ),
       alpha_sigma_g = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       alpha_sigma_s = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       alpha_sigma_e = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       beta_mu_l = rnorm( 1e3 , 0 , 0.4 ),
       beta_mu_c = rnorm( 1e3 , 0 , 0.4 ),
       beta_sigma_l = rtnorm( 1e3 , 0 , 0.4 , 0 ),
       beta_sigma_c = rtnorm( 1e3 , 0 , 0.4 , 0 ),
       mu = exp(
         rnorm( 1e3 , alpha_mu , alpha_sigma_g ) +
           rnorm( 1e3 , 0 , alpha_sigma_s ) +
           rnorm( 1e3 , 0 , alpha_sigma_e ) +
           rnorm( 1e3 , beta_mu_l , beta_sigma_l ) * 0.5 +
           rnorm( 1e3 , beta_mu_c , beta_sigma_c ) * 0.5
       ),
       sigma = rexp( 1e3 , 1 )) %>%
  expand_grid(Day = meta %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(
    Ratio_mu = 1 / ( 1 + exp( 5 / mu * ( Day - mu ) ) ),
    Ratio = rlnorm( n() , log(Ratio_mu) , sigma )
  ) %>%
  pivot_longer(cols = c(Ratio_mu, Ratio),
               names_to = "parameter") %>%
  ggplot(aes(Day, value, group = n)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.2 Stan model ####
require(cmdstanr)
meta_c_model <- here("Meta-analysis", "Stan", "meta_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

meta_nc_model <- here("Meta-analysis", "Stan", "meta_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
meta_c_samples <- meta_c_model$sample(
  data = meta %>%
    select(Day, Ratio, Group,
           Light_effect, Chl_effect, 
           Species, Experiment) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4
)

meta_nc_samples <- meta_nc_model$sample(
  data = meta %>%
    select(Day, Ratio, Group,
           Light_effect, Chl_effect, 
           Species, Experiment) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4
)

# Save draws
meta_c_samples$draws() %>%
  write_rds(here("Meta-analysis", "RDS", "meta_c_samples.rds"))
meta_c_samples$draws(format = "df") %>%
  write_rds(here("Meta-analysis", "RDS", "meta_c_samples_df.rds"))

meta_nc_samples$draws() %>%
  write_rds(here("Meta-analysis", "RDS", "meta_nc_samples.rds"))
meta_nc_samples$draws(format = "df") %>%
  write_rds(here("Meta-analysis", "RDS", "meta_nc_samples_df.rds"))

# 2.3 Model checks ####
# Rhat
meta_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 0.3% of rhat above 1.001. rhat = 1.00 ± 0.000120.

meta_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000117.
# Models are practically identical, but the non-centred 
# model is slightly better.

# Chains
require(bayesplot)
meta_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(units = "cm", height = 100, width = 100, 
         device = cairo_pdf, filename = "meta_c_chains.pdf",
         path = here("Meta-analysis", "Plots"))

meta_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(units = "cm", height = 100, width = 100, 
         device = cairo_pdf, filename = "meta_nc_chains.pdf",
         path = here("Meta-analysis", "Plots"))
# No chains lost, chains look good.

# Pairs
meta_c_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c(
      "alpha_mu", "alpha_g[1]",
      "alpha_s[1]", "alpha_e[1]",
      "beta_mu_l", "beta_l[1]", 
      "beta_mu_c", "beta_c[1]"
    )
  ) %>%
  ggsave(units = "cm", height = 50, width = 50, 
         device = cairo_pdf, filename = "meta_c_pairs.pdf",
         path = here("Meta-analysis", "Plots"))

meta_nc_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c(
      "alpha_mu", "alpha_g[1]",
      "alpha_s[1]", "alpha_e[1]",
      "beta_mu_l", "beta_l[1]", 
      "beta_mu_c", "beta_c[1]"
    )
  ) %>%
  ggsave(units = "cm", height = 50, width = 50, 
         device = cairo_pdf, filename = "meta_nc_pairs.pdf",
         path = here("Meta-analysis", "Plots"))
# No correlations. All looks well.

# 2.4 Prior-posterior comparison ####
# Sample prior
source("functions.R")
meta_prior <- prior_samples(
  model = meta_nc_model, # priors only work well non-centred
  data = meta %>%
    select(Day, Ratio, Group,
           Light_effect, Chl_effect, 
           Species, Experiment) %>%
    compose_data()
  )

# Groups
meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_c_samples,
    group = meta %>% select(Group),
    parameters = c("alpha_mu", "alpha_sigma_g", "alpha_g[Group]",
                   "alpha_sigma_s", "alpha_sigma_e", "sigma",
                   "beta_mu_l", "beta_sigma_l", "beta_l[Group]",
                   "beta_mu_c", "beta_sigma_c", "beta_c[Group]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Group"
  ) %>%
  ggsave(units = "cm", height = 20, width = 30, 
         device = cairo_pdf, filename = "meta_c_prior_posterior_group.pdf",
         path = here("Meta-analysis", "Plots"))

meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_nc_samples,
    group = meta %>% select(Group),
    parameters = c("alpha_mu", "alpha_sigma_g", "alpha_g[Group]",
                   "alpha_sigma_s", "alpha_sigma_e", "sigma",
                   "beta_mu_l", "beta_sigma_l", "beta_l[Group]",
                   "beta_mu_c", "beta_sigma_c", "beta_c[Group]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Group"
  ) %>%
  ggsave(units = "cm", height = 20, width = 30, 
         device = cairo_pdf, filename = "meta_nc_prior_posterior_group.pdf",
         path = here("Meta-analysis", "Plots"))

# Species
meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_c_samples,
    group = meta %>% select(Species),
    parameters = c("alpha_s[Species]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species",
  ) %>%
  ggsave(units = "cm", height = 40, width = 60, 
         device = cairo_pdf, filename = "meta_c_prior_posterior_species.pdf",
         path = here("Meta-analysis", "Plots"))

meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_nc_samples,
    group = meta %>% select(Species),
    parameters = c("alpha_s[Species]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species",
  ) %>%
  ggsave(units = "cm", height = 40, width = 60, 
         device = cairo_pdf, filename = "meta_nc_prior_posterior_species.pdf",
         path = here("Meta-analysis", "Plots"))

# Experiments
meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_c_samples,
    group = meta %>% select(Experiment),
    parameters = c("alpha_e[Experiment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Experiment",
  ) %>%
  ggsave(units = "cm", height = 80, width = 120, 
         device = cairo_pdf, filename = "meta_c_prior_posterior_experiment.pdf",
         path = here("Meta-analysis", "Plots"))

meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_nc_samples,
    group = meta %>% select(Experiment),
    parameters = c("alpha_e[Experiment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Experiment",
  ) %>%
  ggsave(units = "cm", height = 80, width = 120, 
         device = cairo_pdf, filename = "meta_nc_prior_posterior_experiment.pdf",
         path = here("Meta-analysis", "Plots"))
# Both models are practically identical.
# Choose non-centred model because of marginally 
# better performance.

# 2.5 Parameters ####

# 2.6 Continuous prediction ####

# 3. Half-life ####
# 3.1 Data ####
# 3.1.1 Load data ####

# 3.1.2 Check redundancy ####

# 3.2 Model ####
# 3.2.1 Prior simulation ####

# 3.2.2 Stan model ####

# 3.2.3 Model checks ####

# 3.2.4 Prior-posterior comparison ####

# 3.2.5 Parameters ####


# 4. Figures ####
# 4.1 Figure 2 ####

# 4.2 Figure 3 ####

# 5. Tables ####

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