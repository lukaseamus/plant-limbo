#### Plant limbo: detrital photosynthesis meta-analysis ####
#### Luka Seamus Wright                                 ####

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
# median of 3 experiments per study

meta %>% group_by(Reference, Series) %>% n_groups()
# 551 experiments in total

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

meta %>% 
  group_by(Unit %>% str_detect("Proportion") |
             Unit %>% str_detect("%")) %>% 
  summarise(References = n_distinct(Reference)) %>%
  arrange(References)
# 29 out of 129 studies use relative measures
# 105 out of 129 studies use absolute measures
# so 5 studies use both.

# 1.4 Export table of species ####
Table_S1 <- meta %>% 
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
  arrange(as.character(Phylum), as.character(Order), 
          as.character(Family), as.character(Species)) %T>%
  print(n = 100)

Table_S1 %>%
  write_csv(here("Tables", "Table_S1.csv"))

require(officer)
read_docx() %>%
  body_add_table(value = Table_S1) %>%
  print(target = here("Tables", "Table_S1.docx"))

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
         filename = "meta_c_pairs.png",
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
         filename = "meta_nc_pairs.png",
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
    group_name = "Species"
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
    group_name = "Species"
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
    group_name = "Experiment"
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
    group_name = "Experiment"
  ) %>%
  ggsave(units = "cm", height = 80, width = 120, 
         device = cairo_pdf, filename = "meta_nc_prior_posterior_experiment.pdf",
         path = here("Meta-analysis", "Plots"))
# Both models are practically identical.
# Choose non-centred model because of marginally 
# better performance.

# 2.5 Parameter distributions ####
# Global parameters
meta_prior_posterior_global <- meta_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_nc_samples,
    parameters = c("alpha_mu", "alpha_sigma_g", "alpha_sigma_s",
                   "alpha_sigma_e", "beta_mu_l", "beta_sigma_l", 
                   "beta_mu_c", "beta_sigma_c", "sigma"),
    format = "short"
  ) %>%
  mutate( # Calculate parameters for new groups, species and experiments
    alpha = rnorm( n() , alpha_mu , alpha_sigma_g ) +
      rnorm( n() , 0 , alpha_sigma_s ) +
      rnorm( n() , 0 , alpha_sigma_e ),
    beta_l = rnorm( n() , beta_mu_l , beta_sigma_l ),
    beta_c = rnorm( n() , beta_mu_c , beta_sigma_c ),
    mu_Light_Chl = exp( alpha + beta_l * 0.5 + beta_c * 0.5 ),
    mu_Light_P = exp( alpha + beta_l * 0.5 + beta_c * -0.5 ),
    mu_Dark_Chl = exp( alpha + beta_l * -0.5 + beta_c * 0.5 ),
    mu_Dark_P = exp( alpha + beta_l * -0.5 + beta_c * -0.5 )
  ) %T>%
  print()

meta_prior_posterior_global %>% # Save
  write_rds(here("Meta-analysis", "RDS", "meta_prior_posterior_global.rds"))

# Group parameters
meta_prior_posterior_group <- meta %>%
  distinct(Group, Light, Response,
           Light_effect, Chl_effect) %>%
  # Groups are not matched to existing
  # Light and Response effects so this 
  # needs to be done manually
  left_join( # Join Group samples
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = meta %>% select(Group),
        parameters = c("alpha_g[Group]", "beta_l[Group]",
                       "beta_c[Group]", "alpha_sigma_s",
                       "alpha_sigma_e", "sigma"),
        format = "short"
      ),
    by = "Group",
    relationship = "many-to-many"
  ) %>%
  mutate( # Calculate group parameters for new species and 
    # experiments but existing light and response effects
    alpha = rnorm( n() , alpha_g , alpha_sigma_s ) +
      rnorm( n() , 0 , alpha_sigma_e ),
    mu = exp( alpha + beta_l * Light_effect + beta_c * Chl_effect )
  ) %>%
  filter(Group == "Terrestrial" & Light == "Yes" & 
           Response == "Photosynthesis" & distribution == "prior" |
           distribution == "posterior") %>% # Remove redundant priors
  mutate( # Embed prior in Group, Light and Response
    Group = if_else(
      distribution == "prior", "Prior", Group
    ) %>% fct(),
    Light = if_else(
      distribution == "prior", "Prior", Light
    ) %>% fct(),
    Response = if_else(
      distribution == "prior", "Prior", Response
    ) %>% fct()
  ) %>%
  select(starts_with("."), Group, Light, Response,
         alpha, beta_l, beta_c, mu, sigma) %T>%
  print()

meta_prior_posterior_group %>% # Save
  write_rds(here("Meta-analysis", "RDS", "meta_prior_posterior_group.rds"))

# Combine global and group parameters
meta_mu <- meta_prior_posterior_group %>%
  # alpha is identical for all levels of Light and Response
  filter(Group == "Prior" |
           Light == "Yes" & Response == "Chlorophyll") %>%
  select(starts_with("."), Group, alpha, sigma) %>%
  mutate(mu = exp(alpha)) %>% # the group mean of mu is just exp(alpha)
  bind_rows(
    meta_prior_posterior_global %>%
      filter(distribution == "posterior") %>%
      select(starts_with("."), alpha, sigma) %>%
      mutate(Group = "Global" %>% fct(),
             mu = exp(alpha))
  ) %T>%
  print()

# Species parameters
meta_prior_posterior_species <- meta %>%
  distinct(Group, Species, Light, Response,
           Light_effect, Chl_effect) %>%
  # Species are not matched to their Group
  # or Light and Response effects so this 
  # needs to be done manually
  left_join( # Join Group samples
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = meta %>% select(Group),
        parameters = c("alpha_g[Group]", 
                       "beta_l[Group]",
                       "beta_c[Group]"),
        format = "short"
      ),
    by = "Group",
    relationship = "many-to-many"
  ) %>%
  left_join( # Join Species samples
    meta_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_nc_samples,
        group = meta %>% select(Species),
        parameters = c("alpha_s[Species]", 
                       "alpha_sigma_e"),
        format = "short"
      ),
    by = c("Species", ".chain", ".iteration", 
           ".draw", "distribution"),
    relationship = "many-to-many"
  ) %>%
  mutate( # Calculate species parameters for new experiments
    # but existing light and response effects
    alpha = rnorm( n() , alpha_g + alpha_s , alpha_sigma_e ),
    mu = exp( alpha + beta_l * Light_effect + beta_c * Chl_effect )
  ) %>%
  filter(Group == "Terrestrial" & Species == "Acer saccharinum" &
           Light == "Yes" & Response == "Photosynthesis" &
           distribution == "prior" | # Remove redundant priors
           distribution == "posterior") %>%
  mutate( # Embed prior in Group, Species, Light and Response
    Group = if_else(
      distribution == "prior", "Prior", Group
    ) %>% fct(),
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct(),
    Light = if_else(
      distribution == "prior", "Prior", Light
    ) %>% fct(),
    Response = if_else(
      distribution == "prior", "Prior", Response
    ) %>% fct()
  ) %>%
  select(starts_with("."), Group, Species, 
         Light, Response, alpha, mu) %T>%
  print()

meta_prior_posterior_species %>% # Save
  write_rds(here("Meta-analysis", "RDS", "meta_prior_posterior_species.rds"))

# Experiment parameters
# Prior is already fully contained in the other prior_posterior datasets
# and is too memory intensive to compute here too. Only calculate posterior.
meta_posterior_experiment <- meta %>%
  distinct(Group, Species, Experiment, Light, 
           Response, Light_effect, Chl_effect) %>%
  # For this I can use original tidybayes functions
  left_join( # Join Group samples
    meta_nc_samples %>%
      recover_types(meta %>% select(Group)) %>%
      spread_draws(alpha_g[Group], beta_l[Group], beta_c[Group]) %>%
      ungroup(),
    by = "Group",
    relationship = "many-to-many"
  ) %>%
  left_join( # Join Species samples
    meta_nc_samples %>%
      recover_types(meta %>% select(Species)) %>%
      spread_draws(alpha_s[Species]) %>%
      ungroup(),
    by = c("Species", ".chain", ".iteration", ".draw"),
    relationship = "many-to-many"
  ) %>%
  left_join( # Join Experiment samples
    meta_nc_samples %>%
      recover_types(meta %>% select(Experiment)) %>%
      spread_draws(alpha_e[Experiment]) %>%
      ungroup(),
    by = c("Experiment", ".chain", ".iteration", ".draw"),
    relationship = "many-to-many"
  ) %>%
  mutate( # Calculate parameters for existing experiments
    alpha = alpha_g + alpha_s + alpha_e, # No added variance
    mu = exp( alpha + beta_l * Light_effect + beta_c * Chl_effect )
  ) %>%
  select(starts_with("."), Group, Species, Experiment,
         Light, Response, alpha, mu) %T>%
  print()

meta_posterior_experiment %>% # Save
  write_rds(here("Meta-analysis", "RDS", "meta_posterior_experiment.rds"))

# Clean up
rm(meta_c_model, meta_nc_model, meta_prior,
   meta_c_samples, meta_nc_samples)
gc()

# 2.6 Parameter estimates ####
require(glue)
# Global
meta_hyperparameters <- meta_prior_posterior_global %>%
  filter(distribution == "posterior") %>%
  select(-c(distribution, starts_with("."), starts_with("mu"))) %>%
  mutate(ratio_l = exp(beta_l), # betas are log ratios
         ratio_c = exp(beta_c)) %>%
  pivot_longer(cols = everything(),
               names_to = "Parameter",
               values_to = "Sample") %>%
  group_by(Parameter) %>%
  summarise(
    mean = mean(Sample),
    sd = sd(Sample),
    median = median(Sample),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    Summary = glue(
      "{signif(mean, 2)} ± {signif(sd, 2)} ({signif(median, 2)})"
    )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()
# Medians and means match well, so distributions are all
# fairly normal, but this won't be the case for mu, so
# I'll summarise that as median (log mean ± log sd).

meta_mu_summary <- meta_mu %>%
  group_by(Group) %>%
  summarise(
    across(
      everything(), 
      list(
        mean = mean, 
        sd = sd, 
        median = median
      )
    ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    mu = glue(
      "{signif(mu_median, 2)} ({signif(alpha_mean, 2)} ± {signif(alpha_sd, 2)})"
    )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

meta_mu_effect_global <- meta_prior_posterior_global %>%
  filter(distribution == "posterior") %>%
  select(starts_with("mu")) %>%
  pivot_longer(cols = everything(),
               names_to = "Light_Response",
               values_to = "mu",
               names_prefix = "mu_") %>%
  mutate(log_mu = log(mu)) %>%
  group_by(Light_Response) %>%
  summarise(
    across(
      everything(), 
      list(
        mean = mean, 
        sd = sd, 
        median = median
      )
    ),
    n = n()
  ) %>%
  ungroup() %>%
  separate(Light_Response, into = c("Light", "Response"),
           sep = "_") %>%
  mutate(
    mu = glue(
      "{signif(mu_median, 2)} ({signif(log_mu_mean, 2)} ± {signif(log_mu_sd, 2)})"
    )
  ) %>%
  arrange(mu_median) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Group
meta_mu_effect_group <- meta_prior_posterior_group %>%
  filter(Group != "Prior") %>%
  droplevels() %>%
  select(Group, Light, Response, mu) %>%
  mutate(log_mu = log(mu)) %>%
  group_by(Group, Light, Response) %>%
  summarise(
    across(
      everything(), 
      list(
        mean = mean, 
        sd = sd, 
        median = median
      )
    ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    mu = glue(
      "{signif(mu_median, 2)} ({signif(log_mu_mean, 2)} ± {signif(log_mu_sd, 2)})"
    )
  ) %>%
  arrange(mu_median) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Species
meta_mu_effect_species <- meta_prior_posterior_species %>%
  filter(Group != "Prior") %>%
  droplevels() %>%
  select(Group, Species, Light, Response, mu) %>%
  mutate(log_mu = log(mu)) %>%
  group_by(Group, Species, Light, Response) %>%
  summarise(
    across(
      everything(), 
      list(
        mean = mean, 
        sd = sd, 
        median = median
      )
    ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    mu_median_rounded = case_when(
      mu_median < 100 ~ signif(mu_median, 2),
      mu_median < 1e3 ~ signif(mu_median, 3)
    ),
    mu = glue(
      "{mu_median_rounded} ({signif(log_mu_mean, 2)} ± {signif(log_mu_sd, 2)})"
    )
  ) %>%
  arrange(mu_median) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Experiments
meta_mu_effect_experiment <- meta_posterior_experiment %>%
  select(Group, Species, Experiment, Light, Response, mu) %>%
  mutate(log_mu = log(mu)) %>%
  group_by(Group, Species, Experiment, Light, Response) %>%
  summarise(
    across(
      everything(), 
      list(
        mean = mean, 
        sd = sd, 
        median = median
      )
    ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    mu_median_rounded = case_when(
      mu_median < 100 ~ signif(mu_median, 2),
      mu_median < 1e3 ~ signif(mu_median, 3)
    ),
    mu = glue(
      "{mu_median_rounded} ({signif(log_mu_mean, 2)} ± {signif(log_mu_sd, 2)})"
    )
  ) %>%
  arrange(mu_median) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# 2.7 Contrast distributions ####
# Effects
meta_beta <- meta_prior_posterior_group %>%
  # Remove redundant posteriors (repeated effects) by filtering
  # for Light and Response that is present in all groups
  filter(Group == "Prior" |
           Light == "Yes" & Response == "Chlorophyll") %>%
  select(starts_with("."), Group, starts_with("beta")) %>%
  pivot_longer(cols = starts_with("beta"), 
               names_to = "Effect", 
               values_to = "beta",
               names_prefix = "beta_") %>%
  bind_rows(
    meta_prior_posterior_global %>%
      filter(distribution == "posterior") %>%
      select(starts_with("."), beta_l, beta_c) %>%
      pivot_longer(cols = starts_with("beta"), 
                   names_to = "Effect", 
                   values_to = "beta",
                   names_prefix = "beta_") %>%
      mutate(Group = "Global" %>% fct())
  ) %>%
  mutate(ratio = exp(beta), # beta is the difference of logs, so a log ratio
         Effect = Effect %>% fct()) %T>% 
  print()
 
# Group means
# The light and temperature effects are already expressed as contrasts
# but the group means marginalised over these effects also need contrasting.
meta_contrast <- meta_mu %>%
  filter(!Group %in% c("Prior", "Global")) %>%
  select(-c(alpha, sigma)) %>%
  pivot_wider(names_from = Group,
              values_from = mu) %>%
  mutate(
    # Calculate differences
    Seaweed_Seagrass_diff = Seaweed - Seagrass,
    Seaweed_Freshwater_diff = Seaweed - Freshwater,
    Seaweed_Terrestrial_diff = Seaweed - Terrestrial,
    Seagrass_Freshwater_diff = Seagrass - Freshwater,
    Seagrass_Terrestrial_diff = Seagrass - Terrestrial,
    Freshwater_Terrestrial_diff = Freshwater - Terrestrial,
    # Calculate ratios
    Seaweed_Seagrass_ratio = Seaweed / Seagrass,
    Seaweed_Freshwater_ratio = Seaweed / Freshwater,
    Seaweed_Terrestrial_ratio = Seaweed / Terrestrial,
    Seagrass_Freshwater_ratio = Seagrass / Freshwater,
    Seagrass_Terrestrial_ratio = Seagrass / Terrestrial,
    Freshwater_Terrestrial_ratio = Freshwater / Terrestrial
  ) %>%
  select(-c(Terrestrial, Freshwater, Seagrass, Seaweed)) %>%
  pivot_longer(cols = -c(.chain, .iteration, .draw)) %>%
  # This takes long but str_extract/str_remove is similar
  separate(name, into = c("Contrast", "Type"),
           sep = "_(?=[^_]+$)") %>% 
  pivot_wider(names_from = Type,
              values_from = value) %>%
  mutate(Contrast = Contrast %>% fct()) %T>%
  print()
  
# 2.8 Contrast estimates ####
# Effects
meta_beta_summary <- meta_beta %>%
  group_by(Group, Effect) %>%
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
    P = mean( beta > 0 ) %>% signif(2)
  ) %>%
  ungroup() %>%
  mutate(
    ratio = glue(
      "{signif(ratio_median, 2)} ({signif(beta_mean, 2)} ± {signif(beta_sd, 2)})"
    )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Group means
meta_contrast_summary <- meta_contrast %>%
  mutate(log_ratio = log(ratio)) %>%
  group_by(Contrast) %>%
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
    P = mean( diff > 0 ) %>% signif(2)
  ) %>%
  ungroup() %>%
  mutate(
    diff = glue(
      "{signif(diff_mean, 2)} ± {signif(diff_sd, 2)} ({signif(diff_median, 2)})"
    ),
    ratio = glue(
      "{signif(ratio_median, 2)} ({signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 2)})"
    )
  ) %>%
  arrange(desc(ratio_median)) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()
# Difference median and mean are far from matching and it cannot be
# nicely expressed in log space, Ratio is much better in that regard.

# 2.9 Continuous prediction ####
# Group means
# Remove redundant posteriors (no effects) by filtering
# for Light and Response that is present in all groups
meta_prediction_group <- meta_mu %>% 
  filter(Group != "Global") %>%
  droplevels() %>%
  spread_continuous(data = meta,
                    group_name = "Group",
                    predictor_name = "Day",
                    length = 200) %>%
  mutate(
    Ratio_mu = 1 / ( 1 + exp( 5 / mu * ( Day - mu ) ) ),
    Ratio = rlnorm( n() , log(Ratio_mu) , sigma )
  ) %>% # Summarise predictions
  group_by(Day, Group) %>%
  median_qi(Ratio_mu, Ratio, .width = c(.5, .8, .9)) %T>%
  print()

meta_prediction_group %>% # Save
  write_rds(here("Meta-analysis", "RDS", "meta_prediction_group.rds"))

# Groups with effects
# Too many groupings makes spread_continuous() slow and awkward,
# so I'll nest and calculate sequentially.
meta_prediction_group_effect <- meta_prior_posterior_group %>%
  group_by(Group, Light, Response) %>% 
  nest(.key = "prior_posterior") %>%
  ungroup() %>% # Nesting causes grouping which should be undone
  mutate(
    predictor = case_when(
        Group == "Terrestrial" ~ list( seq(0, 40, length.out = 100) ),
        Group == "Freshwater" ~ list( seq(0, 60, length.out = 100) ),
        Group == "Seagrass" ~ list( seq(0, 100, length.out = 100) ),
        Group %in% c("Seaweed", "Prior") ~ list( seq(0, 300, length.out = 100) ),
      ),
    prediction = map2(
      prior_posterior, predictor,
      ~.x %>% 
        slice( rep( 1:n() , each = length(.y) ) ) %>%
        mutate(
          Day = rep( .y , times = nrow(.x) ),
          Ratio_mu = 1 / ( 1 + exp( 5 / mu * ( Day - mu ) ) ),
          Ratio = rlnorm( n() , log(Ratio_mu) , sigma )
        ) %>%
        group_by(Day) %>%
        median_qi(Ratio_mu, Ratio, .width = c(.5, .8, .9)) %T>%
        print() # printing helps keep track of progress
    )
  ) %>% 
  select(Group, Light, Response, prediction) %>%
  unnest(prediction) %T>%
  print()

meta_prediction_group_effect %>% # Save
  write_rds(here("Meta-analysis", "RDS", "meta_prediction_group_effect.rds"))

# Experiments
meta_prediction_experiment <- meta_posterior_experiment %>%
  group_by(Group, Experiment) %>% # Too many experiments to run in one go
  nest(.key = "posterior") %>% # -> nesting and sequential calculation
  ungroup() %>% # Nesting causes grouping which should be undone
  left_join(
    meta %>%
      group_by(Group, Experiment) %>%
      nest() %>%
      ungroup(),
    by = c("Group", "Experiment")
  ) %>%
  mutate(
    predictor = data %>% map(
      # Ensure that min(Day) = 0 for prediction
      ~seq( 0 , max(.x$Day) , length.out = 100 )
    ),
    prediction = map2(
      posterior, predictor,
      ~.x %>% 
        slice( rep( 1:n() , each = length(.y) ) ) %>%
        mutate(
          Day = rep( .y , times = nrow(.x) ),
          Ratio_mu = 1 / ( 1 + exp( 5 / mu * ( Day - mu ) ) )
        ) %>%
        group_by(Day) %>%
        median_qi(Ratio_mu, .width = c(.5, .8, .9)) %T>%
        print() # printing helps keep track of progress
    )
  ) %>% 
  select(Group, Experiment, prediction) %>%
  unnest(prediction) %>% 
  ungroup() %T>%
  print()

meta_prediction_experiment %>% # Save
  write_rds(here("Meta-analysis", "RDS", "meta_prediction_experiment.rds"))

# 3. Half-life ####
# 3.1 Data ####
# 3.1.1 Load ForC data ####
# These data can be downloaded from github.com/forc-db/ForC/tree/master/data.
# The data publication is Anderson-Teixeira 2018 (doi: 10.1002/ecy.2229).
ForC <- here("Meta-analysis", "ForC_measurements.csv") %>% 
  read_csv() %T>%
  print()
# Some variables are problematic
problems(ForC) %>% print(n = 30)
ForC[,c(37, 40)] # Not important

# 3.1.2 Select relevant variables #### 
# Consult ForC_variables.csv from github.com/forc-db/ForC/tree/master/data
ForC %<>%
  select(citation.ID, measurement.ID, sites.sitename, plot.name, 
         stand.age, dominant.veg, dominant.life.form, scientific.name, 
         variable.name, date, mean, required.citations) %>%
  filter(
    variable.name %in% 
      c("ANPP_litterfall_1_C", "ANPP_litterfall_0_C", "R_het_litter_C", 
        "litter_C", "organic.layer_C", "O.horizon_C")
  ) %>%
  droplevels() %>%
  mutate(
    Reference = citation.ID %>% str_remove("_[^_]*$"),
    Site = sites.sitename,
    Plot = plot.name,
    Review = required.citations %>% str_remove("_[^_]*$")
  ) %>%
  select(-c(citation.ID, sites.sitename, plot.name, required.citations)) %T>%
  print()

# 3.1.3 Tidy data #### 
ForC %<>%
  select(-measurement.ID) %>% # Measurment ID is unique per variable
  pivot_wider(names_from = variable.name,
              values_from = mean,
              # Summarise duplicates as mean per site
              values_fn = mean) %T>%
  print()

# 3.1.4 Coalesce and calculate #### 
ForC %<>%
  mutate(
    litterfall = coalesce(ANPP_litterfall_1_C, ANPP_litterfall_0_C),
    litter = coalesce(litter_C, organic.layer_C, O.horizon_C),
    # All fluxed are given in Mg C ha^-1 y^-1 and all stocks in Mg C ha^-1,
    # so the following rates have units y^-1. Divide by 365 to get d^-1.
    PB = litterfall / litter, # detrital turnover from detrital production
    DB = R_het_litter_C / litter, # detrital turnover from decomposition
    # DB directly measures decomposition, so is preferred over PB
    k = coalesce(DB, PB) / 365 # (d^-1) 
  ) %>%
  select(-c(contains("litter"), ends_with("_C"))) %>%
  filter(is.finite(k)) %>% # drop NAs and Infs
  droplevels() %T>% # drop unused factor levels
  print()

ForC %>% nrow() # 103 observations
ForC %>% group_by(Reference) %>% n_groups() # 46 references
ForC %>% distinct(Review) 
# Anderson-Teixeira 2018 (doi: 10.1002/ecy.2229) is the reference to cite

# 3.1.5 Load Just's data #### 
# See github.com/lukaseamus/plant-carbon for details
Just <- here("Meta-analysis", "Carbon.csv") %>% 
  read_csv() %T>%
  print()

# 3.1.6 Select and calculate relevant data ####
Just %<>%
  mutate( # Create plant group to match Group in meta
    Group = case_when( # There is some disagreement between Biome and Community
      Biome == "Macroalgal forest" ~ "Seaweed",
      Biome == "Seagrass meadow" ~ "Seagrass",
      Biome %in% 
        c("Woodland", "Grassland", "Tundra", 
          "Saltmarsh", "Mangrove forest") | # Wetland contains freshwater plants
        Community == "freshwater and marine marshes" ~ "Terrestrial",
      # Only submerged freshwater plants are considered for Freshwater
      Community == "meadows of freshwater submerged macrophytes" ~ "Freshwater"
    )
  ) %>%
  drop_na(Group) %>% # drop all others (microphytobenthos and phytoplankton)
  mutate(
    # Same as before, detrital turnover from detrital production (y^-1)
    PB = `Detrital production (g C m^-2 y^-1)` / `Detrital biomass (g C m^-2)`,
    # Same as before, detrital turnover from decomposition (y^-1)
    DB = `Decomposition (g C m^-2 y^-1)` / `Detrital biomass (g C m^-2)`,
    # k (d^-1) is also directly reported and obviously preferred
    k = coalesce(`k (d^-1)`, DB/365, PB/365) # (d^-1)
  ) %>%
  drop_na(k) %>% # Filter out NAs
  distinct(Reference, Group, k) %T>% # Filter out duplicates within reference
  print()

Just %>% nrow() # 212 observations
Just %>% group_by(Reference) %>% n_groups() # 85 references

# 3.1.7 Check redundancy and merge ####
ForC %>%
  distinct(Reference) %>%
  mutate(Reference = Reference %>% str_extract("^[^_ ]+")) %>%
  inner_join(
    Just %>%
      distinct(Reference) %>%
      mutate(Reference = Reference %>% str_extract("^[^_ ]+"))
  )
# Two potentially redundant references

ForC %>% 
  filter(Reference %>% str_detect("O'Connell")) %>%
  select(Reference, k) %>%
  mutate(data = "ForC") %>%
  bind_rows(
    Just %>% 
      filter(Reference %>% str_detect("O'Connell")) %>%
      select(Reference, k) %>%
      mutate(data = "Just")
  )
# O'Connell is not redundant

ForC %>% 
  filter(Reference %>% str_detect("Fahey")) %>%
  select(Reference, k) %>%
  mutate(data = "ForC") %>%
  bind_rows(
    Just %>% 
      filter(Reference %>% str_detect("Fahey")) %>%
      select(Reference, k) %>%
      mutate(data = "Just")
  )
# Fahey is not redundant
# No redundancy after all

# Merge
deco <- ForC %>%
  select(Reference, k) %>%
  # Add plant group to ForC
  mutate(Group = "Terrestrial") %>%
  bind_rows(Just) %T>%
  print()

deco %>% nrow() # 315 observations
deco %>% group_by(Reference) %>% n_groups() # 131 references

deco %>%
  group_by(Group) %>%
  count()
# Seaweeds are severely underrepresented. Fortunately I compiled
# a meta-dataset of k values for seaweeds.

# 3.1.8 Load macroalgal k data #### 
# See github.com/lukaseamus/Ma-k-roalgen for details
Makro <- here("Meta-analysis", "k.csv") %>% 
  read_csv() %T>%
  print()
# This dataset goes into far more detail than required for my
# purpose here, so I will simplify by summarising k for each
# species and reference as the mean across reviews. This is
# the only way if I want to compare these data to the others.
# I will also filter out small-scale studies because I believe 
# this is more representative of the other data. Finally, this
# dataset contains negative k values because seaweed detritus
# can grow. These need to be removed to calculate half-life.

# 3.1.9 Merge with deco #### 
# Filter and average across reviews
Makro %<>%
  mutate(
    # I am counting Filbee-Dexter et al. 2024 and
    # Moten Foldager Pedersen unpublished as one
    # review since they were not conducted independently.
    Review = if_else(
      Review %>% str_detect("Filbee-Dexter et al. 2024"),
      "Filbee-Dexter et al. 2024", Review
    ),
    # Similarly, I am counting these references as one.
    Reference = if_else(
      Reference == "Filbee-Dexter et al. unpublished",
      "Filbee-Dexter et al. 2022", Reference
    )
  ) %>%
  filter(!Dissolved & !Smallscale) %>%
  group_by(Reference, Species, Wrack, Prekilled, Buried, 
           Anoxic, Lab, Darkness) %>%
  summarise(
    k = mean(k),
    n_review = n_distinct(Review),
    n = n()
  ) %>%
  ungroup() %>%
  filter(k > 0) %T>%
  print()
# n varies quite a bit, but I can think, mostly because of 
# varying number of reviews, but sometimes because data
# are not correctly stratified because not all levels are
# encoded in the group_by call. Nonetheless, this is the
# best I can think of to represent all reviews.

# Data in Just are accounted for in Makro, so all I need 
# to do is replace the seaweed observations in deco with
# those in Makro.
deco %<>%
  filter(Group != "Seaweed") %>%
  bind_rows(
    Makro %>% 
      select(Reference, k) %>%
      mutate(Group = "Seaweed")
  ) %>%
  mutate(
    Reference = Reference %>% fct(),
    Group = Group %>% 
      fct_relevel("Terrestrial", "Freshwater", "Seagrass")
  ) %>%
  select(Reference, Group, k) %T>%
  print()

# 3.1.10 Calculate half-life #### 
deco %<>%
  mutate(
    Halflife = log(2) / k # half-life in days
  ) %T>%
  print()

deco %>% nrow() # 475 observations
deco %>% group_by(Reference) %>% n_groups() # 207 studies
deco %>%
  group_by(Group) %>%
  count()
# More balanced but most observations are on
# terrestrial plants and seaweeds, similar
# to the imbalance in meta.

# Clean up
rm(ForC, Just, Makro)

# 3.2 Model ####
# 3.2.1 Stan model ####
# I am picking the same prior as for photosynthetic half-life, 
# lognormal on 14 days, albeit with more uncertainty.
halflife_c_model <- here("Meta-analysis", "Stan", "halflife_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

halflife_nc_model <- here("Meta-analysis", "Stan", "halflife_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

halflife_c_samples <- halflife_c_model$sample(
  data = deco %>%
    select(-k) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4
)

halflife_nc_samples <- halflife_nc_model$sample(
  data = deco %>%
    select(-k) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4
)

# Save draws
halflife_c_samples$draws() %>%
  write_rds(here("Meta-analysis", "RDS", "halflife_c_samples.rds"))
halflife_c_samples$draws(format = "df") %>%
  write_rds(here("Meta-analysis", "RDS", "halflife_c_samples_df.rds"))

halflife_nc_samples$draws() %>%
  write_rds(here("Meta-analysis", "RDS", "halflife_nc_samples.rds"))
halflife_nc_samples$draws(format = "df") %>%
  write_rds(here("Meta-analysis", "RDS", "halflife_nc_samples_df.rds"))

# 3.2.2 Model checks ####
# Rhat
halflife_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000917.

halflife_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000792.
# Models are identical, but the centred model 
# ran faster and had no divergent transitions.

# Chains
halflife_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(units = "cm", height = 40, width = 40, 
         device = cairo_pdf, filename = "halflife_c_chains.pdf",
         path = here("Meta-analysis", "Plots"))

halflife_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay() %>%
  ggsave(units = "cm", height = 60, width = 60, 
         device = cairo_pdf, filename = "halflife_nc_chains.pdf",
         path = here("Meta-analysis", "Plots"))
# No chains lost, chains look good.

# Pairs
halflife_c_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c(
      "alpha_mu", "alpha_g[1]", "alpha_g[2]",
      "alpha_g[3]", "alpha_g[4]", "alpha_r[1]",
      "alpha_sigma_g", "alpha_sigma_r"
    )
  ) %>%
  ggsave(units = "cm", height = 50, width = 50, 
         filename = "halflife_c_pairs.png",
         path = here("Meta-analysis", "Plots"))

halflife_nc_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c(
      "alpha_mu", "alpha_g[1]", "alpha_g[2]",
      "alpha_g[3]", "alpha_g[4]", "alpha_r[1]",
      "alpha_sigma_g", "alpha_sigma_r"
    )
  ) %>%
  ggsave(units = "cm", height = 50, width = 50, 
         filename = "halflife_nc_pairs.png",
         path = here("Meta-analysis", "Plots"))
# No correlations. All looks well.

# 3.2.3 Prior-posterior comparison ####
# Sample prior
halflife_prior <- prior_samples(
  model = halflife_nc_model, # priors only work well non-centred
  data = deco %>%
    select(-k) %>%
    compose_data()
  )

# Groups
halflife_prior %>% 
  prior_posterior_draws(
    posterior_samples = halflife_c_samples,
    group = deco %>% select(Group),
    parameters = c("alpha_mu", "alpha_sigma_g", "alpha_g[Group]",
                   "alpha_sigma_r", "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Group"
  ) %>%
  ggsave(units = "cm", height = 20, width = 30, 
         device = cairo_pdf, filename = "halflife_c_prior_posterior_group.pdf",
         path = here("Meta-analysis", "Plots"))

halflife_prior %>% 
  prior_posterior_draws(
    posterior_samples = halflife_nc_samples,
    group = deco %>% select(Group),
    parameters = c("alpha_mu", "alpha_sigma_g", "alpha_g[Group]",
                   "alpha_sigma_r", "sigma"),
    format = "long"
  ) %>%
  prior_posterior_plot(
    group_name = "Group"
  ) %>%
  ggsave(units = "cm", height = 20, width = 30, 
         device = cairo_pdf, filename = "halflife_nc_prior_posterior_group.pdf",
         path = here("Meta-analysis", "Plots"))

# References
halflife_prior %>% 
  prior_posterior_draws(
    posterior_samples = halflife_c_samples,
    group = deco %>% select(Reference),
    parameters = c("alpha_r[Reference]"),
    format = "long"
  ) %>%
  prior_posterior_plot(
    group_name = "Reference"
  ) %>%
  ggsave(units = "cm", height = 80, width = 120, 
         device = cairo_pdf, filename = "halflife_c_prior_posterior_reference.pdf",
         path = here("Meta-analysis", "Plots"))

halflife_prior %>% 
  prior_posterior_draws(
    posterior_samples = halflife_nc_samples,
    group = deco %>% select(Reference),
    parameters = c("alpha_r[Reference]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Reference"
  ) %>%
  ggsave(units = "cm", height = 80, width = 120, 
         device = cairo_pdf, filename = "halflife_nc_prior_posterior_reference.pdf",
         path = here("Meta-analysis", "Plots"))
# Predictions are identical. Proceed with centred model.

# 3.2.4 Parameter distributions ####
# Global
halflife_prior_posterior_global <- halflife_prior %>% 
  prior_posterior_draws(
    posterior_samples = halflife_c_samples,
    parameters = c("alpha_mu", "alpha_sigma_g", 
                   "alpha_sigma_r", "sigma"),
    format = "short"
  ) %>%
  mutate( # Calculate halflife for new groups and studies
    mu = rnorm( n() , alpha_mu , alpha_sigma_g ) +
          rnorm( n() , 0 , alpha_sigma_r ),
    Halflife = rlnorm( n() , mu , sigma )
  ) %>%
  select(starts_with("."), distribution, mu, Halflife) %T>%
  print()

halflife_prior_posterior_global %>% # Save
  write_rds(here("Meta-analysis", "RDS", "halflife_prior_posterior_global.rds"))

# Group
halflife_prior_posterior_group <- halflife_prior %>% 
  prior_posterior_draws(
    posterior_samples = halflife_c_samples,
    group = deco %>% select(Group),
    parameters = c("alpha_g[Group]", "alpha_sigma_r", "sigma"),
    format = "short"
  ) %>%
  mutate( # Calculate group parameters for new studies
    mu = rnorm( n() , alpha_g , alpha_sigma_r ),
    Halflife = rlnorm( n() , mu , sigma )
  ) %>%
  filter(Group == "Terrestrial" & distribution == "prior" |
           distribution == "posterior") %>% # Remove redundant priors
  mutate( # Embed prior in Group, Light and Response
    Group = if_else(
      distribution == "prior", "Prior", Group
    ) %>% fct()
  ) %>%
  select(starts_with("."), Group, mu, Halflife) %T>%
  print()

halflife_prior_posterior_group %>% # Save
  write_rds(here("Meta-analysis", "RDS", "halflife_prior_posterior_group.rds"))

# Combine global and group parameters
halflife_prior_posterior <- halflife_prior_posterior_group %>%
  bind_rows(
    halflife_prior_posterior_global %>%
      filter(distribution == "posterior") %>%
      select(-distribution) %>%
      mutate(Group = "Global" %>% fct())
  ) %T>%
  print()
  
# 3.2.5 Parameter estimates ####
halflife_estimates <- halflife_prior_posterior %>%
  # I want to express Halflife on the same log scale as mu
  mutate(log_Halflife = log(Halflife),
         exp_mu = exp(mu)) %>%
  group_by(Group) %>%
  summarise(
    across(
      everything(), 
      list(
        mean = mean, 
        sd = sd, 
        median = median
      )
    ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    exp_mu_median_rounded = case_when(
      exp_mu_median < 100 ~ signif(exp_mu_median, 2),
      exp_mu_median < 1e3 ~ signif(exp_mu_median, 3)
    ),
    mu = glue(
      "{exp_mu_median_rounded} ({signif(mu_mean, 2)} ± {signif(mu_sd, 2)})"
    ),
    Halflife_median_rounded = case_when(
      Halflife_median < 100 ~ signif(Halflife_median, 2),
      Halflife_median < 1e3 ~ signif(Halflife_median, 3)
    ),
    Halflife = glue(
      "{Halflife_median_rounded} ({signif(log_Halflife_mean, 2)} ± {signif(log_Halflife_sd, 2)})"
    )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Clean up
rm(halflife_c_model, halflife_nc_model, 
   halflife_c_samples, halflife_nc_samples,
   halflife_prior)
gc()

# 3.2.6 Contrast distributions ####
halflife_contrast <- halflife_prior_posterior_group %>%
  filter(Group != "Prior") %>%
  select(-mu) %>%
  pivot_wider(names_from = Group,
              values_from = Halflife) %>%
  mutate(
    # Calculate differences
    Seaweed_Seagrass_diff = Seaweed - Seagrass,
    Seaweed_Freshwater_diff = Seaweed - Freshwater,
    Seaweed_Terrestrial_diff = Seaweed - Terrestrial,
    Seagrass_Freshwater_diff = Seagrass - Freshwater,
    Seagrass_Terrestrial_diff = Seagrass - Terrestrial,
    Freshwater_Terrestrial_diff = Freshwater - Terrestrial,
    # Calculate ratios
    Seaweed_Seagrass_ratio = Seaweed / Seagrass,
    Seaweed_Freshwater_ratio = Seaweed / Freshwater,
    Seaweed_Terrestrial_ratio = Seaweed / Terrestrial,
    Seagrass_Freshwater_ratio = Seagrass / Freshwater,
    Seagrass_Terrestrial_ratio = Seagrass / Terrestrial,
    Freshwater_Terrestrial_ratio = Freshwater / Terrestrial
  ) %>%
  select(-c(Terrestrial, Freshwater, Seagrass, Seaweed)) %>%
  pivot_longer(cols = -c(.chain, .iteration, .draw)) %>%
  # This takes long but str_extract/str_remove is similar
  separate(name, into = c("Contrast", "Type"),
           sep = "_(?=[^_]+$)") %>% 
  pivot_wider(names_from = Type,
              values_from = value) %>%
  mutate(Contrast = Contrast %>% fct()) %T>%
  print()

# Contrast with mu
# To express the relative physiological determinism of decomposition, 
# I propose the ratio of the photosynthetic half-life, that is mu,
# to detrital half-life, that is log(2) / k. I need to use an average
# for mu and the easiest way is to marginalise over the light and
# chlorophyll effects, that is use exp(alpha) instead of mu.
mu_vs_halflife <- meta_mu %>% 
  select(-c(alpha, sigma)) %>%
  full_join(
    halflife_prior_posterior %>% select(-mu),
    by = c(".chain", ".iteration", ".draw", "Group")
  ) %>%
  mutate( # Calculate difference and ratio
    diff = mu - Halflife,
    ratio = mu / Halflife
  ) %T>%
  print()

# 3.2.7 Contrast estimates ####
halflife_contrast_summary <- halflife_contrast %>%
  mutate(log_ratio = log(ratio)) %>%
  group_by(Contrast) %>%
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
    P = mean( diff > 0 ) %>% signif(2)
  ) %>%
  ungroup() %>%
  mutate(
    diff = glue(
      "{signif(diff_mean, 2)} ± {signif(diff_sd, 2)} ({signif(diff_median, 2)})"
    ),
    ratio = glue(
      "{signif(ratio_median, 2)} ({signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 2)})"
    )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Contrast with mu
mu_vs_halflife_summary <- mu_vs_halflife %>%
  mutate(log_ratio = log(ratio)) %>%
  group_by(Group) %>%
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
    P = mean( diff > 0 ) %>% signif(2),
    P_max = pmax( mean( diff < 0 ) , mean( diff > 0 ) ) %>% 
      signif(2)
  ) %>%
  ungroup() %>%
  mutate(
    diff = glue(
      "{signif(diff_mean, 2)} ± {signif(diff_sd, 2)} ({signif(diff_median, 2)})"
    ),
    ratio = glue(
      "{signif(ratio_median, 2)} ({signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 2)})"
    )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()
# Difference is incredible skewed. The only way to represent this is
# as a log ratio.

# 4. Figures ####
# 4.1 Figure 2 ####
# 4.1.1 Figure 2a ####
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
                 panel.spacing.x = unit(0.8, "cm"),
                 panel.spacing.y = unit(0.6, "cm"),
                 text = element_text(family = "Futura"))

require(ggh4x)
Fig_2a <- meta_prediction_group_effect %>%
  filter(Group != "Prior") %>%
  ggplot() +
    geom_hline(yintercept = 1) +
    geom_point(data = meta %>% # Filtering looks better than clipping
                 filter(
                   Ratio <= 1.2 & (
                     Group == "Terrestrial" & Day <= 40 |
                       Group == "Freshwater" & Day <= 60 |
                       Group == "Seagrass" & Day <= 100 |
                       Group == "Seaweed" & Day <= 300
                   )
                 ), 
               aes(Day, Ratio, colour = Group),
               shape = 16, alpha = 0.1, size = 1.2) +
    geom_ribbon(aes(Day, ymin = Ratio_mu.lower, ymax = Ratio_mu.upper,
                    alpha = factor(.width), fill = Group)) +
    geom_line(aes(Day, Ratio_mu, colour = Group)) +
    geom_hline(yintercept = 0) +
    scale_colour_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), guide = "none") +
    scale_fill_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_y_continuous(breaks = seq(0, 1.2, 0.4),
                       labels = scales::label_number(accuracy = c(1, rep(0.1, 3)))) +
    facet_nested(Response ~ Group + Light, 
                 switch = "y", scale = "free",
                 nest_line = TRUE, 
                 # strip = strip_nested(
                 #   text_y = element_text(angle = 0, hjust = 0, vjust = 1)
                 # ),
                 labeller = labeller(
                   Light = as_labeller(
                     c("Yes" = "Light",
                       "No" = "Dark")
                   ),
                   Group = as_labeller(
                     c("Terrestrial" = "Terrestrial plants",
                       "Freshwater" = "Freshwater plants",
                       "Seagrass" = "Seagrasses",
                       "Seaweed" = "Seaweeds")
                   )
                 )) +
    facetted_pos_scales(
      x = list(
        Group == "Terrestrial" ~
          scale_x_continuous(limits = c(0, 40),
                             breaks = seq(0, 40, by = 20)),
        Group == "Freshwater" ~
          scale_x_continuous(limits = c(0, 60),
                             breaks = seq(0, 60, by = 30)),
        Group == "Seagrass" ~
          scale_x_continuous(limits = c(0, 100),
                             breaks = seq(0, 100, by = 50)),
        Group == "Seaweed" ~
          scale_x_continuous(limits = c(0, 300),
                             breaks = seq(0, 300, by = 150))
      )
    ) +
    coord_cartesian(ylim = c(0, 1.2), expand = F, clip = "off") +
    labs(x = "Detrital age (days)",
         y = expression("Relative photosynthesis ("*italic(p)*"/"*italic(p)[0]*")")) +
    mytheme +
    theme(axis.line.x = element_blank(),
          plot.margin = margin(0, 0.5, 0.2, 0.2, unit = "cm"))

Fig_2a

# 4.1.2 Figure 2b ####
require(ggridges)
Fig_2b <- meta_prior_posterior_group %>%
  filter(Group != "Prior") %>%
  mutate(Response = Response %>% fct_relevel("Chlorophyll")) %>%
  ggplot() +
    stat_density_ridges(aes(mu, y = Response, fill = Group), 
                        colour = NA, n = 2^10, scale = 0.9, 
                        bandwidth = 4*0.02) +
    geom_text(
      data = expand_grid(
        Group = meta %$% levels(Group) %>% fct(),
        Response = meta %$% levels(Response) %>% fct(),
        Light = meta %$% levels(Light) %>% fct()
      ) %>%
        mutate(
          label = case_when(
            Group == "Terrestrial" & Light == "Yes" &
              Response == "Photosynthesis" ~ "Photosynth.",
            Group == "Terrestrial" & Light == "Yes" &
              Response == "Chlorophyll" ~ "Chlorophyll"
          )
        ),
      aes(x = 10^-7.1, y = Response, label = label),
      family = "Futura", size = 12, size.unit = "pt",
      hjust = 0, vjust = 0
    ) +
    scale_fill_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), guide = "none") +
    scale_x_log10(breaks = 10^seq(0, 4, 2),
                  labels = scales::label_log()) +
    facet_grid(~ Group + Light) +
    labs(x = expression("Photosynthetic half-life ("*italic(μ)*", days)")) +
    coord_cartesian(xlim = 10^c(0, 4), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_2b
# Safely ignore warning, which is due to intentional NAs in geom_text.

# 4.1.3 Figure 2c ####
# Build densities manually
meta_beta_dens <- meta_beta %>%
  group_by(Group, Effect) %>%
  reframe(x = c(-2.3, density(beta, n = 2^10, from = -2.3, to = 2.8, bw = 5.1 * 0.02)$x, 2.8),
          y = c(0, density(beta, n = 2^10, from = -2.3, to = 2.8, bw = 5.1 * 0.02)$y, 0)) %>%
  group_by(Group, Effect) %>% # Standardise area with Riemann sum (avoid manually added x[1]).
  mutate(y = y * 0.42 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

# Add labels to meta_beta
meta_beta %<>%
  left_join(
    meta_beta_summary %>%
      select(Group, Effect, P),
    by = c("Group", "Effect")
  ) %>%
  mutate(label = ( P * 100 ) %>% str_c("%"))

# Make annotation
annotation <- expand_grid(
    Group = levels(meta$Group) %>% fct(),
    Effect = c("l", "c") %>% fct()
  ) %>%
  mutate(
    label = case_when(
      Group == "Terrestrial" &
        Effect == "l" ~ "frac(Light, Dark)",
      Group == "Terrestrial" &
        Effect == "c" ~ "frac(Chlorophyll, Photosynth.)"
    ),
    y = if_else(Effect == "l", 1, 0)
  ) %T>%
  print()

require(geomtextpath)
Fig_2c <-  meta_beta_dens %>%
  filter(!Group %in% c("Prior", "Global")) %>%
  ggplot() +
    geom_polygon(data = . %>% filter(Effect == "l"),
                 aes(exp(x), y + 1, fill = Group)) +
    geom_polygon(data = . %>% filter(Effect == "c"),
                 aes(exp(x), y, fill = Group)) +
    geom_textdensity(data = meta_beta %>% 
                       filter(!Group %in% c("Prior", "Global") &
                                Effect == "l") %>%
                       mutate(hjust = if_else(Group == "Freshwater", 0.65, 0.7)),
                     aes(ratio, after_stat(density) * 0.18 + 1,
                         label = label, colour = Group, hjust = hjust),
                     family = "Futura", size = 3.5, vjust = 0, 
                     n = 2^10, bw = 2.2*0.02, # log10(exp(2.8)) - log10(exp(-2.3))
                     text_only = TRUE) +
    geom_textdensity(data = meta_beta %>% 
                       filter(!Group %in% c("Prior", "Global") &
                                Effect == "c"),
                     aes(ratio, after_stat(density) * 0.18,
                         label = label, colour = Group),
                     family = "Futura", size = 3.5, hjust = 0.65,
                     vjust = 0, n = 2^10, bw = 2.2*0.02,
                     text_only = TRUE) +
    geom_text(data = annotation, 
              aes(x = 10^-2.372, y = y, label = label),
              family = "Futura", size = 12, size.unit = "pt",
              hjust = 0, vjust = 0, parse = TRUE) +
    geom_vline(aes(xintercept = 10^0)) +
    scale_colour_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), 
                        guide = "none") +
    scale_fill_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), 
                      guide = "none") +
    scale_x_log10(limits = 10^c(-1, 1.2),
                  breaks = 10^(-1:1),
                  labels = scales::label_log(),
                  oob = scales::oob_keep) +
    facet_grid(~ Group) +
    labs(x = "Relative photosynthetic half-life (ratio)") +
    coord_cartesian(xlim = 10^c(-1, 1), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_2c

# 4.1.4 Combine panels ####
require(patchwork)
Fig_2 <- ( Fig_2a / Fig_2b / Fig_2c ) +
  plot_layout(heights = c(1, 0.4, 0.5))
Fig_2
# On second thought, Fig_2b doesn't add much

Fig_2 <- ( Fig_2a / Fig_2c ) +
  plot_layout(heights = c(1, 0.475))

Fig_2 %>%
  ggsave(filename = "Fig_2.pdf", path = "Figures",
         device = cairo_pdf, height = 15, width = 20, units = "cm")

# 4.2 Figure 3 ####
# 4.2.1 Figure 3a ####
Fig_3a <- meta_prediction_experiment %>%
  filter(.width == 0.5) %>% # I only want medians, not intervals
  ggplot() +
    geom_line(aes(Day, Ratio_mu, colour = Group,
                  group = Experiment),
              alpha = 0.2) +
    geom_line(data = meta_prediction_group %>%
                filter(Group != "Prior" & .width == 0.5),
              aes(Day, Ratio_mu), colour = "black") +
    # geom_ribbon(data = meta_prediction_group %>%
    #             filter(Group != "Prior" & .width == 0.9),
    #             aes(Day, ymin = Ratio_mu.lower, ymax = Ratio_mu.upper), 
    #             fill = NA, colour = "black") +
    scale_colour_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), 
                        guide = "none") +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = scales::label_number(accuracy = c(1, rep(0.1, 4), 1))
    ) +
    facet_grid(~ Group, scales = "free",
               labeller = labeller(
                 Group = as_labeller(
                   c("Terrestrial" = "Terrestrial plants",
                     "Freshwater" = "Freshwater plants",
                     "Seagrass" = "Seagrasses",
                     "Seaweed" = "Seaweeds")
                 ))) +
    facetted_pos_scales(
      x = list(
        Group == "Terrestrial" ~
          scale_x_continuous(limits = c(0, 28),
                             breaks = seq(0, 28, by = 7)),
        Group == "Freshwater" ~
          scale_x_continuous(limits = c(0, 28),
                             breaks = seq(0, 28, by = 7)),
        Group == "Seagrass" ~
          scale_x_continuous(limits = c(0, 56),
                             breaks = seq(0, 56, by = 14)),
        Group == "Seaweed" ~
          scale_x_continuous(limits = c(0, 120),
                             breaks = seq(0, 120, by = 30))
      )
    ) +
    coord_cartesian(expand = F) +
    labs(x = "Detrital age (days)",
         y = expression("Photosynthesis ("*italic(p)*"/"*italic(p)[0]*")")) +
    mytheme +
    theme(plot.margin = margin(0, 0.5, 0, 0.2, unit = "cm"))

Fig_3a

# 4.2.2 Figure 3b ####
# Build densities manually
meta_experiment_dens <- meta_posterior_experiment %>%
  mutate(
    min = if_else(Group == "Terrestrial", -0.5, 0),
    max = if_else(Group == "Seaweed", 7.4, log(1e3)),
    range = max - min
  ) %>%
  group_by(Group, Experiment) %>%
  reframe(
    x = c(
      min[1], 
      density(
        log(mu), n = 2^10, 
        from = min[1], to = max[1], 
        bw = range[1] * 0.02
      )$x,
      max[1]
    ),
    y = c(
      0, 
      density(
        log(mu), n = 2^10, 
        from = min[1], to = max[1], 
        bw = range[1] * 0.02
      )$y, 
      0
    )
  ) %>%
  group_by(Group, Experiment) %>% # Standardise area with Riemann sum (avoid manually added x[1]).
  mutate(y = y * 0.5 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

meta_group_dens <- meta_prior_posterior_group %>%
  filter(Light == "Yes" & Response == "Chlorophyll") %>%
  select(starts_with("."), Group, alpha) %>% # alpha is just mean log(mu)
  group_by(Group) %>%
  reframe(
    x = c(
      0, 
      density(
        alpha, n = 2^10, 
        from = 0, to = log(1e3), 
        bw = log(1e3) * 0.02
      )$x,
      log(1e3)
    ),
    y = c(
      0, 
      density(
        alpha, n = 2^10, 
        from = 0, to = log(1e3), 
        bw = log(1e3) * 0.02
      )$y, 
      0
    )
  ) %>%
  group_by(Group) %>% # Standardise area with Riemann sum (avoid manually added x[1]).
  mutate(y = y * 0.5 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %T>%
  print()

Fig_3b <- meta_experiment_dens %>%
  ggplot() +
    geom_line(aes(exp(x), y, colour = Group,
                  group = Experiment),
              alpha = 0.2) +
    geom_line(data = meta_group_dens, aes(exp(x), y),
              colour = "black") +
    scale_colour_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), 
                        guide = "none") +
    scale_x_log10(breaks = 10^(0:3),
                  label = scales::label_log()) +
    facet_grid(~ Group) +
    labs(x = expression("Photosynthetic half-life ("*italic(μ)*", days)")) +
    coord_cartesian(xlim = 10^c(0, 3), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_3b
# Safely ignore warning, which is due to intentional NAs in geom_text.

# 4.2.3 Figure 3c ####
Fig_3c <- halflife_prior_posterior_group %>%
  filter(Group != "Prior") %>%
  ggplot() +
    geom_point(data = deco,
               aes(Halflife, 0.5, colour = Group),
               shape = 16, alpha = 0.2, size = 2.4,
               position = position_jitter(height = 0.3)) +
    stat_density_ridges(aes(Halflife, 1, fill = Group),
                        colour = NA, n = 2^10,
                        bandwidth = 6*0.02, scale = 1.2) +
    scale_colour_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), 
                        guide = "none") +
    scale_fill_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), 
                      guide = "none") +
    scale_x_log10(breaks = 10^seq(-1, 5, 2),
                  labels = scales::label_log()) +
    facet_grid(~ Group) +
    labs(x = expression("Detrital half-life ("*italic(t)["½"]*", days)")) +
    coord_cartesian(xlim = 10^c(-1, 5), ylim = c(0, NA), 
                    expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_3c

# 4.2.4 Figure 3d ####
Fig_3d <- mu_vs_halflife %>% 
  left_join( # Add labels to mu_vs_halflife
    mu_vs_halflife_summary %>%
      mutate(label = ( P_max * 100 ) %>% str_c("%"),
             hjust = case_when(
               Group == "Seaweed" ~ 0.8,
               Group == "Terrestrial" ~ 0.2,
               TRUE ~ 0.35
             )) %>%
      select(Group, label, hjust),
    by = "Group"
  ) %>%
  filter(!Group %in% c("Prior", "Global")) %>%
  ggplot() +
    geom_density(aes(ratio, fill = Group), 
                 n = 2^10, bw = 8*0.02, colour = NA) +
    geom_textdensity(aes(ratio, label = label, hjust = hjust,
                         colour = Group),
                     family = "Futura", size = 3.5, text_only = TRUE, 
                     vjust = 0, n = 2^10, bw = 8*0.02) +
    scale_fill_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), 
                      guide = "none") +
    scale_colour_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), 
                        guide = "none") +
    geom_vline(aes(xintercept = 10^0)) +
    scale_x_log10(breaks = 10^seq(-4, 4, 2),
                  labels = scales::label_log(),
                  oob = scales::oob_keep) +
    facet_grid(~ Group) +
    labs(x = expression("Relative physiological influence on decomposition ("*italic(μ)*"/"*italic(t)["½"]*")")) +
    coord_cartesian(xlim = 10^c(-4, 4), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_3d

# 4.2.5 Combine panels ####
Fig_3 <- ( Fig_3a / Fig_3b / Fig_3c / Fig_3d ) +
  plot_layout(heights = c(1, 0.6, 0.4, 0.2))
Fig_3

Fig_3 %>%
  ggsave(filename = "Fig_3.pdf", path = "Figures",
         device = cairo_pdf, height = 18, width = 20, units = "cm")

# 4.3 Figure 4 ####
# This figure is a schematic which is mostly illustrated in Affinity Designer.
# The following plots make up the simple data component of the schematic.

# 4.3.1 Figure 4a ####
meta_prediction_group <- here("Meta-analysis", "RDS", "meta_prediction_group.rds") %>%
  read_rds() %T>%
  print()

Fig_4a <- meta_prediction_group %>%
  filter(Group != "Prior" & .width == 0.5) %>%
  ggplot() +
    geom_line(aes(Day, Ratio, colour = Group)) +
    scale_x_continuous(position = "top", breaks = seq(0, 100, 10)) +
    scale_colour_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400")) +
    labs(x = "Detrital age (days)",
         y = "Photosynthesis") +
    coord_cartesian(xlim = c(0, 100), expand = F) +
    mytheme +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

Fig_4a

Fig_4a_simple <- meta_prediction_group %>%
  filter(Group != "Prior" & .width == 0.5) %>%
  ggplot() +
    geom_slab(aes(x = Day, y = Group %>% fct_rev(), slab_alpha = Ratio, fill = Group), 
              thickness = 1, fill_type = "gradient") +
    annotate("text", x = 100, y = 4.5:1.5,
             label = c("Terrestrial plants", "Freshwater plants", 
                       "Seagrasses", "Seaweeds"),
             hjust = 1, size.unit = "pt", size = 12, family = "Futura") +
    scale_x_continuous(position = "top", breaks = seq(0, 100, 10)) +
    scale_fill_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"),
                      guide = "none") +
    guides(slab_alpha = "none") +
    labs(x = "Detrital age (days)",
         y = "Photosynthesis") +
    coord_cartesian(xlim = c(0, 100), expand = F) +
    mytheme +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

Fig_4a_simple

# 4.3.2 Figure 4b ####
halflife_prior_posterior_group <- here("Meta-analysis", "RDS", "halflife_prior_posterior_group.rds") %>%
  read_rds() %T>%
  print()

decay <- halflife_prior_posterior_group %>%
  mutate(k = log(2) / Halflife) %>%
  summarise(k = median(k), .by = Group) %>%
  expand_grid(Day = seq(0, 100, length.out = 200)) %>%
  mutate(Ratio = exp(-k * Day)) %T>% # Simple exponential decay of m/m0
  print()

Fig_4b <- decay %>%
  filter(Group != "Prior") %>%
  ggplot() +
    geom_line(aes(Day, Ratio, colour = Group)) +
    scale_x_continuous(position = "top", breaks = seq(0, 100, 10)) +
    scale_colour_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400")) +
    labs(x = "Detrital age (days)",
         y = "Detrital mass") +
    coord_cartesian(xlim = c(0, 100), expand = F) +
    mytheme +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

Fig_4b

Fig_4b_simple <- decay %>%
  filter(Group != "Prior") %>%
  ggplot() +
    geom_slab(aes(x = Day, y = Group %>% fct_rev(), slab_alpha = Ratio, fill = Group), 
              thickness = 1, fill_type = "gradient") +
    # annotate("text", x = 100, y = 4.5:1.5,
    #          label = c("Terrestrial plants", "Freshwater plants", 
    #                    "Seagrasses", "Seaweeds"),
    #          hjust = 1, vjust = 0.5, size.unit = "pt", size = 12, family = "Futura") +
    scale_x_continuous(position = "top", breaks = seq(0, 100, 10)) +
    scale_fill_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"),
                      guide = "none") +
    guides(slab_alpha = "none") +
    labs(x = "Detrital age (days)",
         y = "Detrital mass") +
    coord_cartesian(xlim = c(0, 100), expand = F) +
    mytheme +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

Fig_4b_simple

# 4.3.3 Combine panels ####
Fig_4 <- Fig_4a / 
  Fig_4b + theme(legend.position = "none",
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.line.x = element_blank(),
                 axis.title.x = element_blank())
Fig_4

Fig_4_simple <- Fig_4a_simple /
  Fig_4b_simple + theme(legend.position = "none",
                        axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.line.x = element_blank(),
                        axis.title.x = element_blank())
Fig_4_simple

Fig_4_simple %>%
  ggsave(filename = "Fig_4_data.pdf", path = "Figures",
         device = cairo_pdf, height = 8.5, width = 20, units = "cm")


# 4.4 Figure S2 ####
# 4.4.1 Figure S2a ####
# Present contrasts as fractions
fracs <- c(
  "Seaweed_Seagrass" = "frac(Seaweeds, Seagrasses)",
  "Seaweed_Freshwater" = "frac(Seaweeds, Freshwater~plants)",
  "Seaweed_Terrestrial" = "frac(Seaweeds, Terrestrial~plants)",
  "Seagrass_Freshwater" = "frac(Seagrasses, Freshwater~plants)",
  "Seagrass_Terrestrial" = "frac(Seagrasses, Terrestrial~plants)",
  "Freshwater_Terrestrial" = "frac(Freshwater~plants, Terrestrial~plants)"
)

Fig_S2a <- meta_contrast %>%
  ggplot() +
    stat_density_ridges(
      aes(
        ratio, 
        Contrast %>% fct_rev(), 
        fill = case_when(
          after_stat(x) > 10^0 & after_stat(y) %in% c(4,5,6) ~ "Seaweed",
          after_stat(x) < 10^0 & after_stat(y) == 6 |
            after_stat(x) > 10^0 & after_stat(y) %in% c(2,3) ~ "Seagrass",
          after_stat(x) < 10^0 & after_stat(y) %in% c(3,5) |
            after_stat(x) > 10^0 & after_stat(y) == 1 ~ "Freshwater",
          after_stat(x) < 10^0 & after_stat(y) %in% c(1,2,4) ~ "Terrestrial",
        ) %>%
          fct_relevel("Terrestrial", "Freshwater", "Seagrass")
      ),
      geom = "density_ridges_gradient", colour = NA, n = 2^10,
      bandwidth = 4*0.02, scale = 0.9, rel_min_height = 0.002
    ) +
    geom_vline(xintercept = 10^0) +
    scale_fill_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"),
                      labels = c(
                        "Terrestrial" = "Terrestrial plants",
                        "Freshwater" = "Freshwater plants",
                        "Seagrass" = "Seagrasses",
                        "Seaweed" = "Seaweeds"
                      )) +
    scale_x_log10(breaks = 10^(-2:2),
                  labels = scales::label_log()) +
    scale_y_discrete(labels = function(x) scales::label_parse()(fracs[x])) +
    labs(x = "Relative photosynthetic half-life (ratio)") +
    coord_cartesian(xlim = 10^c(-2, 2), ylim = c(1, 7), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12, vjust = 0, hjust = 0.5),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

Fig_S2a

# 4.4.2 Figure S2b ####
Fig_S2b <- halflife_contrast %>%
  ggplot() +
    stat_density_ridges(
      aes(
        ratio, 
        Contrast %>% fct_rev(), 
        fill = case_when(
          after_stat(x) > 10^0 & after_stat(y) %in% c(4,5,6) ~ "Seaweed",
          after_stat(x) < 10^0 & after_stat(y) == 6 |
            after_stat(x) > 10^0 & after_stat(y) %in% c(2,3) ~ "Seagrass",
          after_stat(x) < 10^0 & after_stat(y) %in% c(3,5) |
            after_stat(x) > 10^0 & after_stat(y) == 1 ~ "Freshwater",
          after_stat(x) < 10^0 & after_stat(y) %in% c(1,2,4) ~ "Terrestrial",
        ) %>%
          fct_relevel("Terrestrial", "Freshwater", "Seagrass")
      ),
      geom = "density_ridges_gradient", colour = NA, n = 2^10,
      bandwidth = 8*0.02, scale = 0.9, rel_min_height = 0.002
    ) +
    geom_vline(xintercept = 10^0) +
    scale_fill_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"),
                      labels = c(
                        "Terrestrial" = "Terrestrial plants",
                        "Freshwater" = "Freshwater plants",
                        "Seagrass" = "Seagrasses",
                        "Seaweed" = "Seaweeds"
                      )) +
    scale_x_log10(breaks = 10^seq(-4, 4, 2),
                  labels = scales::label_log()) +
    scale_y_discrete(labels = function(x) scales::label_parse()(fracs[x])) +
    labs(x = "Relative detrital half-life (ratio)") +
    coord_cartesian(xlim = 10^c(-4, 4), ylim = c(1, 7), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12, vjust = 0, hjust = 0.5),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

Fig_S2b

# 4.4.3 Combine panels ####
Fig_S2 <- ( Fig_S2a | plot_spacer() | Fig_S2b + 
              theme(axis.text.y = element_blank(),
                    legend.position = "none") ) +
  plot_layout(widths = c(1, 0.1, 1))
Fig_S2

Fig_S2 %>%
  ggsave(filename = "Fig_S2.pdf", path = "Figures",
         device = cairo_pdf, height = 12, width = 20, units = "cm")

# 4.5 Graphical abstract ####
Fig_0a <- mu_vs_halflife %>%
  filter(!Group %in% c("Prior", "Global")) %>%
  ggplot() +
    stat_gradientinterval(aes(mu, Group %>% fct_rev(), fill = Group),
                          n = 2^10, show_point = FALSE,
                          show_interval = FALSE) +
    geom_crossbar(
      data = . %>%
        summarise(median = median(mu), .by = Group),
      aes(median, Group %>% fct_rev()), xmin = NA, xmax = NA,
      colour = "black", middle.linewidth = 0.5
    ) +
    scale_fill_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), 
                      guide = "none") +
    scale_x_log10(breaks = c(1, 7, 30, 365),
                  labels = c("Day", "Week", "Month", "Year")) +
    scale_y_discrete(
      labels = c(
      "Terrestrial" = "Terrestrial plants",
      "Freshwater" = "Freshwater plants",
      "Seagrass" = "Seagrasses",
      "Seaweed" = "Seaweeds"
      )
    ) +
    coord_cartesian(xlim = c(1, 365), clip = "off") +
    mytheme +
    theme(axis.text = element_text(size = 12),
          axis.text.y = element_text(hjust = 0, margin = margin(r = -30)),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.text = element_blank(),
          plot.margin = margin(0, 0.5, 0, 0, unit = "cm"))
    
Fig_0a

Fig_0b <- mu_vs_halflife %>%
  filter(!Group %in% c("Prior", "Global")) %>%
  ggplot() +
    stat_gradientinterval(aes(ratio, Group %>% fct_rev(), fill = Group),
                          n = 2^10, show_point = FALSE,
                          show_interval = FALSE) +
    geom_crossbar(
      data = . %>%
        summarise(median = median(ratio), .by = Group),
      aes(median, Group %>% fct_rev()), xmin = NA, xmax = NA,
      colour = "black", middle.linewidth = 0.5
    ) +
    # geom_vline(xintercept = 10^0) +
    scale_fill_manual(values = c("#004237", "#81a512", "#00a88f", "#a29400"), 
                      guide = "none") +
    # ratio = 1 means that detrital and photosynthetic half-lives are equally long.
    # I set qualitative breaks "scarcely" and "partly" for ratio < 1 and "fully" for
    # ratio > 1. Breaks are adjusted for readability.
    scale_x_log10(breaks = 10^c(-1.8, -0.5, 0.5),
                  labels = c("Scarcely", "Partly", "Fully")) +
    coord_cartesian(xlim = 10^c(-2, 2), clip = "off") +
    mytheme +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.text = element_blank(),
          plot.margin = margin(0, 0.5, 0, 0, unit = "cm"))

Fig_0b

Fig_0 <- ( Fig_0a | plot_spacer() | Fig_0b ) +
  plot_layout(widths = c(1, 0.15, 1)) +
  plot_annotation(
    title = "How long can detritus photosynthesise and how affected is its decomposition?",
  ) &
  theme(
    plot.title = element_text(family = "Futura", size = 12, face = "bold",
                              margin = margin(b = -0.2, unit = "cm")),
    plot.margin = margin(0.2, 0.3, 0.07, 0.07, unit = "cm")
  )
Fig_0

Fig_0 %>%
  ggsave(filename = "Fig_0.pdf", path = "Figures",
         device = cairo_pdf, height = 6, width = 18, units = "cm")

# 5. Tables ####
# 5.1 Table 2 ####
Table_2 <- meta_mu_summary %>%
  full_join(
    meta_beta_summary %>%
      pivot_wider(names_from = Effect,
                  values_from = c(ratio, P))
  ) %>%
  full_join(
    halflife_estimates %>%
      select(-mu)
  ) %>%
  full_join(
    mu_vs_halflife_summary %>%
      select(-c(P_max, diff))
  ) %>%
  select(-n) %>%
  mutate(across(-Group, as.character)) %>%
  pivot_longer(cols = -Group,
               names_to = "Parameter") %>%
  pivot_wider(names_from = Group,
              values_from = value) %T>%
  print()

Table_2 %>%
  write_csv(here("Tables", "Table_2.csv"))

read_docx() %>%
  body_add_table(value = Table_2) %>%
  print(target = here("Tables", "Table_2.docx"))

# 5.2 Table 3 ####
Table_3 <- meta_contrast_summary %>%
  select(-diff) %>%
  rename(ratio_mu = ratio, P_mu = P) %>%
  full_join(
    halflife_contrast_summary %>%
      select(-diff) %>%
      rename(ratio_t0.5 = ratio, P_t0.5 = P)
  ) %>%
  select(Contrast, ratio_mu, P_mu, ratio_t0.5, P_t0.5) %T>%
  print()

Table_3 %>%
  write_csv(here("Tables", "Table_3.csv"))

read_docx() %>%
  body_add_table(value = Table_3) %>%
  print(target = here("Tables", "Table_3.docx"))

# 5.3 Table S2 ####
Table_S2 <- meta_mu_effect_group %>% 
  bind_rows(
    meta_mu_effect_global %>%
      mutate(Group = "Global" %>% fct(),
             Light = Light %>% fct_recode(
               Yes = "Light", No = "Dark"
             ),
             Response = Response %>% fct_recode(
               Photosynthesis = "P", Chlorophyll = "Chl"
             ))
  ) %>%
  select(-n) %>%
  mutate(Light = Light %>% fct_recode(
    Light = "Yes", Dark = "No"
  )) %>%
  arrange(Group, Light, Response) %T>%
  print()

Table_S2 %>%
  write_csv(here("Tables", "Table_S2.csv"))

read_docx() %>%
  body_add_table(value = Table_S2) %>%
  print(target = here("Tables", "Table_S2.docx"))

# 5.4 Table S3 ####
Table_S3 <- meta_mu_effect_species %>%
  mutate(Light = Light %>% fct_recode(
    Light = "Yes", Dark = "No"
  )) %>%
  arrange(desc(Group), as.character(Species), Light, Response) %T>%
  print(n = 200)

Table_S3 %>%
  write_csv(here("Tables", "Table_S3.csv"))

read_docx() %>%
  body_add_table(value = Table_S3) %>%
  print(target = here("Tables", "Table_S3.docx"))

# 5.5 Experiments ####
meta_experiments <- meta_mu_effect_experiment %>%
  mutate(Light = Light %>% fct_recode(
    Light = "Yes", Dark = "No"
  )) %>%
  arrange(desc(Group), as.character(Species), 
          Light, Response) %T>%
  print(n = 551)

meta_experiments %>%
  write_csv(here("Tables", "meta_experiments.csv"))

read_docx() %>%
  body_add_table(value = meta_experiments) %>%
  print(target = here("Tables", "meta_experiments.docx"))