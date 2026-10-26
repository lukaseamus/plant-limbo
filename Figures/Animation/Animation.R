# 1. Load data ####
require(tidyverse)
require(magrittr)
require(here)
require(extraDistr) # R doesn't have a built-in truncated normal

set.seed(100) # Ensure reproducibility of simulation (see Meta.R)
meta <- here("Meta-analysis", "Meta.csv") %>% 
  read_csv(
    col_types = list(
      "f", "c", "f", "f", "f", "f",
      "f", "f", "f", "f", "d", "d",
      "d", "d", "f", "f", "c", "c"
    )
  ) %>%
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
  unnest(Observation) %>%
  group_by(Reference, Series) %>%
  mutate(
    Ratio = if_else( 
      is.na(SEM) & N == 1,
      Observation / mean( Observation[ Day == min(Day) ] ),
      Observation / first( Mean[ Day == min(Day) ] )
    ),
    Experiment = cur_group_id() %>% 
      as_factor()
  ) %>%
  ungroup() %>%
  mutate(
    Group = Group %>%
      fct_relevel("Terrestrial", "Freshwater", "Seagrass"),
    Light = Light %>% fct_relevel("Yes"),
    Water = Water %>% fct_relevel("Yes"),
    Response = Response %>% fct_relevel("Photosynthesis"),
    Ratio = if_else(
      Ratio <= 0, 
      min( Ratio[Ratio > 0] ), 
      Ratio 
    ),
    Light_effect = if_else(
      Light == "Yes", 0.5, -0.5 
    ),
    Chl_effect = if_else(
      Response == "Chlorophyll", 0.5, -0.5
    )
  ) %T>%
  print()

meta_prior_posterior <- here("Meta-analysis", "RDS", "meta_prior_posterior_group.rds") %>%
  read_rds() %T>%
  print()

halflife_prior_posterior <- here("Meta-analysis", "RDS", "halflife_prior_posterior_group.rds") %>%
  read_rds() %T>%
  print()

# Predict across full predictor range
require(tidybayes)
meta_prediction <- meta_prior_posterior %>%
  # alpha is identical for all levels of Light and Response
  filter(Group == "Prior" |
           Light == "Yes" & Response == "Chlorophyll") %>%
  select(starts_with("."), Group, alpha, sigma) %>%
  mutate(mu = exp(alpha)) %>% # the group mean of mu is just exp(alpha)
  expand_grid(
    Day = exp(
      seq( log(1/24/60/60) , log(36500) , length.out = 200 )
    )
  ) %>%
  mutate(
    Ratio_mu = plogis( -5 / mu * ( Day - mu ) ),
    Ratio = rlnorm( n() , log(Ratio_mu) , sigma )
  ) %>%
  # below is a longer but much much (!) more efficient version of 
  # group_by(Day, Treatment) %>% 
  # tidybayes::median_qi(Ratio_mu, Ratio, .width = c(.5, .8, .9))
  summarise(
    across(
      c(Ratio_mu, Ratio),
      list(median = median, 
           lower_0.5 = ~ qi(.x, .width = .5)[1],
           upper_0.5 = ~ qi(.x, .width = .5)[2],
           lower_0.8 = ~ qi(.x, .width = .8)[1],
           upper_0.8 = ~ qi(.x, .width = .8)[2],
           lower_0.9 = ~ qi(.x, .width = .9)[1],
           upper_0.9 = ~ qi(.x, .width = .9)[2]),
      .names = "{.col}.{.fn}"
    ),
    .by = c(Day, Group)
  ) %>%
  rename(Ratio_mu = Ratio_mu.median, Ratio = Ratio.median) %>%
  pivot_longer(cols = contains("lower") | contains("upper")) %>%
  separate(col = name, into = c("name", ".width"), sep = "_(?=[^_]*$)") %>%
  pivot_wider(names_from = name, values_from = value) %T>%
  print()

# 2. Aesthetics ####
# Data points
meta_points <- meta %>% 
  filter(Day > 0 & Ratio < 1.04) %>% # remove day 0 and ratios outside plot
  mutate(
    alpha = .1, # Low baseline transparency because there are many points
    colour = case_when( # Colours
      Group == "Terrestrial" ~ "#004237",
      Group == "Freshwater" ~ "#81a512",
      Group == "Seagrass" ~ "#00a88f",
      Group == "Seaweed" ~ "#a29400"
    )
  ) %>%
  rename(x = Day, y = Ratio) %>% # Rename to x and y for consistency
  select(Group, x, y, alpha, colour) %T>% # Remove irrelevant variables
  print()

# Lines
meta_line <- meta_prediction %>%
  select(Group, Day, Ratio_mu) %>%
  rename(x = Day, y = Ratio_mu) %>%
  mutate(
    colour = case_when(
      Group == "Prior" ~ "#b5b8ba", # here I also have Prior
      Group == "Terrestrial" ~ "#004237",
      Group == "Freshwater" ~ "#81a512",
      Group == "Seagrass" ~ "#00a88f",
      Group == "Seaweed" ~ "#a29400"
    )
  ) %T>%
  print()

# Ribbons
meta_ribbon <- meta_prediction %>%
  group_by(Group, .width) %>%
  reframe(
    x = c(Day, rev(Day)),
    y = c(Ratio_mu.upper, rev(Ratio_mu.lower))
  ) %>%
  mutate(
    alpha = case_when(
      .width == .9 ~ .5, 
      .width == .8 ~ .4, 
      .width == .5 ~ .3
    ),
    fill = case_when(
      Group == "Prior" ~ "#b5b8ba",
      Group == "Terrestrial" ~ "#004237",
      Group == "Freshwater" ~ "#81a512",
      Group == "Seagrass" ~ "#00a88f",
      Group == "Seaweed" ~ "#a29400"
    )
  ) %>%
  select(-.width) %T>%
  print()

# Densities
mu_dens <- meta_prior_posterior %>%
  group_by(Group) %>% 
  reframe( 
    # Close polygons at zero and max with c(0, density, max) for x and
    # at zero and zero with c(0, density, 0) for y.
    x = c(
      log(1/24/60/60), 
      density(
        alpha, n = 2^10, 
        from = log(1/24/60/60), to = log(36500), 
        bw = diff( c( log(1/24/60/60), log(36500) ) ) * 0.01
      )$x, 
      log(36500)
    ),
    y = c(
      0, 
      density(
        alpha, n = 2^10, 
        from = log(1/24/60/60), to = log(36500), 
        bw = diff( c( log(1/24/60/60), log(36500) ) ) * 0.01
      )$y, 
      0
    )
  ) %>%
  group_by(Group) %>% # Standardise area with Riemann sum (avoid manually added x[1]).
  mutate( y = y / ( sum(y) * ( x[3] - x[2] ) ) ) %>%
  ungroup() %>%
  mutate(
    fill = case_when(
      Group == "Prior" ~ "#b5b8ba",
      Group == "Terrestrial" ~ "#004237",
      Group == "Freshwater" ~ "#81a512",
      Group == "Seagrass" ~ "#00a88f",
      Group == "Seaweed" ~ "#a29400"
    )
  ) %T>%
  print()

halflife_dens <- halflife_prior_posterior %>%
  group_by(Group) %>% 
  reframe(
    x = c(
      log(1/24/60/60), 
      density(
        log(Halflife), n = 2^10, 
        from = log(1/24/60/60), to = log(36500), 
        bw = diff( c( log(1/24/60/60), log(36500) ) ) * 0.01
      )$x, 
      log(36500)
    ),
    y = c(
      0, 
      density(
        log(Halflife), n = 2^10, 
        from = log(1/24/60/60), to = log(36500), 
        bw = diff( c( log(1/24/60/60), log(36500) ) ) * 0.01
      )$y, 
      0
    )
  ) %>%
  group_by(Group) %>%  # note I use same multiplier (3) as for mu_dens 
  mutate( y = y / ( sum(y) * ( x[3] - x[2] ) ) ) %>%
  ungroup() %>%
  mutate(
    fill = case_when(
      Group == "Prior" ~ "#b5b8ba",
      Group == "Terrestrial" ~ "#004237",
      Group == "Freshwater" ~ "#81a512",
      Group == "Seagrass" ~ "#00a88f",
      Group == "Seaweed" ~ "#a29400"
    )
  ) %T>%
  print()

# Text
text <- meta_prior_posterior %>%
  distinct(Group) %>%
  mutate(
    colour = case_when(
      Group == "Prior" ~ "#b5b8ba",
      Group == "Terrestrial" ~ "#004237",
      Group == "Freshwater" ~ "#81a512",
      Group == "Seagrass" ~ "#00a88f",
      Group == "Seaweed" ~ "#a29400"
    ),
    label = case_when(
      Group == "Prior" ~ "Prior",
      Group == "Terrestrial" ~ "Terrestrial plants",
      Group == "Freshwater" ~ "Freshwater plants",
      Group == "Seagrass" ~ "Seagrasses",
      Group == "Seaweed" ~ "Seaweeds"
    ),
    x = 1.36e-5, y = 0.05
  ) %T>%
  print()

# Illustrations
image <- meta %>%
  distinct(Group) %>%
  mutate(
    image = case_when(
      Group == "Terrestrial" ~ here("Figures", "Illustrations", "Quercus.pdf"),
      Group == "Freshwater" ~ here("Figures", "Illustrations", "Elodea.pdf"),
      Group == "Seagrass" ~ here("Figures", "Illustrations", "Enhalus.pdf"),
      Group == "Seaweed" ~ here("Figures", "Illustrations", "Saccharina.pdf")
    ),
    alpha = 1, # baseline opacity
    x = 1e4, y = .66
  ) %T>%
  print()

# 3. Static figures ####
# Define custom theme
mytheme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  plot.margin = margin(0.2, 0.7, 0.2, 0.2, unit = "cm"),
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
  panel.spacing = unit(1, "cm"),
  text = element_text(family = "Futura")
)

# Top panel
require(ggimage)
ggplot() +
  geom_point(
    data = meta_points,
    aes(x, y, colour = colour, alpha = alpha),
    shape = 16, size = 2.35 # this size matches default size of point range
  ) +
  geom_line(
    data = meta_line,
    aes(x, y, colour = colour)
  ) +
  geom_polygon(
    data = meta_ribbon,
    aes(x, y, fill = fill, alpha = alpha, group = alpha)
  ) +
  geom_text(
    data = text,
    aes(x, y, label = label, colour = colour),
    hjust = 0, vjust = 0, size.unit = "pt", size = 12, family = "Futura"
  ) +
  geom_image(
    data = image,
    aes(x, y, image = image, alpha = alpha),
    size = .4, by = "height"
  ) +
  scale_alpha_identity() +
  scale_colour_identity() +
  scale_fill_identity() +
  scale_x_log10(
    breaks = c(1/24/60/60, 1/24/60, 1/24, 1, 7, 30, 365, 3650, 36500),
    label = c("Second", "Minute", "Hour", "Day", "Week", "Month", 
              "Year", "Decade", "Century")
  ) +
  scale_y_continuous(labels = scales::label_number(accuracy = c(1, 1e-2, 0.1, 1e-2, 1))) +
  facet_grid(~ Group, space = "free", scales = "free") +
  labs(x = "Detrital age", 
       y = expression("Photosynthesis ("*italic(p)*"/"*italic(p)[0]*")")) +
  coord_cartesian(xlim = c(1/24/60/60, 36500), ylim = c(0, 1), expand = F, clip = "off") +
  mytheme

# Bottom panel v1
ggplot() +
  geom_polygon(
    data = mu_dens,
    aes(exp(x), y, fill = fill)
  ) +
  scale_fill_identity() +
  scale_x_log10(
    breaks = c(1/24/60/60, 1/24/60, 1/24, 1, 7, 30, 365, 3650, 36500),
    label = c("Second", "Minute", "Hour", "Day", "Week", "Month", 
              "Year", "Decade", "Century")
  ) +
  # Ensure same spacing as top panel
  scale_y_continuous(labels = scales::label_number(accuracy = 1e-2)) +
  facet_grid(~ Group, space = "free", scales = "free") +
  labs(x = "Half-life",
       y = expression("Photosynthesis ("*italic(p)*"/"*italic(p)[0]*")")) +
  coord_cartesian(
    xlim = c(1/24/60/60, 36500), ylim = c(0, .7), expand = F, clip = "off"
  ) +
  mytheme +
  # These will be hidden (transparent) but ensure identical spacing
  theme(axis.title.y = element_text(colour = "transparent"),
        axis.text.y = element_text(colour = "transparent"),
        axis.ticks.y = element_line(colour = "transparent"),
        axis.line.y = element_line(colour = "transparent"))

# Bottom panel v2
ggplot() +
  geom_polygon(
    data = mu_dens,
    aes(exp(x), y, fill = fill), alpha = 0.8
  ) +
  geom_polygon(
    data = halflife_dens,
    aes(exp(x), y, fill = fill), alpha = 0.8
  ) +
  scale_fill_identity() +
  scale_x_log10(
    breaks = c(1/24/60/60, 1/24/60, 1/24, 1, 7, 30, 365, 3650, 36500),
    label = c("Second", "Minute", "Hour", "Day", "Week", "Month", 
              "Year", "Decade", "Century")
  ) +
  # Ensure same spacing as top panel
  scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
  facet_grid(~ Group, space = "free", scales = "free") +
  labs(x = "Half-life",
       y = expression("Photosynthesis ("*italic(p)*"/"*italic(p)[0]*")")) +
  coord_cartesian(
    xlim = c(1/24/60/60, 36500), ylim = c(0, .7), expand = F, clip = "off"
  ) +
  mytheme +
  # These will be hidden (transparent) but ensure identical spacing
  theme(axis.title.y = element_text(colour = "transparent"),
        axis.text.y = element_text(colour = "transparent"),
        axis.ticks.y = element_line(colour = "transparent"),
        axis.line.y = element_line(colour = "transparent"))

# 4. Animation ####
# Define enter/exit functions
enter <- function(data){
  data$alpha <- 0
  data
}

exit <- function(data){
  data$alpha <- 0
  data
}

# Tween points
require(tweenr)
meta_points_ani <- bind_rows( 
  tween_state( # Points are arbitrarily paired
    meta_points %>% filter(Group == "Prior"), 
    meta_points %>% filter(Group == "Seaweed"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50),
  tween_state(
    meta_points %>% filter(Group == "Seaweed"), 
    meta_points %>% filter(Group == "Seagrass"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_state(
    meta_points %>% filter(Group == "Seagrass"), 
    meta_points %>% filter(Group == "Freshwater"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_state(
    meta_points %>% filter(Group == "Freshwater"), 
    meta_points %>% filter(Group == "Terrestrial"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_state(
    meta_points %>% filter(Group == "Terrestrial"), 
    meta_points %>% filter(Group == "Prior"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# Tween lines
# tween_state also works for some lines, but 
# I believe it is more stable to use tween_path.
require(transformr)
meta_line_ani <- bind_rows(
  tween_path(
    meta_line %>% filter(Group == "Prior"), 
    meta_line %>% filter(Group == "Seaweed"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_path(
    meta_line %>% filter(Group == "Seaweed"), 
    meta_line %>% filter(Group == "Seagrass"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_path(
    meta_line %>% filter(Group == "Seagrass"), 
    meta_line %>% filter(Group == "Freshwater"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_path(
    meta_line %>% filter(Group == "Freshwater"), 
    meta_line %>% filter(Group == "Terrestrial"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_path(
    meta_line %>% filter(Group == "Terrestrial"), 
    meta_line %>% filter(Group == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# Tween ribbons
# tween_path is more stable than tween_polygon for ribbons.
meta_ribbon_ani <- bind_rows(
  tween_path(
    meta_ribbon %>% filter(Group == "Prior"),
    meta_ribbon %>% filter(Group == "Seaweed"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_path(
    meta_ribbon %>% filter(Group == "Seaweed"),
    meta_ribbon %>% filter(Group == "Seagrass"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_path(
    meta_ribbon %>% filter(Group == "Seagrass"),
    meta_ribbon %>% filter(Group == "Freshwater"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_path(
    meta_ribbon %>% filter(Group == "Freshwater"),
    meta_ribbon %>% filter(Group == "Terrestrial"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_path(
    meta_ribbon %>% filter(Group == "Terrestrial"),
    meta_ribbon %>% filter(Group == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# Tween densities
mu_dens_ani <- bind_rows(
  tween_polygon(
    mu_dens %>% filter(Group == "Prior"),
    mu_dens %>% filter(Group == "Seaweed"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_polygon(
    mu_dens %>% filter(Group == "Seaweed"),
    mu_dens %>% filter(Group == "Seagrass"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_polygon(
    mu_dens %>% filter(Group == "Seagrass"),
    mu_dens %>% filter(Group == "Freshwater"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_polygon(
    mu_dens %>% filter(Group == "Freshwater"),
    mu_dens %>% filter(Group == "Terrestrial"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_polygon(
    mu_dens %>% filter(Group == "Terrestrial"),
    mu_dens %>% filter(Group == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

halflife_dens_ani <- bind_rows(
  tween_polygon(
    halflife_dens %>% filter(Group == "Prior"),
    halflife_dens %>% filter(Group == "Seaweed"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_polygon(
    halflife_dens %>% filter(Group == "Seaweed"),
    halflife_dens %>% filter(Group == "Seagrass"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_polygon(
    halflife_dens %>% filter(Group == "Seagrass"),
    halflife_dens %>% filter(Group == "Freshwater"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_polygon(
    halflife_dens %>% filter(Group == "Freshwater"),
    halflife_dens %>% filter(Group == "Terrestrial"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_polygon(
    halflife_dens %>% filter(Group == "Terrestrial"),
    halflife_dens %>% filter(Group == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# Tween text
text_ani <- bind_rows( 
  tween_state(
    text %>% filter(Group == "Prior"), 
    text %>% filter(Group == "Seaweed"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_state(
    text %>% filter(Group == "Seaweed"), 
    text %>% filter(Group == "Seagrass"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_state(
    text %>% filter(Group == "Seagrass"), 
    text %>% filter(Group == "Freshwater"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_state(
    text %>% filter(Group == "Freshwater"), 
    text %>% filter(Group == "Terrestrial"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_state(
    text %>% filter(Group == "Terrestrial"), 
    text %>% filter(Group == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# Tween illustrations
image_ani <- bind_rows( 
  tween_state(
    image %>% filter(Group == "Prior"), 
    image %>% filter(Group == "Seaweed"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50),
  tween_state(
    image %>% filter(Group == "Seaweed"), 
    image %>% filter(Group == "Seagrass"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_state(
    image %>% filter(Group == "Seagrass"), 
    image %>% filter(Group == "Freshwater"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_state(
    image %>% filter(Group == "Freshwater"), 
    image %>% filter(Group == "Terrestrial"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_state(
    image %>% filter(Group == "Terrestrial"), 
    image %>% filter(Group == "Prior"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# 5. Animated figures ####
require(gganimate) # also install gifski
require(ggtext) # for html syntax in gif

# Top panel
( ggplot() +
    geom_point(
      data = meta_points_ani,
      aes(x, y, colour = colour, alpha = alpha),
      shape = 16, size = 2.35
    ) +
    geom_line(
      data = meta_line_ani,
      aes(x, y, colour = colour)
    ) +
    geom_polygon(
      data = meta_ribbon_ani,
      aes(x, y, fill = fill, alpha = alpha, group = alpha)
    ) +
    geom_text(
      data = text_ani,
      aes(x, y, label = label, colour = colour),
      hjust = 0, vjust = 0, size.unit = "pt", size = 12, family = "Futura"
    ) +
    geom_image(
      data = image_ani,
      aes(x, y, image = image, alpha = alpha), 
      size = .6, by = "height"
    ) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_x_log10(
      breaks = c(1/24/60/60, 1/24/60, 1/24, 1, 7, 30, 365, 3650, 36500),
      label = c("Second", "Minute", "Hour", "Day", "Week", "Month", 
                "Year", "Decade", "Century")
    ) +
    scale_y_continuous(
      labels = scales::label_number(accuracy = c(1, 1e-2, 0.1, 1e-2, 1))
    ) +
    labs(
      x = "Detrital age", 
      y = "Photosynthesis (*p*/*p*<sub><span style='font-size:8.4pt'>0</span></sub>)"
    ) +
    coord_cartesian(
      xlim = c(1/24/60/60, 36500), ylim = c(0, 1), 
      expand = F, clip = "off"
    ) +
    transition_manual(.frame) +
    mytheme +
    theme(
      axis.title.y = element_markdown() # only specify for y axis
    ) ) %>%
  animate(
    fps = 25, # optimal gif frame rate (4 centiseconds per frame) with reasonable file size
    duration = 5 * 3, # number of transitions (5) * time per transition (3 s) to set speed
    detail = 2, # 25*5*3 is 1/2 of all frames (5*150), so the skipped frames can add detail
    width = 20, height = 6,
    units = "cm", res = 300,
    renderer = gifski_renderer(),
    device = "ragg_png" # necessary for italic!
  ) %>%
  anim_save(
    filename = "meta_top.gif", 
    path = here("Figures", "Animation")
  )

# Bottom panel v1
( ggplot() +
    geom_polygon(
      data = mu_dens_ani,
      aes(exp(x), y, fill = fill)
    ) +
    scale_fill_identity() +
    scale_x_log10(
      breaks = c(1/24/60/60, 1/24/60, 1/24, 1, 7, 30, 365, 3650, 36500),
      label = c("Second", "Minute", "Hour", "Day", "Week", "Month", 
                "Year", "Decade", "Century")
    ) +
    # Ensure same spacing as top panel
    scale_y_continuous(labels = scales::label_number(accuracy = 1e-2)) +
    labs(
      x = "Half-life",
      y = "Photosynthesis (*p*/*p*<sub><span style='font-size:8.4pt'>0</span></sub>)"
    ) +
    coord_cartesian(
      xlim = c(1/24/60/60, 36500), ylim = c(0, .7), 
      expand = F, clip = "off"
    ) +
    transition_manual(.frame) +
    mytheme +
    theme(
      # These will be hidden (transparent) but ensure identical spacing
      axis.title.y = element_markdown(colour = "transparent"),
      axis.text.y = element_text(colour = "transparent"),
      axis.ticks.y = element_line(colour = "transparent"),
      axis.line.y = element_line(colour = "transparent")
    ) ) %>%
  animate(
    fps = 25, duration = 5 * 3, detail = 2,
    width = 20, height = 4, units = "cm", res = 300, 
    renderer = gifski_renderer(), device = "ragg_png"
  ) %>%
  anim_save(
    filename = "meta_bottom_1.gif", 
    path = here("Figures", "Animation")
  )

# Bottom panel v2
( ggplot() +
    geom_polygon(
      data = mu_dens_ani,
      aes(exp(x), y, fill = fill), alpha = .8
    ) +
    geom_polygon(
      data = halflife_dens_ani,
      aes(exp(x), y, fill = fill), alpha = .8
    ) +
    scale_fill_identity() +
    scale_x_log10(
      breaks = c(1/24/60/60, 1/24/60, 1/24, 1, 7, 30, 365, 3650, 36500),
      label = c("Second", "Minute", "Hour", "Day", "Week", "Month", 
                "Year", "Decade", "Century")
    ) +
    scale_y_continuous(labels = scales::label_number(accuracy = 1e-2)) +
    labs(
      x = "Half-life",
      y = "Photosynthesis (*p*/*p*<sub><span style='font-size:8.4pt'>0</span></sub>)"
    ) +
    coord_cartesian(
      xlim = c(1/24/60/60, 36500), ylim = c(0, .7), 
      expand = F, clip = "off"
    ) +
    transition_manual(.frame) +
    mytheme +
    theme(
      axis.title.y = element_markdown(colour = "transparent"),
      axis.text.y = element_text(colour = "transparent"),
      axis.ticks.y = element_line(colour = "transparent"),
      axis.line.y = element_line(colour = "transparent")
    ) ) %>%
  animate(
    fps = 25, duration = 5 * 3, detail = 2,
    width = 20, height = 4, units = "cm", res = 300, 
    renderer = gifski_renderer(), device = "ragg_png"
  ) %>%
  anim_save(
    filename = "meta_bottom_2.gif", 
    path = here("Figures", "Animation")
  )