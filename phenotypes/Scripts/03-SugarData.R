# Erika M Hansson 2025
# erika.hansson@manchester.ac.uk

#### ---------------------------------------------------- ####
# 1) Plotting sugar data 
# 2) Fitting a model to sugar data to predict release and release amount
# 3) Plots with model fit
# Then: repeating all of the above with the 86 strains in the phylogeny for robustness
#### ---------------------------------------------------- ####

if(!require("tidyverse")) {install.packages("tidyverse")} # wrangling, plotting
if(!require("emmeans")) {install.packages("emmeans")} # for fitting mixed effects models
if(!require("janitor")) {install.packages("janitor")} # for wrangling
if(!require("brms")) {install.packages("brms")} # for multivariate response modeling
if(!require("bayesplot")) {install.packages("bayesplot")} # for plotting Bayesian models

# Load data frame
df_tidy_sugar <- read_csv("Data/df_tidy_sugar.csv")

# Look at it
glimpse(df_tidy_sugar)

# Sugar release for each strain by pH
df_tidy_sugar %>%
  group_by(ph, lifestyle, strain) %>%
  summarise(
    mean = mean(release),
    se = sd(release)/sqrt(n())
  ) %>%
  ggplot(aes(x = factor(ph), y = mean, colour = lifestyle, group = lifestyle)) +
  geom_hline(yintercept = 0) +
  geom_point(data = df_tidy_sugar %>% mutate(mean = release), alpha = 0.3) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
  theme_bw() +
  labs(y = "Mean total sugar release (ppm)",
       colour = "Strain lifestyle",
       x = "pH") +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  facet_wrap(~strain)

### Fitting a model
# The data is both right skewed and has a huge amount of 0s
df_tidy_sugar %>%
  ggplot(aes(x = release, fill = lifestyle)) +
  geom_histogram() +
  facet_wrap(~sugar)

# Need to fit a model to account for
# 1. Do they release anything at all?
# 2. If so, how much do they release

## Fitting Bayesian regression model with brms()

# Pivot the data frame to allow fitting each sugar as a separate response variable
brms_df <- df_tidy_sugar %>%
  pivot_wider(names_from = "sugar", values_from = "release") %>%
  clean_names()

# Look at it
glimpse(brms_df)

# Define the model formula for each sugar
# Random effect excluded from hurdle to improve fit
f_glucose <- bf(glucose ~ ph * lifestyle + (1 | strain), hu ~ ph * lifestyle)
f_maltose <- bf(maltose ~ ph * lifestyle + (1 | strain), hu ~ ph * lifestyle)
f_trehalose <- bf(trehalose ~ ph * lifestyle + (1 | strain), hu ~ ph * lifestyle)

# The model is pretty big so switching to cmdstanr backend
# Run the following hashed lines if it is first time using this backend

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# Then install CmdStan itself
#cmdstanr::install_cmdstan()

library(cmdstanr)

# Setting weakly defined priors to improve performance
priors <- c(
  # Fixed effects for glucose
  set_prior("normal(0, 2)", class = "b", resp = "glucose"),
  set_prior("normal(0, 2)", class = "b", dpar = "hu", resp = "glucose"),
  set_prior("normal(-2, 1.5)", class = "Intercept", resp = "glucose"),
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "hu", resp = "glucose"),
  set_prior("exponential(1)", class = "sd", group = "strain", resp = "glucose"),
  set_prior("exponential(1)", class = "sd", group = "strain", coef = "Intercept", resp = "glucose"),
  set_prior("exponential(1)", class = "shape", resp = "glucose"),
  
  # Repeat for maltose
  set_prior("normal(0, 2)", class = "b", resp = "maltose"),
  set_prior("normal(0, 2)", class = "b", dpar = "hu", resp = "maltose"),
  set_prior("normal(-2, 1.5)", class = "Intercept", resp = "maltose"),
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "hu", resp = "maltose"),
  set_prior("exponential(1)", class = "sd", group = "strain", resp = "maltose"),
  set_prior("exponential(1)", class = "sd", group = "strain", coef = "Intercept", resp = "maltose"),
  set_prior("exponential(1)", class = "shape", resp = "maltose"),
  
  # Repeat for trehalose
  set_prior("normal(0, 2)", class = "b", resp = "trehalose"),
  set_prior("normal(0, 2)", class = "b", dpar = "hu", resp = "trehalose"),
  set_prior("normal(-2, 1.5)", class = "Intercept", resp = "trehalose"),
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "hu", resp = "trehalose"),
  set_prior("exponential(1)", class = "sd", group = "strain", resp = "trehalose"),
  set_prior("exponential(1)", class = "sd", group = "strain", coef = "Intercept", resp = "trehalose"),
  set_prior("exponential(1)", class = "shape", resp = "trehalose")
)

#Why these priors?
#normal(0, 2): Allows moderate fixed effects on the log scale — enough flexibility without letting the model go wild.
#Intercept = normal(-2, 1.5): On log scale, this reflects prior expectation of small release amounts when nothing else influences it.
#logistic(0, 1) for hurdle Intercept: Centers on 0.5 probability of zero, still flexible, aligns with data (heavy zero-inflation).
#exponential(1) for all SD and shape parameters: Regularizes these to be small but doesn’t rule out larger variance.

# Fit the model and run (or load from saved file using readRDS below to save time)
#brms_model <- brm(
#  f_glucose + f_maltose + f_trehalose,
#  data = brms_df,
#  prior = priors,
#  sample_prior = "no",
#  family = hurdle_gamma(link = "log"),
#  chains = 4, cores = 4, iter = 4000,
#  control = list(adapt_delta = 0.95, max_treedepth = 15),
#  backend = "cmdstanr"
#)
# adapt_delta = 0.95 and max_treedepth = 15 are common tweaks for tricky models 
# with strong zero-inflation and correlations.

# saveRDS(brms_model, "Data/brms_model.rds")

# Import the model from file
brms_model <- readRDS("Data/brms_model.rds")
 
# Model summary
summary(brms_model)
 
# Divergent transitions:
 rstan::get_num_divergent(brms_model$fit)
# Maximum tree depth:
 rstan::get_num_max_treedepth(brms_model$fit)
# Both are 0, which is ideal -- it means we're exploring the parameter space well
 
# Do posterior predictive checks, this is looking at the simulated data, the distribution should look like the observed data
pp_check(brms_model, resp = "glucose", ndraws = 200)
pp_check(brms_model, resp = "maltose", ndraws = 200)
pp_check(brms_model, resp = "trehalose", ndraws = 200)
 
## Plotting hurdle and conditional simulated vs observed values in a loop
## pp_check should be able to do this but is not cooperating so using bayesplot instead
 
 # Set number of draws for plotting
 ndraws <- 500
 
 # List of response variable names
 response_vars <- c("glucose", "maltose", "trehalose")
 
 # Loop through each sugar
 for (resp in response_vars) {
   cat("\n==== Plotting for:", resp, "====\n")
   
   # Extract posterior predictive draws
   yrep <- posterior_predict(brms_model, resp = resp, ndraws = ndraws)
   
   # Extract observed values
   y <- brms_model$data[[resp]]
   
   ## === HURDLE COMPONENT === ##
   y_binary <- as.numeric(y == 0)
   yrep_binary <- apply(yrep, 1:2, function(x) as.numeric(x == 0))
   
   print(
     ppc_bars(y = y_binary, yrep = yrep_binary) +
       ggtitle(paste(resp, "- Hurdle (Zero vs Positive)")) +
       theme_minimal()
   )
   
   ## === CONDITIONAL COMPONENT === ##
   pos_idx <- which(y > 0)
   y_pos <- y[pos_idx]
   yrep_pos <- yrep[, pos_idx, drop = FALSE]
   
   if (length(pos_idx) > 5) {
     print(
       ppc_dens_overlay(y = y_pos, yrep = yrep_pos) +
         ggtitle(paste(resp, "- Conditional (Positive Only)")) +
         theme_minimal()
     )
   } else {
     message(paste("Too few positive values for conditional plot of", resp))
   }
 }
 
 # Some overestimation on the low end in conditional model for glucose and maltose, 
 # but otherwise a good fit for both hurdle and conditional
 
 # Residuals vs. fitted values -- these look fine
 pp_check(brms_model, lifestyle = "error_binned", resp = "glucose")
 pp_check(brms_model, lifestyle = "error_binned", resp = "maltose")
 pp_check(brms_model, lifestyle = "error_binned", resp = "trehalose")
 
 # Leave‑One‑Out Cross‑Validation (LOO)
 # Compare predictive accuracy or check for problematic observations.
 loo_fit <- loo(brms_model)
 print(loo_fit)
 plot(loo_fit)
 # These values fine, the bad ones are not frequent enough to worry about
 
 # Bayesian R2 -- how much variation is explained by model
 bayes_R2(brms_model)
 
 # Strong predictive power for glucose and maltose, very low predctive power for trehalose, just as expected
 
 # Conditional effects pull out the fitted lines and points with 95% confidence intervals
 conditional_effects(brms_model, resp = "glucose")
 conditional_effects(brms_model, resp = "maltose")
 conditional_effects(brms_model, resp = "trehalose")
 # These look very sensible
 
# Plotting the observed data vs estimates for publication
# Posterior predictions for the hurdle (probability of > 0)
hurdle_df <- bind_rows(data.frame(conditional_effects(brms_model, dpar = "hu", effects = "lifestyle:ph")[[1]], sugar = "Glucose"),
                        data.frame(conditional_effects(brms_model, dpar = "hu", effects = "lifestyle:ph")[[2]], sugar = "Maltose"),
                        data.frame(conditional_effects(brms_model, dpar = "hu", effects = "lifestyle:ph")[[3]], sugar = "Trehalose")) %>%
   filter(ph != 6.550000) %>%
   dplyr::select(lifestyle, effect2__, estimate__, se__, lower__, upper__, sugar) %>%
   dplyr::rename("ph" = effect2__,
                 "prop_zero" = estimate__,
                 "se" = se__,
                 "lower" = lower__,
                 "upper" = upper__) %>%
  mutate(lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                        lifestyle == "FL" ~ "Free-living"))
 
glimpse(hurdle_df)
 
prop_zero_df <- df_tidy_sugar %>%
   mutate(zero = release == 0,
          lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                                lifestyle == "FL" ~ "Free-living")) %>%
   group_by(lifestyle, ph, sugar) %>%
   summarise(prop_zero = mean(zero), n = n(), .groups = "drop")
 
# Observed vs. estimates for hurdle

 ggplot(prop_zero_df, aes(x = factor(ph), y = 1-prop_zero)) +
   geom_col(aes(fill = lifestyle), width = 0.6, alpha = 0.5) +
   geom_point(data = hurdle_df, aes(colour = lifestyle), size = 2) +
   geom_errorbar(data = hurdle_df, aes(ymin = 1-lower, ymax = 1-upper, 
                 colour = lifestyle), width = 0.2) +
   labs(
     y = "Proportion of non-zero values",
     x = "pH"
   ) +
   theme_bw () +
   theme(legend.position = "none") +
   scale_fill_manual(values =
                       c("Free-living" = "green3",
                         "Symbiotic" = "cornflowerblue")) +
   scale_colour_manual(values =
                         c("Free-living" = "green3",
                           "Symbiotic" = "cornflowerblue")) +
   facet_grid(lifestyle~sugar)
 
# Posterior predictions for the conditional (amount)
cond_df <- bind_rows(data.frame(conditional_effects(brms_model, effects = "lifestyle:ph")[[1]], sugar = "Glucose"), 
                      data.frame(conditional_effects(brms_model, effects = "lifestyle:ph")[[2]], sugar = "Maltose"),
                      data.frame(conditional_effects(brms_model, effects = "lifestyle:ph")[[3]], sugar = "Trehalose")) %>%
   filter(ph != 6.550000) %>%
   dplyr::select(lifestyle, effect2__, estimate__, se__, lower__, upper__, sugar) %>%
   dplyr::rename("ph" = effect2__,
                 "estimate" = estimate__,
                 "se" = se__,
                 "lower" = lower__,
                 "upper" = upper__) %>%
   mutate(ph = as.numeric(as.character(ph)),
          lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                                lifestyle == "FL" ~ "Free-living"))
 
glimpse(cond_df)
 
nonzero_df <- df_tidy_sugar %>%
   filter(release > 0) %>%
  mutate(lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                               lifestyle == "FL" ~ "Free-living"))
 
cond_df %>%
   ggplot(aes(x = factor(ph), y = estimate, colour = lifestyle)) +
   geom_point() +
   geom_jitter(data = nonzero_df %>% mutate(estimate = release), alpha = 0.3) +
   geom_errorbar(aes(ymin = lower, ymax = upper)) +
   facet_wrap(~sugar) +
  scale_colour_manual(values =
                        c("Free-living" = "green3",
                          "Symbiotic" = "cornflowerblue")) +
  theme_bw() +
  labs(y = "Sugar release (ppm)",
       x = "pH",
       colour = "Lifestyle") +
  theme(legend.position = "bottom")
 
# Plot strain level effects
# Strains that don't overlap with 0 meaningfully diverge from the average
ranef_array <- ranef(brms_model)$strain
glucose_ranef <- as.data.frame(ranef_array[, , "glucose_Intercept"]) %>%
  mutate(sugar = "Glucose")
glucose_ranef$strain <- rownames(ranef_array)
maltose_ranef <- as.data.frame(ranef_array[, , "maltose_Intercept"]) %>%
  mutate(sugar = "Maltose")
maltose_ranef$strain <- rownames(ranef_array)
trehalose_ranef <- as.data.frame(ranef_array[, , "trehalose_Intercept"]) %>%
  mutate(sugar = "Trehalose")
trehalose_ranef$strain <- rownames(ranef_array)

ranef_df <- bind_rows(glucose_ranef, maltose_ranef, trehalose_ranef)

ranef_df %>%
  mutate(lifestyle = ifelse(str_detect(strain, "s") == TRUE, "Symbiotic", "Free-living")) %>%
  ggplot(aes(x = strain, y = Estimate, colour = lifestyle)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5)) +
  scale_colour_manual(values =
                        c("Free-living" = "green3",
                          "Symbiotic" = "cornflowerblue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  labs(x = "Strain",
       y = "Strain-level intercept deviation",
       colour = "Lifestyle") +
  facet_wrap(~ sugar, ncol = 1)


### Pairwise contrasts
# Extract slope for conditional and hurdle part of model for each sugar
mu_trend_glucose <- emtrends(brms_model, var = "ph", specs = "lifestyle", resp = "glucose")
hu_trend_glucose <- emtrends(brms_model, specs = "lifestyle", var = "ph", dpar = "hu", resp = "glucose")
mu_trend_maltose <- emtrends(brms_model, var = "ph", specs = "lifestyle", resp = "maltose")
hu_trend_maltose <- emtrends(brms_model, specs = "lifestyle", var = "ph", dpar = "hu", resp = "maltose")
mu_trend_trehalose <- emtrends(brms_model, var = "ph", specs = "lifestyle", resp = "trehalose")
hu_trend_trehalose <- emtrends(brms_model, specs = "lifestyle", var = "ph", dpar = "hu", resp = "trehalose")

# View estimated marginal slopes
summary(mu_trend_glucose)
summary(hu_trend_glucose)
summary(mu_trend_maltose)
summary(hu_trend_maltose)
summary(mu_trend_trehalose)
summary(hu_trend_trehalose)

# Direct contrasts, if the HPD interval overlaps 0 there's no difference
bind_rows(data.frame(pairs(mu_trend_glucose),
           sugar = "glucose",
           model = "mu"),
          data.frame(pairs(hu_trend_glucose),
                     sugar = "glucose",
                     model = "hu"),
          data.frame(pairs(mu_trend_maltose),
                     sugar = "maltose",
                     model = "mu"),
          data.frame(pairs(hu_trend_maltose),
                     sugar = "maltose",
                     model = "hu"),
          data.frame(pairs(mu_trend_trehalose),
                     sugar = "trehalose",
                     model = "mu"),
          data.frame(pairs(hu_trend_trehalose),
                     sugar = "trehalose",
                     model = "hu")) %>%
  mutate(`overlap?` = ifelse(lower.HPD < 0 & upper.HPD > 0, "yes", "no"))

# Reanalysis using only the 86 strains from the phylogeny to check robustness of results ###########################
# Sugar release for each strain by pH

df_tidy_sugar_86 <- df_tidy_sugar %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354")


df_tidy_sugar %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  group_by(ph, lifestyle, strain) %>%
  summarise(
    mean = mean(release),
    se = sd(release)/sqrt(n())
  ) %>%
  ggplot(aes(x = factor(ph), y = mean, colour = lifestyle, group = lifestyle)) +
  geom_hline(yintercept = 0) +
  geom_point(data = df_tidy_sugar %>% mutate(mean = release) %>%  filter(strain != "155",
                                                                         strain != "192",
                                                                         strain != "242",
                                                                         strain != "252",
                                                                         strain != "258",
                                                                         strain != "s354"), alpha = 0.3) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
  theme_bw() +
  labs(y = "Mean total sugar release (ppm)",
       colour = "Strain lifestyle",
       x = "pH") +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  facet_wrap(~strain)

### Fitting a model
# The data is both right skewed and has a huge amount of 0s
df_tidy_sugar %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  ggplot(aes(x = release, fill = lifestyle)) +
  geom_histogram() +
  facet_wrap(~sugar)

# Need to fit a model to account for
# 1. Do they release anything at all?
# 2. If so, how much do they release

## Fitting Bayesian regression model with brms()

# Pivot the data frame to allow fitting each sugar as a separate response variable
brms_df_86 <- df_tidy_sugar %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  pivot_wider(names_from = "sugar", values_from = "release") %>%
  clean_names()

# Look at it
glimpse(brms_df_86)

# Define the model formula for each sugar
# Random effect excluded from hurdle to improve fit
f_glucose <- bf(glucose ~ ph * lifestyle + (1 | strain), hu ~ ph * lifestyle)
f_maltose <- bf(maltose ~ ph * lifestyle + (1 | strain), hu ~ ph * lifestyle)
f_trehalose <- bf(trehalose ~ ph * lifestyle + (1 | strain), hu ~ ph * lifestyle)

# The model is pretty big so switching to cmdstanr backend
# Run the following hashed lines if it is first time using this backend

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# Then install CmdStan itself
#cmdstanr::install_cmdstan()

library(cmdstanr)

# Setting weakly defined priors to improve performance
priors <- c(
  # Fixed effects for glucose
  set_prior("normal(0, 2)", class = "b", resp = "glucose"),
  set_prior("normal(0, 2)", class = "b", dpar = "hu", resp = "glucose"),
  set_prior("normal(-2, 1.5)", class = "Intercept", resp = "glucose"),
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "hu", resp = "glucose"),
  set_prior("exponential(1)", class = "sd", group = "strain", resp = "glucose"),
  set_prior("exponential(1)", class = "sd", group = "strain", coef = "Intercept", resp = "glucose"),
  set_prior("exponential(1)", class = "shape", resp = "glucose"),
  
  # Repeat for maltose
  set_prior("normal(0, 2)", class = "b", resp = "maltose"),
  set_prior("normal(0, 2)", class = "b", dpar = "hu", resp = "maltose"),
  set_prior("normal(-2, 1.5)", class = "Intercept", resp = "maltose"),
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "hu", resp = "maltose"),
  set_prior("exponential(1)", class = "sd", group = "strain", resp = "maltose"),
  set_prior("exponential(1)", class = "sd", group = "strain", coef = "Intercept", resp = "maltose"),
  set_prior("exponential(1)", class = "shape", resp = "maltose"),
  
  # Repeat for trehalose
  set_prior("normal(0, 2)", class = "b", resp = "trehalose"),
  set_prior("normal(0, 2)", class = "b", dpar = "hu", resp = "trehalose"),
  set_prior("normal(-2, 1.5)", class = "Intercept", resp = "trehalose"),
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "hu", resp = "trehalose"),
  set_prior("exponential(1)", class = "sd", group = "strain", resp = "trehalose"),
  set_prior("exponential(1)", class = "sd", group = "strain", coef = "Intercept", resp = "trehalose"),
  set_prior("exponential(1)", class = "shape", resp = "trehalose")
)

#Why these priors?
#normal(0, 2): Allows moderate fixed effects on the log scale — enough flexibility without letting the model go wild.
#Intercept = normal(-2, 1.5): On log scale, this reflects prior expectation of small release amounts when nothing else influences it.
#logistic(0, 1) for hurdle Intercept: Centers on 0.5 probability of zero, still flexible, aligns with data (heavy zero-inflation).
#exponential(1) for all SD and shape parameters: Regularizes these to be small but doesn’t rule out larger variance.

# Fit the model and run (or load from saved file using readRDS below to save time)
#brms_model_86 <- brm(
#  f_glucose + f_maltose + f_trehalose,
#  data = brms_df_86,
#  prior = priors,
#  sample_prior = "no",
#  family = hurdle_gamma(link = "log"),
#  chains = 4, cores = 4, iter = 4000,
#  control = list(adapt_delta = 0.95, max_treedepth = 15),
#  backend = "cmdstanr"
#)
# adapt_delta = 0.95 and max_treedepth = 15 are common tweaks for tricky models 
# with strong zero-inflation and correlations.

#saveRDS(brms_model, "Data/brms_model_86.rds")

# Import the model from file
brms_model_86 <- readRDS("Data/brms_model_86.rds")

# Model summary
summary(brms_model_86)

# Divergent transitions:
rstan::get_num_divergent(brms_model_86$fit)
# Maximum tree depth:
rstan::get_num_max_treedepth(brms_model_86$fit)
# Both are 0, which is ideal -- it means we're exploring the parameter space well

# Do posterior predictive checks, this is looking at the simulated data, the distribution should look like the observed data
pp_check(brms_model_86, resp = "glucose", ndraws = 200)
pp_check(brms_model_86, resp = "maltose", ndraws = 200)
pp_check(brms_model_86, resp = "trehalose", ndraws = 200)

## Plotting hurdle and conditional simulated vs observed values in a loop
## pp_check should be able to do this but is not cooperating so using bayesplot instead

# Set number of draws for plotting
ndraws <- 500

# List of response variable names
response_vars <- c("glucose", "maltose", "trehalose")

# Loop through each sugar
for (resp in response_vars) {
  cat("\n==== Plotting for:", resp, "====\n")
  
  # Extract posterior predictive draws
  yrep <- posterior_predict(brms_model_86, resp = resp, ndraws = ndraws)
  
  # Extract observed values
  y <- brms_model_86$data[[resp]]
  
  ## === HURDLE COMPONENT === ##
  y_binary <- as.numeric(y == 0)
  yrep_binary <- apply(yrep, 1:2, function(x) as.numeric(x == 0))
  
  print(
    ppc_bars(y = y_binary, yrep = yrep_binary) +
      ggtitle(paste(resp, "- Hurdle (Zero vs Positive)")) +
      theme_minimal()
  )
  
  ## === CONDITIONAL COMPONENT === ##
  pos_idx <- which(y > 0)
  y_pos <- y[pos_idx]
  yrep_pos <- yrep[, pos_idx, drop = FALSE]
  
  if (length(pos_idx) > 5) {
    print(
      ppc_dens_overlay(y = y_pos, yrep = yrep_pos) +
        ggtitle(paste(resp, "- Conditional (Positive Only)")) +
        theme_minimal()
    )
  } else {
    message(paste("Too few positive values for conditional plot of", resp))
  }
}

# Some overestimation on the low end in conditional model for glucose and maltose, 
# but otherwise a good fit for both hurdle and conditional

# Residuals vs. fitted values -- these look fine
pp_check(brms_model_86, lifestyle = "error_binned", resp = "glucose")
pp_check(brms_model_86, lifestyle = "error_binned", resp = "maltose")
pp_check(brms_model_86, lifestyle = "error_binned", resp = "trehalose")

# Leave‑One‑Out Cross‑Validation (LOO)
# Compare predictive accuracy or check for problematic observations.
loo_fit_86 <- loo(brms_model_86)
print(loo_fit_86)
plot(loo_fit_86)
# These values fine, the bad ones are not frequent enough to worry about

# Bayesian R2 -- how much variation is explained by model
bayes_R2(brms_model_86)

# Strong predictive power for glucose and maltose, very low predctive power for trehalose, just as expected

# Conditional effects pull out the fitted lines and points with 95% confidence intervals
conditional_effects(brms_model_86, resp = "glucose")
conditional_effects(brms_model_86, resp = "maltose")
conditional_effects(brms_model_86, resp = "trehalose")
# These look very sensible

# Plotting the observed data vs estimates for publication
# Posterior predictions for the hurdle (probability of > 0)
hurdle_df_86 <- bind_rows(data.frame(conditional_effects(brms_model_86, dpar = "hu", effects = "lifestyle:ph")[[1]], sugar = "Glucose"),
                          data.frame(conditional_effects(brms_model_86, dpar = "hu", effects = "lifestyle:ph")[[2]], sugar = "Maltose"),
                          data.frame(conditional_effects(brms_model_86, dpar = "hu", effects = "lifestyle:ph")[[3]], sugar = "Trehalose")) %>%
  filter(ph != 6.550000) %>%
  dplyr::select(lifestyle, effect2__, estimate__, se__, lower__, upper__, sugar) %>%
  dplyr::rename("ph" = effect2__,
                "prop_zero" = estimate__,
                "se" = se__,
                "lower" = lower__,
                "upper" = upper__) %>%
  mutate(lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                               lifestyle == "FL" ~ "Free-living"))

glimpse(hurdle_df_86)

prop_zero_df_86 <- df_tidy_sugar_86 %>%
  mutate(zero = release == 0,
         lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                               lifestyle == "FL" ~ "Free-living")) %>%
  group_by(lifestyle, ph, sugar) %>%
  summarise(prop_zero = mean(zero), n = n(), .groups = "drop")

# Observed vs. estimates for hurdle

ggplot(prop_zero_df_86, aes(x = factor(ph), y = 1-prop_zero)) +
  geom_col(aes(fill = lifestyle), width = 0.6, alpha = 0.5) +
  geom_point(data = hurdle_df_86, aes(colour = lifestyle), size = 2) +
  geom_errorbar(data = hurdle_df, aes(ymin = 1-lower, ymax = 1-upper, 
                                      colour = lifestyle), width = 0.2) +
  labs(
    y = "Proportion of non-zero values",
    x = "pH"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values =
                      c("Free-living" = "green3",
                        "Symbiotic" = "cornflowerblue")) +
  scale_colour_manual(values =
                        c("Free-living" = "green3",
                          "Symbiotic" = "cornflowerblue")) +
  facet_grid(lifestyle~sugar)

# Posterior predictions for the conditional (amount)
cond_df_86 <- bind_rows(data.frame(conditional_effects(brms_model_86, effects = "lifestyle:ph")[[1]], sugar = "Glucose"), 
                        data.frame(conditional_effects(brms_model_86, effects = "lifestyle:ph")[[2]], sugar = "Maltose"),
                        data.frame(conditional_effects(brms_model_86, effects = "lifestyle:ph")[[3]], sugar = "Trehalose")) %>%
  filter(ph != 6.550000) %>%
  dplyr::select(lifestyle, effect2__, estimate__, se__, lower__, upper__, sugar) %>%
  dplyr::rename("ph" = effect2__,
                "estimate" = estimate__,
                "se" = se__,
                "lower" = lower__,
                "upper" = upper__) %>%
  mutate(ph = as.numeric(as.character(ph)),
         lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                               lifestyle == "FL" ~ "Free-living"))

glimpse(cond_df_86)

nonzero_df_86 <- df_tidy_sugar_86 %>%
  filter(release > 0) %>%
  mutate(lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                               lifestyle == "FL" ~ "Free-living"))

cond_df_86 %>%
  ggplot(aes(x = factor(ph), y = estimate, colour = lifestyle)) +
  geom_point() +
  geom_jitter(data = nonzero_df_86 %>% mutate(estimate = release), alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  facet_wrap(~sugar) +
  scale_colour_manual(values =
                        c("Free-living" = "green3",
                          "Symbiotic" = "cornflowerblue")) +
  theme_bw() +
  labs(y = "Sugar release (ppm)",
       x = "pH",
       colour = "Lifestyle") +
  theme(legend.position = "bottom")

# Plot strain level effects
# Strains that don't overlap with 0 meaningfully diverge from the average
ranef_array_86 <- ranef(brms_model_86)$strain
glucose_ranef_86 <- as.data.frame(ranef_array_86[, , "glucose_Intercept"]) %>%
  mutate(sugar = "Glucose")
glucose_ranef_86$strain <- rownames(ranef_array_86)
maltose_ranef_86 <- as.data.frame(ranef_array_86[, , "maltose_Intercept"]) %>%
  mutate(sugar = "Maltose")
maltose_ranef_86$strain <- rownames(ranef_array_86)
trehalose_ranef_86 <- as.data.frame(ranef_array_86[, , "trehalose_Intercept"]) %>%
  mutate(sugar = "Trehalose")
trehalose_ranef_86$strain <- rownames(ranef_array_86)

ranef_df_86 <- bind_rows(glucose_ranef_86, maltose_ranef_86, trehalose_ranef_86)

ranef_df_86 %>%
  mutate(lifestyle = ifelse(str_detect(strain, "s") == TRUE, "Symbiotic", "Free-living")) %>%
  ggplot(aes(x = strain, y = Estimate, colour = lifestyle)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5)) +
  scale_colour_manual(values =
                        c("Free-living" = "green3",
                          "Symbiotic" = "cornflowerblue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  labs(x = "Strain",
       y = "Strain-level intercept deviation",
       colour = "Lifestyle") +
  facet_wrap(~ sugar, ncol = 1)


### Pairwise contrasts
# Extract slope for conditional and hurdle part of model for each sugar
mu_trend_glucose_86 <- emtrends(brms_model_86, var = "ph", specs = "lifestyle", resp = "glucose")
hu_trend_glucose_86 <- emtrends(brms_model_86, specs = "lifestyle", var = "ph", dpar = "hu", resp = "glucose")
mu_trend_maltose_86 <- emtrends(brms_model_86, var = "ph", specs = "lifestyle", resp = "maltose")
hu_trend_maltose_86 <- emtrends(brms_model_86, specs = "lifestyle", var = "ph", dpar = "hu", resp = "maltose")
mu_trend_trehalose_86 <- emtrends(brms_model_86, var = "ph", specs = "lifestyle", resp = "trehalose")
hu_trend_trehalose_86 <- emtrends(brms_model_86, specs = "lifestyle", var = "ph", dpar = "hu", resp = "trehalose")

# View estimated marginal slopes
summary(mu_trend_glucose_86)
summary(hu_trend_glucose_86)
summary(mu_trend_maltose_86)
summary(hu_trend_maltose_86)
summary(mu_trend_trehalose_86)
summary(hu_trend_trehalose_86)

# Direct contrasts, if the HPD interval overlaps 0 there's no difference
bind_rows(data.frame(pairs(mu_trend_glucose_86),
                     sugar = "glucose",
                     model = "mu"),
          data.frame(pairs(hu_trend_glucose_86),
                     sugar = "glucose",
                     model = "hu"),
          data.frame(pairs(mu_trend_maltose_86),
                     sugar = "maltose",
                     model = "mu"),
          data.frame(pairs(hu_trend_maltose_86),
                     sugar = "maltose",
                     model = "hu"),
          data.frame(pairs(mu_trend_trehalose_86),
                     sugar = "trehalose",
                     model = "mu"),
          data.frame(pairs(hu_trend_trehalose_86),
                     sugar = "trehalose",
                     model = "hu")) %>%
  mutate(`overlap?` = ifelse(lower.HPD < 0 & upper.HPD > 0, "yes", "no"))

# No changed conclusions


# Reanalysis using only the 92 strains from the phylogeny to check robustness of results ###########################
# Sugar release for each strain by pH

df_tidy_sugar_92 <- df_tidy_sugar %>%
  filter(strain != "192")


df_tidy_sugar %>%
  filter(strain != "192") %>%
  group_by(ph, lifestyle, strain) %>%
  summarise(
    mean = mean(release),
    se = sd(release)/sqrt(n())
  ) %>%
  ggplot(aes(x = factor(ph), y = mean, colour = lifestyle, group = lifestyle)) +
  geom_hline(yintercept = 0) +
  geom_point(data = df_tidy_sugar %>% mutate(mean = release) %>%  filter(strain != "192"), alpha = 0.3) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
  theme_bw() +
  labs(y = "Mean total sugar release (ppm)",
       colour = "Strain lifestyle",
       x = "pH") +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  facet_wrap(~strain)

### Fitting a model
# The data is both right skewed and has a huge amount of 0s
df_tidy_sugar %>%
  filter(strain != "192") %>%
  ggplot(aes(x = release, fill = lifestyle)) +
  geom_histogram() +
  facet_wrap(~sugar)

# Need to fit a model to account for
# 1. Do they release anything at all?
# 2. If so, how much do they release

## Fitting Bayesian regression model with brms()

# Pivot the data frame to allow fitting each sugar as a separate response variable
brms_df_92 <- df_tidy_sugar %>%
  filter(strain != "192") %>%
  pivot_wider(names_from = "sugar", values_from = "release") %>%
  clean_names()

# Look at it
glimpse(brms_df_92)

# Define the model formula for each sugar
# Random effect excluded from hurdle to improve fit
f_glucose <- bf(glucose ~ ph * lifestyle + (1 | strain), hu ~ ph * lifestyle)
f_maltose <- bf(maltose ~ ph * lifestyle + (1 | strain), hu ~ ph * lifestyle)
f_trehalose <- bf(trehalose ~ ph * lifestyle + (1 | strain), hu ~ ph * lifestyle)

# The model is pretty big so switching to cmdstanr backend
# Run the following hashed lines if it is first time using this backend

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# Then install CmdStan itself
#cmdstanr::install_cmdstan()

library(cmdstanr)

# Setting weakly defined priors to improve performance
priors <- c(
  # Fixed effects for glucose
  set_prior("normal(0, 2)", class = "b", resp = "glucose"),
  set_prior("normal(0, 2)", class = "b", dpar = "hu", resp = "glucose"),
  set_prior("normal(-2, 1.5)", class = "Intercept", resp = "glucose"),
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "hu", resp = "glucose"),
  set_prior("exponential(1)", class = "sd", group = "strain", resp = "glucose"),
  set_prior("exponential(1)", class = "sd", group = "strain", coef = "Intercept", resp = "glucose"),
  set_prior("exponential(1)", class = "shape", resp = "glucose"),
  
  # Repeat for maltose
  set_prior("normal(0, 2)", class = "b", resp = "maltose"),
  set_prior("normal(0, 2)", class = "b", dpar = "hu", resp = "maltose"),
  set_prior("normal(-2, 1.5)", class = "Intercept", resp = "maltose"),
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "hu", resp = "maltose"),
  set_prior("exponential(1)", class = "sd", group = "strain", resp = "maltose"),
  set_prior("exponential(1)", class = "sd", group = "strain", coef = "Intercept", resp = "maltose"),
  set_prior("exponential(1)", class = "shape", resp = "maltose"),
  
  # Repeat for trehalose
  set_prior("normal(0, 2)", class = "b", resp = "trehalose"),
  set_prior("normal(0, 2)", class = "b", dpar = "hu", resp = "trehalose"),
  set_prior("normal(-2, 1.5)", class = "Intercept", resp = "trehalose"),
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "hu", resp = "trehalose"),
  set_prior("exponential(1)", class = "sd", group = "strain", resp = "trehalose"),
  set_prior("exponential(1)", class = "sd", group = "strain", coef = "Intercept", resp = "trehalose"),
  set_prior("exponential(1)", class = "shape", resp = "trehalose")
)

#Why these priors?
#normal(0, 2): Allows moderate fixed effects on the log scale — enough flexibility without letting the model go wild.
#Intercept = normal(-2, 1.5): On log scale, this reflects prior expectation of small release amounts when nothing else influences it.
#logistic(0, 1) for hurdle Intercept: Centers on 0.5 probability of zero, still flexible, aligns with data (heavy zero-inflation).
#exponential(1) for all SD and shape parameters: Regularizes these to be small but doesn’t rule out larger variance.

# Fit the model and run (or load from saved file using readRDS below to save time)
#brms_model_92 <- brm(
#  f_glucose + f_maltose + f_trehalose,
#  data = brms_df_92,
#  prior = priors,
#  sample_prior = "no",
#  family = hurdle_gamma(link = "log"),
#  chains = 4, cores = 4, iter = 4000,
#  control = list(adapt_delta = 0.95, max_treedepth = 15),
#  backend = "cmdstanr"
#)
# adapt_delta = 0.95 and max_treedepth = 15 are common tweaks for tricky models 
# with strong zero-inflation and correlations.

 #saveRDS(brms_model, "Data/brms_model_92.rds")

# Import the model from file
brms_model_92 <- readRDS("Data/brms_model_92.rds")

# Model summary
summary(brms_model_92)

# Divergent transitions:
rstan::get_num_divergent(brms_model_92$fit)
# Maximum tree depth:
rstan::get_num_max_treedepth(brms_model_92$fit)
# Both are 0, which is ideal -- it means we're exploring the parameter space well

# Do posterior predictive checks, this is looking at the simulated data, the distribution should look like the observed data
pp_check(brms_model_92, resp = "glucose", ndraws = 200)
pp_check(brms_model_92, resp = "maltose", ndraws = 200)
pp_check(brms_model_92, resp = "trehalose", ndraws = 200)

## Plotting hurdle and conditional simulated vs observed values in a loop
## pp_check should be able to do this but is not cooperating so using bayesplot instead

# Set number of draws for plotting
ndraws <- 500

# List of response variable names
response_vars <- c("glucose", "maltose", "trehalose")

# Loop through each sugar
for (resp in response_vars) {
  cat("\n==== Plotting for:", resp, "====\n")
  
  # Extract posterior predictive draws
  yrep <- posterior_predict(brms_model_92, resp = resp, ndraws = ndraws)
  
  # Extract observed values
  y <- brms_model_92$data[[resp]]
  
  ## === HURDLE COMPONENT === ##
  y_binary <- as.numeric(y == 0)
  yrep_binary <- apply(yrep, 1:2, function(x) as.numeric(x == 0))
  
  print(
    ppc_bars(y = y_binary, yrep = yrep_binary) +
      ggtitle(paste(resp, "- Hurdle (Zero vs Positive)")) +
      theme_minimal()
  )
  
  ## === CONDITIONAL COMPONENT === ##
  pos_idx <- which(y > 0)
  y_pos <- y[pos_idx]
  yrep_pos <- yrep[, pos_idx, drop = FALSE]
  
  if (length(pos_idx) > 5) {
    print(
      ppc_dens_overlay(y = y_pos, yrep = yrep_pos) +
        ggtitle(paste(resp, "- Conditional (Positive Only)")) +
        theme_minimal()
    )
  } else {
    message(paste("Too few positive values for conditional plot of", resp))
  }
}

# Some overestimation on the low end in conditional model for glucose and maltose, 
# but otherwise a good fit for both hurdle and conditional

# Residuals vs. fitted values -- these look fine
pp_check(brms_model_92, lifestyle = "error_binned", resp = "glucose")
pp_check(brms_model_92, lifestyle = "error_binned", resp = "maltose")
pp_check(brms_model_92, lifestyle = "error_binned", resp = "trehalose")

# Leave‑One‑Out Cross‑Validation (LOO)
# Compare predictive accuracy or check for problematic observations.
loo_fit_92 <- loo(brms_model_92)
print(loo_fit_92)
plot(loo_fit_92)
# These values fine, the bad ones are not frequent enough to worry about

# Bayesian R2 -- how much variation is explained by model
bayes_R2(brms_model_92)

# Strong predictive power for glucose and maltose, very low predctive power for trehalose, just as expected

# Conditional effects pull out the fitted lines and points with 95% confidence intervals
conditional_effects(brms_model_92, resp = "glucose")
conditional_effects(brms_model_92, resp = "maltose")
conditional_effects(brms_model_92, resp = "trehalose")
# These look very sensible

# Plotting the observed data vs estimates for publication
# Posterior predictions for the hurdle (probability of > 0)
hurdle_df_92 <- bind_rows(data.frame(conditional_effects(brms_model_92, dpar = "hu", effects = "lifestyle:ph")[[1]], sugar = "Glucose"),
                       data.frame(conditional_effects(brms_model_92, dpar = "hu", effects = "lifestyle:ph")[[2]], sugar = "Maltose"),
                       data.frame(conditional_effects(brms_model_92, dpar = "hu", effects = "lifestyle:ph")[[3]], sugar = "Trehalose")) %>%
  filter(ph != 6.550000) %>%
  dplyr::select(lifestyle, effect2__, estimate__, se__, lower__, upper__, sugar) %>%
  dplyr::rename("ph" = effect2__,
                "prop_zero" = estimate__,
                "se" = se__,
                "lower" = lower__,
                "upper" = upper__) %>%
  mutate(lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                               lifestyle == "FL" ~ "Free-living"))

glimpse(hurdle_df_92)

prop_zero_df_92 <- df_tidy_sugar_92 %>%
  mutate(zero = release == 0,
         lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                               lifestyle == "FL" ~ "Free-living")) %>%
  group_by(lifestyle, ph, sugar) %>%
  summarise(prop_zero = mean(zero), n = n(), .groups = "drop")

# Observed vs. estimates for hurdle

ggplot(prop_zero_df_92, aes(x = factor(ph), y = 1-prop_zero)) +
  geom_col(aes(fill = lifestyle), width = 0.6, alpha = 0.5) +
  geom_point(data = hurdle_df_92, aes(colour = lifestyle), size = 2) +
  geom_errorbar(data = hurdle_df, aes(ymin = 1-lower, ymax = 1-upper, 
                                      colour = lifestyle), width = 0.2) +
  labs(
    y = "Proportion of non-zero values",
    x = "pH"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values =
                      c("Free-living" = "green3",
                        "Symbiotic" = "cornflowerblue")) +
  scale_colour_manual(values =
                        c("Free-living" = "green3",
                          "Symbiotic" = "cornflowerblue")) +
  facet_grid(lifestyle~sugar)

# Posterior predictions for the conditional (amount)
cond_df_92 <- bind_rows(data.frame(conditional_effects(brms_model_92, effects = "lifestyle:ph")[[1]], sugar = "Glucose"), 
                     data.frame(conditional_effects(brms_model_92, effects = "lifestyle:ph")[[2]], sugar = "Maltose"),
                     data.frame(conditional_effects(brms_model_92, effects = "lifestyle:ph")[[3]], sugar = "Trehalose")) %>%
  filter(ph != 6.550000) %>%
  dplyr::select(lifestyle, effect2__, estimate__, se__, lower__, upper__, sugar) %>%
  dplyr::rename("ph" = effect2__,
                "estimate" = estimate__,
                "se" = se__,
                "lower" = lower__,
                "upper" = upper__) %>%
  mutate(ph = as.numeric(as.character(ph)),
         lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                               lifestyle == "FL" ~ "Free-living"))

glimpse(cond_df_92)

nonzero_df_92 <- df_tidy_sugar_92 %>%
  filter(release > 0) %>%
  mutate(lifestyle = case_when(lifestyle == "S" ~ "Symbiotic",
                               lifestyle == "FL" ~ "Free-living"))

cond_df_92 %>%
  ggplot(aes(x = factor(ph), y = estimate, colour = lifestyle)) +
  geom_point() +
  geom_jitter(data = nonzero_df_92 %>% mutate(estimate = release), alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  facet_wrap(~sugar) +
  scale_colour_manual(values =
                        c("Free-living" = "green3",
                          "Symbiotic" = "cornflowerblue")) +
  theme_bw() +
  labs(y = "Sugar release (ppm)",
       x = "pH",
       colour = "Lifestyle") +
  theme(legend.position = "bottom")

# Plot strain level effects
# Strains that don't overlap with 0 meaningfully diverge from the average
ranef_array_92 <- ranef(brms_model_92)$strain
glucose_ranef_92 <- as.data.frame(ranef_array_92[, , "glucose_Intercept"]) %>%
  mutate(sugar = "Glucose")
glucose_ranef_92$strain <- rownames(ranef_array_92)
maltose_ranef_92 <- as.data.frame(ranef_array_92[, , "maltose_Intercept"]) %>%
  mutate(sugar = "Maltose")
maltose_ranef_92$strain <- rownames(ranef_array_92)
trehalose_ranef_92 <- as.data.frame(ranef_array_92[, , "trehalose_Intercept"]) %>%
  mutate(sugar = "Trehalose")
trehalose_ranef_92$strain <- rownames(ranef_array_92)

ranef_df_92 <- bind_rows(glucose_ranef_92, maltose_ranef_92, trehalose_ranef_92)

ranef_df_92 %>%
  mutate(lifestyle = ifelse(str_detect(strain, "s") == TRUE, "Symbiotic", "Free-living")) %>%
  ggplot(aes(x = strain, y = Estimate, colour = lifestyle)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5)) +
  scale_colour_manual(values =
                        c("Free-living" = "green3",
                          "Symbiotic" = "cornflowerblue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom") +
  labs(x = "Strain",
       y = "Strain-level intercept deviation",
       colour = "Lifestyle") +
  facet_wrap(~ sugar, ncol = 1)


### Pairwise contrasts
# Extract slope for conditional and hurdle part of model for each sugar
mu_trend_glucose_92 <- emtrends(brms_model_92, var = "ph", specs = "lifestyle", resp = "glucose")
hu_trend_glucose_92 <- emtrends(brms_model_92, specs = "lifestyle", var = "ph", dpar = "hu", resp = "glucose")
mu_trend_maltose_92 <- emtrends(brms_model_92, var = "ph", specs = "lifestyle", resp = "maltose")
hu_trend_maltose_92 <- emtrends(brms_model_92, specs = "lifestyle", var = "ph", dpar = "hu", resp = "maltose")
mu_trend_trehalose_92 <- emtrends(brms_model_92, var = "ph", specs = "lifestyle", resp = "trehalose")
hu_trend_trehalose_92 <- emtrends(brms_model_92, specs = "lifestyle", var = "ph", dpar = "hu", resp = "trehalose")

# View estimated marginal slopes
summary(mu_trend_glucose_92)
summary(hu_trend_glucose_92)
summary(mu_trend_maltose_92)
summary(hu_trend_maltose_92)
summary(mu_trend_trehalose_92)
summary(hu_trend_trehalose_92)

# Direct contrasts, if the HPD interval overlaps 0 there's no difference
bind_rows(data.frame(pairs(mu_trend_glucose_92),
                     sugar = "glucose",
                     model = "mu"),
          data.frame(pairs(hu_trend_glucose_92),
                     sugar = "glucose",
                     model = "hu"),
          data.frame(pairs(mu_trend_maltose_92),
                     sugar = "maltose",
                     model = "mu"),
          data.frame(pairs(hu_trend_maltose_92),
                     sugar = "maltose",
                     model = "hu"),
          data.frame(pairs(mu_trend_trehalose_92),
                     sugar = "trehalose",
                     model = "mu"),
          data.frame(pairs(hu_trend_trehalose_92),
                     sugar = "trehalose",
                     model = "hu")) %>%
  mutate(`overlap?` = ifelse(lower.HPD < 0 & upper.HPD > 0, "yes", "no"))

# No changed conclusions
  