# Erika M Hansson 2025
# erika.hansson@manchester.ac.uk

#### ---------------------------------------------------- ####
# 1) Plotting fluorimetry data
# 2) Fitting a model
# Then: repeating all of the above with the 86 strains in the phylogeny for robustness
#### ---------------------------------------------------- ####

# Load packages
if(!require("tidyverse")) {install.packages("tidyverse")} # wrangling, plotting
if(!require("broom")) {install.packages("broom")} # for wrangling
if(!require("nlme")) {install.packages("nlme")} # for mixed effects models
if(!require("car")) {install.packages("car")} # stats tools

# Load data frame
df_tidy_fl <- read_csv("Data/df_tidy_fluorimetry.csv")

# Plot the means
df_tidy_fl %>%
  group_by(lifestyle) %>%
  summarise(
    mean = mean(mean_FvFm),
    se = sd(mean_FvFm)/sqrt(n())
  ) %>%
  ggplot(aes(x = lifestyle, y = mean, colour = lifestyle)) +
  geom_jitter(data = df_tidy_fl %>% mutate(mean = mean_FvFm), alpha = 0.2, width = 0.05) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, linewidth = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Lifestyle", y = "Mean Fv/Fm") +
  guides(colour = "none")

# And each individual strain
df_tidy_fl %>%
  group_by(strain, lifestyle) %>%
  summarise(
    mean = mean(mean_FvFm),
    se = sd(mean_FvFm)/sqrt(n())
  ) %>%
  ggplot(aes(x = strain, y = mean, colour = lifestyle)) +
  geom_point(data = df_tidy_fl %>% mutate(mean = mean_FvFm), alpha = 0.4) +
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 1, linewidth = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Strain", y = "Mean Fv/Fm", colour = "Lifestyle")

# Distribution looks fine
df_tidy_fl %>%
  ggplot(aes(x = log(mean_FvFm))) +
  geom_histogram() +
  facet_wrap(~lifestyle)

# Is there a difference between means?
# As it looks like we have more variance in Fl than S, let's model it explicitly
# (so using nlme rather than lme4 for our ME model)

fl_mod <- lme(
  mean_FvFm ~ lifestyle,
  random = ~1 | strain,
  weights = varIdent(form = ~1 | lifestyle),
  data = df_tidy_fl
)

# Make a data frame for diagnostics
diagnostics_fl <- broom.mixed::augment(fl_mod)

# Here's a Levene's test to check that variance actually isn't homogeneous between groups
leveneTest(.resid ~ lifestyle, data = diagnostics_fl)

# But let's do our due diligence and compare our model to one with homogeneous variance
fl_mod_hom <- lme(
  mean_FvFm ~ lifestyle,
  random = ~1 | strain,
  data = df_tidy_fl
)

anova(fl_mod_hom, fl_mod)  # LRT to see if varIdent improves fit

# It's marginal, so let's keep it in

# QQ plot -- fine
ggplot(diagnostics_fl, aes(sample = .resid)) +
  geom_qq() + 
  geom_qq_line() +
  theme_bw()

# Residuals vs. fitted values -- fine
ggplot(diagnostics_fl, aes(x = .fitted, y = .resid, colour = lifestyle)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw()

# Now let's look at what the model tells us
summary(fl_mod)
Anova(fl_mod)

# This suggets marginal significance but let's do this properly and LRT a p-value

fl_mod_null_ml <- lme(
  mean_FvFm ~ 1,                     # intercept only
  random = ~1 | strain,              
  weights = varIdent(~1 | lifestyle), # will be ignored
  data = df_tidy_fl,
  method = "ML" # needed for LRT
)

# Updating our model to be "ML"
fl_mod_ml <- update(fl_mod, method = "ML")

# LRT
anova(fl_mod_null_ml, fl_mod_ml)

# That actually shifts our interpretation, there is a difference in means

# How big is the difference in spread?

# SD and variance per lifestyle
df_tidy_fl %>%
  group_by(lifestyle) %>%
  summarise(
    sd_FvFm = sd(mean_FvFm), 
            var_FvFm = var(mean_FvFm)
            )

# Variance ratio (FL/S)
var_ratio <- var(df_tidy_fl$mean_FvFm[df_tidy_fl$lifestyle == "FL"]) /
  var(df_tidy_fl$mean_FvFm[df_tidy_fl$lifestyle == "S"])

var_ratio
# Variance is 3.3x in FL compared to S

# Let's bootstrap a CI
set.seed(123)
B <- 5000
boot_var_ratio <- replicate(B, {
  var(sample(df_tidy_fl$mean_FvFm[df_tidy_fl$lifestyle == "FL"], replace = TRUE)) /
    var(sample(df_tidy_fl$mean_FvFm[df_tidy_fl$lifestyle == "S"], replace = TRUE))
})
quantile(boot_var_ratio, c(0.025, 0.5, 0.975))

# So 95% CI: 2.22--5.01



# Reanalysis using only the 86 strains from the phylogeny to check robustness of results ###########################
# Plot the means
df_tidy_fl %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  group_by(lifestyle) %>%
  summarise(
    mean = mean(mean_FvFm),
    se = sd(mean_FvFm)/sqrt(n())
  ) %>%
  ggplot(aes(x = lifestyle, y = mean, colour = lifestyle)) +
  geom_jitter(data = df_tidy_fl %>% filter(strain != "155",
                                           strain != "192",
                                           strain != "242",
                                           strain != "252",
                                           strain != "258",
                                           strain != "s354") %>% mutate(mean = mean_FvFm), alpha = 0.2, width = 0.05) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, linewidth = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Lifestyle", y = "Mean Fv/Fm") +
  guides(colour = "none")

# And each individual strain
df_tidy_fl %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  group_by(strain, lifestyle) %>%
  summarise(
    mean = mean(mean_FvFm),
    se = sd(mean_FvFm)/sqrt(n())
  ) %>%
  ggplot(aes(x = strain, y = mean, colour = lifestyle)) +
  geom_point(data = df_tidy_fl %>% filter(strain != "155",
                                          strain != "192",
                                          strain != "242",
                                          strain != "252",
                                          strain != "258",
                                          strain != "s354") %>% mutate(mean = mean_FvFm), alpha = 0.4) +
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 1, linewidth = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Strain", y = "Mean Fv/Fm", colour = "Lifestyle")

# Distribution looks fine
df_tidy_fl %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  ggplot(aes(x = log(mean_FvFm))) +
  geom_histogram() +
  facet_wrap(~lifestyle)

# Is there a difference between means?
# As it looks like we have more variance in Fl than S, let's model it explicitly
# (so using nlme rather than lme4 for our ME model)

df_tidy_fl_86 <- df_tidy_fl %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354")

fl_mod_86 <- lme(
  mean_FvFm ~ lifestyle,
  random = ~1 | strain,
  weights = varIdent(form = ~1 | lifestyle),
  data = df_tidy_fl_86
)

# Make a data frame for diagnostics
diagnostics_fl_86 <- broom.mixed::augment(fl_mod_86)

# Here's a Levene's test to check that variance actually isn't homogeneous between groups
leveneTest(.resid ~ lifestyle, data = diagnostics_fl_86)

# But let's do our due diligence and compare our model to one with homogeneous variance
fl_mod_hom_86 <- lme(
  mean_FvFm ~ lifestyle,
  random = ~1 | strain,
  data = df_tidy_fl_86
)

anova(fl_mod_hom_86, fl_mod_86)  # LRT to see if varIdent improves fit

# It's marginal, so let's keep it in

# QQ plot -- fine
ggplot(diagnostics_fl_86, aes(sample = .resid)) +
  geom_qq() + 
  geom_qq_line() +
  theme_bw()

# Residuals vs. fitted values -- fine
ggplot(diagnostics_fl_86, aes(x = .fitted, y = .resid, colour = lifestyle)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw()

# Now let's look at what the model tells us
summary(fl_mod_86)
Anova(fl_mod_86)

# This suggets marginal significance but let's do this properly and LRT a p-value

fl_mod_null_ml_86 <- lme(
  mean_FvFm ~ 1,                     # intercept only
  random = ~1 | strain,              
  weights = varIdent(~1 | lifestyle), # will be ignored
  data = df_tidy_fl_86,
  method = "ML" # needed for LRT
)

# Updating our model to be "ML"
fl_mod_ml_86 <- update(fl_mod_86, method = "ML")

# LRT
anova(fl_mod_null_ml_86, fl_mod_ml_86)

# That actually shifts our interpretation, there is a difference in means

# How big is the difference in spread?

# SD and variance per lifestyle
df_tidy_fl_86 %>%
  group_by(lifestyle) %>%
  summarise(
    sd_FvFm = sd(mean_FvFm), 
    var_FvFm = var(mean_FvFm)
  )

# Variance ratio (FL/S)
var_ratio_86 <- var(df_tidy_fl_86$mean_FvFm[df_tidy_fl_86$lifestyle == "FL"]) /
  var(df_tidy_fl_86$mean_FvFm[df_tidy_fl_86$lifestyle == "S"])

var_ratio_86
# Variance is 3.3x in FL compared to S

# Let's bootstrap a CI
set.seed(123)
B <- 5000
boot_var_ratio_86 <- replicate(B, {
  var(sample(df_tidy_fl_86$mean_FvFm[df_tidy_fl_86$lifestyle == "FL"], replace = TRUE)) /
    var(sample(df_tidy_fl_86$mean_FvFm[df_tidy_fl_86$lifestyle == "S"], replace = TRUE))
})
quantile(boot_var_ratio_86, c(0.025, 0.5, 0.975))

# So 95% CI: 2.22--5.01

# So with reanalysis we lose the difference in mean Fv/Fm by lifestyle, and the model no longer needs
# split variance. Variance ratio difference still considerable though.

# Reanalysis using only the 92 strains from the phylogeny to check robustness of results ###########################
# Plot the means
df_tidy_fl %>%
  filter(strain != "192") %>%
  group_by(lifestyle) %>%
  summarise(
    mean = mean(mean_FvFm),
    se = sd(mean_FvFm)/sqrt(n())
  ) %>%
  ggplot(aes(x = lifestyle, y = mean, colour = lifestyle)) +
  geom_jitter(data = df_tidy_fl %>% filter(strain != "192") %>% mutate(mean = mean_FvFm), alpha = 0.2, width = 0.05) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, linewidth = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Lifestyle", y = "Mean Fv/Fm") +
  guides(colour = "none")

# And each individual strain
df_tidy_fl %>%
  filter(strain != "192") %>%
  group_by(strain, lifestyle) %>%
  summarise(
    mean = mean(mean_FvFm),
    se = sd(mean_FvFm)/sqrt(n())
  ) %>%
  ggplot(aes(x = strain, y = mean, colour = lifestyle)) +
  geom_point(data = df_tidy_fl %>% filter(strain != "192") %>% mutate(mean = mean_FvFm), alpha = 0.4) +
  geom_point(size = 1) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 1, linewidth = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Strain", y = "Mean Fv/Fm", colour = "Lifestyle")

# Distribution looks fine
df_tidy_fl %>%
  filter(strain != "192") %>%
  ggplot(aes(x = log(mean_FvFm))) +
  geom_histogram() +
  facet_wrap(~lifestyle)

# Is there a difference between means?
# As it looks like we have more variance in Fl than S, let's model it explicitly
# (so using nlme rather than lme4 for our ME model)

df_tidy_fl_92 <- df_tidy_fl %>%
  filter(strain != "192")

fl_mod_92 <- lme(
  mean_FvFm ~ lifestyle,
  random = ~1 | strain,
  weights = varIdent(form = ~1 | lifestyle),
  data = df_tidy_fl_92
)

# Make a data frame for diagnostics
diagnostics_fl_92 <- broom.mixed::augment(fl_mod_92)

# Here's a Levene's test to check that variance actually isn't homogeneous between groups
leveneTest(.resid ~ lifestyle, data = diagnostics_fl_92)

# But let's do our due diligence and compare our model to one with homogeneous variance
fl_mod_hom_92 <- lme(
  mean_FvFm ~ lifestyle,
  random = ~1 | strain,
  data = df_tidy_fl_92
)

anova(fl_mod_hom_92, fl_mod_92)  # LRT to see if varIdent improves fit

# It's marginal, so let's keep it in

# QQ plot -- fine
ggplot(diagnostics_fl_92, aes(sample = .resid)) +
  geom_qq() + 
  geom_qq_line() +
  theme_bw()

# Residuals vs. fitted values -- fine
ggplot(diagnostics_fl_92, aes(x = .fitted, y = .resid, colour = lifestyle)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw()

# Now let's look at what the model tells us
summary(fl_mod_92)
Anova(fl_mod_92)

# No significance, but let's do this properly and LRT a p-value

fl_mod_null_ml_92 <- lme(
  mean_FvFm ~ 1,                     # intercept only
  random = ~1 | strain,              
  weights = varIdent(~1 | lifestyle), # will be ignored
  data = df_tidy_fl_92,
  method = "ML" # needed for LRT
)

# Updating our model to be "ML"
fl_mod_ml_92 <- update(fl_mod_92, method = "ML")

# LRT
anova(fl_mod_null_ml_92, fl_mod_ml_92)

# Still not significant but it's now marginal

# How big is the difference in spread?

# SD and variance per lifestyle
df_tidy_fl_92 %>%
  group_by(lifestyle) %>%
  summarise(
    sd_FvFm = sd(mean_FvFm), 
    var_FvFm = var(mean_FvFm)
  )

# Variance ratio (FL/S)
var_ratio_92 <- var(df_tidy_fl_92$mean_FvFm[df_tidy_fl_92$lifestyle == "FL"]) /
  var(df_tidy_fl_92$mean_FvFm[df_tidy_fl_92$lifestyle == "S"])

var_ratio_92
# Variance is 3.3x in FL compared to S

# Let's bootstrap a CI
set.seed(123)
B <- 5000
boot_var_ratio_92 <- replicate(B, {
  var(sample(df_tidy_fl_92$mean_FvFm[df_tidy_fl_92$lifestyle == "FL"], replace = TRUE)) /
    var(sample(df_tidy_fl_92$mean_FvFm[df_tidy_fl_92$lifestyle == "S"], replace = TRUE))
})
quantile(boot_var_ratio_92, c(0.025, 0.5, 0.975))

# So 95% CI: 2.22--5.01

# So with reanalysis we lose the difference in mean Fv/Fm by lifestyle, and the model no longer needs
# split variance. Variance ratio difference still considerable though.