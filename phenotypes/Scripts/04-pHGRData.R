# Erika M Hansson 2025
# erika.hansson@manchester.ac.uk

#### ---------------------------------------------------- ####
# 1) Plotting pH growth rate data
# 2) Fitting a model to pH growth rate data
# 3) Plotting pH growth rate data against sugar release data
# 4) Testing for a correlation between pH growth rate and sugar release
# Then: repeating all of the above with the 86 strains in the phylogeny for robustness
#### ---------------------------------------------------- ####

# Load packages
if(!require("tidyverse")) {install.packages("tidyverse")} # wrangling, plotting
if(!require("emmeans")) {install.packages("emmeans")} # for fitting mixed effects models
if(!require("janitor")) {install.packages("janitor")} # for wrangling
if(!require("car")) {install.packages("car")} # stats tools
if(!require("broom")) {install.packages("broom")} # for wrangling
if(!require("lme4")) {install.packages("lme4")} # for fitting mixed effects models
if(!require("pbkrtest")) {install.packages("pbkrtest")} # for parametric bootstrapping of mixed effects models

# Load data frame
df_tidy_ph <- read_csv("Data/df_tidy_pH.csv") 

glimpse(df_tidy_ph)

# Plot
df_tidy_ph %>%
  group_by(strain, day, lifestyle, ph) %>%
  summarise(mean = mean(od750),
            se = sd(od750)/sqrt(n())) %>%
  ggplot(aes(x = day, y = mean, colour = factor(ph), groups = factor(ph))) +
  geom_point(data = df_tidy_ph %>% mutate(mean = od750), alpha = 0.2) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
  geom_line() +
  facet_wrap(~ strain)

# Calculate growth rate
# A function for growth rate
gr_func <- function(x, y, t) {
  gr <- (log(x/y))/t
  return(gr)
}

# Add to data frame
df_tidy_phgr <- df_tidy_ph %>%
  pivot_wider(names_from = day, names_prefix = "day", values_from = od750) %>%
  mutate(gr = gr_func(day7, day0, 7))
  
# Plot
df_tidy_phgr %>%
  group_by(lifestyle, ph) %>%
  summarise(mean = mean(gr),
            se = sd(gr)/sqrt(n())) %>%
  ggplot(aes(x = factor(ph), y = mean, colour = lifestyle, group = lifestyle)) +
  geom_jitter(data = df_tidy_phgr %>% mutate(mean = gr), alpha = 0.2, width = 0.125) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, linewidth = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "pH", y = "Mean daily growth rate", colour = "Lifestyle")

# Distribution -- we have a bit of biomodality going on in one group
df_tidy_phgr %>%
  ggplot(aes(x = gr, fill = lifestyle)) +
  geom_histogram() +
  facet_grid(lifestyle~ph)


# Fit ME model
mod_phgr <- lmer(gr ~ ph*lifestyle + (1 | strain), data = df_tidy_phgr)

# Make a data frame for diagnostics
diagnostics_phgr <- broom.mixed::augment(mod_phgr)

# QQ plot -- fine
ggplot(diagnostics_phgr, aes(sample = .resid)) +
  geom_qq() + 
  geom_qq_line() +
  theme_bw()

# Residuals vs. fitted values -- some patterns here, but model should be robust given QQ plot is fine and sample size is large
ggplot(diagnostics_phgr, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw() +
  facet_grid(lifestyle~ph)

# Signficance testing
Anova(mod_phgr)

# Fit reduced models for parametric bootstrapping to verify Anova() output
mod_phgr_red_lifestyle <- lmer(gr ~ lifestyle + (1 |strain), data = df_tidy_phgr)
mod_phgr_red_ph <- lmer(gr ~ ph + (1 | strain), data = df_tidy_phgr)
mod_phgr_red_int <- lmer(gr ~ lifestyle + ph + (1 | strain), data = df_tidy_phgr)

# Run simulations
PBmodcomp(mod_phgr, mod_phgr_red_lifestyle, nsim = 1000)
PBmodcomp(mod_phgr, mod_phgr_red_ph, nsim = 1000)
PBmodcomp(mod_phgr, mod_phgr_red_int, nsim = 1000)

# Emmeans
# Comparison of means
phgr_means <- emmeans(mod_phgr , specs = ~ ph * lifestyle)

# This gives you the estimated trends, their SEs and confidence limits
phgr_means

# Here's a helpful plot to visualise
plot(phgr_means)

# Pairwise comparisons
pairs(phgr_means)

### Correlate to sugar
cor_df <- df_tidy_phgr %>%
  group_by(strain, ph, lifestyle) %>%
  summarise(mean_gr = mean(gr),
            se_gr = sd(gr)/sqrt(n())) %>%
  left_join(read_csv("Data/df_tidy_sugar.csv") %>%
              pivot_wider(names_from = "sugar", values_from = "release") %>%
              clean_names() %>%
              group_by(strain, ph, lifestyle) %>%
              summarise(mean_glu = mean(glucose),
                        se_glu = sd(glucose)/sqrt(n()),
                        mean_mal = mean(maltose),
                        se_mal = sd(maltose)/sqrt(n()),
                        mean_tre = mean(trehalose),
                        se_tre = sd(trehalose)/sqrt(n())), by = c("strain", "ph", "lifestyle"))

cor_df_5.5 <- cor_df %>% filter(ph == 5.5)
cor_df_7.6 <- cor_df %>% filter(ph == 7.6)

# Make a plot df, add rank-transformed variables within each facet
plot_df <- cor_df %>%
  pivot_longer(cols = matches("^(mean|se)_(glu|mal|tre)$"),  # only sugar columns starting with mean_ or se_
               names_to = c(".value", "sugar"),       # .value says part of name becomes the column name
               names_sep = "_"                         # split at the underscore
  ) %>%
  dplyr::filter(sugar != "tre") %>% # Remove trehalose
  mutate(sugar = case_when(sugar == "mal" ~ "Maltose",
                           sugar == "glu" ~ "Glucose")) 

# Plot
plot_df %>%
ggplot(aes(x = mean_gr, y = mean)) +
  geom_point(aes(color = lifestyle, shape = factor(ph))) +
  # geom_smooth(method = "gam", se = FALSE, colour = "black") + # Lazy fit across all points
  geom_smooth(aes(linetype = factor(ph)), method = "lm", se = FALSE, colour = "black", linewidth = 0.5) + # Lazy fit by ph
  facet_wrap(~sugar) +
  theme_bw() +
  labs(x = "Mean growth rate", y = "Mean sugar release", colour = "Lifestyle", shape = "pH", linetype = "pH") +
  #scale_y_continuous(limits = c(0, max(cor_df$mean_mal))) +
    scale_colour_manual(values = c(
      "FL" = "green3",
      "S" = "cornflowerblue"),
      labels = c(
        "FL" = "Free-living",
        "S" = "Symbiotic")) +
  scale_shape_manual(values = c(
    "5.5" = 16,
    "7.6" = 1)) +
    theme(legend.position = "bottom")

# Spearman tests for the correlation (don't mind the warning, the p is fine)
cor.test(cor_df_5.5$mean_gr, cor_df_5.5$mean_mal, method = "spearman") #sig
cor.test(cor_df_5.5$mean_gr, cor_df_5.5$mean_glu, method = "spearman") #sig
cor.test(cor_df_5.5$mean_gr, cor_df_5.5$mean_tre, method = "spearman")  #ns

cor.test(cor_df_7.6$mean_gr, cor_df_7.6$mean_mal, method = "spearman") #sig
cor.test(cor_df_7.6$mean_gr, cor_df_7.6$mean_glu, method = "spearman") #sig
cor.test(cor_df_7.6$mean_gr, cor_df_7.6$mean_tre, method = "spearman") #ns

# Reanalysis using only the 86 strains from the phylogeny to check robustness of results ###########################

# Plot
df_tidy_ph %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  group_by(strain, day, lifestyle, ph) %>%
  summarise(mean = mean(od750),
            se = sd(od750)/sqrt(n())) %>%
  ggplot(aes(x = day, y = mean, colour = factor(ph), groups = factor(ph))) +
  geom_point(data = df_tidy_ph %>%   filter(strain != "155",
                                            strain != "192",
                                            strain != "242",
                                            strain != "252",
                                            strain != "258",
                                            strain != "s354") %>% mutate(mean = od750), alpha = 0.2) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
  geom_line() +
  facet_wrap(~ strain)

# Calculate growth rate
# A function for growth rate
gr_func <- function(x, y, t) {
  gr <- (log(x/y))/t
  return(gr)
}

# Add to data frame
df_tidy_phgr_86 <- df_tidy_ph %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  pivot_wider(names_from = day, names_prefix = "day", values_from = od750) %>%
  mutate(gr = gr_func(day7, day0, 7))

# Plot
df_tidy_phgr_86 %>%
  group_by(lifestyle, ph) %>%
  summarise(mean = mean(gr),
            se = sd(gr)/sqrt(n())) %>%
  ggplot(aes(x = factor(ph), y = mean, colour = lifestyle, group = lifestyle)) +
  geom_jitter(data = df_tidy_phgr_86 %>% mutate(mean = gr), alpha = 0.2, width = 0.125) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, linewidth = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "pH", y = "Mean daily growth rate", colour = "Lifestyle")

# Distribution -- we have a bit of biomodality going on in one group
df_tidy_phgr_86 %>%
  ggplot(aes(x = gr, fill = lifestyle)) +
  geom_histogram() +
  facet_grid(lifestyle~ph)


# Fit ME model
mod_phgr_86 <- lmer(gr ~ ph*lifestyle + (1 | strain), data = df_tidy_phgr_86)

# Make a data frame for diagnostics
diagnostics_phgr_86 <- broom.mixed::augment(mod_phgr_86)

# QQ plot -- fine
ggplot(diagnostics_phgr_86, aes(sample = .resid)) +
  geom_qq() + 
  geom_qq_line() +
  theme_bw()

# Residuals vs. fitted values -- some patterns here, but model should be robust given QQ plot is fine and sample size is large
ggplot(diagnostics_phgr_86, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw() +
  facet_grid(lifestyle~ph)

# Signficance testing
Anova(mod_phgr_86)

# Fit reduced models for parametric bootstrapping to verify Anova() output
mod_phgr_red_lifestyle_86 <- lmer(gr ~ lifestyle + (1 |strain), data = df_tidy_phgr_86)
mod_phgr_red_ph_86 <- lmer(gr ~ ph + (1 | strain), data = df_tidy_phgr_86)
mod_phgr_red_int_86 <- lmer(gr ~ lifestyle + ph + (1 | strain), data = df_tidy_phgr_86)

# Run simulations
PBmodcomp(mod_phgr_86, mod_phgr_red_lifestyle_86, nsim = 1000)
PBmodcomp(mod_phgr_86, mod_phgr_red_ph_86, nsim = 1000)
PBmodcomp(mod_phgr_86, mod_phgr_red_int_86, nsim = 1000)

# Emmeans
# Comparison of means
phgr_means_86 <- emmeans(mod_phgr_86, specs = ~ ph * lifestyle)

# This gives you the estimated trends, their SEs and confidence limits
phgr_means_86

# Here's a helpful plot to visualise
plot(phgr_means_86)

# Pairwise comparisons
pairs(phgr_means_86)

### Correlate to sugar
cor_df_86 <- df_tidy_phgr_86 %>%
  group_by(strain, ph, lifestyle) %>%
  summarise(mean_gr = mean(gr),
            se_gr = sd(gr)/sqrt(n())) %>%
  left_join(read_csv("Data/df_tidy_sugar.csv") %>%
              pivot_wider(names_from = "sugar", values_from = "release") %>%
              clean_names() %>%
              group_by(strain, ph, lifestyle) %>%
              summarise(mean_glu = mean(glucose),
                        se_glu = sd(glucose)/sqrt(n()),
                        mean_mal = mean(maltose),
                        se_mal = sd(maltose)/sqrt(n()),
                        mean_tre = mean(trehalose),
                        se_tre = sd(trehalose)/sqrt(n())), by = c("strain", "ph", "lifestyle"))

cor_df_5.5_86 <- cor_df_86 %>% filter(ph == 5.5)
cor_df_7.6_86 <- cor_df_86 %>% filter(ph == 7.6)



# Make a plot df, add rank-transformed variables within each facet
plot_df_86 <- cor_df_86 %>%
  pivot_longer(cols = matches("^(mean|se)_(glu|mal|tre)$"),  # only sugar columns starting with mean_ or se_
               names_to = c(".value", "sugar"),       # .value says part of name becomes the column name
               names_sep = "_"                         # split at the underscore
  ) %>%
  dplyr::filter(sugar != "tre") %>% # Remove trehalose
  mutate(sugar = case_when(sugar == "mal" ~ "Maltose",
                           sugar == "glu" ~ "Glucose")) 

# Plot
plot_df_86 %>%
  ggplot(aes(x = mean_gr, y = mean)) +
  geom_point(aes(color = lifestyle, shape = factor(ph))) +
  # geom_smooth(method = "gam", se = FALSE, colour = "black") + # Lazy fit across all points
  geom_smooth(aes(linetype = factor(ph)), method = "lm", se = FALSE, colour = "black", linewidth = 0.5) + # Lazy fit by ph
  facet_wrap(~sugar) +
  theme_bw() +
  labs(x = "Mean growth rate", y = "Mean sugar release", colour = "Lifestyle", shape = "pH", linetype = "pH") +
  #scale_y_continuous(limits = c(0, max(cor_df$mean_mal))) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  scale_shape_manual(values = c(
    "5.5" = 16,
    "7.6" = 1)) +
  theme(legend.position = "bottom")

# Spearman tests for the correlation (don't mind the warning, the p is fine)
cor.test(cor_df_5.5_86$mean_gr, cor_df_5.5_86$mean_mal, method = "spearman") #sig
cor.test(cor_df_5.5_86$mean_gr, cor_df_5.5_86$mean_glu, method = "spearman") #sig
cor.test(cor_df_5.5_86$mean_gr, cor_df_5.5_86$mean_tre, method = "spearman")  #ns

cor.test(cor_df_7.6_86$mean_gr, cor_df_7.6_86$mean_mal, method = "spearman") #sig
cor.test(cor_df_7.6_86$mean_gr, cor_df_7.6_86$mean_glu, method = "spearman") #sig
cor.test(cor_df_7.6_86$mean_gr, cor_df_7.6_86$mean_tre, method = "spearman") #ns

# No changed conclusions

# Reanalysis using only the 92 strains from the phylogeny to check robustness of results ###########################

# Plot
df_tidy_ph %>%
  filter(strain != "192") %>%
  group_by(strain, day, lifestyle, ph) %>%
  summarise(mean = mean(od750),
            se = sd(od750)/sqrt(n())) %>%
  ggplot(aes(x = day, y = mean, colour = factor(ph), groups = factor(ph))) +
  geom_point(data = df_tidy_ph %>% filter(strain != "192") %>% mutate(mean = od750), alpha = 0.2) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
  geom_line() +
  facet_wrap(~ strain)

# Calculate growth rate

# A function for growth rate
gr_func <- function(x, y, t) {
  gr <- (log(x/y))/t
  return(gr)
}

# Add to data frame
df_tidy_phgr_92 <- df_tidy_ph %>%
  filter(strain != "192") %>%
  pivot_wider(names_from = day, names_prefix = "day", values_from = od750) %>%
  mutate(gr = gr_func(day7, day0, 7))

# Plot
df_tidy_phgr_92 %>%
  group_by(lifestyle, ph) %>%
  summarise(mean = mean(gr),
            se = sd(gr)/sqrt(n())) %>%
  ggplot(aes(x = factor(ph), y = mean, colour = lifestyle, group = lifestyle)) +
  geom_jitter(data = df_tidy_phgr_92 %>% mutate(mean = gr), alpha = 0.2, width = 0.125) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, linewidth = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "pH", y = "Mean daily growth rate", colour = "Lifestyle")

# Distribution -- we have a bit of biomodality going on in one group
df_tidy_phgr_92 %>%
  ggplot(aes(x = gr, fill = lifestyle)) +
  geom_histogram() +
  facet_grid(lifestyle~ph)


# Fit ME model
mod_phgr_92 <- lmer(gr ~ ph*lifestyle + (1 | strain), data = df_tidy_phgr_92)

# Make a data frame for diagnostics
diagnostics_phgr_92 <- broom.mixed::augment(mod_phgr_92)

# QQ plot -- fine
ggplot(diagnostics_phgr_92, aes(sample = .resid)) +
  geom_qq() + 
  geom_qq_line() +
  theme_bw()

# Residuals vs. fitted values -- some patterns here, but model should be robust given QQ plot is fine and sample size is large
ggplot(diagnostics_phgr_92, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw() +
  facet_grid(lifestyle~ph)

# Signficance testing
Anova(mod_phgr_92)

# Fit reduced models for parametric bootstrapping to verify Anova() output
mod_phgr_red_lifestyle_92 <- lmer(gr ~ lifestyle + (1 |strain), data = df_tidy_phgr_92)
mod_phgr_red_ph_92 <- lmer(gr ~ ph + (1 | strain), data = df_tidy_phgr_92)
mod_phgr_red_int_92 <- lmer(gr ~ lifestyle + ph + (1 | strain), data = df_tidy_phgr_92)

# Run simulations
PBmodcomp(mod_phgr_92, mod_phgr_red_lifestyle_92, nsim = 1000)
PBmodcomp(mod_phgr_92, mod_phgr_red_ph_92, nsim = 1000)
PBmodcomp(mod_phgr_92, mod_phgr_red_int_92, nsim = 1000)

# Emmeans
# Comparison of means
phgr_means_92 <- emmeans(mod_phgr_92, specs = ~ ph * lifestyle)

# This gives you the estimated trends, their SEs and confidence limits
phgr_means_92

# Here's a helpful plot to visualise
plot(phgr_means_92)

# Pairwise comparisons
pairs(phgr_means_92)

### Correlate to sugar
cor_df_92 <- df_tidy_phgr_92 %>%
  group_by(strain, ph, lifestyle) %>%
  summarise(mean_gr = mean(gr),
            se_gr = sd(gr)/sqrt(n())) %>%
  left_join(read_csv("Data/df_tidy_sugar.csv") %>%
              pivot_wider(names_from = "sugar", values_from = "release") %>%
              clean_names() %>%
              group_by(strain, ph, lifestyle) %>%
              summarise(mean_glu = mean(glucose),
                        se_glu = sd(glucose)/sqrt(n()),
                        mean_mal = mean(maltose),
                        se_mal = sd(maltose)/sqrt(n()),
                        mean_tre = mean(trehalose),
                        se_tre = sd(trehalose)/sqrt(n())), by = c("strain", "ph", "lifestyle"))

cor_df_5.5_92 <- cor_df_92 %>% filter(ph == 5.5)
cor_df_7.6_92 <- cor_df_92 %>% filter(ph == 7.6)

# Make a plot df, add rank-transformed variables within each facet
plot_df_92 <- cor_df_92 %>%
  pivot_longer(cols = matches("^(mean|se)_(glu|mal|tre)$"),  # only sugar columns starting with mean_ or se_
               names_to = c(".value", "sugar"),       # .value says part of name becomes the column name
               names_sep = "_"                         # split at the underscore
  ) %>%
  dplyr::filter(sugar != "tre") %>% # Remove trehalose
  mutate(sugar = case_when(sugar == "mal" ~ "Maltose",
                           sugar == "glu" ~ "Glucose")) 

# Plot
plot_df_92 %>%
  ggplot(aes(x = mean_gr, y = mean)) +
  geom_point(aes(color = lifestyle, shape = factor(ph))) +
  # geom_smooth(method = "gam", se = FALSE, colour = "black") + # Lazy fit across all points
  geom_smooth(aes(linetype = factor(ph)), method = "lm", se = FALSE, colour = "black", linewidth = 0.5) + # Lazy fit by ph
  facet_wrap(~sugar) +
  theme_bw() +
  labs(x = "Mean growth rate", y = "Mean sugar release", colour = "Lifestyle", shape = "pH", linetype = "pH") +
  #scale_y_continuous(limits = c(0, max(cor_df$mean_mal))) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")) +
  scale_shape_manual(values = c(
    "5.5" = 16,
    "7.6" = 1)) +
  theme(legend.position = "bottom")

# Spearman tests for the correlation (don't mind the warning, the p is fine)
cor.test(cor_df_5.5_92$mean_gr, cor_df_5.5_92$mean_mal, method = "spearman") #sig
cor.test(cor_df_5.5_92$mean_gr, cor_df_5.5_92$mean_glu, method = "spearman") #sig
cor.test(cor_df_5.5_92$mean_gr, cor_df_5.5_92$mean_tre, method = "spearman")  #ns

cor.test(cor_df_7.6_92$mean_gr, cor_df_7.6_92$mean_mal, method = "spearman") #sig
cor.test(cor_df_7.6_92$mean_gr, cor_df_7.6_92$mean_glu, method = "spearman") #sig
cor.test(cor_df_7.6_92$mean_gr, cor_df_7.6_92$mean_tre, method = "spearman") #ns

# No changed conclusions
