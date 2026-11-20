# Erika M Hansson 2025
# erika.hansson@manchester.ac.uk

#### ---------------------------------------------------- ####
# 1) Plotting N data
# 2) Fitting a model to N data for significance testing
# Then: repeating all of the above with the 86 strains in the phylogeny for robustness
#### ---------------------------------------------------- ####

# Packages
if(!require("tidyverse")) {install.packages("tidyverse")} # wrangling, plotting
if(!require("ggpattern")) {install.packages("ggpattern")} # plotting patterns
if(!require("car")) {install.packages("car")} # stats tools
if(!require("broom")) {install.packages("broom")} # for wrangling
if(!require("lme4")) {install.packages("lme4")} # for fitting mixed effects models
if(!require("e1017")) {install.packages("e1017")} # for handy skewness function
if(!require("pbkrtest")) {install.packages("pbkrtest")} # for parametric bootstrapping of mixed effects models
if(!require("emmeans")) {install.packages("emmeans")} # for fitting mixed effects models
if(!require("multcompView")) {install.packages("multcompView")} # for post hoc testing
if(!require("multcomp")) {install.packages("multcomp")} # for post hoc testing
if(!require("ggrepel")) {install.packages("ggrepel")} # plot labelling
if(!require("vegan")) {install.packages("vegan")} # ordination methods etc.
if(!require("randomForest")) {install.packages("randomForest")} # random forest modelling
if(!require("igraph")) {install.packages("igraph")} # network analysis
if(!require("Hmisc")) {install.packages("Hmisc")} # for correlation matrix with p-values
if(!require("corrplot")) {install.packages("corrplot")} # for correlation matrix with p-values

# Load the data
df_tidy_n <- read_csv("Data/df_tidy_n.csv") %>%
  # Pivoting the data frame and calculating growth rate
pivot_wider(values_from = reading, names_from = c(day), names_prefix = "day") %>%
  mutate(gr = log(day7/day0)/7) %>%
  dplyr::select(-day0, -day7)

# Look at it
glimpse(df_tidy_n) 

df_tidy_n %>% filter(n_source == "Arginine", lifestyle == "S") %>%
  print(n = nrow(.))

## A summary plot showing growth in each medium
df_tidy_n %>%
  group_by(n_source, lifestyle, n_type) %>%
  summarise(
    mean = mean(gr),
    se = sd(gr)/sqrt(n())
  ) %>%
  ungroup() %>%
  ggplot(aes(x = lifestyle, y = mean, colour = lifestyle)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_point(data = df_tidy_n %>% mutate(mean = gr), alpha = 0.3) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.5) +
  theme_bw() +
  facet_wrap(n_type~n_source, nrow = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")
  ) +
  labs(
    y = "Mean growth rate per day",
    x = "Lifestyle") +
  guides(colour = "none")

## Corresponding histogram to highlight how symbiotic strains have bimodal distribution for some N sources
df_tidy_n %>%
  group_by(strain, n_source, lifestyle, n_type) %>%
  summarise(
    gr = mean(gr)
  ) %>%
  ungroup() %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                          lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = lifestyle)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram() +
  facet_grid(lifestyle~n_source) +
  scale_fill_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

df_tidy_n$n_source_order <- with(df_tidy_n, factor(n_source, levels = unique(n_source[order(n_type)])))

# Changing the the faceting order -- Irma's preferred plot!
df_tidy_n %>%
  group_by(strain, n_source, lifestyle, n_type, n_source_order) %>%
  summarise(
    gr = mean(gr)
  ) %>%
  ungroup() %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                               lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = n_type)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram() +
  facet_grid(n_source_order ~ lifestyle) +
  scale_fill_manual(values = c(
    "Inorganic" = "mediumorchid",
    "Organic" = "#FF8C00")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        theme(strip.text.y.right = element_blank()))

# With lifestyle colouring and pattern for N type
df_tidy_n %>%
  group_by(strain, n_source, lifestyle, n_type, n_source_order) %>%
  summarise(
    gr = mean(gr)
  ) %>%
  ungroup() %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                               lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = lifestyle, pattern = n_type)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram_pattern(
    bins = 30,
    pattern_density = 0.03,
    colour = "white",
    pattern_fill = "white"
  ) +
  facet_grid(n_type*n_source ~ lifestyle) +
  scale_fill_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none", pattern = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

df_tidy_n %>%
  ggplot(aes(x = strain, y = gr)) +
  geom_point() +
  geom_hline(yintercept = 0)

#### Fitting a model

# A bit of a left-skew
df_tidy_n %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                          lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = lifestyle)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram() +
  scale_fill_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Let's look at the skewness -- mild/moderate, totally fine for this model as long as the diagnostic plots behave
#skewness(df_tidy_n$gr) 

# Fit the model
mod_n <- lmer(gr ~ lifestyle*n_source + (1 | strain), data = df_tidy_n)

# Make a data frame for diagnostics
diagnostics_n <- broom.mixed::augment(mod_n)

# QQ plot -- fine
ggplot(diagnostics_n, aes(sample = .resid)) +
  geom_qq() + 
  geom_qq_line() +
  theme_bw()

# Residuals vs. fitted values -- fine
ggplot(diagnostics_n, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw()

# Summary
summary(mod_n)

# Significance testing
Anova(mod_n)

# Fit reduced models for parametric bootstrapping to verify Anova() output
mod_n_red_lifestyle <- lmer(gr ~ lifestyle + (1 |strain), data = df_tidy_n)
mod_n_red_source <- lmer(gr ~ n_source + (1 | strain), data = df_tidy_n)
mod_n_red_int <- lmer(gr ~ lifestyle + n_source + (1 | strain), data = df_tidy_n)

# Run simulations
PBmodcomp(mod_n, mod_n_red_lifestyle, nsim = 1000)
PBmodcomp(mod_n, mod_n_red_source, nsim = 1000)
PBmodcomp(mod_n, mod_n_red_int, nsim = 1000)

# Same result here

## Post hoc testing

# Let's compare all means to each other in a fully factorial way
all_n_means <- emmeans(mod_n, specs = ~ lifestyle * n_source)

# This gives you the estimated means, their SEs and confidence limits
all_n_means

# Here's a helpful plot to visualise
plot(all_n_means)

# And this essentially applies a Tukey test across all combinations
pairs(all_n_means)

# Let's make letter groupings 
get_intercept_letters <- function(model, model_name, group_vars = c("lifestyle", "n_source")) {
  # Step 1: Estimate marginal means (intercepts)
  emms <- emmeans(model, specs = group_vars)
  
  # Step 2: Pairwise comparisons using Tukey adjustment
  pairs_emms <- contrast(emms, method = "pairwise", adjust = "tukey")
  summary_emms <- summary(pairs_emms)
  
  # Step 3: Build p-value matrix
  pval_vec <- setNames(summary_emms$p.value, summary_emms$contrast)
  split_names <- strsplit(names(pval_vec), " - ")
  levels <- unique(unlist(split_names))
  pmat <- matrix(1, nrow = length(levels), ncol = length(levels),
                 dimnames = list(levels, levels))
  
  for (i in seq_along(pval_vec)) {
    comps <- unlist(split_names[[i]])
    pmat[comps[1], comps[2]] <- pval_vec[i]
    pmat[comps[2], comps[1]] <- pval_vec[i]
  }
  
  # Step 4: Generate letter groupings
  letters <- multcompLetters(pmat, threshold = 0.05)$Letters
  
  # Step 5: Return as named list
  out <- list()
  out[[model_name]] <- letters
  return(out)
}

# List of models
model_list <- list(
  mod_n
)

# Apply get_slope_letters to each model in the list
intercept_letters_list <- mapply(
  FUN = get_intercept_letters,
  model = model_list,
  model_name = paste0("Model_", seq_along(model_list)),
  SIMPLIFY = FALSE
)

# Combine results into a single named list
intercept_letters_all <- do.call(c, intercept_letters_list)

# Inspect results
as.data.frame(intercept_letters_all)

# Reanalysis using only the 86 strains from the phylogeny to check robustness of results ###########################
# Load the data
df_tidy_n_86 <- read_csv("Data/df_tidy_n.csv") %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  # Pivoting the data frame and calculating growth rate
  pivot_wider(values_from = reading, names_from = c(day), names_prefix = "day") %>%
  mutate(gr = log(day7/day0)/7) %>%
  dplyr::select(-day0, -day7)

# Look at it
glimpse(df_tidy_n_86) 

## A summary plot showing growth in each medium
df_tidy_n_86 %>%
  group_by(n_source, lifestyle, n_type) %>%
  summarise(
    mean = mean(gr),
    se = sd(gr)/sqrt(n())
  ) %>%
  ggplot(aes(x = lifestyle, y = mean, colour = lifestyle)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_point(data = df_tidy_n %>% mutate(mean = gr), alpha = 0.3) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.5) +
  theme_bw() +
  facet_wrap(n_type~n_source, nrow = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")
  ) +
  labs(
    y = "Mean growth rate per day",
    x = "Lifestyle") +
  guides(colour = "none")

## Corresponding histogram to highlight how symbiotic strains have bimodal distribution for some N sources
df_tidy_n_86 %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                               lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = lifestyle)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram() +
  facet_grid(lifestyle~n_source) +
  scale_fill_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

df_tidy_n_86$n_source_order <- with(df_tidy_n_86, factor(n_source, levels = unique(n_source[order(n_type)])))

# Changing the the faceting order -- Irma's preferred plot!
df_tidy_n_86 %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                               lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = n_type)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram() +
  facet_grid(n_source_order ~ lifestyle) +
  scale_fill_manual(values = c(
    "Inorganic" = "mediumorchid",
    "Organic" = "#FF8C00")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        theme(strip.text.y.right = element_blank()))

# With lifestyle colouring and pattern for N type
df_tidy_n_86 %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                               lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = lifestyle, pattern = n_type)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram_pattern(
    bins = 30,
    pattern_density = 0.03,
    colour = "white",
    pattern_fill = "white"
  ) +
  facet_grid(n_type*n_source ~ lifestyle) +
  scale_fill_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none", pattern = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#### Fitting a model

# A bit of a left-skew
df_tidy_n_86 %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                               lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = lifestyle)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram() +
  scale_fill_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Let's look at the skewness -- mild/moderate, totally fine for this model as long as the diagnostic plots behave
#skewness(df_tidy_n$gr) 

# Fit the model
mod_n_86 <- lmer(gr ~ lifestyle*n_source + (1 | strain), data = df_tidy_n_86)

# Make a data frame for diagnostics
diagnostics_n_86 <- broom.mixed::augment(mod_n_86)

# QQ plot -- fine
ggplot(diagnostics_n_86, aes(sample = .resid)) +
  geom_qq() + 
  geom_qq_line() +
  theme_bw()

# Residuals vs. fitted values -- fine
ggplot(diagnostics_n_86, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw()

# Summary
summary(mod_n_86)

# Significance testing
Anova(mod_n_86)

# Fit reduced models for parametric bootstrapping to verify Anova() output
mod_n_red_lifestyle_86 <- lmer(gr ~ lifestyle + (1 |strain), data = df_tidy_n_86)
mod_n_red_source_86 <- lmer(gr ~ n_source + (1 | strain), data = df_tidy_n_86)
mod_n_red_int_86 <- lmer(gr ~ lifestyle + n_source + (1 | strain), data = df_tidy_n_86)

# Run simulations
PBmodcomp(mod_n_86, mod_n_red_lifestyle_86, nsim = 1000)
PBmodcomp(mod_n_86, mod_n_red_source_86, nsim = 1000)
PBmodcomp(mod_n_86, mod_n_red_int_86, nsim = 1000)

# Same result here

## Post hoc testing

# Let's compare all means to each other in a fully factorial way
all_n_means_86 <- emmeans(mod_n_86, specs = ~ lifestyle * n_source)

# This gives you the estimated means, their SEs and confidence limits
all_n_means_86

# Here's a helpful plot to visualise
plot(all_n_means_86)

# And this essentially applies a Tukey test across all combinations
pairs(all_n_means_86)

# Let's make letter groupings 
get_intercept_letters <- function(model, model_name, group_vars = c("lifestyle", "n_source")) {
  # Step 1: Estimate marginal means (intercepts)
  emms <- emmeans(model, specs = group_vars)
  
  # Step 2: Pairwise comparisons using Tukey adjustment
  pairs_emms <- contrast(emms, method = "pairwise", adjust = "tukey")
  summary_emms <- summary(pairs_emms)
  
  # Step 3: Build p-value matrix
  pval_vec <- setNames(summary_emms$p.value, summary_emms$contrast)
  split_names <- strsplit(names(pval_vec), " - ")
  levels <- unique(unlist(split_names))
  pmat <- matrix(1, nrow = length(levels), ncol = length(levels),
                 dimnames = list(levels, levels))
  
  for (i in seq_along(pval_vec)) {
    comps <- unlist(split_names[[i]])
    pmat[comps[1], comps[2]] <- pval_vec[i]
    pmat[comps[2], comps[1]] <- pval_vec[i]
  }
  
  # Step 4: Generate letter groupings
  letters <- multcompLetters(pmat, threshold = 0.05)$Letters
  
  # Step 5: Return as named list
  out <- list()
  out[[model_name]] <- letters
  return(out)
}

# List of models
model_list <- list(
  mod_n_86
)

# Apply get_slope_letters to each model in the list
intercept_letters_list_86 <- mapply(
  FUN = get_intercept_letters,
  model = model_list,
  model_name = paste0("Model_", seq_along(model_list)),
  SIMPLIFY = FALSE
)

# Combine results into a single named list
intercept_letters_all_86 <- do.call(c, intercept_letters_list_86)

# Inspect results
as.data.frame(intercept_letters_all_86)

# Conclusions unchanged


# Reanalysis using only the 92 strains from the phylogeny to check robustness of results ###########################
# Load the data
df_tidy_n_92 <- read_csv("Data/df_tidy_n.csv") %>%
  filter(strain != "192") %>%
  # Pivoting the data frame and calculating growth rate
  pivot_wider(values_from = reading, names_from = c(day), names_prefix = "day") %>%
  mutate(gr = log(day7/day0)/7) %>%
  dplyr::select(-day0, -day7)

# Look at it
glimpse(df_tidy_n_92) 

## A summary plot showing growth in each medium
df_tidy_n_92 %>%
  group_by(n_source, lifestyle, n_type) %>%
  summarise(
    mean = mean(gr),
    se = sd(gr)/sqrt(n())
  ) %>%
  ggplot(aes(x = lifestyle, y = mean, colour = lifestyle)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_point(data = df_tidy_n %>% mutate(mean = gr), alpha = 0.3) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.5) +
  theme_bw() +
  facet_wrap(n_type~n_source, nrow = 1) +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic")
  ) +
  labs(
    y = "Mean growth rate per day",
    x = "Lifestyle") +
  guides(colour = "none")

## Corresponding histogram to highlight how symbiotic strains have bimodal distribution for some N sources
df_tidy_n_92 %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                               lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = lifestyle)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram() +
  facet_grid(lifestyle~n_source) +
  scale_fill_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

df_tidy_n_92$n_source_order <- with(df_tidy_n_92, factor(n_source, levels = unique(n_source[order(n_type)])))

# Changing the the faceting order -- Irma's preferred plot!
df_tidy_n_92 %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                               lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = n_type)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram() +
  facet_grid(n_source_order ~ lifestyle) +
  scale_fill_manual(values = c(
    "Inorganic" = "mediumorchid",
    "Organic" = "#FF8C00")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        theme(strip.text.y.right = element_blank()))

# With lifestyle colouring and pattern for N type
df_tidy_n_92 %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                               lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = lifestyle, pattern = n_type)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram_pattern(
    bins = 30,
    pattern_density = 0.03,
    colour = "white",
    pattern_fill = "white"
  ) +
  facet_grid(n_type*n_source ~ lifestyle) +
  scale_fill_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none", pattern = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#### Fitting a model

# A bit of a left-skew
df_tidy_n_92 %>%
  mutate(lifestyle = case_when(lifestyle =="FL" ~ "Free-living",
                               lifestyle == "S" ~ "Symbiotic")) %>%
  ggplot(aes(x = gr, fill = lifestyle)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_histogram() +
  scale_fill_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue")
  ) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Let's look at the skewness -- mild/moderate, totally fine for this model as long as the diagnostic plots behave
#skewness(df_tidy_n$gr) 

# Fit the model
mod_n_92 <- lmer(gr ~ lifestyle*n_source + (1 | strain), data = df_tidy_n_92)

# Make a data frame for diagnostics
diagnostics_n_92 <- broom.mixed::augment(mod_n_92)

# QQ plot -- fine
ggplot(diagnostics_n_92, aes(sample = .resid)) +
  geom_qq() + 
  geom_qq_line() +
  theme_bw()

# Residuals vs. fitted values -- fine
ggplot(diagnostics_n_92, aes(x = .fitted, y = .resid)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw()

# Summary
summary(mod_n_92)

# Significance testing
Anova(mod_n_92)

# Fit reduced models for parametric bootstrapping to verify Anova() output
mod_n_red_lifestyle_92 <- lmer(gr ~ lifestyle + (1 |strain), data = df_tidy_n_92)
mod_n_red_source_92 <- lmer(gr ~ n_source + (1 | strain), data = df_tidy_n_92)
mod_n_red_int_92 <- lmer(gr ~ lifestyle + n_source + (1 | strain), data = df_tidy_n_92)

# Run simulations
PBmodcomp(mod_n_92, mod_n_red_lifestyle_92, nsim = 1000)
PBmodcomp(mod_n_92, mod_n_red_source_92, nsim = 1000)
PBmodcomp(mod_n_92, mod_n_red_int_92, nsim = 1000)

# Same result here

## Post hoc testing

# Let's compare all means to each other in a fully factorial way
all_n_means_92 <- emmeans(mod_n_92, specs = ~ lifestyle * n_source)

# This gives you the estimated means, their SEs and confidence limits
all_n_means_92

# Here's a helpful plot to visualise
plot(all_n_means_92)

# And this essentially applies a Tukey test across all combinations
pairs(all_n_means_92)

# Let's make letter groupings 
get_intercept_letters <- function(model, model_name, group_vars = c("lifestyle", "n_source")) {
  # Step 1: Estimate marginal means (intercepts)
  emms <- emmeans(model, specs = group_vars)
  
  # Step 2: Pairwise comparisons using Tukey adjustment
  pairs_emms <- contrast(emms, method = "pairwise", adjust = "tukey")
  summary_emms <- summary(pairs_emms)
  
  # Step 3: Build p-value matrix
  pval_vec <- setNames(summary_emms$p.value, summary_emms$contrast)
  split_names <- strsplit(names(pval_vec), " - ")
  levels <- unique(unlist(split_names))
  pmat <- matrix(1, nrow = length(levels), ncol = length(levels),
                 dimnames = list(levels, levels))
  
  for (i in seq_along(pval_vec)) {
    comps <- unlist(split_names[[i]])
    pmat[comps[1], comps[2]] <- pval_vec[i]
    pmat[comps[2], comps[1]] <- pval_vec[i]
  }
  
  # Step 4: Generate letter groupings
  letters <- multcompLetters(pmat, threshold = 0.05)$Letters
  
  # Step 5: Return as named list
  out <- list()
  out[[model_name]] <- letters
  return(out)
}

# List of models
model_list <- list(
  mod_n_92
)

# Apply get_slope_letters to each model in the list
intercept_letters_list_92 <- mapply(
  FUN = get_intercept_letters,
  model = model_list,
  model_name = paste0("Model_", seq_along(model_list)),
  SIMPLIFY = FALSE
)

# Combine results into a single named list
intercept_letters_all_92 <- do.call(c, intercept_letters_list_92)

# Inspect results
as.data.frame(intercept_letters_all_92)

# Conclusions unchanged

# More on N -- looking at multivariate space and lost pathways
mult_df_n <- 
  read_csv("Data/df_tidy_n.csv") %>%
  dplyr::select(-n_type) %>%
  pivot_wider(values_from = reading, names_from = c(day), names_prefix = "day") %>%
  mutate(gr = log(day7/day0)/7) %>%
  dplyr::select(-day0, -day7) %>%
  pivot_wider(values_from = "gr", names_from = "n_source") %>%
  group_by(strain, lifestyle) %>%
  summarise(
    Arginine = mean(Arginine),
    Asparagine = mean(Asparagine),
    Glutamine = mean(Glutamine),
    Glutamate = mean(Glutamate),
    Serine = mean(Serine),
    Ammonium = mean(Ammonium),
    Nitrate = mean(Nitrate),
    Urea = mean(Urea)
  )

glimpse(mult_df_n)

# Chaaracterise pathway as lost if average growth is < 0 for a strain
binary_n <- mult_df_n %>%
  mutate(
    across(
      Arginine:Urea,
      ~ if_else(. > 0, 1L, 0L)
    ),
    n_lost = rowSums(across(Arginine:Urea) == 0)) 

# Let's just look at Arginine loss
binary_n %>% filter(Arginine == 0)

mult_df_scaled_n <- mult_df_n %>%
  ungroup() %>%
  dplyr::select(-strain, -lifestyle) %>%
  drop_na() %>%                # remove rows with missing values
  scale()                      # standardize to mean 0, sd 1

# Run PCA
pca_result_n <- prcomp(mult_df_scaled_n, center = FALSE, scale. = FALSE)

# View summary of explained variance
summary(pca_result_n)

# Scree plot to visualize variance explained
plot(pca_result_n, type = "l")

# Tidy version of PCA result
pca_df_n <- as_tibble(pca_result_n$x)

# Add metadata for plotting
pca_df_n <- bind_cols(mult_df_n %>% drop_na() %>% dplyr::select(strain, lifestyle), pca_df_n)

# Plot first two PCs
ggplot(pca_df_n, aes(x = PC1, y = PC2, colour = lifestyle, label = strain)) +
  labs(x = "PC1", y = "PC2") +
  geom_label() +
  theme_bw()

# Get loadings
pca_result_n$rotation

# Get scores for each PC by type
type_scores_n <- augment(pca_result_n) %>%
  bind_cols(mult_df_n) %>%
  dplyr::select(2:4, strain, lifestyle) %>%
  unique() %>%
  left_join(binary_n %>% dplyr::select(strain, n_lost), by = "strain")

# Get variable loadings
loadings_n <- tidy(pca_result_n, matrix = "rotation") %>%
  pivot_wider(names_from = PC, values_from = value, names_prefix = "PC") %>%
  mutate(column = case_when(
    column == "Glutamine" ~ "Glutamine growth",
    column == "Asparagine" ~ "Asparagine growth",
    column == "Urea" ~ "Urea growth",
    column == "Nitrate" ~ "Nitrate growth",
    column == "Serine" ~ "Serine growth",
    column == "Ammonium" ~ "Ammonium growth",
    column == "Glutamate" ~  "Glutamate growth",
    column == "Arginine" ~  "Arginine growth",
  ))

# PC1 vs PC2 plot

hulls <- type_scores_n %>%
  group_by(lifestyle) %>%
  slice(chull(.fittedPC1, .fittedPC2))

n_pca <- ggplot(type_scores_n, aes(.fittedPC1, .fittedPC2)) +
  geom_polygon(
    data = hulls,
    aes(fill = lifestyle, group = lifestyle),
    alpha = 0.4
  ) +
  scale_fill_manual(values = c("FL" = "green3",
                               "S" = "cornflowerblue"),
                    name = "Lifestyle",
                    labels = c("FL" = "Free-living",
                               "S" = "Symbiotic")) +
  geom_segment(data = loadings_n, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_n, aes(x = PC1 * 10.2, y = PC2 * 10.2, label = column), size = 3, alpha = 0.9) +
  ggnewscale::new_scale_fill() +
  geom_label(aes(label = strain, fill = factor(n_lost)), size = 3, colour = "white") +
  labs(x = "PC1 (55%)", y = "PC2 (14%)",
       fill = "Number of lost N pathways") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(fill = "grey50")))

n_pca


# Let's look at whether lifestyle clustering is significant
# Creating a new data frame for statistical testing of the phenotypic separation we see in the PCA
rf_df <- mult_df_n %>%
  ungroup() %>%
  mutate(lifestyle = as.factor(lifestyle)) %>%
  dplyr::select(-strain)

rf_df

# PERMANOVA to test whether lifestyle explains the clustering in the multivariate space
# and whether there's an impact of the environmental variables
X <- scale(rf_df[,-1])

adonis2(X ~ lifestyle, method = "euclidean", data = rf_df, by = "term")

# Checking betadispersion
anova(betadisper(dist(X), rf_df$lifestyle))

# Using envfit to test which trait vectors contribute to this separation
pca <- prcomp(rf_df[,-1], scale. = TRUE, center = TRUE)

# Fit trait vectors to PCA space (note this is basically a linear correlation with the multivariate space)
fit <- envfit(pca, rf_df[,-1], permutations = 999)

# View results
fit

# Random forest model to tell us which traits are most predictive for a given lifestyle
set.seed(123)  # for reproducibility
rf_df$lifestyle <- factor(rf_df$lifestyle)

rf_model <- randomForest(lifestyle ~ ., data = rf_df, 
                         importance = TRUE, 
                         ntree = 1000)

# Inspect the model output
print(rf_model)

# Rank the most important predictors
varImpPlot(rf_model, n.var = 8, main = "Top N Pathway Predictors of Lifestyle")

# Summary table of lost pathways
binary_n %>% 
  group_by(lifestyle) %>%
  summarise(
    Arg_lost = n() - sum(Arginine),
    Asp_lost = n() - sum(Asparagine),
    Gln_lost = n() - sum(Glutamine),
    Glu_lost = n() - sum(Glutamate),
    Ser_lost = n() - sum(Serine),
    NH4_lost = n() - sum(Ammonium),
    NO3_lost = n() - sum(Nitrate),
    Urea_lost = n() - sum(Urea),
  )

# Which strains have losses?
binary_n %>%
  filter(n_lost > 0) %>%
  summarise(
    across(
      Arginine:Urea,
      ~ mean(. == 0)
    )
  ) %>%
  print(n = nrow(.))

# Make a pathway presence matrix
pres_mat <- binary_n %>%
  ungroup() %>%
  dplyr::select(Arginine:Urea) %>%
mutate(across(everything(), ~ . == 0))

# Marginal probability of loss
colMeans(pres_mat)

# Conditional loss matrix
cond_prob <- outer(
  names(pres_mat),
  names(pres_mat),
  Vectorize(function(x, y) {
    mean(pres_mat[[y]][pres_mat[[x]]])
  })
)

dimnames(cond_prob) <- list(
  lost_x = names(pres_mat),
  lost_y = names(pres_mat)
)
cond_prob

# Framing this in terms of losses rather than whether they can use them for NODF (nestedness) analysis using vegan
loss_mat <- 1 - binary_n %>%
  ungroup() %>%
  dplyr::select(Arginine:Urea) %>%
  as.matrix()

loss_mat

# Do the analysis
nestednodf(loss_mat)

# Note matrix fill is only 0.15

# Compare to null distribution
nulls <- oecosimu(
  loss_mat,
  nestednodf,
  method = "swap",
  nsimul = 9999
)

summary(nulls)
nestedtemp(loss_mat)

# High nestedness, but significantly less so than the null distribution
# Probably because ammonium is so commonly lost, but no clear order otherwise

# Let's look at pathway loss patterns
binary_n %>%
  ungroup() %>%
  mutate(
    pattern = apply(
      dplyr::select(., Arginine:Urea) == 0,
      1,
      function(x)
        paste(names(x)[x], collapse = ", ")
    )
  ) %>%
  count(pattern, sort = TRUE)
# By far the most common pattern is to not have lost any, really no clearly dominant pattern otherwise

# Split by lifestyle
binary_n %>%
  ungroup() %>%
  mutate(
    pattern = apply(
      dplyr::select(., Arginine:Urea) == 0,
      1,
      function(x)
        paste(names(x)[x], collapse = ", ")
    )
  ) %>%
  count(pattern, lifestyle, sort = TRUE) %>%
  print(n = nrow(.))

binary_n %>%
  pivot_longer(
    Arginine:Urea,
    names_to = "pathway",
    values_to = "present"
  ) %>%
  filter(present == 0) %>%
  group_by(pathway) %>%
  summarise(
    mean_n_lost = mean(n_lost),
    median_n_lost = median(n_lost),
    n = n()
  ) %>%
  arrange(mean_n_lost)

# Let's look at average number of lost pathways when a given pathway is lost to see if we can estimate if loss of a particular pathway
# tends to happen early or late
loss_long <- binary_n %>%
  pivot_longer(
    Arginine:Urea,
    names_to = "pathway",
    values_to = "lost"
  ) %>%
  filter(lost == 0)

set.seed(123)

boot_means <- loss_long %>%
  group_by(pathway) %>%
  group_modify(~{
    
    boots <- replicate(
      10000,
      mean(sample(.x$n_lost,
                  size = nrow(.x),
                  replace = TRUE))
    )
    
    tibble(
      mean = mean(.x$n_lost),
      lower = quantile(boots, 0.025),
      upper = quantile(boots, 0.975),
      n = nrow(.x)
    )
  }) %>%
  ungroup() %>%
  arrange(mean)

boot_means

ggplot(
  boot_means,
  aes(
    x = reorder(pathway, mean),
    y = mean
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower,
        ymax = upper),
    width = 0.2
  ) +
  coord_flip() +
  labs(
    x = "Pathway",
    y = "Mean number of pathways lost"
  )

# All that is clear from this is that arginine is only ever lost on its own. 
# There's a bit of a trend for ammonium to be generally early and urea to be generally late, but everything else overlaps

# Finally let's look at whether pathways tend to be co-lost

# Co-loss matrix
coloss_mat <- binary_n %>%
  ungroup() %>%
  dplyr::select(Arginine:Urea) %>%
  mutate(across(everything(),
                ~ . == 0)) %>%
  as.matrix()

co_loss <- t(coloss_mat) %*% coloss_mat
co_loss

# Make network on raw values
# Raw values
g_raw <- graph_from_adjacency_matrix(
  co_loss,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

plot(
  g_raw,
  edge.width = E(g_raw)$weight,
  vertex.size = diag(co_loss) * 2
)

cluster_walktrap(g_raw)
strength(g_raw)

# Might be really biased by the omnipresence of ammonium loss though, so let's try some association values instead

# Jaccard distances instead
jacc <- vegdist(
  t(loss_mat),
  method = "jaccard",
  binary = TRUE
)

jacc_mat <- 1 - as.matrix(jacc)

g_jacc <- graph_from_adjacency_matrix(
  jacc_mat,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

plot(
  g_jacc,
  edge.width = E(g_jacc)$weight,
  vertex.size = diag(co_loss) * 2
)

cluster_walktrap(g_jacc)
strength(g_jacc)

# Same result really, but this doesn't account for frequency

# Phi coefficient instead
phi_mat <- rcorr(as.matrix(loss_mat))

g_phi <- graph_from_adjacency_matrix(
  phi_mat$r,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

plot(
  g_phi,
  edge.width = E(g_phi)$weight,
  vertex.size = diag(co_loss) * 2
)

cluster_walktrap(g_phi)

# We no longer get modules when accounting for frequency
# This will be because Arginine is so rare, so it can't really be reliably estimated, see p-values:
round(phi_mat$P, 3)

# Let's visualise
phi_plot <- phi_mat$r[
  c("Arginine","Asparagine","Glutamine",
    "Serine","Ammonium","Nitrate","Urea"),
  c("Arginine","Asparagine","Glutamine",
    "Serine","Ammonium","Nitrate","Urea")
]

corrplot(
  phi_plot,
  method = "color",
  type = "lower",
  addCoef.col = "black",
  tl.col = "black",
  na.label = "",
  mar = c(0, 0, 0, 0))

corr_df <- as.data.frame(as.table(phi_plot)) %>%
  rename(
    Var1 = Var1,
    Var2 = Var2,
    corr = Freq
  )

# Keep lower triangle only
corr_df <- corr_df %>%
  mutate(
    Var1 = factor(Var1, levels = rownames(phi_plot)),
    Var2 = factor(Var2, levels = colnames(phi_plot))
  ) %>%
  filter(as.numeric(Var1) > as.numeric(Var2))

gg_corr <- ggplot(corr_df,
                  aes(Var2, Var1, fill = corr)) +
  geom_tile() +
  geom_text(aes(label = round(corr, 2)),
            size = 3) +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank()
  ) +
  labs(fill = "Correlation")

gg_corr
