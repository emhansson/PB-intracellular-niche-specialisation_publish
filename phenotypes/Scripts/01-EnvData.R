# Erika M Hansson 2025
# erika.hansson@manchester.ac.uk

#### ---------------------------------------------------- ####
# 1) Summary of environmental breadth covered
# Then: repeating all of the above with the 86 strains in the phylogeny for robustness
#### ---------------------------------------------------- ####

# Packages
if(!require("tidyverse")) {install.packages("tidyverse")} # wrangling, plotting

# Load data frame
df_tidy_env <- read_csv("Data/df_tidy_env.csv")

# Look at it
glimpse(df_tidy_env)

df_tidy_env

## What's the breadth of environments represented by these ponds?
# Summary data frame giving us min, max, mean and sd for each vairable
df_tidy_env %>%
  dplyr::select(ph, conductivity, tds, no3) %>%
  summarise(across(c(ph, conductivity, tds, no3), list(
    min = ~min(. , na.rm = TRUE),
    max = ~max(. , na.rm = TRUE),
    mean = ~mean(. , na.rm = TRUE),
    sd   = ~sd(. , na.rm = TRUE)
  ), .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = everything(),
               names_to = c("variable", "statistic"),
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = statistic, values_from = value)

# Violin plot of spread
df_tidy_env %>%
  dplyr::select(ph, conductivity, tds, no3, pond_id) %>%
  pivot_longer(cols = c(ph, conductivity, tds, no3),
               names_to = "variable",
               values_to = "value") %>%
ggplot(aes(x = variable, y = value)) +
  geom_violin(fill = "lightgreen", alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkgreen") +
  theme_minimal() +
  labs(x = "Parameter",
       y = "Value") +
  scale_x_discrete(labels = "", position = "top") +
  facet_wrap(~variable, scales = "free", ncol = 4)

# Reanalysis using only the 86 strains from the phylogeny to check robustness of results ###########################

## What's the breadth of environments represented by these ponds?
# Summary data frame giving us min, max, mean and sd for each vairable
df_tidy_env %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  dplyr::select(ph, conductivity, tds, no3) %>%
  summarise(across(c(ph, conductivity, tds, no3), list(
    min = ~min(. , na.rm = TRUE),
    max = ~max(. , na.rm = TRUE),
    mean = ~mean(. , na.rm = TRUE),
    sd   = ~sd(. , na.rm = TRUE)
  ), .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = everything(),
               names_to = c("variable", "statistic"),
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = statistic, values_from = value)

# Violin plot of spread
df_tidy_env %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  dplyr::select(ph, conductivity, tds, no3, pond_id) %>%
  pivot_longer(cols = c(ph, conductivity, tds, no3),
               names_to = "variable",
               values_to = "value") %>%
  ggplot(aes(x = variable, y = value)) +
  geom_violin(fill = "lightgreen", alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkgreen") +
  theme_minimal() +
  labs(x = "Parameter",
       y = "Value") +
  scale_x_discrete(labels = "", position = "top") +
  facet_wrap(~variable, scales = "free", ncol = 4)

# Slight shift in mean and sd, but same environmental breadth covered

# Reanalysis using only the 92 strains from the phylogeny to check robustness of results ###########################

## What's the breadth of environments represented by these ponds?
# Summary data frame giving us min, max, mean and sd for each vairable
df_tidy_env %>%
  filter(strain != "192") %>%
  dplyr::select(ph, conductivity, tds, no3) %>%
  summarise(across(c(ph, conductivity, tds, no3), list(
    min = ~min(. , na.rm = TRUE),
    max = ~max(. , na.rm = TRUE),
    mean = ~mean(. , na.rm = TRUE),
    sd   = ~sd(. , na.rm = TRUE)
  ), .names = "{.col}_{.fn}")) %>%
  pivot_longer(cols = everything(),
               names_to = c("variable", "statistic"),
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = statistic, values_from = value)

# Violin plot of spread
df_tidy_env %>%
  filter(strain != "192") %>%
  dplyr::select(ph, conductivity, tds, no3, pond_id) %>%
  pivot_longer(cols = c(ph, conductivity, tds, no3),
               names_to = "variable",
               values_to = "value") %>%
  ggplot(aes(x = variable, y = value)) +
  geom_violin(fill = "lightgreen", alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkgreen") +
  theme_minimal() +
  labs(x = "Parameter",
       y = "Value") +
  scale_x_discrete(labels = "", position = "top") +
  facet_wrap(~variable, scales = "free", ncol = 4)

# No material change
