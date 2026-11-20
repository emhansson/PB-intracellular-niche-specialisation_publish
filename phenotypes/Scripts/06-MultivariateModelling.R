# Erika M Hansson 2025
# erika.hansson@manchester.ac.uk

#### ---------------------------------------------------- ####
# 1) Import all data and summarise by strain for multivariate analysis
# 2) Fit PCA
# 3) Create biplot and overlay strain identity and lifestyle
# 4) Overlay environmental factors
# 5) PERMANOVA, envfit and random forest modelling on variables
# Then: repeating all of the above with the 86 strains in the phylogeny for robustness
# 6) FAMD modelling to test whether phylogroup should be included in multivariate analysis
# 7) PCA, PERMANOVA on phylogroup only
#### ---------------------------------------------------- ####

# Load packages
if(!require("tidyverse")) {install.packages("tidyverse")} # wrangling, plotting
if(!require("gridExtra")) {install.packages("gridExtra")} # plotting
if(!require("cowplot")) {install.packages("cowplot")} # plotting
if(!require("ggrepel")) {install.packages("ggrepel")} # plot labelling
if(!require("vegan")) {install.packages("vegan")} # ordination methods etc.
if(!require("randomForest")) {install.packages("randomForest")} # random forest modelling
if(!require("broom.mixed")) {install.packages("broom.mixed")} # wrangling untidy outputs
if(!require("FactoMineR")) {install.packages("FactoMineR")} # FADM
if(!require("factoextra")) {install.packages("factoextra")} # FADM

# Import sugar data and summarise as slopes for each sugar for each strain
mult_df_sug <- 
  read_csv("Data/df_tidy_sugar.csv") %>%
  # Remove trehalose readings
  filter(sugar != "Trehalose") %>%
  # Calculate slope for each strain in each pH
  group_by(strain, ph, sugar, lifestyle) %>%
  summarise(
    mean = mean(release)
  ) %>%
  pivot_wider(names_from = "ph", values_from = "mean") %>%
  mutate(slope = `5.5` - `7.6`) %>%
  dplyr::select(-`7.6`, -`5.5`) %>%
  pivot_wider(names_from = "sugar", values_from = "slope")

# Add in means for release as well (means represent overall release, slope represents resopsiveness to pH change)
mult_df_sug_combo <- 
  read_csv("Data/df_tidy_sugar.csv") %>%
  # Remove trehalose readings
  filter(sugar != "Trehalose") %>%
  group_by(strain, sugar, lifestyle) %>%
  summarise(
    mean = mean(release)
  ) %>%
  pivot_wider(names_from = "sugar", values_from = "mean") %>%
  rename("Glucose_mean" = Glucose,
         "Maltose_mean" = Maltose) %>%
  left_join(mult_df_sug) %>%
  rename("Glucose_slope" = Glucose,
         "Maltose_slope" = Maltose)

glimpse(mult_df_sug_combo)

# Import nitrogen growth data, summarise as growth rates by strain
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


# Import pHGR data, summarise as slopes for growth rate by strain
mult_df_phgr <- read_csv("Data/df_tidy_pH.csv") %>%
  group_by(strain, ph, lifestyle) %>%
  summarise(
    mean = mean(od750)
  ) %>%
  pivot_wider(names_from = "ph", values_from = "mean") %>%
  mutate(pH = `5.5` - `7.6`) %>%
  dplyr::select(-`7.6`, -`5.5`)

# Doing the same thing as for sugars above with using both means and slopes
mult_df_phgr_combo <- read_csv("Data/df_tidy_pH.csv") %>%
  group_by(strain, lifestyle) %>%
  summarise(
    mean = mean(od750)
  ) %>%
  rename("ph_mean" = mean) %>%
  left_join(mult_df_phgr) %>%
  rename("ph_slope" = pH)

glimpse(mult_df_phgr_combo)

# Import fluorimetry data, summarise as mean Fv/Fm by strain
# Has 186b, sample for now, will be removed later
mult_df_fl <- read_csv("Data/df_tidy_fluorimetry.csv") %>%
  group_by(strain, lifestyle) %>%
  summarise(
    FvFm = mean(mean_FvFm)
  )

glimpse(mult_df_fl)

# Put together and do PCA -- with combo data frames using both slopes and means where appropriate
mult_df_combo <- mult_df_sug_combo %>%
  merge(mult_df_n, by = c("strain", "lifestyle")) %>%
  merge(mult_df_phgr_combo, by = c("strain", "lifestyle")) %>%
  merge(mult_df_fl, by = c("strain", "lifestyle"))

glimpse(mult_df_combo)

write_csv(mult_df_combo, "Data/mult_df_combo.csv")

# PCA of strains 
mult_df_scaled_combo <- mult_df_combo %>%
  dplyr::select(-strain, -lifestyle) %>%
  drop_na() %>%                # remove rows with missing values
  scale()                      # standardize to mean 0, sd 1

# Run PCA
pca_result_combo <- prcomp(mult_df_scaled_combo, center = FALSE, scale. = FALSE)

# View summary of explained variance
summary(pca_result_combo)

# Scree plot to visualize variance explained
plot(pca_result_combo, type = "l")

# Tidy version of PCA result
pca_df_combo <- as_tibble(pca_result_combo$x)

pca_df_combo

# Add metadata for plotting
pca_df_combo <- bind_cols(mult_df_combo %>% drop_na() %>% dplyr::select(strain, lifestyle), pca_df_combo)

# Plot first two PCs
ggplot(pca_df_combo, aes(x = PC1, y = PC2, colour = lifestyle, label = strain)) +
  labs(x = "PC1", y = "PC2") +
  geom_label() +
  theme_bw()

# Get loadings
pca_result_combo$rotation

# Get scores for each PC by type
type_scores_combo <- augment(pca_result_combo) %>%
  bind_cols(mult_df_combo) %>%
  dplyr::select(2:4, strain, lifestyle) %>%
  unique()

# Get variable loadings
loadings_combo <- tidy(pca_result_combo, matrix = "rotation") %>%
  pivot_wider(names_from = PC, values_from = value, names_prefix = "PC") %>%
  mutate(column = case_when(
    column == "Maltose_mean" ~ "Average maltose release",
    column == "Maltose_slope" ~ "Maltose release pH sensitivity",
    column == "Glucose_mean" ~ "Average glucose release",
    column == "Glucose_slope" ~ "Glucose release pH sensitivity",
    column == "ph_mean" ~ "Average growth rate",
    column == "ph_slope" ~ "Growth rate pH sensitivity",
    column == "FvFm" ~ "Fv/Fm",
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
ggplot(type_scores_combo, aes(.fittedPC1, .fittedPC2)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings_combo, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo, aes(x = PC1 * 10.2, y = PC2 * 10.2, label = column), size = 3) +
  labs(x = "PC1 (44%)", y = "PC2 (13%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# PC1 v PC3
ggplot(type_scores_combo, aes(.fittedPC1, .fittedPC3)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings_combo, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC3 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo, aes(x = PC1 * 10.2, y = PC3 * 10.2, label = column), size = 3) +
  labs(x = "PC1 (44%)", y = "PC3 (10%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# PC2 v PC3
ggplot(type_scores_combo, aes(.fittedPC2, .fittedPC3)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings_combo, aes(x = 0, y = 0, xend = PC2 * 10, yend = PC3 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo, aes(x = PC2 * 10.2, y = PC3 * 10.2, label = column), size = 3) +
  labs(x = "PC2 (13%)", y = "PC3 (10%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# PC1 v PC3
ggplot(type_scores_combo, aes(.fittedPC1, .fittedPC3)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings_combo, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC3 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo, aes(x = PC1 * 10.2, y = PC3 * 10.2, label = column), size = 3) +
  labs(x = "PC1 (44%)", y = "PC3 (10%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# Highlighting the strains in the top cluster for the phylogeny
# PC1 vs PC2 plot
p.pc1_pc2 <- ggplot(type_scores_combo %>%
         mutate(cluster = case_when(strain == "25" | strain == "63" | strain == "34" | strain == "1" |
                                        strain == "270" | strain == "29" | strain == "318" | strain == "273" |
                                        strain == "27" | strain == "31" | strain == "126" ~ "FL (top cluster)",
                                      strain != "s70" & lifestyle == "S" ~ "S (top cluster)",
                                      TRUE ~ lifestyle)), aes(.fittedPC1, .fittedPC2)) +
  geom_label(aes(label = strain, fill = lifestyle, colour = cluster), size = 3) +
  geom_segment(data = loadings_combo, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo, aes(x = PC1 * 10.2, y = PC2 * 10.2, label = column), size = 3) +
  labs(x = "PC1 (44%)", y = "PC2 (13%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle") +
  scale_colour_manual(values = c(
    "S (top cluster)" = "darkblue",
    "FL (top cluster)" = "darkgreen",
    "S" = "lightblue",
    "FL" = "lightgreen"
  ), guide = "none") +
  theme(legend.position = "bottom")

# PC2 v PC3
p.pc2_pc3 <- ggplot(type_scores_combo %>%
         mutate(cluster = case_when(strain == "25" | strain == "63" | strain == "34" | strain == "1" |
                                      strain == "270" | strain == "29" | strain == "318" | strain == "273" |
                                      strain == "27" | strain == "31" | strain == "126" ~ "FL (top cluster)",
                                    strain != "s70" & lifestyle == "S" ~ "S (top cluster)",
                                    TRUE ~ lifestyle)), aes(.fittedPC2, .fittedPC3)) +
  geom_label(aes(label = strain, fill = lifestyle, colour = cluster), size = 3) +
  geom_segment(data = loadings_combo, aes(x = 0, y = 0, xend = PC2 * 10, yend = PC3 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo, aes(x = PC2 * 10.2, y = PC3 * 10.2, label = column), size = 3) +
  labs(x = "PC2 (13%)", y = "PC3 (10%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  scale_colour_manual(values = c(
    "S (top cluster)" = "darkblue",
    "FL (top cluster)" = "darkgreen",
    "S" = "lightblue",
    "FL" = "lightgreen"
  ), guide = "none") +
  theme(legend.position = "bottom")

# PC1 v PC3
p.pc1_pc3 <- ggplot(type_scores_combo %>%
         mutate(cluster = case_when(strain == "25" | strain == "63" | strain == "34" | strain == "1" |
                                      strain == "270" | strain == "29" | strain == "318" | strain == "273" |
                                      strain == "27" | strain == "31" | strain == "126" ~ "FL (top cluster)",
                                    strain != "s70" & lifestyle == "S" ~ "S (top cluster)",
                                    TRUE ~ lifestyle)), aes(.fittedPC1, .fittedPC3)) +
  geom_label(aes(label = strain, fill = lifestyle, colour = cluster), size = 3) +
  geom_segment(data = loadings_combo, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC3 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo, aes(x = PC1 * 10.2, y = PC3 * 10.2, label = column), size = 3) +
  labs(x = "PC1 (44%)", y = "PC3 (10%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  scale_colour_manual(values = c(
    "S (top cluster)" = "darkblue",
    "FL (top cluster)" = "darkgreen",
    "S" = "lightblue",
    "FL" = "lightgreen"
  ), guide = "none") +
  theme(legend.position = "bottom")

# Extract legend
legend <- cowplot::get_legend(
  p.pc1_pc2 + theme(legend.position = "top")
)

# Remove legends from the individual plots
p1 <- p.pc1_pc2 + theme(legend.position = "none")
p2 <- p.pc1_pc3 + theme(legend.position = "none")
p3 <- p.pc2_pc3 + theme(legend.position = "none")

right_col <- cowplot::plot_grid(
  legend,      # grob from get_legend()
  p2, 
  p3,
  ncol = 1,
  rel_heights = c(0.12, 0.44, 0.44),  # tweak these to taste (must sum to 1-ish)
  align = "v"
)

# final layout: p1 (left) + right_col (right)
final_plot <- cowplot::plot_grid(
  p1,
  right_col,
  ncol = 2,
  rel_widths = c(0.6, 0.4),   # left = 60% width, right = 40% width
  align = "hv",
  axis = "tblr"
)

final_plot

# Joining with environmental data for overlay
type_scores_combo <- type_scores_combo %>%
  left_join(read_csv("Data/df_tidy_env.csv"), by = c("strain", "lifestyle")) %>%
  # Adding phylotypes
  left_join(read_csv("Data/phylotypes.csv"), by = "strain") %>%
  mutate(phylotype = as.character(phylotype))


overlay_biplots2 <- function(overlay, title) {
  
  overlay_var <- deparse(substitute(overlay))
  is_continuous <- is.numeric(type_scores_combo[[overlay_var]])
  
  base_plot <- function(x, y, xlab, ylab) {
    ggplot(type_scores_combo, aes_string(x, y)) +
      geom_label(aes_string(label = "strain", fill = overlay_var),
                 size = 3, colour = "white") +
      labs(x = xlab, y = ylab) +
      theme_bw() +
      theme(legend.position = "bottom")
  }
  
  p1 <- base_plot(".fittedPC1", ".fittedPC2", "PC1 (44%)", "PC2 (13%)")
  p2 <- base_plot(".fittedPC1", ".fittedPC3", "PC1 (44%)", "PC3 (10%)")
  p3 <- base_plot(".fittedPC2", ".fittedPC3", "PC2 (13%)", "PC3 (10%)")
  
  # Apply appropriate legend scaling
  if (is_continuous) {
    p1 <- p1 + guides(fill = guide_colorbar(barwidth = 20, barheight = 1))
  } else {
    # Get number of unique levels in the overlay variable
    n_levels <- length(unique(type_scores_combo[[overlay_var]]))
    
    # Adjust number of columns based on how many rows we want
    ncol_legend <- ceiling(n_levels / 3)
    
    p1 <- p1 + guides(
      fill = guide_legend(
        ncol = ncol_legend,
        byrow = TRUE
      )
    )
  }
  
  # Extract legend
  legend <- cowplot::get_legend(
    p1 + theme(legend.position = "top") + labs(fill = paste(title))
  )
  
  # Remove legends from the individual plots
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- p3 + theme(legend.position = "none")
  
  aligned <- cowplot::align_plots(p1, p2, p3, align = "hv", axis = "tblr")
  
  p <- cowplot::plot_grid(
    cowplot::plot_grid(plotlist = aligned, ncol = 3, align = "v"),
    legend,
    ncol = 1,
    rel_heights = c(1, 0.35)
  )
  
  return(p)
}

# Origin pond
p.pond_id <- overlay_biplots2(pond_id, "Pond ID")

# Pond pH
p.pond_ph <- overlay_biplots2(ph, "Pond pH")

# Pond area
p.pond_area <- overlay_biplots2(pond_area, "Pond Area")

# Pond NO3
p.pond_no3 <- overlay_biplots2(no3, "Pond Nitrate Content")

# Pond TDS
p.pond_tds <- overlay_biplots2(tds, "Pond TDS")

# Pond conductivity
p.pond_conductivity <- overlay_biplots2(conductivity, "Pond Conductivity")

# Put together plots of interest
p.env_overlay <- plot_grid(p.pond_ph, p.pond_no3, p.pond_tds, p.pond_conductivity, p.pond_id,
                           labels = c("a)", "b)", "c)", "d)", "e)"), 
                           label_size = 14,
                           label_fontface = "bold",
                              ncol = 1,
                           scale = 1)

p.env_overlay

# Also checking phylotype
overlay_biplots2(phylotype, "Phylotype")

# Creating a new data frame for statistical testing of the phenotypic separation we see in the PCA
rf_df <- mult_df_combo %>%
  mutate(lifestyle = as.factor(lifestyle)) %>%
  rename("Autonomous_growth_mean" = ph_mean,
         "Autonomous_growth_slope" = ph_slope) %>%
  left_join(read_csv("Data/df_tidy_env.csv"), by = c("strain", "lifestyle")) %>%
  dplyr::select(-strain)

rf_df

# PERMANOVA to test whether lifestyle explains the clustering in the multivariate space
# and whether there's an impact of the environmental variables
X <- scale(rf_df[,-c(1, 17, 18, 19, 20, 21, 22)])

adonis2(X ~ lifestyle, method = "euclidean", data = rf_df, by = "term")
adonis2(X ~ lifestyle*pond_id, method = "euclidean", data = rf_df, by = "term")
adonis2(X ~ lifestyle*no3, method = "euclidean", data = rf_df, by = "term")
adonis2(X ~ lifestyle*tds, method = "euclidean", data = rf_df, by = "term")
adonis2(X ~ lifestyle*ph, method = "euclidean", data = rf_df, by = "term")
adonis2(X ~ lifestyle*conductivity, method = "euclidean", data = rf_df, by = "term")

# Checking betadispersion
anova(betadisper(dist(X), rf_df$lifestyle))
anova(betadisper(dist(X), rf_df$pond_id))
anova(betadisper(dist(X), rf_df$no3))
anova(betadisper(dist(X), rf_df$tds))
anova(betadisper(dist(X), rf_df$ph))
anova(betadisper(dist(X), rf_df$conductivity))

# Using envfit to test which trait vectors contribute to this separation
pca <- prcomp(rf_df[,-c(1, 17, 18, 19, 20, 21, 22)], scale. = TRUE, center = TRUE)

# Fit trait vectors to PCA space (note this is basically a linear correlation with the multivariate space)
fit <- envfit(pca, rf_df[,-c(1, 17, 18, 19, 20, 21, 22)], permutations = 999)

# View results
fit

# Random forest model to tell us which traits are most predictive for a given lifestyle
set.seed(123)  # for reproducibility
rf_df$lifestyle <- factor(rf_df$lifestyle)

rf_model <- randomForest(lifestyle ~ ., data = rf_df[,-c(17, 18, 19, 20, 21, 22)], 
                         importance = TRUE, 
                         ntree = 1000)

# Inspect the model output
print(rf_model)

# Rank the most important predictors
varImpPlot(rf_model, n.var = 15, main = "Top Trait Predictors of Lifestyle")


# Reanalysis using only the 86 strains from the phylogeny to check robustness of results ###########################
# Put together and do PCA -- with combo data frames using both slopes and means where appropriate
mult_df_combo_86 <- mult_df_sug_combo %>%
  filter(strain != "155",
         strain != "192",
         strain != "242",
         strain != "252",
         strain != "258",
         strain != "s354") %>%
  merge(mult_df_n, by = c("strain", "lifestyle")) %>%
  merge(mult_df_phgr_combo, by = c("strain", "lifestyle")) %>%
  merge(mult_df_fl, by = c("strain", "lifestyle"))


glimpse(mult_df_combo_86)

# PCA of strains 
mult_df_scaled_combo_86 <- mult_df_combo_86 %>%
  dplyr::select(-strain, -lifestyle) %>%
  drop_na() %>%                # remove rows with missing values
  scale()                      # standardize to mean 0, sd 1

# Run PCA
pca_result_combo_86 <- prcomp(mult_df_scaled_combo_86, center = FALSE, scale. = FALSE)

# View summary of explained variance
summary(pca_result_combo_86)

# Scree plot to visualize variance explained
plot(pca_result_combo_86, type = "l")

# Tidy version of PCA result
pca_df_combo_86 <- as_tibble(pca_result_combo_86$x)

# Add metadata for plotting
pca_df_combo_86 <- bind_cols(mult_df_combo_86 %>% drop_na() %>% dplyr::select(strain, lifestyle), pca_df_combo_86)

# Plot first two PCs
ggplot(pca_df_combo_86, aes(x = PC1, y = PC2, colour = lifestyle, label = strain)) +
  labs(x = "PC1", y = "PC2") +
  geom_label() +
  theme_bw()

# Get loadings
pca_result_combo_86$rotation

# Get scores for each PC by type
type_scores_combo_86 <- augment(pca_result_combo_86) %>%
  bind_cols(mult_df_combo_86) %>%
  dplyr::select(2:4, strain, lifestyle) %>%
  unique()

# Get variable loadings
loadings_combo_86 <- tidy(pca_result_combo_86, matrix = "rotation") %>%
  pivot_wider(names_from = PC, values_from = value, names_prefix = "PC") %>%
  mutate(column = case_when(
    column == "Maltose_mean" ~ "Average maltose release",
    column == "Maltose_slope" ~ "Maltose release pH sensitivity",
    column == "Glucose_mean" ~ "Average glucose release",
    column == "Glucose_slope" ~ "Glucose release pH sensitivity",
    column == "ph_mean" ~ "Average growth rate",
    column == "ph_slope" ~ "Growth rate pH sensitivity",
    column == "FvFm" ~ "Fv/Fm",
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
ggplot(type_scores_combo_86, aes(.fittedPC1, .fittedPC2)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings_combo_86, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo_86, aes(x = PC1 * 10.2, y = PC2 * 10.2, label = column), size = 3) +
  labs(x = "PC1 (44%)", y = "PC2 (13%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# PC1 v PC3
ggplot(type_scores_combo_86, aes(.fittedPC1, .fittedPC3)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings_combo_86, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC3 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo_86, aes(x = PC1 * 10.2, y = PC3 * 10.2, label = column), size = 3) +
  labs(x = "PC1 (44%)", y = "PC3 (10%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# PC2 v PC3
ggplot(type_scores_combo_86, aes(.fittedPC2, .fittedPC3)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings_combo_86, aes(x = 0, y = 0, xend = PC2 * 10, yend = PC3 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo_86, aes(x = PC2 * 10.2, y = PC3 * 10.2, label = column), size = 3) +
  labs(x = "PC2 (13%)", y = "PC3 (10%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# Joining with environmental data for overlay
type_scores_combo_86 <- type_scores_combo_86 %>%
  left_join(read_csv("Data/df_tidy_env.csv"), by = c("strain", "lifestyle")) %>%
  # Adding phylotypes
  left_join(read_csv("Data/phylotypes.csv"), by = "strain") %>%
  mutate(phylotype = as.character(phylotype))


overlay_biplots2_86 <- function(overlay, title) {
  
  overlay_var <- deparse(substitute(overlay))
  is_continuous <- is.numeric(type_scores_combo[[overlay_var]])
  
  base_plot <- function(x, y, xlab, ylab) {
    ggplot(type_scores_combo_86, aes_string(x, y)) +
      geom_label(aes_string(label = "strain", fill = overlay_var),
                 size = 3, colour = "white") +
      labs(x = xlab, y = ylab) +
      theme_bw() +
      theme(legend.position = "bottom")
  }
  
  p1 <- base_plot(".fittedPC1", ".fittedPC2", "PC1 (44%)", "PC2 (13%)")
  p2 <- base_plot(".fittedPC1", ".fittedPC3", "PC1 (44%)", "PC3 (10%)")
  p3 <- base_plot(".fittedPC2", ".fittedPC3", "PC2 (13%)", "PC3 (10%)")
  
  # Apply appropriate legend scaling
  if (is_continuous) {
    p1 <- p1 + guides(fill = guide_colorbar(barwidth = 20, barheight = 1))
  } else {
    # Get number of unique levels in the overlay variable
    n_levels <- length(unique(type_scores_combo[[overlay_var]]))
    
    # Adjust number of columns based on how many rows we want
    ncol_legend <- ceiling(n_levels / 3)
    
    p1 <- p1 + guides(
      fill = guide_legend(
        ncol = ncol_legend,
        byrow = TRUE
      )
    )
  }
  
  # Extract legend
  legend <- cowplot::get_legend(
    p1 + theme(legend.position = "top") + labs(fill = paste(title))
  )
  
  # Remove legends from the individual plots
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- p3 + theme(legend.position = "none")
  
  aligned <- cowplot::align_plots(p1, p2, p3, align = "hv", axis = "tblr")
  
  p <- cowplot::plot_grid(
    cowplot::plot_grid(plotlist = aligned, ncol = 3, align = "v"),
    legend,
    ncol = 1,
    rel_heights = c(1, 0.35)
  )
  
  return(p)
}

# Origin pond
p.pond_id_86 <- overlay_biplots2_86(pond_id, "Pond ID")

# Pond pH
p.pond_ph_86 <- overlay_biplots2_86(ph, "Pond pH")

# Pond area
p.pond_area_86 <- overlay_biplots2_86(pond_area, "Pond Area")

# Pond NO3
p.pond_no3_86 <- overlay_biplots2_86(no3, "Pond Nitrate Content")

# Pond TDS
p.pond_tds_86 <- overlay_biplots2_86(tds, "Pond TDS")

# Pond conductivity
p.pond_conductivity_86 <- overlay_biplots2_86(conductivity, "Pond Conductivity")


overlay_biplots2_86(phylotype, "Phylotypes")

# Put together plots of interest

p.env_overlay_86 <- plot_grid(p.pond_ph_86, p.pond_no3_86, p.pond_tds_86, p.pond_conductivity_86, p.pond_id_86,
                              labels = c("a)", "b)", "c)", "d)", "e)"), 
                              label_size = 14,
                              label_fontface = "bold",
                              ncol = 1,
                              scale = 1)

p.env_overlay_86


# Creating a new data frame for statistical testing of the phenotypic separation we see in the PCA
rf_df_86 <- mult_df_combo_86 %>%
  mutate(lifestyle = as.factor(lifestyle)) %>%
  rename("Autonomous_growth_mean" = ph_mean,
         "Autonomous_growth_slope" = ph_slope) %>%
  left_join(read_csv("Data/df_tidy_env.csv"), by = c("strain", "lifestyle")) %>%
  dplyr::select(-strain)

rf_df_86

# PERMANOVA to test whether lifestyle explains the clustering in the multivariate space
# and whether there's an impact of the environmental variables
X_86 <- scale(rf_df_86[,-c(1, 17, 18, 19, 20, 21, 22)])
adonis2(X_86 ~ lifestyle, method = "euclidean", data = rf_df_86, by = "term")
adonis2(X_86 ~ lifestyle*pond_id, method = "euclidean", data = rf_df_86, by = "term")
adonis2(X_86 ~ lifestyle*no3, method = "euclidean", data = rf_df_86, by = "term")
adonis2(X_86 ~ lifestyle*tds, method = "euclidean", data = rf_df_86, by = "term")
adonis2(X_86 ~ lifestyle*ph, method = "euclidean", data = rf_df_86, by = "term")
adonis2(X_86 ~ lifestyle*conductivity, method = "euclidean", data = rf_df_86, by = "term")

# Checking betadispersion
anova(betadisper(dist(X_86), rf_df_86$lifestyle))
anova(betadisper(dist(X_86), rf_df_86$pond_id))
anova(betadisper(dist(X_86), rf_df_86$no3))
anova(betadisper(dist(X_86), rf_df_86$tds))
anova(betadisper(dist(X_86), rf_df_86$ph))
anova(betadisper(dist(X_86), rf_df_86$conductivity))

# Using envfit to test which trait vectors contribute to this separation
pca <- prcomp(rf_df_86[,-c(1, 17, 18, 19, 20, 21, 22)], scale. = TRUE, center = TRUE)

# Fit trait vectors to PCA space (note this is basically a linear correlation with the multivariate space)
fit <- envfit(pca, rf_df_86[,-c(1, 17, 18, 19, 20, 21, 22)], permutations = 999)

# View results
fit

# Random forest model to tell us which traits are most predictive for a given lifestyle
set.seed(123)  # for reproducibility
rf_df_86$lifestyle <- factor(rf_df_86$lifestyle)

rf_model_86 <- randomForest(lifestyle ~ ., data = rf_df_86[,-c(17, 18, 19, 20, 21, 22)], 
                            importance = TRUE, 
                            ntree = 1000)

# Inspect the model output
print(rf_model_86)

# Rank the most important predictors
varImpPlot(rf_model_86, n.var = 15, main = "Top Trait Predictors of Lifestyle")

# Conclusions unchanged


## FAMD to test inclusion of phylotypes in the multivariate model
# Make a new data frame
famd_df <- mult_df_combo %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.)))) %>%
  left_join(read_csv("Data/phylotypes.csv") %>% filter(strain != 126), by = "strain") %>%
  drop_na(phylotype) %>% 
  column_to_rownames("strain") %>%
  dplyr::select(-lifestyle, -strain_long) %>%
  mutate(phylotype = as.factor(phylotype))

glimpse(famd_df)

# Fit the model, get some plots
res <- FAMD(
  famd_df,
  ncp = 16             # number of dimensions to keep
)

# Look at important variables from model
# Eigenvalues/variance
res$eig
# These are the coordinates for our strains
res$ind$coord
plot(res, axes = c(1, 2))

# Contribution of each variable
res$var$contrib

# factoextra plot
fviz_famd_ind(
  res,
  habillage = "phylotype",
  addEllipses = TRUE
)

# Recreating the PCA plot
famd_scores <- as.data.frame(res$ind$coord[, 1:2]) %>%
  mutate(strain = row.names(res$ind$coord),
         lifestyle = ifelse(str_detect(strain, "s"), "S", "FL")) %>%
  left_join(read_csv("Data/phylotypes.csv") %>% filter(strain != 126), by = "strain")
famd_scores

famd_loadings <- res$var$coord # I know they're not actually loadings but it's close as we can get
famd_loadings

ggplot(famd_scores, aes(-Dim.1, Dim.2)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3) +
  geom_segment(data = famd_loadings, aes(x = 0, y = 0, xend = -Dim.1 * 10, yend = Dim.2 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = famd_loadings, aes(x = -Dim.1 * 10.2, y = Dim.2 * 10.2, label = row.names(famd_loadings)), size = 3) +
  stat_ellipse(aes(group = phylotype, colour = factor(phylotype)), type = "norm", linewidth = 0.8) +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")


# PERMANOVA for phylotype
famd_df_perm <- famd_df %>%
  mutate(strain = row.names(famd_df),
         lifestyle = ifelse(str_detect(strain, "s"), "S", "FL"))

famd_df_perm

X_famd <- scale(famd_df_perm[,-c(16:18)])

adonis2(X_famd ~ phylotype + lifestyle, method = "euclidean", data = famd_df_perm, by = "term")

anova(betadisper(dist(X_famd), famd_df_perm$phylotype))

pca_famd <- prcomp(famd_df[, -16], center = TRUE, scale. = TRUE)

axis_df <- data.frame(
  PCA1  = pca_famd$x[,1],
  FAMD1 = res$ind$coord[,1],
  PCA2  = pca_famd$x[,2],
  FAMD2 = res$ind$coord[,2],
  phylogroup = famd_df$phylotype
)

# How phylogroup fits to axis 1
anova(lm(FAMD1 ~ phylogroup, data = axis_df))
anova(lm(PCA1 ~ phylogroup, data = axis_df))

# Correlation between axes in the two models
cor(axis_df$PCA1, axis_df$FAMD1)
cor(axis_df$PCA2, axis_df$FAMD2)


# PCA on just phylogroup 1
famd_df1 <- famd_df %>%
  filter(phylotype == "Phylogroup 1") %>%
  dplyr::select(-phylotype)

# Run PCA
pca_result1 <- prcomp(famd_df1, center = TRUE, scale. = TRUE)

# View summary of explained variance
summary(pca_result1)

# Scree plot to visualize variance explained
plot(pca_result1, type = "l")

# Tidy version of PCA result
pca_df1 <- as_tibble(pca_result1$x)

# Add metadata for plotting
pca_df1 <- bind_cols(data.frame(strain = row.names(famd_df1)) %>% mutate(lifestyle = ifelse(str_detect(strain, "s"), "S", "FL")), pca_df1)

# Plot first two PCs
ggplot(pca_df1, aes(x = PC1, y = PC2, colour = lifestyle, label = strain)) +
  labs(x = "PC1", y = "PC2") +
  geom_label() +
  theme_bw()

# Get loadings
pca_result1$rotation

# Get scores for each PC by type
type_scores1 <- augment(pca_result1) %>%
  bind_cols(data.frame(strain = row.names(famd_df1)) %>% mutate(lifestyle = ifelse(str_detect(strain, "s"), "S", "FL"))) %>%
  dplyr::select(2:4, strain, lifestyle) %>%
  unique()

# Get variable loadings
loadings1 <- tidy(pca_result1, matrix = "rotation") %>%
  pivot_wider(names_from = PC, values_from = value, names_prefix = "PC") %>%
  mutate(column = case_when(
    column == "Maltose_mean" ~ "Average maltose release",
    column == "Maltose_slope" ~ "Maltose release pH sensitivity",
    column == "Glucose_mean" ~ "Average glucose release",
    column == "Glucose_slope" ~ "Glucose release pH sensitivity",
    column == "ph_mean" ~ "Average growth rate",
    column == "ph_slope" ~ "Growth rate pH sensitivity",
    column == "FvFm" ~ "Fv/Fm",
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
ggplot(type_scores1, aes(.fittedPC1, .fittedPC2)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings1, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings1, aes(x = PC1 * 10.2, y = PC2 * 10.2, label = column), size = 3) +
  labs(x = "PC1 (34%)", y = "PC2 (19%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# PERMANOVA to test whether lifestyle explains the clustering in the multivariate space
famd_df1_perm <- famd_df1 %>%
  mutate(strain = row.names(famd_df1),
         lifestyle = ifelse(str_detect(strain, "s"), "S", "FL"))

famd_df1_perm

X1 <- scale(famd_df1_perm[,-(16:17)])

adonis2(X1 ~ lifestyle, method = "euclidean", data = famd_df1_perm)


# Reanalysis using only the 92 strains from the phylogeny to check robustness of results ###########################
# Put together and do PCA -- with combo data frames using both slopes and means where appropriate
mult_df_combo_92 <- mult_df_sug_combo %>%
  filter(strain != "192") %>%
  merge(mult_df_n, by = c("strain", "lifestyle")) %>%
  merge(mult_df_phgr_combo, by = c("strain", "lifestyle")) %>%
  merge(mult_df_fl, by = c("strain", "lifestyle"))


glimpse(mult_df_combo_92)

# PCA of strains 
mult_df_scaled_combo_92 <- mult_df_combo_92 %>%
  dplyr::select(-strain, -lifestyle) %>%
  drop_na() %>%                # remove rows with missing values
  scale()                      # standardize to mean 0, sd 1

# Run PCA
pca_result_combo_92 <- prcomp(mult_df_scaled_combo_92, center = FALSE, scale. = FALSE)

# View summary of explained variance
summary(pca_result_combo_92)

# Scree plot to visualize variance explained
plot(pca_result_combo_92, type = "l")

# Tidy version of PCA result
pca_df_combo_92 <- as_tibble(pca_result_combo_92$x)

# Add metadata for plotting
pca_df_combo_92 <- bind_cols(mult_df_combo_92 %>% drop_na() %>% dplyr::select(strain, lifestyle), pca_df_combo_92)

# Plot first two PCs
ggplot(pca_df_combo_92, aes(x = PC1, y = PC2, colour = lifestyle, label = strain)) +
  labs(x = "PC1", y = "PC2") +
  geom_label() +
  theme_bw()

# Get loadings
pca_result_combo_92$rotation

# Get scores for each PC by type
type_scores_combo_92 <- augment(pca_result_combo_92) %>%
  bind_cols(mult_df_combo_92) %>%
  dplyr::select(2:4, strain, lifestyle) %>%
  unique()

# Get variable loadings
loadings_combo_92 <- tidy(pca_result_combo_92, matrix = "rotation") %>%
  pivot_wider(names_from = PC, values_from = value, names_prefix = "PC") %>%
  mutate(column = case_when(
    column == "Maltose_mean" ~ "Average maltose release",
    column == "Maltose_slope" ~ "Maltose release pH sensitivity",
    column == "Glucose_mean" ~ "Average glucose release",
    column == "Glucose_slope" ~ "Glucose release pH sensitivity",
    column == "ph_mean" ~ "Average growth rate",
    column == "ph_slope" ~ "Growth rate pH sensitivity",
    column == "FvFm" ~ "Fv/Fm",
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
ggplot(type_scores_combo_92, aes(.fittedPC1, .fittedPC2)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings_combo_92, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo_92, aes(x = PC1 * 10.2, y = PC2 * 10.2, label = column), size = 3) +
  labs(x = "PC1 (44%)", y = "PC2 (13%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# PC1 v PC3
ggplot(type_scores_combo_92, aes(.fittedPC1, .fittedPC3)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings_combo_92, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC3 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo_92, aes(x = PC1 * 10.2, y = PC3 * 10.2, label = column), size = 3) +
  labs(x = "PC1 (44%)", y = "PC3 (10%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# PC2 v PC3
ggplot(type_scores_combo_92, aes(.fittedPC2, .fittedPC3)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings_combo_92, aes(x = 0, y = 0, xend = PC2 * 10, yend = PC3 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo_92, aes(x = PC2 * 10.2, y = PC3 * 10.2, label = column), size = 3) +
  labs(x = "PC2 (13%)", y = "PC3 (10%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# Joining with environmental data for overlay
type_scores_combo_92 <- type_scores_combo_92 %>%
  left_join(read_csv("Data/df_tidy_env.csv"), by = c("strain", "lifestyle")) %>%
  # Adding phylotypes
  left_join(read_csv("Data/phylotypes.csv"), by = "strain") %>%
  mutate(phylotype = as.character(phylotype))


overlay_biplots2_92 <- function(overlay, title) {
  
  overlay_var <- deparse(substitute(overlay))
  is_continuous <- is.numeric(type_scores_combo[[overlay_var]])
  
  base_plot <- function(x, y, xlab, ylab) {
    ggplot(type_scores_combo_92, aes_string(x, y)) +
      geom_label(aes_string(label = "strain", fill = overlay_var),
                 size = 3, colour = "white") +
      labs(x = xlab, y = ylab) +
      theme_bw() +
      theme(legend.position = "bottom")
  }
  
  p1 <- base_plot(".fittedPC1", ".fittedPC2", "PC1 (44%)", "PC2 (13%)")
  p2 <- base_plot(".fittedPC1", ".fittedPC3", "PC1 (44%)", "PC3 (10%)")
  p3 <- base_plot(".fittedPC2", ".fittedPC3", "PC2 (13%)", "PC3 (10%)")
  
  # Apply appropriate legend scaling
  if (is_continuous) {
    p1 <- p1 + guides(fill = guide_colorbar(barwidth = 20, barheight = 1))
  } else {
    # Get number of unique levels in the overlay variable
    n_levels <- length(unique(type_scores_combo[[overlay_var]]))
    
    # Adjust number of columns based on how many rows we want
    ncol_legend <- ceiling(n_levels / 3)
    
    p1 <- p1 + guides(
      fill = guide_legend(
        ncol = ncol_legend,
        byrow = TRUE
      )
    )
  }
  
  # Extract legend
  legend <- cowplot::get_legend(
    p1 + theme(legend.position = "top") + labs(fill = paste(title))
  )
  
  # Remove legends from the individual plots
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- p3 + theme(legend.position = "none")
  
  aligned <- cowplot::align_plots(p1, p2, p3, align = "hv", axis = "tblr")
  
  p <- cowplot::plot_grid(
    cowplot::plot_grid(plotlist = aligned, ncol = 3, align = "v"),
    legend,
    ncol = 1,
    rel_heights = c(1, 0.35)
  )
  
  return(p)
}

# Origin pond
p.pond_id_92 <- overlay_biplots2_92(pond_id, "Pond ID")

# Pond pH
p.pond_ph_92 <- overlay_biplots2_92(ph, "Pond pH")

# Pond area
p.pond_area_92 <- overlay_biplots2_92(pond_area, "Pond Area")

# Pond NO3
p.pond_no3_92 <- overlay_biplots2_92(no3, "Pond Nitrate Content")

# Pond TDS
p.pond_tds_92 <- overlay_biplots2_92(tds, "Pond TDS")

# Pond conductivity
p.pond_conductivity_92 <- overlay_biplots2_92(conductivity, "Pond Conductivity")

# Phylotypes
overlay_biplots2_92(phylotype, "Phylotypes")

# Put together plots of interest

p.env_overlay_92 <- plot_grid(p.pond_ph_92, p.pond_no3_92, p.pond_tds_92, p.pond_conductivity_92, p.pond_id_92,
                           labels = c("a)", "b)", "c)", "d)", "e)"), 
                           label_size = 14,
                           label_fontface = "bold",
                           ncol = 1,
                           scale = 1)

p.env_overlay_92


# Creating a new data frame for statistical testing of the phenotypic separation we see in the PCA
rf_df_92 <- mult_df_combo_92 %>%
  mutate(lifestyle = as.factor(lifestyle)) %>%
  rename("Autonomous_growth_mean" = ph_mean,
         "Autonomous_growth_slope" = ph_slope) %>%
  left_join(read_csv("Data/df_tidy_env.csv"), by = c("strain", "lifestyle")) %>%
  dplyr::select(-strain)

rf_df_92

# PERMANOVA to test whether lifestyle explains the clustering in the multivariate space
# and whether there's an impact of the environmental variables
X_92 <- scale(rf_df_92[,-c(1, 17, 18, 19, 20, 21, 22)])
adonis2(X_92 ~ lifestyle, method = "euclidean", data = rf_df_92, by = "term")
adonis2(X_92 ~ lifestyle*pond_id, method = "euclidean", data = rf_df_92, by = "term")
adonis2(X_92 ~ lifestyle*no3, method = "euclidean", data = rf_df_92, by = "term")
adonis2(X_92 ~ lifestyle*tds, method = "euclidean", data = rf_df_92, by = "term")
adonis2(X_92 ~ lifestyle*ph, method = "euclidean", data = rf_df_92, by = "term")
adonis2(X_92 ~ lifestyle*conductivity, method = "euclidean", data = rf_df_92, by = "term")

# Checking betadispersion
anova(betadisper(dist(X_92), rf_df_92$lifestyle))
anova(betadisper(dist(X_92), rf_df_92$pond_id))
anova(betadisper(dist(X_92), rf_df_92$no3))
anova(betadisper(dist(X_92), rf_df_92$tds))
anova(betadisper(dist(X_92), rf_df_92$ph))
anova(betadisper(dist(X_92), rf_df_92$conductivity))

# Using envfit to test which trait vectors contribute to this separation
pca <- prcomp(rf_df_92[,-c(1, 17, 18, 19, 20, 21, 22)], scale. = TRUE, center = TRUE)

# Fit trait vectors to PCA space (note this is basically a linear correlation with the multivariate space)
fit <- envfit(pca, rf_df_92[,-c(1, 17, 18, 19, 20, 21, 22)], permutations = 999)

# View results
fit

# Random forest model to tell us which traits are most predictive for a given lifestyle
set.seed(123)  # for reproducibility
rf_df_92$lifestyle <- factor(rf_df_92$lifestyle)

rf_model_92 <- randomForest(lifestyle ~ ., data = rf_df_92[,-c(17, 18, 19, 20, 21, 22)], 
                         importance = TRUE, 
                         ntree = 1000)

# Inspect the model output
print(rf_model_92)

# Rank the most important predictors
varImpPlot(rf_model_92, n.var = 15, main = "Top Trait Predictors of Lifestyle")

# Conclusions unchanged


## FAMD to test inclusion of phylotypes in the multivariate model
# Make a new data frame
famd_df <- mult_df_combo %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.)))) %>%
  left_join(read_csv("Data/phylotypes.csv") %>% filter(strain != 126), by = "strain") %>%
  drop_na(phylotype) %>% 
  column_to_rownames("strain") %>%
  dplyr::select(-lifestyle, -strain_long) %>%
  mutate(phylotype = as.factor(phylotype))
  
glimpse(famd_df)

# Fit the model, get some plots
res <- FAMD(
  famd_df,
  ncp = 16             # number of dimensions to keep
)

# Look at important variables from model
# Eigenvalues/variance
res$eig
# These are the coordinates for our strains
res$ind$coord
plot(res, axes = c(1, 2))

# Contribution of each variable
res$var$contrib

# factoextra plot
fviz_famd_ind(
  res,
  habillage = "phylotype",
  addEllipses = TRUE
)

# Recreating the PCA plot
famd_scores <- as.data.frame(res$ind$coord[, 1:2]) %>%
  mutate(strain = row.names(res$ind$coord),
         lifestyle = ifelse(str_detect(strain, "s"), "S", "FL")) %>%
  left_join(read_csv("Data/phylotypes.csv") %>% filter(strain != 126), by = "strain")
famd_scores

famd_loadings <- res$var$coord # I know they're not actually loadings but it's close as we can get
famd_loadings

ggplot(famd_scores, aes(-Dim.1, Dim.2)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3) +
  geom_segment(data = famd_loadings, aes(x = 0, y = 0, xend = -Dim.1 * 10, yend = Dim.2 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = famd_loadings, aes(x = -Dim.1 * 10.2, y = Dim.2 * 10.2, label = row.names(famd_loadings)), size = 3) +
  stat_ellipse(aes(group = phylotype, colour = factor(phylotype)), type = "norm", linewidth = 0.8) +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")


# PERMANOVA for phylotype
famd_df_perm <- famd_df %>%
  mutate(strain = row.names(famd_df),
         lifestyle = ifelse(str_detect(strain, "s"), "S", "FL"))

famd_df_perm

X_famd <- scale(famd_df_perm[,-c(16:18)])

adonis2(X_famd ~ phylotype + lifestyle, method = "euclidean", data = famd_df_perm, by = "term")

anova(betadisper(dist(X_famd), famd_df_perm$phylotype))

pca_famd <- prcomp(famd_df[, -16], center = TRUE, scale. = TRUE)

axis_df <- data.frame(
  PCA1  = pca_famd$x[,1],
  FAMD1 = res$ind$coord[,1],
  PCA2  = pca_famd$x[,2],
  FAMD2 = res$ind$coord[,2],
  phylogroup = famd_df$phylotype
)

# How phylogroup fits to axis 1
anova(lm(FAMD1 ~ phylogroup, data = axis_df))
anova(lm(PCA1 ~ phylogroup, data = axis_df))

# Correlation between axes in the two models
cor(axis_df$PCA1, axis_df$FAMD1)
cor(axis_df$PCA2, axis_df$FAMD2)


# PCA on just phylogroup 1
famd_df1 <- famd_df %>%
  filter(phylotype == "Phylogroup 1") %>%
  dplyr::select(-phylotype)

# Run PCA
pca_result1 <- prcomp(famd_df1, center = TRUE, scale. = TRUE)

# View summary of explained variance
summary(pca_result1)

# Scree plot to visualize variance explained
plot(pca_result1, type = "l")

# Tidy version of PCA result
pca_df1 <- as_tibble(pca_result1$x)

# Add metadata for plotting
pca_df1 <- bind_cols(data.frame(strain = row.names(famd_df1)) %>% mutate(lifestyle = ifelse(str_detect(strain, "s"), "S", "FL")), pca_df1)

# Plot first two PCs
ggplot(pca_df1, aes(x = PC1, y = PC2, colour = lifestyle, label = strain)) +
  labs(x = "PC1", y = "PC2") +
  geom_label() +
  theme_bw()

# Get loadings
pca_result1$rotation

# Get scores for each PC by type
type_scores1 <- augment(pca_result1) %>%
  bind_cols(data.frame(strain = row.names(famd_df1)) %>% mutate(lifestyle = ifelse(str_detect(strain, "s"), "S", "FL"))) %>%
  dplyr::select(2:4, strain, lifestyle) %>%
  unique()

# Get variable loadings
loadings1 <- tidy(pca_result1, matrix = "rotation") %>%
  pivot_wider(names_from = PC, values_from = value, names_prefix = "PC") %>%
  mutate(column = case_when(
    column == "Maltose_mean" ~ "Average maltose release",
    column == "Maltose_slope" ~ "Maltose release pH sensitivity",
    column == "Glucose_mean" ~ "Average glucose release",
    column == "Glucose_slope" ~ "Glucose release pH sensitivity",
    column == "ph_mean" ~ "Average growth rate",
    column == "ph_slope" ~ "Growth rate pH sensitivity",
    column == "FvFm" ~ "Fv/Fm",
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
ggplot(type_scores1, aes(.fittedPC1, .fittedPC2)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings1, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings1, aes(x = PC1 * 10.2, y = PC2 * 10.2, label = column), size = 3) +
  labs(x = "PC1 (34%)", y = "PC2 (19%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle"
  ) +
  theme(legend.position = "bottom")

# PERMANOVA to test whether lifestyle explains the clustering in the multivariate space
famd_df1_perm <- famd_df1 %>%
  mutate(strain = row.names(famd_df1),
         lifestyle = ifelse(str_detect(strain, "s"), "S", "FL"))

famd_df1_perm

X1 <- scale(famd_df1_perm[,-(16:17)])

adonis2(X1 ~ lifestyle, method = "euclidean", data = famd_df1_perm)


