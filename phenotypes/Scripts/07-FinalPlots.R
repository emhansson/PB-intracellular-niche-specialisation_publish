# Erika M Hansson 2025
# erika.hansson@manchester.ac.uk

# All the plots in one place + code to export as TIFF
# NB this script only works when all the objects from the related scripts are loaded
library(ragg)
library(cowplot)

# Plot grid of FvFm, GR, sugar, sugar v gr, N data + N-loss PCA
p.fvfm <- df_tidy_fl %>%
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
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) +
  labs(x = "Strain", y = "Mean Fv/Fm", colour = "Lifestyle")

p.gr <- df_tidy_phgr %>%
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

p.sug_hu <- ggplot(prop_zero_df, aes(x = factor(ph), y = 1-prop_zero)) +
  geom_col(aes(fill = lifestyle), width = 0.6, alpha = 0.5) +
  geom_point(data = hurdle_df, aes(colour = lifestyle), size = 2) +
  geom_errorbar(data = hurdle_df, aes(ymin = 1-lower, ymax = 1-upper, 
                                      colour = lifestyle), width = 0.2, linewidth = 1) +
  labs(
    y = "Proportion of non-zero values",
    x = "pH"
  ) +
  theme_bw () +
  theme(legend.position = "none", strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values =
                      c("Free-living" = "green3",
                        "Symbiotic" = "cornflowerblue")) +
  scale_colour_manual(values =
                        c("Free-living" = "green3",
                          "Symbiotic" = "cornflowerblue")) +
  facet_grid(lifestyle~sugar)

p.sug_mu <- cond_df %>%
  ggplot(aes(x = factor(ph), y = estimate, colour = lifestyle)) +
  geom_point() +
  geom_jitter(data = nonzero_df %>% mutate(estimate = release), alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), linewidth = 1) +
  facet_wrap(~sugar) +
  scale_colour_manual(values =
                        c("Free-living" = "green3",
                          "Symbiotic" = "cornflowerblue")) +
  theme_bw() +
  labs(y = "Sugar release (ppm)",
       x = "pH",
       colour = "Lifestyle") +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "white"))

p.gr_v_sug <- plot_df %>%
  ggplot(aes(x = mean_gr, y = mean)) +
  geom_point(aes(color = lifestyle, shape = factor(ph))) +
  geom_smooth(aes(linetype = factor(ph)), method = "lm", se = FALSE, colour = "black", linewidth = 0.5) + # Lazy fit by ph
  facet_wrap(~sugar) +
  theme_bw() +
  labs(x = "Mean growth rate", y = "Mean sugar release", colour = "Lifestyle", shape = "pH", linetype = "pH") +
  scale_colour_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c(
      "FL" = "Free-living",
      "S" = "Symbiotic"),
    guide = "none") +
  scale_shape_manual(values = c(
    "5.5" = 16,
    "7.6" = 1)) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "white"))

# Nitrogen use plot
library(ggtext)

n_plot <- df_tidy_n %>%
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
    pattern_density = 0.01,
    pattern_fill = "white",
    colour = "white",
    size = 0.2
  ) +
  facet_grid(
    n_type*n_source ~ lifestyle,
    labeller = labeller(
      n_type = as_labeller(c(
        "Organic" = "<span style='color:#A33D00; font-weight:bold;'>Organic</span>",
        "Inorganic" = "<span style='color:#54204F; font-weight:bold;'>Inorganic</span>"
      ))
    )
  ) +
  scale_fill_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue")
  ) +
  scale_pattern_manual(values = c("Inorganic" = "stripe", "Organic" = "none")) +
  labs(y = "Count", x = "Mean growth rate per day") +
  guides(fill = "none", pattern = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        strip.background = element_rect(fill = "white"),
        strip.text = element_markdown())

n_plot

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
  geom_label_repel(
    data = loadings_n %>%
      mutate(
        names = str_extract(column, "^[^ ]+"),
        n_colour = ifelse(
          grepl("Ammonium|Nitrate|Urea", names),
          "Inorganic",
          "Organic"
        )
      ),
    aes(x = PC1 * 10.2, y = PC2 * 10.2,
        label = names, colour = n_colour),
    size = 3,
    alpha = 0.9
  ) +
  scale_colour_manual(
    values = c(
      "Organic" = "#A33D00",
      "Inorganic" = "#54204F"
    ),
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_colour() +
  geom_label(aes(label = strain, fill = factor(n_lost), colour = lifestyle), size = 3) +
  scale_colour_manual(values = c("FL" = "white",
                                 "S" = "black"),
                      name = "Lifestyle",
                      labels = c("FL" = "Free-living",
                                 "S" = "Symbiotic")) +
  labs(x = "PC1 (55%)", y = "PC2 (14%)",
       fill = "# of N sources without growth") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(fill = "grey70")))


n_pca
n_plots <- plot_grid(
  n_plot,
  n_pca,
  ncol = 2,
  rel_widths = c(0.5, 1),
  labels = c("F", "G")
)

n_plots

legend <- get_legend(p.gr + theme(legend.position = "top"))

# Remove legends from the individual plots
p1 <- p.gr + theme(legend.position = "none")
p2 <- p.fvfm + theme(legend.position = "none")
p3 <- p.sug_hu + theme(legend.position = "none")
p4 <- p.sug_mu + theme(legend.position = "none")
p5 <- p.gr_v_sug #+ theme(legend.position = "none")

# Align plots (ensures axes & facets line up neatly)
aligned <- align_plots(p1, p2, p3, p4, p5, align = "hv", axis = "tblr")

# Stack the legend, p1 and p2 horizontally, then p3, 4, 5 horizontally
gr_sug_fl_n_grid <- plot_grid(legend,
  plot_grid(p1, p2, #aligned[[1]], aligned[[2]], 
            ncol = 2, 
            #align = "h", 
            rel_widths = c(1, 2),
            labels = c("A", "B")),
  plot_grid(p3, p4, p5, #aligned[[3]], aligned[[4]], aligned[[5]], 
            ncol = 3, 
            #align = "h",
            labels = c("C", "D", "E")),
  ncol = 1,
  rel_heights = c(0.15, 1, 1, 2.2),
  n_plots
)

gr_sug_fl_n_grid

# PCA plot

loadings_combo

loadings_combo <- loadings_combo %>%
  mutate(labels = case_when(column == "Average growth rate" |
                              column == "Average glucose release" |
                              column == "Glucose release pH sensitivity" |
                              column == "Growth rate pH sensitivity" |
                              column == "Glutamate growth" |
                              column == "Fv/Fm" ~ "difficult",
                            TRUE ~ "easy"))

loadings_combo_easy <- loadings_combo %>% filter(labels == "easy")

loadings_combo_difficult <- loadings_combo %>% filter(labels == "difficult")
loadings_combo_difficult$PC1[loadings_combo_difficult$column == "Average growth rate"] <- -0.274 - 0.08
loadings_combo_difficult$PC2[loadings_combo_difficult$column == "Average growth rate"] <- -6.42e-3 - 0.03
loadings_combo_difficult$PC1[loadings_combo_difficult$column == "Fv/Fm"] <- 0.0934 + 0.02
loadings_combo_difficult$PC1[loadings_combo_difficult$column == "Glutamate growth"] <- 0.141 + 0.035
loadings_combo_difficult$PC2[loadings_combo_difficult$column == "Glutamate growth"] <- 8.19e-2 - 0.01
loadings_combo_difficult$PC1[loadings_combo_difficult$column == "Average glucose release"] <- 0.324
loadings_combo_difficult$PC2[loadings_combo_difficult$column == "Average glucose release"] <- 2.47e-1 + 0.005
loadings_combo_difficult$PC2[loadings_combo_difficult$column == "Glucose release pH sensitivity"] <- 2.39e-1 - 0.02
loadings_combo_difficult$PC1[loadings_combo_difficult$column == "Glucose release pH sensitivity"] <- 0.318 - 0.02
loadings_combo_difficult$PC1[loadings_combo_difficult$column == "Growth rate pH sensitivity"] <- 0.210 + 0.08
loadings_combo_difficult$PC2[loadings_combo_difficult$column == "Growth rate pH sensitivity"] <- 1.31e-1

set.seed(123)
p.pca <- ggplot(type_scores_combo %>%
                  mutate(cluster = case_when(strain == "25" | strain == "63" | strain == "34" | strain == "1" |
                                               strain == "270" | strain == "29" | strain == "318" | strain == "273" |
                                               strain == "27" | strain == "31" | strain == "126" ~ "FL (top cluster)",
                                             strain != "s70" & lifestyle == "S" ~ "S (top cluster)",
                                             TRUE ~ lifestyle)), aes(.fittedPC1, .fittedPC2)) +
  geom_segment(data = loadings_combo, aes(x = 0, y = 0, xend = PC1 * 12, yend = PC2 * 12),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label(aes(label = strain, fill = lifestyle, colour = cluster), size = 3) +
  geom_label_repel(data = loadings_combo_easy, aes(x = PC1 * 12, y = PC2 * 12, label = column), size = 3, alpha = 0.8) +
  geom_label_repel(data = loadings_combo_difficult, aes(x = PC1 * 12, y = PC2 * 12, label = column), size = 3, alpha = 0.8) +
  labs(x = "PC1 (44%)", y = "PC2 (13%)") +
  theme_bw() +
  scale_fill_manual(values = c(
    "FL" = "green3",
    "S" = "cornflowerblue"),
    labels = c("FL" = "Free-living",
               "S" = "Symbiotic"),
    name = "Lifestyle") +
  scale_colour_manual(values = c(
    "S (top cluster)" = "#000033",
    "FL (top cluster)" = "#002B1A",
    "S" = "lightcyan",
    "FL" = "darkseagreen1"
  ), guide = "none") +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

p.pca


# Saving plots function
save_plot_custom <- function(plot, width, height,
                             format = c("png", "pdf", "tiff"), units = "px", res = 300, path = ".") {
  format <- match.arg(format)
  
  
  # Get the name of the plot object
  plot_name <- deparse(substitute(plot))
  filename <- file.path(path, paste0("Plots/Final/", plot_name, ".", format))
  
  
  if (format == "png") {
    agg_png(filename, width = width, height = height, units = units, res = res)
    print(plot)
    dev.off()
  } else if (format == "pdf") {
    agg_pdf(filename, width = width, height = height, units = units)
    print(plot)
    dev.off()
  } else if (format == "tiff") {
    agg_tiff(filename, width = width, height = height, units = units, res = res)
    print(plot)
    dev.off()
  }
  
  message("Saved plot as: ", filename)
}


#Submission guidelines
#Width: 789 – 2250 pixels (at 300 dpi). Height maximum: 2625 pixels (at 300 dpi).

scale <- 1.5
save_plot_custom(gr_sug_fl_n_grid, 2250*scale, 2625*scale, "tiff")
save_plot_custom(p.pca, 2200*scale, 2000*scale, "tiff")

# For supp mats, PC2 and PC3
# PC2 v PC3
p.pca2 <- ggplot(type_scores_combo %>%
                   mutate(cluster = case_when(strain == "25" | strain == "63" | strain == "34" | strain == "1" |
                                                strain == "270" | strain == "29" | strain == "318" | strain == "273" |
                                                strain == "27" | strain == "31" | strain == "126" ~ "FL (top cluster)",
                                              strain != "s70" & lifestyle == "S" ~ "S (top cluster)",
                                              TRUE ~ lifestyle)), aes(.fittedPC2, .fittedPC3)) +
  geom_label(aes(label = strain, fill = lifestyle, colour = cluster), size = 3) +
  geom_segment(data = loadings_combo, aes(x = 0, y = 0, xend = PC2 * 10, yend = PC3 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo, aes(x = PC2 * 10.2, y = PC3 * 10.2, label = column), size = 3, alpha = 0.8) +
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
    "S (top cluster)" = "#000033",
    "FL (top cluster)" = "#002B1A",
    "S" = "lightcyan",
    "FL" = "darkseagreen1"
  ), guide = "none") +
  theme(legend.position = "bottom")

# PC1 v PC3
p.pca3 <- ggplot(type_scores_combo %>%
                   mutate(cluster = case_when(strain == "25" | strain == "63" | strain == "34" | strain == "1" |
                                                strain == "270" | strain == "29" | strain == "318" | strain == "273" |
                                                strain == "27" | strain == "31" | strain == "126" ~ "FL (top cluster)",
                                              strain != "s70" & lifestyle == "S" ~ "S (top cluster)",
                                              TRUE ~ lifestyle)), aes(.fittedPC1, .fittedPC3)) +
  geom_label(aes(label = strain, fill = lifestyle, colour = cluster), size = 3) +
  geom_segment(data = loadings_combo, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC3 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings_combo, aes(x = PC1 * 10.2, y = PC3 * 10.2, label = column), size = 3, alpha = 0.8) +
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
    "S (top cluster)" = "#000033",
    "FL (top cluster)" = "#002B1A",
    "S" = "lightcyan",
    "FL" = "darkseagreen1"
  ), guide = "none") +
  theme(legend.position = "bottom")


legend <- get_legend(p.pca + theme(legend.position = "top"))

# Remove legends from the individual plots
p2 <- p.pca2 + theme(legend.position = "none")
p3 <- p.pca3 + theme(legend.position = "none")

# Align plots (ensures axes & facets line up neatly)
aligned <- align_plots(p2, 
                       p3, 
                       align = "hv", axis = "tblr")

# Stack the legend, p2, and p3 vertically for the right column
bottom_row <- plot_grid(aligned[[1]], aligned[[2]], ncol = 2, align = "hv")

just_p2_p3 <-  pca_plot <- plot_grid(
  bottom_row, 
  ncol = 1,
  legend,
  rel_heights = c(1, 0.05)
)

just_p2_p3

save_plot_custom(just_p2_p3, 2200*scale, 2000*scale, "tiff")


# Supp mats -- Correlation of N growth loss
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

save_plot_custom(gg_corr, 1000*scale, 1000*scale, "tiff")

# Supp mats -- environmental overlay PCA
overlay_biplots2 <- function(overlay, title) {
  
  overlay_var <- deparse(substitute(overlay))
  is_continuous <- is.numeric(type_scores_combo[[overlay_var]])
  
  base_plot <- function(x, y, xlab, ylab) {
    ggplot(type_scores_combo, aes_string(x, y)) +
      geom_point(aes_string(colour = overlay_var),
                 size = 3, alpha = 0.7) +
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
    p1 + theme(legend.position = "top") + labs(colour = paste(title))
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
                           labels = "AUTO", 
                           label_size = 14,
                           label_fontface = "bold",
                           ncol = 1,
                           scale = 1)

p.env_overlay

save_plot_custom(p.env_overlay, 2250*scale, 3200*scale, "tiff")


# RF imp plot

imp_df <- as.data.frame(importance(rf_model))

imp_df$variable <- row.names(imp_df)
imp_df$variable <- gsub("_", " ", imp_df$variable)

imp_dfA <- imp_df %>%
  dplyr::select(variable, MeanDecreaseAccuracy) %>%
  mutate(variable = as.factor(variable))

imp_dfA$variable <- reorder(imp_dfA$variable, imp_dfA$MeanDecreaseAccuracy)

imp_dfG <- imp_df %>%
  dplyr::select(variable, MeanDecreaseGini) %>%
  mutate(variable = as.factor(variable))

imp_dfG$variable <- reorder(imp_dfG$variable, imp_dfG$MeanDecreaseGini)


rf_imp <- plot_grid(
  imp_dfA %>%
    ggplot(aes(x = MeanDecreaseAccuracy, y = variable)) +
    geom_point(size = 3) +
    theme_bw() +
    labs(y = ""),
  imp_dfG %>%
    ggplot(aes(x = MeanDecreaseGini, y = variable)) +
    geom_point(size = 3) +
    theme_bw() +
    labs(y = ""),
  ncol = 2)

rf_imp

save_plot_custom(rf_imp, 2250*scale, 1000*scale, "tiff")


# PCA for just phylogroup 1
pca1 <- ggplot(type_scores1, aes(.fittedPC1, .fittedPC2)) +
  geom_label(aes(label = strain, fill = lifestyle), size = 3, colour = "white") +
  geom_segment(data = loadings1, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black", linewidth = 1) +
  geom_label_repel(data = loadings1, aes(x = PC1 * 10.2, y = PC2 * 10.2, label = column), size = 3, alpha = 0.8) +
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


pca1

save_plot_custom(pca1, 2250*scale, 1000*scale, "tiff")