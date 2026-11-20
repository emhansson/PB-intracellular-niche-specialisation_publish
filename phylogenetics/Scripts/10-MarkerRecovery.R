library(tidyverse)
library(gridExtra)

# As the extraction test was really successful, markers were extracted for all samples regardless of assembly quality

# Importing qc output 
df <- read_tsv("Data/Input/master_assembly_qc_table.tsv")

# Importing marker summary (this is for field strains only, markers for references were extracted separately)
markers <- read_tsv("Data/Input/blast_marker_summary.tsv") %>%
  # Adding metadata and tidying up names
  mutate(strain_no = if_else(
    str_detect(sample_id, "-s-"),
    paste0(
      "s",
      str_extract(sample_id, "(?<=-s-)\\d+(?=_getorganelle)")
    ),
    str_extract(sample_id, "(?<=-FL-)\\d+(?=_getorganelle)")
  ),
  lifestyle = ifelse(str_detect(sample_id, "-s-") == TRUE, "s", "FL"),
  sample_id = str_remove(sample_id, "_getorganelle")) 

glimpse(markers) 
glimpse(df)

strains <- unique(markers$strain_no)
strains

# Inspecting marker recovery by sample
markers %>%
  ggplot(aes(y = top_identity, x = marker, colour = top_contig)) +
  facet_wrap(~sample_id) +
  geom_point() +
  geom_hline(aes(yintercept = 90)) +
  theme(legend.position = "none")

markers %>%
  filter(lifestyle == "FL") %>%
  ggplot(aes(y = top_identity, x = marker, colour = top_contig)) +
  facet_wrap(~strain_no*sample_id) +
  geom_point() +
  geom_hline(aes(yintercept = 90)) +
  theme(legend.position = "none")

markers %>%
  filter(lifestyle == "s") %>%
  ggplot(aes(y = top_identity, x = marker, colour = top_contig)) +
  facet_wrap(~strain_no*sample_id) +
  geom_point() +
  geom_hline(aes(yintercept = 90)) +
  theme(legend.position = "none")

markers %>%
  filter(strain_no == strains[1]) %>%
  ggplot(aes(y = top_identity, x = marker, colour = top_contig)) +
  facet_wrap(~sample_id) +
  geom_point() +
  geom_hline(aes(yintercept = 90))

markers %>%
  filter(strain_no == strains[2]) %>%
  ggplot(aes(y = top_identity, x = marker, colour = top_contig)) +
  facet_wrap(~sample_id) +
  geom_point() +
  geom_hline(aes(yintercept = 90))

markers %>%
  filter(strain_no == strains[3]) %>%
  ggplot(aes(y = top_identity, x = marker, colour = top_contig)) +
  facet_wrap(~sample_id) +
  geom_point() +
  geom_hline(aes(yintercept = 90))

markers %>%
  filter(strain_no == "192") %>%
  ggplot(aes(y = top_identity, x = marker, colour = top_contig)) +
  facet_wrap(~sample_id) +
  geom_point() +
  geom_hline(aes(yintercept = 90))

markers %>%
  filter(strain_no == "126") %>%
  ggplot(aes(y = top_identity, x = marker, colour = top_contig)) +
  facet_wrap(~sample_id) +
  geom_point() +
  geom_hline(aes(yintercept = 90))

# Recovery by marker -- based on this I decided to exclude rpoA, rpoC1, rps11
markers %>%
  drop_na() %>%
  group_by(marker) %>%
  summarise(
    recovery = n()
  ) %>%
  print(n = nrow(.))

# Summary of marker recovery by sample
marker_summary <- markers %>%
  drop_na() %>%
  filter(marker != "rpoA",
         marker != "rpoC1",
         marker != "rps11") %>%
  group_by(sample_id, strain_no, assembly_type) %>%
  summarise(
    markers_recovered = n(),
    mean_id = mean(top_identity),
    min_id = min(top_identity),
    sd_id = sd(top_identity),
    n_unique_contigs = n_distinct(top_contig)
  ) %>%
  arrange(n_unique_contigs) %>%
  print(n = nrow(.))

marker_summary

# Merging the marker hits with assembly qc data frame
df2 <- df %>% 
  merge(markers, by = c("sample_id", "strain_no", "lifestyle")) 

df2

# Best samples determined based on marker recovery, number of unique contigs, mean_id, N50, assembly type and qc batch priority
best_samples <- marker_summary %>%
  left_join(df, by = c("sample_id", "strain_no")) %>%
  mutate(
    batch_priority = ifelse(str_detect(batch, "_x$"), 1, 0)
  ) %>%
  arrange(
    strain_no,
    desc(markers_recovered),
    n_unique_contigs,
    desc(mean_id),
    desc(N50),
    assembly_type,
    batch_priority
  ) %>%
  group_by(strain_no) %>%
  slice(1) %>%
  ungroup()

best_samples

# Making a final data frame highlighting which are the best samples
final_df <- df2 %>%
  mutate(best_sample = ifelse(sample_id %in% best_samples$sample_id, "yes", "no"))

# Saving this
write_csv(final_df, "Data/Output/final_sample_list.csv")

# Marker summary for best samples
marker_summary2 <- final_df %>%
  filter(best_sample == "yes",
         n_hits > 0) %>%
  group_by(marker) %>%
  summarise(
    n_recovered = n(),
    prop_recovered = n() / n_distinct(final_df$strain_no),
    mean_identity = mean(top_identity),
    min_identity = min(top_identity),
    mean_length = mean(top_length),
    sd_identity = sd(top_identity)
  ) %>%
  arrange(desc(prop_recovered))

marker_summary2 %>%
  print(n = nrow(.))

# top_identity for each marker by lifestyle
final_df %>%
  filter(best_sample == "yes",
         n_hits > 0,
         marker != "rpoA",
         marker != "rpoC1",
         marker != "rps11") %>%
  ggplot(aes(marker, top_identity, colour = lifestyle)) +
  geom_jitter(alpha = 0.5)

# Marker recovery by lifestyle
final_df %>%
  filter(best_sample == "yes",
         n_hits > 0,
         marker != "rpoA",
         marker != "rpoC1",
         marker != "rps11") %>%
  count(marker, lifestyle)

# Marker recovery by strain
strain_summary <- final_df %>%
  filter(best_sample == "yes",
         n_hits > 0,
         marker != "rpoA",
         marker != "rpoC1",
         marker != "rps11") %>%
  group_by(strain_no) %>%
  summarise(
    recovered = n(),
    mean_identity = mean(top_identity)
  )

strain_summary %>%
  print(n = nrow(.))

# Final list of samples to use (excluding 192 for poor overall recovery 
# and suspiciously low identity for some of the markers that were recovered even if others matched the average pattern)
sample_list <- final_df %>%
  filter(best_sample == "yes",
         strain_no != "192") %>%
  mutate(file_name = paste0(sample_id, "_getorganelle")) %>%
  select(file_name) %>%
  unique()

sample_list

# Final list of markers to use
marker_list <- final_df %>%
  filter(best_sample == "yes",
         n_hits > 0,
         marker != "rpoA",
         marker != "rpoC1",
         marker != "rps11") %>%
  select(marker) %>%
  unique()

# Saving for use in scripts
write.table(sample_list, "Data/Output/best_samples.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(marker_list, "Data/Output/core_markers.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)