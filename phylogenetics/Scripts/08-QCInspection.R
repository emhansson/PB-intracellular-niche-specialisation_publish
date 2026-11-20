library(tidyverse)
library(gridExtra)

# Inspecting qc output 
df <- read_tsv("Data/Input/master_assembly_qc_table.tsv")

# Fixing the pond typo
df$pond[df$pond == "RX10"] <- "RIX10"

glimpse(df)

# Some plots to get a sense of distribution
p.total_length <- df %>%
ggplot(aes(x = total_length)) +
  geom_density()

p.n_contig <- df %>%
  ggplot(aes(x = n_contigs)) +
  geom_density()

p.largest_contig <- df %>%
  ggplot(aes(x = largest_contig)) +
  geom_density()

grid.arrange(p.total_length, p.n_contig, p.largest_contig, ncol = 3)

# Inspecting duplicates
dupes <- df %>%
  group_by(strain_no) %>%
  summarise(
    n = n()
  ) %>% 
  filter(n > 1)

df %>%
  filter(strain_no %in% dupes$strain_no) %>%
  select(strain_no, batch, n_contigs, total_length, largest_contig, N50) %>%
  arrange(strain_no) %>%
  print(n = nrow(.))

# Ranking based on n_contigs, N50 and largest_contig
# Classifying potential issues such as suspicious length, fragmentation and low continutity
# Classifying assembly qualities to get a sense of how many are useful for plastid phylogeny/marker recovery
ranked <- df %>%
  group_by(strain_no) %>%
  arrange(
    n_contigs,
    desc(N50),
    desc(largest_contig),
    .by_group = TRUE
  ) %>%
  mutate(
    assembly_rank = row_number()
  ) %>%
  ungroup() %>%
  mutate(
    suspicious_length =
      total_length > 200000 |
      total_length < 50000,
    
    fragmented =
      n_contigs > 100,
    
    low_continuity =
      largest_contig < 10000
  ) %>%
  mutate(
    assembly_quality = case_when(
      
      largest_contig >= 80000 &
        n_contigs <= 10 ~ "excellent",
      
      largest_contig >= 80000 &
        n_contigs <= 500 ~ "likely_usable",
      
      largest_contig >= 50000 ~ "possible_marker_recovery",
      
      TRUE ~ "poor"
    )
  )

ranked %>%
  select(strain_no, batch, n_contigs, total_length, largest_contig, N50, assembly_rank) %>%
  arrange(strain_no) %>%
  print(n = length(df$strain_no))

# Some more plotting to see how the distibutions change if we just use the top ranked duplicates
p.total_length_ranked <- ranked %>%
  filter(assembly_rank == 1) %>%
  ggplot(aes(x = total_length)) +
  geom_density()

p.n_contig_ranked <- ranked %>%
  filter(assembly_rank == 1) %>%
  ggplot(aes(x = n_contigs)) +
  geom_density()

p.largest_contig_ranked <- ranked %>%
  filter(assembly_rank == 1) %>%
  ggplot(aes(x = largest_contig)) +
  geom_density()


grid.arrange(p.total_length, p.n_contig, p.largest_contig, 
             p.total_length_ranked, p.n_contig_ranked, p.largest_contig_ranked, ncol = 3)

# How does largest contig size and number of contigs correlate?
ranked %>%
  ggplot(aes(
    x = n_contigs,
    y = largest_contig,
    colour = assembly_rank == 1
  )) +
  geom_point()

# Summary table of assembly quality
ranked %>%
  filter(assembly_rank == 1) %>%
  count(assembly_quality) %>%
  mutate(
    proportion = n / sum(n)
  )

# Let's look at the poor one
ranked %>%
  filter(
    assembly_rank == 1,
    assembly_quality == "poor"
  ) %>%
  select(
    strain_no,
    batch,
    n_contigs,
    total_length,
    largest_contig,
    N50
  ) %>%
  arrange(desc(n_contigs)) %>%
  print(n = nrow(.))

# n_contigs by assembly_quality
ranked %>%
  filter(assembly_rank == 1) %>%
  ggplot(aes(
    x = assembly_quality,
    y = n_contigs
  )) +
  geom_boxplot()

# total_length by assembly quality
ranked %>%
  filter(assembly_rank == 1) %>%
  ggplot(aes(
    x = total_length,
    fill = assembly_quality
  )) +
  geom_histogram(bins = 50)

ranked %>%
  select(sample_id, strain_no, batch, assembly_quality, suspicious_length, fragmented, low_continuity) %>%
  print(n = nrow(.))

# Test set for marker extraction:
#148 LIMS317701_45-s-73    s73       LIMS317701  possible_marker_recovery TRUE              TRUE       FALSE   
#13 LIMS317701_22-FL-172  172       LIMS317701  excellent                FALSE             FALSE      FALSE         
#14 LIMS32034_x_35-FL-172 172       LIMS32034_x likely_usable            TRUE              TRUE       FALSE          
#146 LIMS32034_69-s-70     s70       LIMS32034   excellent                FALSE             FALSE      FALSE         
#98 LIMS317701_61-s-193   s193      LIMS317701  excellent                FALSE             FALSE      FALSE         
#99 LIMS32034_x_75-s-193  s193      LIMS32034_x possible_marker_recovery TRUE              TRUE       FALSE         
#6 LIMS32034_x_26-FL-126 126       LIMS32034_x excellent                FALSE             FALSE      FALSE         
#7 LIMS317701_19-FL-126  126       LIMS317701  possible_marker_recovery TRUE              TRUE       FALSE         
#8 LIMS32034_x_31-FL-153 153       LIMS32034_x likely_usable            FALSE             FALSE      FALSE         
#9 LIMS32034_31-FL-153   153       LIMS32034   poor                     FALSE             FALSE      FALSE   





