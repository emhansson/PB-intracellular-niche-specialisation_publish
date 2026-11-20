library(tidyverse)
library(ape)
library(pheatmap)
library(ggtree)
library(ggnewscale)

# Tree4 -- with outgroups for rooting and 186b reference but untrimmed alignments
df <- read_tsv("Data/Input/master_assembly_qc_table.tsv")
tree4 <- read.tree("Data/Input/plastid_tree4.treefile")

# Metadata for figure

metadata4 <- data.frame(strain = tree4$tip.label) %>%
  mutate(name = case_when(
    # Free-living strains
    str_detect(strain, "-FL-") ~
      str_extract(strain, "(?<=-FL-)\\d+"),
    # Symbiotic strains
    str_detect(strain, "-s-") ~
      paste0(
        "s",
        str_extract(strain, "(?<=-s-)\\d+")
      ),
    # Reference strains
    !str_detect(strain, "^LIMS") ~ {
      
      # split at underscores
      parts <- str_split(strain, "_", simplify = TRUE)
      
      # accession = final two elements joined
      accession <- paste0(parts[, ncol(parts)-1],
                          parts[, ncol(parts)])
      
      # taxon name = everything except accession
      taxon <- apply(
        parts[, 1:(ncol(parts)-2), drop = FALSE],
        1,
        paste,
        collapse = " "
      )
      
      paste0(taxon, " [", accession, "]")
    },
    
    TRUE ~ strain),
    lifestyle = case_when(str_detect(strain, "^LIMS") == TRUE & str_detect(strain, "FL") == TRUE ~ "Free-living",
                          str_detect(strain, "^LIMS") == TRUE & str_detect(strain, "s") == TRUE ~ "Symbiotic",
                          TRUE ~ "Reference")) %>%
  left_join(df %>% dplyr::select(strain_no, pond) %>% rename("name" = strain_no), by = "name")

metadata4

# Fixing 186b name and RX10 typo
metadata4$name[
  metadata4$name == "m conductrix [186bchloroplast]"
] <- "Micractinium conductrix [symb_186b]"

metadata4$pond[
  metadata4$pond == "RX10"
] <- "RIX10"


# Ordering the ponds so it makes sense in plots later on
pond_order <- c("Reference",
                "CHOR8",
                "CHOR9",
                "COX4",
                "COX7",
                "COX8",
                "DAR4",
                "DAR9",
                "DUN2",
                "DUN4",
                "N6",
                "N12",
                "QB1",
                "RIX4",
                "RIX9",
                "RIX10",
                "SM3",
                "SM6",
                "SM8",
                "SM12",
                "SM10",
                "SP2") 

metadata4$pond <- factor(
  metadata4$pond,
  levels = pond_order
)

metadata4$lifestyle <- factor(
  metadata4$lifestyle,
  levels = c("Free-living", "Symbiotic", "Reference")
)


# Rooting to outgroups
tree_rooted <- root(
  tree4,
  outgroup = c(
    "Chlamydomonas_reinhardtii_NC_005353.1",
    "Dunaliella_salina_NC_016732.1",
    "Choricystis_sp_OY253666.1"
  ),
  resolve.root = TRUE
)


# Rooted circular tree, with ponds
# Making references not have pond as NA
metadata4$pond[is.na(metadata4$pond)] <- "Reference"

# Creating just the tree
p <- ggtree(tree_rooted,
            layout = "circular") %<+% metadata4
p

# Creating pond annotation
ann <- metadata4 %>%
  distinct() %>%
  select(strain, pond)

rownames(ann) <- ann$strain
ann$strain <- NULL

# Making the pond annotation ring
p2 <- gheatmap(
  p,
  ann,
  offset = 0.1,
  width = 0.03,
  colnames = FALSE
)

p2

# Creating colour vector for the ponds
pond_cols <- setNames(
  rainbow(length(pond_order), end = 0.8, s = 0.6, v = 0.9),
  pond_order
)

pond_cols

pond_cols["Reference"] <- "grey95"

# Adding everything else to the figure
p3 <- p2 + geom_treescale(x = -0.1,
                          y = 0,
                          width = 0.05,
                          offset = 10) +
  geom_tippoint(aes(colour = lifestyle),
                size = 1) +
  geom_tiplab(aes(label = name,
                  colour = lifestyle),
              size = 3,
              offset = 0.005) +
  scale_color_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue",
    "Reference" = "black"
  )) +
  scale_fill_manual(values = pond_cols, breaks = pond_order) +
  labs(fill = "Pond", colour = "Lifestyle") +
  guides(
    colour = guide_legend(
      override.aes = list(
        label = "",
        shape = 16,
        size = 5
      )
    )
  ) +
  theme(legend.position = "left") +
  geom_point2(
    aes(
      subset =
        !isTip &
        branch.length > 0.005 &
        as.numeric(label) < 95,
    ),
    size = 1.5, shape = 21, fill = "white"
  )

p3


# Adding phylogroups
metadata4 <- metadata4 %>%
  left_join(read_csv("Data/Input/phylogroups.csv"), by = "strain")

# Creating just the tree
p <- ggtree(tree_rooted,
            layout = "circular") %<+% metadata4
p

# Creating pond annotation
ann_pond <- metadata4 %>%
  distinct() %>%
  select(strain, pond)

rownames(ann_pond) <- ann_pond$strain
ann_pond$strain <- NULL

# Making the pond annotation ring
# Creating colour vector for the ponds
pond_cols <- setNames(
  rainbow(length(pond_order), end = 0.8, s = 0.6, v = 0.9),
  pond_order
)

pond_cols

pond_cols["Reference"] <- "grey95"

p2 <- gheatmap(
  p,
  ann_pond,
  offset = 0.1,
  width = 0.03,
  colnames = FALSE
) +
  scale_fill_manual(values = pond_cols, breaks = pond_order, name = "Pond")

p2

# Making the phylogroup annotation ring

phylogroup_order <- c("Reference", "Phylogroup 1", "Phylogroup 2", "Phylogroup 3", "Phylogroup 4", "Phylogroup 5", 
                      "Phylogroup 6", "Phylogroup 7", "Phylogroup 8", "Phylogroup 9", "Phylogroup 10", "Phylogroup 11", "Phylogroup 12", 
                      "Phylogroup 13")

unique(metadata4$phylogroups)

metadata4$phylogroups <- factor(
  metadata4$phylogroups,
  levels = phylogroup_order
)


ann_phylogroup <- metadata4 %>%
  distinct() %>%
  select(strain, phylogroups)

rownames(ann_phylogroup) <- ann_phylogroup$strain
ann_phylogroup$strain <- NULL

p3 <- p2 +
  ggnewscale::new_scale_fill()

phylogroup_cols <- setNames(
  c("white", "grey50", rep("white", 12)),
  phylogroup_order
)

p3 <- gheatmap(
  p3,
  ann_phylogroup,
  offset = 0.135,
  width = 0.03,
  colnames = FALSE
) +
  scale_fill_manual(values = phylogroup_cols,
                    name = "Phylogroup", guide = "none")

p3


# Adding everything else to the figure
p4 <- p3 + geom_treescale(x = -0.1,
                          y = 0,
                          width = 0.05,
                          offset = 10) +
  geom_tippoint(aes(colour = lifestyle),
                size = 1) +
  geom_tiplab(aes(label = name,
                  colour = lifestyle),
              size = 3,
              offset = 0.005) +
  scale_color_manual(values = c(
    "Free-living" = "green3",
    "Symbiotic" = "cornflowerblue",
    "Reference" = "black"
  )) +
  labs(fill = "Pond", colour = "Lifestyle") +
  guides(
    colour = guide_legend(
      override.aes = list(
        label = "",
        shape = 16,
        size = 5
      )
    )
  ) +
  theme(legend.position = "left") +
  geom_point2(
    aes(
      subset =
        !isTip &
        branch.length > 0.005 &
        as.numeric(label) < 95,
    ),
    size = 1.5, shape = 21, fill = "white"
  )

p4

ggsave("Plots/tree4_rooted_ponds_phylogroups.png", p4, units = "in", width = 12, height = 13)
ggsave("Plots/tree4_rooted_ponds_phylogroups.tiff", p4, units = "in", width = 12, height = 13)


## Estimating marker evolutionary rates from marker trees

files <- list.files(
  "Data/Input/marker_trees",
  pattern = "\\.treefile$",
  full.names = TRUE
)

tree_lengths <- sapply(files, function(f) {
  sum(read.tree(f)$edge.length)
})

marker_rates <- data.frame(
  marker = sub("\\.treefile$", "", basename(files)),
  tree_length = tree_lengths
)

marker_rates %>%
  ggplot(aes(x = marker, y = tree_length)) +
  geom_point()
