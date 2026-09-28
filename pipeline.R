# ============================================================================
# Skin microbiota across zebrafish development
# PIPELINE
# ============================================================================
# This script performs data import, preprocessing and statistical analyses.
# It does not generate figures.
#
# Expected input files in the local working directory:
#   otu_table_merged.txt
#   OTU_assigned_taxonomies.txt
#   samples_table.txt
#   OTU_assigned_taxonomies_metacoder_2.txt
#
# The script saves the complete workspace to pipeline_workspace.RData.
# This file is local and should NOT be uploaded to GitHub.
# ============================================================================

# Packages -----------------------------------------------------------------
library(phyloseq)
library(remotes)
library(psadd)
library(tidyverse)
library(ggplot2)
library(metacoder)
library(RColorBrewer)
library(data.table)
library(qiime2R)
library(vegan)
library(agricolae)

# Function -----------------------------------------------------------------
phyloseq_ntaxa_by_tax <- function(physeq, TaxRank = "Phylum", relative = F,
                                  add_meta_data = T){
  mm <- phyloseq::psmelt(physeq)
  
  count_otus <- function(z, TaxRank = TaxRank, relative = relative){
    zero_otu <- z$Abundance > 0
    if(any(zero_otu)){ z <- z[zero_otu, ] }
    
    rez <- as.data.frame(table(z[, TaxRank]), stringsAsFactors = F)
    colnames(rez) <- c(TaxRank, "N.OTU")
    
    if(relative == TRUE){
      rez$N.OTU <- with(rez, N.OTU / sum(N.OTU))
    }
    
    return(rez)
  }
  
  res <- plyr::ddply(.data = mm, .variables = "Sample",
                     .fun = count_otus,
                     TaxRank = TaxRank,
                     relative = relative)
  
  if(add_meta_data == TRUE){
    if(!is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE))){
      metad <- data.frame(Sample = phyloseq::sample_names(physeq),
                          phyloseq::sample_data(physeq))
      main_cols <- c("Sample")
      meta_cols <- colnames(metad)[!colnames(metad) %in% main_cols]
      res <- cbind(res,
                   metad[match(x = res$Sample,
                               table = metad$Sample), meta_cols])
    }
  }
  
  rownames(res) <- NULL
  return(res)
}

# Import data --------------------------------------------------------------
otu_mat <- read.table("otu_table_merged.txt", header = T, sep = "\t")
otu_mat <- otu_mat %>% tibble::column_to_rownames("OTU.ID")
otu_mat <- as.matrix(otu_mat)

tax_mat <- read.table("OTU_assigned_taxonomies.txt", header = T, sep = "\t")
tax_mat <- tax_mat %>% tibble::column_to_rownames("OTU_ID")
tax_mat <- as.matrix(tax_mat)

samples_df <- read.table("samples_table.txt", header = T, sep = "\t")
samples_df$Sample <- as.character(samples_df$Sample)
samples_df <- samples_df %>% tibble::column_to_rownames("Sample")

# Create phyloseq object ---------------------------------------------------
OTU <- otu_table(otu_mat, taxa_are_rows = T)
TAX <- tax_table(tax_mat)
samples <- sample_data(samples_df)

carbom <- phyloseq(OTU, TAX, samples)

# Normalize by sequencing depth -------------------------------------------
total <- median(sample_sums(carbom))
standf <- function(x, t = total) round(t * (x / sum(x)))
carbom <- transform_sample_counts(carbom, standf)

# Keep sample order --------------------------------------------------------
sample_order <- c("H2O_1", "H2O_2", "H2O_3",
                  "L_3dpf_1", "L_3dpf_2", "L_3dpf_3",
                  "L_10dpf_1", "L_10dpf_2", "L_10dpf_3",
                  "A_1", "A_2", "A_3")

sample_data(carbom)$NewID <- factor(sample_names(carbom),
                                     levels = sample_order)

# Bacterial community -----------------------------------------------------
bacteria <- subset_taxa(carbom, Domain == "Bacteria")
sampleOrder <- unique(sample_names(bacteria))

# Taxonomic counts --------------------------------------------------------
PhylumCounts <- phyloseq_ntaxa_by_tax(bacteria, TaxRank = "Phylum")
ClassCounts <- phyloseq_ntaxa_by_tax(bacteria, TaxRank = "Class")
OrderCounts <- phyloseq_ntaxa_by_tax(bacteria, TaxRank = "Order")
FamilyCounts <- phyloseq_ntaxa_by_tax(bacteria, TaxRank = "Family")
GenusCounts <- phyloseq_ntaxa_by_tax(bacteria, TaxRank = "Genus")

# Alpha diversity ---------------------------------------------------------
alpha_phyloseq <- estimate_richness(bacteria,
                                     measures = c("Shannon", "Observed"))
alpha_phyloseq$Sample <- rownames(alpha_phyloseq)

# NMDS and Bray-Curtis distance -------------------------------------------
carbom.ord <- ordinate(bacteria, "NMDS", "bray")
distance_matrix <- distance(bacteria, method = "bray")

# PERMANOVA ----------------------------------------------------------------
metadata <- as(sample_data(bacteria), "data.frame")

dist_result <- as_tibble(adonis2(distance_matrix ~ Source,
                                  data = metadata))

# Homogeneity of dispersion ----------------------------------------------
despersion_result <- betadisper(distance_matrix, metadata$Source)

# Taxonomic counts per sample --------------------------------------------
phyl_count <- PhylumCounts %>%
  filter(N.OTU > 0) %>%
  count(Sample, name = "Phylum_count")

Ord_count <- OrderCounts %>%
  filter(N.OTU > 0) %>%
  count(Sample, name = "Order_count")

samples_data <- metadata %>%
  rownames_to_column("Sample") %>%
  left_join(phyl_count, by = "Sample") %>%
  left_join(Ord_count, by = "Sample")

samples_data$Source <- factor(samples_data$Source,
                              levels = c("water", "larvae", "juvenile", "adult"))

# ANOVA: taxonomic counts -------------------------------------------------
anova_result_2 <- aov(Phylum_count ~ Source, samples_data)
anova_result_3 <- aov(Order_count ~ Source, samples_data)

summary(anova_result_2)
summary(anova_result_3)

# Tukey tests -------------------------------------------------------------
tukey_result_2 <- HSD.test(anova_result_2, "Source", group = T)
print(tukey_result_2)

tukey_result_3 <- HSD.test(anova_result_3, "Source", group = T)
print(tukey_result_3)

group_data_2 <- tukey_result_2$groups[
  order(rownames(tukey_result_2$groups)),
]

group_data_3 <- tukey_result_3$groups[
  order(rownames(tukey_result_3$groups)),
]

# Metacoder ---------------------------------------------------------------
lineage <- read.table("OTU_assigned_taxonomies_metacoder_2.txt",
                      sep = "\t", header = T)

uniqTax <- lineage %>% distinct(Feature.ID, .keep_all = T)
rownames(uniqTax) <- uniqTax$Feature.ID

tax_metacoder <- parse_taxonomy(uniqTax, tax_sep = ";")

OTU_data <- read.table("OTU_assigned_taxonomies_metacoder_2.txt",
                       header = T, sep = "\t")

metacoder_samples_data <- read.table("samples_table.txt",
                                     header = T, sep = "\t")

obj <- parse_tax_data(OTU_data,
                      class_cols = "Feature.ID",
                      class_sep = ",",
                      class_key = c(tax_rank = "info",
                                    tax_name = "taxon_name"),
                      class_regex = "^(.+)__(.+)$")

obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
                                       cols = metacoder_samples_data$Sample)

obj$data$diff_table <- compare_groups(obj,
                                      data = "tax_abund",
                                      cols = metacoder_samples_data$Sample,
                                      groups = metacoder_samples_data$Source)

print(obj$data$diff_table)

# Alpha diversity from metacoder taxonomic table --------------------------
metacoder_samples_data$Source <- factor(
  metacoder_samples_data$Source,
  levels = c("water", "larvae", "juvenile", "adult")
)

metacoder_samples_data$inv_simp <- diversity(
  obj$data$tax_data[, metacoder_samples_data$Sample],
  index = "invsimpson",
  MARGIN = 2
)

metacoder_samples_data$Shannon <- diversity(
  obj$data$tax_data[, metacoder_samples_data$Sample],
  index = "shannon",
  MARGIN = 2
)

# ANOVA: alpha diversity --------------------------------------------------
anova_result <- aov(inv_simp ~ Source, metacoder_samples_data)
anova_shannon <- aov(Shannon ~ Source, metacoder_samples_data)

summary(anova_result)
summary(anova_shannon)

# Tukey tests -------------------------------------------------------------
tukey_result <- HSD.test(anova_result, "Source", group = T)
print(tukey_result)

tukey_shannon <- HSD.test(anova_shannon, "Source", group = T)
print(tukey_shannon)

group_data <- tukey_result$groups[
  order(rownames(tukey_result$groups)),
]

group_data_shannon <- tukey_shannon$groups[
  order(rownames(tukey_shannon$groups)),
]

# Save workspace for figures ---------------------------------------------
save.image("pipeline_workspace.RData")

message("Pipeline completed successfully.")
message("Workspace saved as pipeline_workspace.RData")
