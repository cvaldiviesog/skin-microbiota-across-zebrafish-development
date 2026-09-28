# ============================================================================
# Skin microbiota across zebrafish development
# FIGURES
# ============================================================================
# Run pipeline.R first. This script loads the resulting workspace and creates
# all figures. Generated PDF files should remain local.
# ============================================================================

library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(metacoder)

load("pipeline_workspace.RData")

# Taxonomic composition ---------------------------------------------------
pdf("OTUs_per_Phylum.pdf")
plot_bar(bacteria, x = "NewID", fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum),
           stat = "identity", position = "stack")
dev.off()

pdf("OTUs_per_Class.pdf", width = 12, height = 10)
plot_bar(bacteria, x = "NewID", fill = "Class") +
  geom_bar(aes(color = Class, fill = Class),
           stat = "identity", position = "stack")
dev.off()

pdf("OTUs_per_Order.pdf", width = 12, height = 10)
plot_bar(bacteria, x = "NewID", fill = "Order") +
  geom_bar(aes(color = Order, fill = Order),
           stat = "identity", position = "stack")
dev.off()

pdf("OTUs_per_Family.pdf", width = 22, height = 16)
plot_bar(bacteria, x = "NewID", fill = "Family") +
  geom_bar(aes(color = Family, fill = Family),
           stat = "identity", position = "stack")
dev.off()

pdf("OTUs_per_Genus.pdf", width = 42, height = 22)
plot_bar(bacteria, x = "NewID", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus),
           stat = "identity", position = "stack")
dev.off()

# Heatmaps ----------------------------------------------------------------
pdf("OTU_heatmap.pdf", width = 42, height = 22)
plot_heatmap(bacteria,
             sample.order = sampleOrder,
             method = "NMDS",
             distance = "bray",
             low = "#66CCFF",
             high = "#000033",
             na.value = "white")
dev.off()

pdf("OTU_heatmap_v2.pdf", width = 35, height = 22)
heatmap(otu_mat, Colv = NA)
legend(x = "right",
       legend = c(1, 16, 256, 4096),
       fill = colorRampPalette(brewer.pal(8, "Oranges"))(4),
       title = "Abundance",
       cex = 4)
dev.off()

# Alpha diversity ---------------------------------------------------------
pdf("Alpha_diversity.pdf", width = 14)
plot_richness(bacteria,
              sortby = sampleOrder,
              measures = c("Shannon", "Observed"))
dev.off()

pdf("Alpha diversity Inverse Simpson Index.pdf")
ggplot(metacoder_samples_data,
       aes(x = Source, y = inv_simp)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data),
                y = max(metacoder_samples_data$inv_simp) + 1,
                label = group_data$groups),
            col = "black", size = 10) +
  geom_boxplot() +
  ggtitle("Alpha diversity of samples") +
  xlab("Samples") +
  ylab("Inverse Simpson Index")
dev.off()

group_data_shannon <- tukey_shannon$groups[
  order(rownames(tukey_shannon$groups)),
]

pdf("Alpha diversity Shannon Index.pdf")
ggplot(metacoder_samples_data,
       aes(x = Source, y = Shannon)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_shannon),
                y = max(metacoder_samples_data$Shannon) + 1,
                label = group_data_shannon$groups),
            col = "black", size = 10) +
  geom_boxplot() +
  ggtitle("Alpha diversity of samples") +
  xlab("Samples") +
  ylab("Shannon Index")
dev.off()

# NMDS --------------------------------------------------------------------
pdf("OTU_NMDS.pdf", width = 14)
plot_ordination(bacteria,
                carbom.ord,
                type = "taxa",
                color = "Phylum",
                title = "OTUs")
dev.off()

pdf("OTU_NMDS_per_Phylum.pdf", width = 14)
plot_ordination(bacteria,
                carbom.ord,
                type = "taxa",
                color = "Class",
                title = "OTUs") +
  facet_wrap(~Phylum, 3)
dev.off()

pdf("Samples_NMDS.pdf", width = 8)
plot_ordination(bacteria,
                carbom.ord,
                type = "samples",
                color = "Source",
                title = "Samples") +
  geom_text(mapping = aes(label = sample_names(bacteria)), size = 2)
dev.off()

# Taxonomic counts per condition -----------------------------------------
pdf("Phylum count per condition.pdf")
ggplot(samples_data,
       aes(x = Source, y = Phylum_count)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_2),
                y = max(samples_data$Phylum_count) + 1,
                label = group_data_2$groups),
            col = "black", size = 10) +
  geom_boxplot() +
  ggtitle("Phylum counts per condition") +
  xlab("Samples") +
  ylab("Phylum count")
dev.off()

pdf("Order count per condition.pdf")
ggplot(samples_data,
       aes(x = Source, y = Order_count)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_3),
                y = max(samples_data$Order_count) + 1,
                label = group_data_3$groups),
            col = "black", size = 10) +
  geom_boxplot() +
  ggtitle("Order counts per condition") +
  xlab("Samples") +
  ylab("Order count")
dev.off()

# Metacoder heat tree -----------------------------------------------------
pdf("Heat_tree_all_samples.pdf")
heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford")
dev.off()

message("All figures generated successfully.")
