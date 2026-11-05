###############################################################################
# Zebrafish skin microbiota: genus-level heatmap
#
# Purpose
#   Generate a heatmap showing the mean relative abundance of selected bacterial
#   genera across developmental/sample groups.
#
# Input requirements
#   Supply a phyloseq object that:
#     1. Contains the microbial community to be analysed (e.g. bacterial taxa).
#     2. Includes a taxonomic table with a "Genus" rank.
#     3. Includes sample metadata with a grouping variable (default: "Source").
#     4. Uses the following group labels by default:
#        "water", "larvae", "juvenile", and "adult".
#
#   The input object does not need to have a particular name. It is passed to
#   create_genus_heatmap() through the `physeq` argument.
#
# Example
#   heatmap_results <- create_genus_heatmap(
#     physeq = your_phyloseq_object,
#     group_var = "Source",
#     output_dir = "results"
#   )
#
#   print(heatmap_results$plot)
#
# Notes
#   - ASV abundances are converted to within-sample relative abundances before
#     being aggregated at genus level.
#   - The heatmap displays the mean relative abundance across samples within
#     each group. Individual sample values are retained in the returned data.
#   - Missing genus-by-sample combinations are treated as zero abundance.
#   - The ecological categories are descriptive groupings for visualization;
#     they are not inferred statistically by this script.
###############################################################################


# 1. Load required packages ---------------------------------------------------

required_packages <- c(
  "phyloseq",
  "tidyverse",
  "ggplot2",
  "viridis",
  "ggh4x"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following required packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(ggplot2)
  library(viridis)
  library(ggh4x)
})


# 2. Define the heatmap function ----------------------------------------------

create_genus_heatmap <- function(
    physeq,
    group_var = "Source",
    group_order = c("water", "larvae", "juvenile", "adult"),
    output_dir = ".",
    file_prefix = "zebrafish_skin_microbiota_heatmap"
) {

  # 2.1. Validate input -------------------------------------------------------

  if (!inherits(physeq, "phyloseq")) {
    stop("`physeq` must be a valid phyloseq object.", call. = FALSE)
  }

  if (is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    stop("The phyloseq object must contain sample metadata.", call. = FALSE)
  }

  if (!"Genus" %in% rank_names(physeq)) {
    stop(
      'The taxonomic table must contain a rank named "Genus".',
      call. = FALSE
    )
  }

  metadata <- data.frame(
    Sample = sample_names(physeq),
    as(sample_data(physeq), "data.frame"),
    check.names = FALSE
  )

  if (!group_var %in% colnames(metadata)) {
    stop(
      paste0(
        'The sample metadata must contain a column named "',
        group_var,
        '".'
      ),
      call. = FALSE
    )
  }

  metadata[[group_var]] <- as.character(metadata[[group_var]])

  missing_groups <- setdiff(group_order, unique(metadata[[group_var]]))

  if (length(missing_groups) > 0) {
    stop(
      "The following requested groups were not found in the metadata: ",
      paste(missing_groups, collapse = ", "),
      call. = FALSE
    )
  }

  unexpected_groups <- setdiff(unique(metadata[[group_var]]), group_order)

  if (length(unexpected_groups) > 0) {
    stop(
      "The grouping variable contains values not included in `group_order`: ",
      paste(unexpected_groups, collapse = ", "),
      ". Update `group_order` or filter the phyloseq object.",
      call. = FALSE
    )
  }


  # 2.2. Select genera and ecological categories -----------------------------

  genus_groups <- tibble(
    Genus = c(
      # Early colonization
      "Pseudomonas",
      "Aeromonas",
      "Sphingobium",
      "Rhizobium",
      "Chryseobacterium",

      # Developmental reorganization
      "Sphingomonas",
      "Variovorax",
      "Shewanella",
      "Stenotrophomonas",
      "Rheinheimera",

      # Adult-associated community
      "Escherichia/Shigella",
      "Micrococcus",
      "Streptococcus",
      "Nocardioides",
      "Mycobacterium"
    ),
    Category = rep(
      c(
        "Early colonization",
        "Developmental reorganization",
        "Adult-associated community"
      ),
      each = 5
    )
  )

  genus_order <- genus_groups$Genus

  tax_matrix <- as(tax_table(physeq), "matrix")
  available_genera <- unique(as.character(tax_matrix[, "Genus"]))

  missing_genera <- setdiff(genus_order, available_genera)

  if (length(missing_genera) > 0) {
    stop(
      "The following selected genera were not found in the taxonomic table: ",
      paste(missing_genera, collapse = ", "),
      ". Check the taxonomic assignments or edit `genus_groups`.",
      call. = FALSE
    )
  }


  # 2.3. Convert abundances to relative abundances ----------------------------

  physeq_relative <- transform_sample_counts(
    physeq,
    function(x) {
      if (sum(x) == 0) {
        x
      } else {
        x / sum(x)
      }
    }
  )


  # 2.4. Aggregate ASVs at genus level ---------------------------------------

  physeq_genus <- tax_glom(
    physeq_relative,
    taxrank = "Genus",
    NArm = FALSE
  )


  # 2.5. Extract selected genera ----------------------------------------------

  melted_data <- psmelt(physeq_genus)

  genus_abundances <- melted_data %>%
    filter(as.character(Genus) %in% genus_order) %>%
    transmute(
      Genus = as.character(Genus),
      Sample = as.character(Sample),
      Relative_abundance = Abundance * 100
    )

  sample_groups <- metadata %>%
    transmute(
      Sample = as.character(Sample),
      Group = as.character(.data[[group_var]])
    )


  # 2.6. Add zero-abundance combinations -------------------------------------

  # Explicitly add every selected genus to every sample. This ensures that
  # genera not detected in a sample are represented by zero rather than being
  # omitted from the calculation of group means.

  sample_abundances <- expand_grid(
    Genus = genus_order,
    Sample = sample_names(physeq)
  ) %>%
    left_join(sample_groups, by = "Sample") %>%
    left_join(
      genus_abundances,
      by = c("Genus", "Sample")
    ) %>%
    mutate(
      Relative_abundance = replace_na(Relative_abundance, 0),
      Group = factor(Group, levels = group_order)
    ) %>%
    left_join(genus_groups, by = "Genus")


  # 2.7. Calculate mean abundance across samples ------------------------------

  heatmap_summary <- sample_abundances %>%
    group_by(Genus, Category, Group) %>%
    summarise(
      Mean_abundance = mean(Relative_abundance, na.rm = TRUE),
      SD = sd(Relative_abundance, na.rm = TRUE),
      N = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Log_mean_abundance = log10(Mean_abundance + 0.001),
      Category = factor(
        Category,
        levels = c(
          "Early colonization",
          "Developmental reorganization",
          "Adult-associated community"
        )
      ),
      Genus = factor(
        Genus,
        levels = rev(genus_order)
      )
    )


  # 3. Build the heatmap -------------------------------------------------------

  # Subtle grey shades distinguish the ecological categories while keeping
  # the abundance scale as the main visual focus.

  category_colors <- c(
    "Early colonization" = "#F2F2F2",
    "Developmental reorganization" = "#E9E9E9",
    "Adult-associated community" = "#E0E0E0"
  )

  group_labels <- c(
    water = "Water",
    larvae = "Larvae",
    juvenile = "Juvenile",
    adult = "Adult"
  )

  heatmap_plot <- ggplot(
    heatmap_summary,
    aes(
      x = Group,
      y = Genus,
      fill = Log_mean_abundance
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.6
    ) +
    ggh4x::facet_grid2(
      Category ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y",
      strip = ggh4x::strip_themed(
        background_y = ggh4x::elem_list_rect(
          fill = unname(category_colors),
          colour = "black"
        ),
        text_y = ggh4x::elem_list_text(
          face = "bold",
          colour = "#263238",
          size = 10
        )
      )
    ) +
    scale_x_discrete(
      labels = group_labels,
      drop = FALSE
    ) +
    scale_fill_viridis_c(
      option = "C",
      direction = 1,
      name = expression(
        log[10] ~ "(Mean relative abundance, %)"
      )
    ) +
    labs(
      title = "Skin microbiota assembly during zebrafish development",
      subtitle = "Mean relative abundance of selected bacterial genera",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 16
      ),
      plot.subtitle = element_text(
        size = 11,
        color = "grey40"
      ),
      axis.text.x = element_text(
        size = 11,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(
        size = 10,
        face = "italic",
        color = "black"
      ),
      strip.text.y.left = element_text(
        angle = 0,
        hjust = 0.5
      ),
      strip.placement = "outside",
      panel.grid = element_blank(),
      panel.spacing.y = grid::unit(0.15, "cm"),
      legend.title = element_text(
        size = 10,
        face = "bold"
      ),
      legend.text = element_text(size = 9),
      plot.margin = margin(10, 15, 10, 10)
    )


  # 4. Export results ----------------------------------------------------------

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  ggsave(
    filename = file.path(
      output_dir,
      paste0(file_prefix, ".pdf")
    ),
    plot = heatmap_plot,
    width = 10,
    height = 7,
    units = "in"
  )

  ggsave(
    filename = file.path(
      output_dir,
      paste0(file_prefix, ".png")
    ),
    plot = heatmap_plot,
    width = 10,
    height = 7,
    units = "in",
    dpi = 600,
    bg = "white"
  )


  # 5. Return plot and data ----------------------------------------------------

  invisible(
    list(
      plot = heatmap_plot,
      summary = heatmap_summary,
      sample_abundances = sample_abundances
    )
  )
}


###############################################################################
# Run the analysis
#
# Source this script, then pass the phyloseq object to the function:

 heatmap_results <- create_genus_heatmap(
   physeq = your_phyloseq_object,
   group_var = "Source",
   group_order = c("water", "larvae", "juvenile", "adult"),
   output_dir = "results"
 )

# Display the heatmap:
 print(heatmap_results$plot)

# Access the mean-abundance table:
 heatmap_results$summary

# Access the sample-level abundance table:
 heatmap_results$sample_abundances
###############################################################################
