# Skin microbiota across zebrafish development

This repository contains the computational scripts used to process 16S rRNA amplicon sequencing data and perform the downstream analyses reported in:

> **Valdivieso, C., Aburto, J., Kallens, V., Rojas, D., Varas, M., Chávez, F. P., & Allende, M. L. (2026). A tradeoff for zebrafish skin microbiota colonization: changes in bacterial diversity throughout development. _Fisheries Science_.**

**DOI:** https://doi.org/10.1007/s12562-026-01977-8

## Repository contents

The repository is intentionally limited to analysis scripts. Raw sequencing data, intermediate files, and generated figures are not included.

```text
skin-microbiota-across-zebrafish-development/
│
├── README.md
├── read_processing_pipeline.sh
├── pipeline.R
└── figures.R
```

### `read_processing_pipeline.sh`

Bash pipeline for processing the paired-end 16S rRNA amplicon reads.

The workflow includes:

1. Initial read quality assessment with **FastQC**.
2. Quality and length filtering with **fastp v0.22.0**.
3. Paired-end read merging with **PEAR v0.9.6**.
4. Taxonomic assignment with **USEARCH v11.0.667_win32** against the RDP 16S training set.
5. Consolidation of the per-sample OTU tables using `usearch -otutab_merge`.
6. Extraction of the taxonomy corresponding to the OTUs retained in the merged abundance table.

### `pipeline.R`

R script containing the downstream microbiota analyses.

The pipeline includes:

- construction and preprocessing of a `phyloseq` object;
- sequencing-depth normalization;
- filtering of bacterial taxa;
- taxonomic composition summaries;
- alpha-diversity analyses;
- Bray–Curtis dissimilarity and NMDS;
- PERMANOVA and multivariate dispersion;
- statistical comparisons of taxonomic richness;
- taxonomic analyses using `metacoder`;
- ANOVA and Tukey post-hoc tests.

The script saves the objects generated during the analysis to a local `pipeline_workspace.RData` file. This intermediate file is not part of the repository.

### `figures.R`

R script used to generate the figures from the objects produced by `pipeline.R`.

The script includes figures representing:

- taxonomic composition across samples;
- OTU heatmaps;
- alpha diversity;
- NMDS ordinations;
- taxonomic richness across developmental stages;
- metacoder heat trees.

Generated PDF figures are not included in the repository.

## Reproducibility

The intended workflow is:

```text
Raw paired-end reads
        │
        ▼
read_processing_pipeline.sh
        │
        ▼
Processed reads / OTU tables / taxonomy
        │
        ▼
pipeline.R
        │
        ▼
Statistical and ecological analyses
        │
        ▼
figures.R
        │
        ▼
Figures
```

The scripts should be executed in this order:

```bash
bash read_processing_pipeline.sh
```

followed by:

```r
source("pipeline.R")
source("figures.R")
```

The required input files and external databases are not stored in this repository. They should be placed locally in the working directory before running the corresponding scripts.

## Software

The sequencing-processing workflow used the following software versions:

- **fastp v0.22.0**
- **PEAR v0.9.6**
- **USEARCH v11.0.667_win32**
- **R** and the packages required by `pipeline.R` and `figures.R`

The original sequencing-processing commands are preserved in `read_processing_pipeline.sh`.

## Citation

If you use these scripts, please cite the associated publication:

Valdivieso, C., Aburto, J., Kallens, V., Rojas, D., Varas, M., Chávez, F. P., & Allende, M. L. (2026). A tradeoff for zebrafish skin microbiota colonization: changes in bacterial diversity throughout development. _Fisheries Science_.

https://doi.org/10.1007/s12562-026-01977-8
