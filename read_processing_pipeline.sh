#!/usr/bin/env bash

# ==============================================================================
# Skin microbiota across zebrafish development
# 16S rRNA amplicon processing pipeline
# ==============================================================================

# This script consolidates the Bash commands used to process the raw reads,
# merge paired-end reads, assign taxonomy, merge per-sample OTU tables, and
# extract the corresponding taxonomy table.
#
# Software reported in the manuscript:
#   fastp v0.22.0
#   PEAR v0.9.6
#   USEARCH v11.0.667_win32
#
# Expected files/programs in the working directory:
#   *.fastq / *.fastq.gz
#   usearch11.0.667_win32.exe
#   trainset19_072023.rdp.fasta
#   OTU_taxonomies.txt
#
# NOTE:
# The commands below reproduce the commands supplied for this project.
# No parameters have been added beyond those present in the original .sh files.

# ------------------------------------------------------------------------------
# 1. Quality control of raw reads (FastQC)
# ------------------------------------------------------------------------------
# This step was present in the original command history but is not described
# in the Methods paragraph supplied with the manuscript.

mkdir -p 16S_fastqc

fastqc --extract -t 5 -o 16S_fastqc *.fastq.gz


# ------------------------------------------------------------------------------
# 2. Quality and length filtering (fastp v0.22.0)
# ------------------------------------------------------------------------------
# Default fastp options were used.

fastp -i H2O_1_1.fastq -I H2O_1_2.fastq \
      -o out.H2O_1_1.fastq -O out.H2O_1_2.fastq

fastp -i H2O_2_1.fastq -I H2O_2_2.fastq \
      -o out.H2O_2_1.fastq -O out.H2O_2_2.fastq

fastp -i H2O_3_1.fastq -I H2O_3_2.fastq \
      -o out.H2O_3_1.fastq -O out.H2O_3_2.fastq

fastp -i 3dpf_1_1.fastq -I 3dpf_1_2.fastq \
      -o out.3dpf_1_1.fastq -O out.3dpf_1_2.fastq

fastp -i 3dpf_2_1.fastq -I 3dpf_2_2.fastq \
      -o out.3dpf_2_1.fastq -O out.3dpf_2_2.fastq

fastp -i 3dpf_3_1.fastq -I 3dpf_3_2.fastq \
      -o out.3dpf_3_1.fastq -O out.3dpf_3_2.fastq

fastp -i 10dpf_1_1.fastq -I 10dpf_1_2.fastq \
      -o out.10dpf_1_1.fastq -O out.10dpf_1_2.fastq

fastp -i 10dpf_2_1.fastq -I 10dpf_2_2.fastq \
      -o out.10dpf_2_1.fastq -O out.10dpf_2_2.fastq

fastp -i 10dpf_3_1.fastq -I 10dpf_3_2.fastq \
      -o out.10dpf_3_1.fastq -O out.10dpf_3_2.fastq

fastp -i A_1_1.fastq -I A_1_2.fastq \
      -o out.A_1_1.fastq -O out.A_1_2.fastq

fastp -i A_2_1.fastq -I A_2_2.fastq \
      -o out.A_2_1.fastq -O out.A_2_2.fastq

fastp -i A_3_1.fastq -I A_3_2.fastq \
      -o out.A_3_1.fastq -O out.A_3_2.fastq


# ------------------------------------------------------------------------------
# 3. Paired-end read merging (PEAR v0.9.6)
# ------------------------------------------------------------------------------

pear -j5 \
     -f out.H2O_1_1.fastq \
     -r out.H2O_1_2.fastq \
     -o H2O_1.fastq

pear -j5 \
     -f out.H2O_2_1.fastq \
     -r out.H2O_2_2.fastq \
     -o H2O_2.fastq

pear -j5 \
     -f out.H2O_3_1.fastq \
     -r out.H2O_3_2.fastq \
     -o H2O_3.fastq

pear -j5 \
     -f out.3dpf_1_1.fastq \
     -r out.3dpf_1_2.fastq \
     -o 3dpf_1.fastq

pear -j5 \
     -f out.3dpf_2_1.fastq \
     -r out.3dpf_2_2.fastq \
     -o 3dpf_2.fastq

pear -j5 \
     -f out.3dpf_3_1.fastq \
     -r out.3dpf_3_2.fastq \
     -o 3dpf_3.fastq

pear -j5 \
     -f out.10dpf_1_1.fastq \
     -r out.10dpf_1_2.fastq \
     -o 10dpf_1.fastq

pear -j5 \
     -f out.10dpf_2_1.fastq \
     -r out.10dpf_2_2.fastq \
     -o 10dpf_2.fastq

pear -j5 \
     -f out.10dpf_3_1.fastq \
     -r out.10dpf_3_2.fastq \
     -o 10dpf_3.fastq

pear -j5 \
     -f out.A_1_1.fastq \
     -r out.A_1_2.fastq \
     -o A_1.fastq

pear -j5 \
     -f out.A_2_1.fastq \
     -r out.A_2_2.fastq \
     -o A_2.fastq

pear -j5 \
     -f out.A_3_1.fastq \
     -r out.A_3_2.fastq \
     -o A_3.fastq


# ------------------------------------------------------------------------------
# 4. Taxonomic assignment with USEARCH v11.0.667
# ------------------------------------------------------------------------------
# Database:
#   trainset19_072023.rdp.fasta
#
# Identity threshold:
#   97%
#
# The original commands use USEARCH -usearch_global with -otutabout.

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/H2O_1.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout H2O_1_otu_table.txt

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/H2O_2.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout H2O_2_otu_table.txt

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/H2O_3.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout H2O_3_otu_table.txt

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/3dpf_1.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout 3dpf_1_otu_table.txt

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/3dpf_2.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout 3dpf_2_otu_table.txt

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/3dpf_3.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout 3dpf_3_otu_table.txt

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/10dpf_1.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout 10dpf_1_otu_table.txt

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/10dpf_2.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout 10dpf_2_otu_table.txt

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/10dpf_4.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout 10dpf_3_otu_table.txt

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/A_1.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout A_1_otu_table.txt

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/A_3.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout A_2_otu_table.txt

./usearch11.0.667_win32.exe \
    -usearch_global 16S_merged_reads/A_4.fastq.assembled.fastq \
    -db trainset19_072023.rdp.fasta \
    -strand plus \
    -id 0.97 \
    -otutabout A_3_otu_table.txt


# ------------------------------------------------------------------------------
# 5. Merge per-sample OTU tables
# ------------------------------------------------------------------------------

./usearch11.0.667_win32.exe \
    -otutab_merge \
    H2O_1_otu_table.txt,H2O_2_otu_table.txt,H2O_3_otu_table.txt,3dpf_1_otu_table.txt,3dpf_2_otu_table.txt,3dpf_3_otu_table.txt,10dpf_1_otu_table.txt,10dpf_2_otu_table.txt,10dpf_3_otu_table.txt,A_1_otu_table.txt,A_2_otu_table.txt,A_3_otu_table.txt \
    -output otu_table_merged.txt


# ------------------------------------------------------------------------------
# 6. Extract taxonomy corresponding to the OTUs in the merged table
# ------------------------------------------------------------------------------

# Extract OTU IDs from the merged abundance table.
# The original command retains the last 1397 IDs.

less otu_table_merged.txt | \
    cut -f1 | \
    tail -n 1397 > OTU_IDs.txt

# Retrieve taxonomy entries corresponding to those OTU IDs.

grep -f OTU_IDs.txt \
    OTU_taxonomies.txt > OTU_assigned_taxonomies.txt


# ==============================================================================
# End of pipeline
# ==============================================================================
