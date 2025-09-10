# Herbicide Soil Microbiome Analysis Pipeline
Pipeline for analyzing 16S rRNA sequencing data from soil samples to evaluate diversity metrics, community composition, and functional and taxonomic profiles.
information on the scripts used, what order to run them in, how to setup the file directory, and software needed

## Prerequisites
- Windows Subsystem for Linux 2 (WSL2) with Ubuntu
- Conda installed in WSL2 with PICRUsT2 environment (`conda create -n picrust2_env -c conda-forge -c bioconda picrust2 -y`)
- R (version 4.0+) and RStudio
- RDP taxonomy database: rdp_19_toSpecies_trainset.fa.gz

## Required Directory Structure
```
project_folder/
├── data/
│   ├── demultiplexed_16s/          # Raw fastq files
│   ├── MapFileHerb(in).csv         # Metadata file
│   ├── rdp_19_toSpecies_trainset.fa.gz
│   └── CO2_Resp_RAW_Data.xlsx      # For CO2 analysis
├── output/                         # Created by scripts
└── Figures/                        # Created by scripts
```

## Critical File Requirements
**Metadata format (MapFileHerb(in).csv):**
- Columns: `Herbicide`, `Concentration`
- Herbicide values: `Control`, `2-4 D`, `Glyphosate`
- Concentration values: `0`, `0.125`, `0.25`, `0.5`, `1`
- Sample IDs must match fastq file names

## Scripts (Run in Order)

**1. herbicide_beginning_soil_taxonomy.R**
This script runs the complete DADA2 pipeline for 16S rRNA gene sequencing analysis. It processes raw fastq files through quality filtering, error learning, denoising, merging, chimera removal, and taxonomy assignment using the RDP database. Generates ASV tables, diversity metrics, ordination plots (NMDS/MDS), relative abundance plots, and PICRUsT2 input files. Outputs comprehensive alpha diversity analysis with boxplots for observed richness, Shannon diversity, and Simpson diversity across treatments.

**2. bash_script_for_running_picrust2.txt**
This script executes the PICRUsT2 pipeline in WSL2 to predict functional gene content from 16S data. Takes ASV sequences and abundance tables from DADA2, runs functional predictions, and adds MetaCyc pathway descriptions. First-time users must uncomment the conda environment creation lines. Outputs functional pathway abundance tables and enzyme predictions.

**3. heatmaps_picstrust2_FINAL.R**
This script generates publication-quality heatmaps of carbon cycling pathways from PICRUsT2 output. Focuses on five key carbon pathways, calculates statistical significance between treatments and controls using t-tests, and creates separate heatmaps for 2,4-D and Glyphosate treatments. Heatmaps display fold-change relative to controls with significance indicators (* p<0.05, ** p<0.01, *** p<0.001).

**4. CO2_Respiration_figures.R**
This script analyzes soil respiration data from Excel files across multiple sampling timepoints. Processes CO2 measurements for 2,4-D and Glyphosate treatments, calculates cumulative CO2 production, and generates bar plots with error bars showing treatment effects on soil microbial respiration. Creates separate publication-ready figures for each herbicide treatment.

## File Path Customization Required
- Line 15 in herbicide_beginning_soil_taxonomy.R: `setwd("your_path_here")`
- Lines 45-47 in bash_script_for_running_picrust2.txt: Update file paths to your project location
- Lines 10, 170 in heatmaps_picstrust2_FINAL.R: Update working directory and PICRUsT2 results path
- Lines 1, 8 in CO2_Respiration_figures.R: Update working directory and data file path

## Expected Outputs
- ASV abundance tables (.csv)
- Taxonomy assignments (.csv)
- Alpha diversity metrics (.csv)
- Ordination plots (.png)
- Relative abundance plots (.png)
- Carbon pathway heatmaps (.png)
- CO2 respiration plots (.png)
- PICRUsT2 functional predictions (folder)

