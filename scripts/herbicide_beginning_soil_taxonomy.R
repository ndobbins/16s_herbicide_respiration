# ============================================================================
# DADA2 Pipeline for Herbicide Soil Taxonomy Analysis
# Author: Nate Dobbins
# Date: 2025-05-29
# Description: Complete DADA2 pipeline for 16S rRNA gene sequencing analysis
# ============================================================================

# Sources:
# https://benjjneb.github.io/dada2/tutorial_1_8.html
# DECIPHER SILVA species SILVA assignment: https://www2.decipher.codes/Downloads.html

# ============================================================================
# SET WORKING DIRECTORY
# ============================================================================

# Replace quotes in parenthesis with the path to your project folder
# Also replace any backslahes "\" with a forward slash "/"
# Folder SHOULD CONTAIN your demultiplexed read files and map file
# See README file for how your working directory folder shuold be setup
setwd("C:/Users/Nathan/Documents/Practicum/herbicide_soil_DOBBINS")

# ============================================================================
# INSTALL AND LOAD PACKAGES
# ============================================================================

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required CRAN packages if missing
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggsci", quietly = TRUE)) {
  install.packages("ggsci")
}
if (!requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}

# Install required Bioconductor packages if missing
if (!requireNamespace("dada2", quietly = TRUE)) {
  BiocManager::install("dada2", ask = FALSE)
}
if (!requireNamespace("DECIPHER", quietly = TRUE)) {
  BiocManager::install("DECIPHER", ask = FALSE)
}
if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq", ask = FALSE)
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings", ask = FALSE)
}

cran_packages <- c("tidyverse", "ggplot2")
bioc_packages <- c("dada2", "DECIPHER", "phyloseq", "Biostrings")

lapply(c(cran_packages, bioc_packages), library, character.only = TRUE)

# Install RColorBrewer
install.packages("RColorBrewer")

# Load all required libraries
library(dada2);      packageVersion("dada2") 
library(DECIPHER);   packageVersion("DECIPHER")
library(phyloseq);   packageVersion("phyloseq")   
library(Biostrings); packageVersion("Biostrings")  
library(tidyverse);  packageVersion("tidyverse") 
library(ggplot2);    packageVersion("ggplot2")
library(RColorBrewer)
library(cowplot)

# ============================================================================
# SET FILE PATHS
# ============================================================================

# Set as file path to the folder containing the demultiplexed fastQ files
path <- "data/demultiplexed_16s"

# List all files in the path to verify correct directory access
# list.files(path)  

# Get the full paths for all forward (R1) and reverse (R2)
# These are your raw sequencing reads  
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names from file names by joining the first 2 elements separated by underscores
# This helps keep track of which reads belong to which sample
# Extract sample names (single-lane data: Sample_S#_L001_R1_001.fastq -> Sample_S#)
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[1:2], collapse="_"))

# Verify sample naming
head(sample.names)

# ============================================================================
# INSPECT QUALITY OF THE READS
# ============================================================================

# NOTE: The task was to produce one R script for this whole DADA2 pipeline, but it needs to be split in two.
# The R script will be split into two 2 scripts. One that goes all the way until these read quality 
# profile plots are produced, and a second that you will execute after viewing the read quality profile plots.
# You must stop the running the code once the quality profile plots have been produced in order to 
# assess them and choose trimming parameters in the next step (Filter and Trim the Reads)

# Plot quality profiles for the first two Forward reads
# Helps visualize per-base quality scores to potentially choose trimming parameters
png("forward_quality_plot.png", width = 800, height = 600)
plotQualityProfile(fnFs[1:2])
dev.off()

# Plot quality profiles for the first two Reverse reads
# Helps visualize per-base quality scores to potentially choose trimming parameters
png("reverse_quality_plot.png", width = 800, height = 600)
plotQualityProfile(fnRs[1:2])
dev.off()

# ============================================================================
# FILTER AND TRIM THE READS
# ============================================================================

# Create new names for filtered read files, put in new subdirectory
# Create names for the filtered files and save them in a new subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter the reads
out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = c(150, 150),        # (set each as no greater than read length)
  # be mindful of amplicon length
  # both trunc lengths added together minus their overlap should be longer than length of amplicon
  maxN = 0,                      # Discards any reads with ambiguous bases (Ns)
  maxEE = c(2, 2),               # Max expected errors: lower is stricter (usually 2-3 for F, 5-7 R)                   
  # adjust as necessary if none/not many reads pass the filter
  # format : c(maxEE for forward, maxEE for reverse)  
  truncQ = 2,                    # Trim where quality score drops below 2 (increase if quality is poor)
  rm.phix = TRUE,                # (keep TRUE)
  compress = TRUE,               # Save output as compressed .gz files
  multithread = FALSE)           # Set to TRUE for Linux/Mac; FALSE for Windows 
head(out)

# ============================================================================
# LEARN THE ERROR RATES
# ============================================================================

# Learn the error rates from the filtered forward reads
# These are essential for the DADA2 algorithm to function
# Can take a while to run
errF <- learnErrors(filtFs, multithread=TRUE)

# Learn the error rates from the filtered reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)

# Visualize how well the estimated error rates (black line) fit the observed errors (points)
# This helps verify that error learning succeeded
plotErrors(errF, nominalQ=TRUE)

# Estimated error rates (black line) are a good fit for the observed rates (the points).

# ============================================================================
# DE-REPLICATION
# ============================================================================

# Combine all identical reads into unique sequences (dereplication)
derepFs <- derepFastq(filtFs, verbose=TRUE) 
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name each dereplicated object using the sample names
# This ensures sample tracking is preserved downstream
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# ============================================================================
# SAMPLE INFERENCE (DENOISING)
# ============================================================================

# Use DADA2 to determine true biological sequences from the forward reads
# This step corrects errors in sequencing data based on the learned error models that we calculated
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

# Repeat de-noising for reverse reads using their error model
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# View the output from the first sample to inspect inference results
dadaFs[[1]]

# ============================================================================
# MERGE PAIRED READS
# ============================================================================

# Merge forward and reverse reads to get full, de-noised sequences
# This step aligns overlapping regions to reconstruct the original amplicons
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# ============================================================================
# CONSTRUCTING SEQUENCE TABLE
# ============================================================================

# Build the sequence table (ASV table) from merged reads
# Rows = samples, Columns = unique sequences (ASVs)
seqtab <- makeSequenceTable(mergers)

# Check dimensions of the sequence table
# Useful to see number of samples (rows) and ASVs (columns)
dim(seqtab)

# Examine the distribution of sequence lengths
# Helps verify consistent amplicon length and detect any outliers
table(nchar(getSequences(seqtab)))

# ============================================================================
# REMOVE CHIMERAS
# ============================================================================

# Remove chimeric sequences
# Chimeras are artifacts from PCR that can be misidentified as true variants
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Calculate the proportion of reads retained after chimera removal
# Helps determine how much data was lost to chimeras
sum(seqtab.nochim)/sum(seqtab)

# Export seqtab.nochim as CSV
# Convert the non-chimeric sequence table to a data frame to be exported
seqtab_nochim_df <- as.data.frame(seqtab.nochim)

# Write the table to a CSV file (row.names=TRUE keeps the sample IDs)
write.csv(seqtab_nochim_df, file="seqtab_nochim_herbicide_soil_NATHAN.csv", row.names=TRUE)

# ============================================================================
# TRACK READS THROUGH PIPELINE
# ============================================================================

# Function for counting unique reads at each step
getN <- function(x) sum(getUniques(x))

# Combine read counts at key pipeline steps into one table:
# input → filtered → denoised → merged → non-chimeric
track <- cbind(
  out,                        # input & filtered reads
  sapply(dadaFs, getN),       # forward denoised reads
  sapply(dadaRs, getN),       # reverse denoised reads
  sapply(mergers, getN),      # merged reads
  rowSums(seqtab.nochim)      # non-chimeric reads
)
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

# Add column headers
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

# Assign sample names to rows
rownames(track) <- sample.names

# Preview the table
head(track)

# ============================================================================
# ASSIGN TAXONOMY (assignTaxonomy)
# ============================================================================

# Why use assignTaxonomy/RDP database?
# Use this with the RDP database to produce assignTaxonomy outputs compatible with downstream analysis such as PICRUsT2.

taxa <- assignTaxonomy(seqtab.nochim, "data/rdp_19_toSpecies_trainset.fa.gz", multithread=TRUE, minBoot =80)

# View the taxonomic assignments from assignTaxonomy
taxa.print <- taxa 
rownames(taxa.print) <- NULL 
head(taxa.print)

# Export assignTaxonomy Assignments as CSV
# Convert the taxonomy assignment matrix to a data frame
taxa_df <- as.data.frame(taxa)

# Write it to a CSV file
write.csv(taxa_df, file="taxonomy_assignments_minBoot80_for_PICRUsT.csv", row.names=TRUE)

head(taxa_df, 20)

# ============================================================================
# ASSIGN TAXONOMY (DECIPHER) - OPTIONAL
# ============================================================================

# Why use DECIPHER/SILVA?
# - SILVA database is more comprehensive than RDP database (more 16s sequences)
# - Validation studies have found it more accurate than using RDP
# - DECIPHER often runs faster than assignTaxonomy

# NOTE: The following DECIPHER code is commented out as it was set to eval=FALSE in the original
# To use DECIPHER, uncomment the following code and ensure you have the SILVA training set

# Assign taxonomy using DECIPHER
# dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
# load("data/SILVA_SSU_r138_2_2024.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
# ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
# ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
# taxid <- t(sapply(ids, function(x) {
#         m <- match(ranks, x$rank)
#         taxa <- x$taxon[m]
#         taxa[startsWith(taxa, "unclassified_")] <- NA
#         taxa
# }))
# colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

# taxid.print <- taxid # Removing sequence rownames for display only
# rownames(taxid.print) <- NULL
# head(taxid.print)

# Export DECIPHER Assignments as CSV
# taxid_df <- as.data.frame(taxid)
# write.csv(taxid_df, file="taxonomy_assignments_DECIPHER.csv", row.names=TRUE)

# ============================================================================
# ALPHA DIVERSITY METRICS
# ============================================================================

theme_set(theme_bw())

# Using MapFile as data frame for sample information
map <- read.csv("data/MapFileHerb(in).csv", row.names = 1, check.names = FALSE)

colnames(map)

# Renaming the Total volume column to contain a valid micro symbol
colnames(map)[colnames(map) == "Total Volume (\xb5L)"] <- "Total Volume (µL)"

# Write the map file as a TSV file for further downstream analysis
write.table(map, file = "output/herb_metadata.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

# Get Rid of _L001 at end of sample names in seqtab.nochim
#rownames(seqtab.nochim) <- gsub("_L001$", "", rownames(seqtab.nochim))

# Construct phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(map), 
               tax_table(as.matrix(taxa))   # use 'taxa' instead
               # or taxa (if using assignTaxonomy from DADA2, not DECIPHER)
)

# Filter out Mitochondria and Chloroplast DNA
# Filter out mitochondrial and chloroplast sequences using Phyloseq
ps_filtered <- subset_taxa(ps, Family != "mitochondria" & Class != "Chloroplast")

# ============================================================================
# ORDINATION PLOTS
# ============================================================================

# Make MDS plot using phyloseq

# All Treatments NMDS using polygons
GP.ord <- ordinate(ps_filtered, "NMDS", "bray")
# Base plot – color & shape by Herbicide
p2 <- plot_ordination(ps_filtered, GP.ord, type = "samples",
                      color = "Herbicide", shape = "Herbicide")
# Add polygons by Concentration
p2 + geom_polygon(aes(fill = factor(Concentration, 
                                    levels = c(0, 0.125, 0.25, 0.5, 1),
                                    labels = c("Control", "12.5%", "25%", "50%", "100%")), 
                      group = interaction(Herbicide, Concentration)),
                  alpha = 0.7) +
  geom_point(size = 3) +
  labs(
    title = "nMDS by Treatment & Concentration",
    x = "nMDS1",
    y = "nMDS2"
  ) + 
  scale_fill_brewer(type = "qual", palette = "Set1", name = "% Recommended Dose")

# 2-4D Treatments NMDS with polygons
# Filter phyloseq object for 2-4 D and Control samples only
ps_24d <- subset_samples(ps_filtered, Herbicide %in% c("Control", "2-4 D"))
# Ordinate the filtered dataset
GP.ord.24d <- ordinate(ps_24d, "NMDS", "bray")
# Base plot – color & shape by Herbicide
p2_24d <- plot_ordination(ps_24d, GP.ord.24d, type = "samples",
                          color = "Herbicide", shape = "Herbicide")
# Add polygons by Concentration
p2_24d + geom_polygon(aes(fill = factor(Concentration, 
                                        levels = c(0, 0.125, 0.25, 0.5, 1),
                                        labels = c("Control", "12.5%", "25%", "50%", "100%")), 
                          group = interaction(Herbicide, Concentration)),
                      alpha = 0.7) +
  geom_point(size = 3) +
  labs(
    title = "nMDS 2-4 D",
    x = "nMDS1",
    y = "nMDS2"
  ) + 
  scale_fill_brewer(type = "qual", palette = "Set1", name = "% Recommended Dose")

# Glyphosate NMDS Plot with polygons
# Filter phyloseq object for Glyphosate and Control samples only
ps_glyph <- subset_samples(ps_filtered, Herbicide %in% c("Control", "Glyphosate"))
# Ordinate the filtered dataset
GP.ord.glyph <- ordinate(ps_glyph, "NMDS", "bray")
# Base plot – color & shape by Herbicide
p2_glyph <- plot_ordination(ps_glyph, GP.ord.glyph, type = "samples",
                            color = "Herbicide", shape = "Herbicide")
# Add polygons by Concentration
p_nmds_glyph <- p2_glyph + geom_polygon(aes(fill = factor(Concentration, 
                                                          levels = c(0, 0.125, 0.25, 0.5, 1),
                                                          labels = c("Control", "12.5%", "25%", "50%", "100%")), 
                                            group = interaction(Herbicide, Concentration)),
                                        alpha = 0.7) +
  geom_point(size = 3) +
  labs(
    title = "nMDS Glyphosate",
    x = "nMDS1",
    y = "nMDS2"
  ) + 
  scale_fill_brewer(type = "qual", palette = "Set1", name = "% Recommended Dose")

print(p_nmds_glyph)
ggsave("Figures/NMDS_glyph_polygons.png", plot = p_nmds_glyph, width = 8, height = 6, dpi = 300)

# 2-4 D NMDS Centroids
# Filter phyloseq object for 2-4 D and Control samples only
ps_24d <- subset_samples(ps_filtered, Herbicide %in% c("Control", "2-4 D"))
# Ordinate the filtered dataset
GP.ord.24d <- ordinate(ps_24d, "NMDS", "bray")
# Base plot – color & shape by Herbicide
p3_24d <- plot_ordination(ps_24d, GP.ord.24d, type = "samples",
                          color = "Herbicide")

# Calculate centroids manually
centroid_data <- p2_24d$data %>%
  mutate(Concentration_Label = factor(Concentration, 
                                      levels = c(0, 0.125, 0.25, 0.5, 1),
                                      labels = c("Control", "12.5%", "25%", "50%", "100%"))) %>%
  group_by(Concentration_Label) %>%
  summarise(
    NMDS1_mean = mean(NMDS1),
    NMDS2_mean = mean(NMDS2),
    .groups = "drop"
  )

# Create plot
p_nmds_24d <- p3_24d + 
  geom_point(aes(color = factor(Concentration, 
                                levels = c(0, 0.125, 0.25, 0.5, 1),
                                labels = c("Control", "12.5%", "25%", "50%", "100%"))),
             size = 3, alpha = 0.6) +  # Individual points
  geom_point(data = centroid_data, 
             aes(x = NMDS1_mean, y = NMDS2_mean, fill = Concentration_Label),
             size = 8, 
             shape = 21, 
             color = "black", stroke = 2, inherit.aes = FALSE) +  # Centroids
  labs(
    title = "2,4-D",
    x = "nMDS1",
    y = "nMDS2"
  ) + 
  scale_fill_brewer(type = "qual", palette = "Set1", name = "% Recommended Dose") +
  scale_color_brewer(type = "qual", palette = "Set1", guide = "none")+
  theme(
    text = element_text(size = 18),            # General text size
    axis.title = element_text(size = 20),      # Axis titles
    axis.text = element_text(size = 16),       # Axis tick labels
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5), # Plot title
    legend.title = element_text(size = 18),    # Legend title
    legend.text = element_text(size = 16)      # Legend labels
  )

print(p_nmds_24d)
ggsave("Figures/NMDS_24D_centroids.png", plot = p_nmds_24d, width = 8, height = 5, dpi = 300)

# Glyphosate NMDS Centroids
# Filter phyloseq object for Glyphosate and Control samples only
ps_glyph <- subset_samples(ps_filtered, Herbicide %in% c("Control", "Glyphosate"))
# Ordinate the filtered dataset
GP.ord.glyph <- ordinate(ps_glyph, "NMDS", "bray")
# Base plot – color & shape by Herbicide
p3_glyph <- plot_ordination(ps_glyph, GP.ord.glyph, type = "samples",
                            color = "Herbicide")

# Calculate centroids manually
centroid_data <- p3_glyph$data %>%
  mutate(Concentration_Label = factor(Concentration, 
                                      levels = c(0, 0.125, 0.25, 0.5, 1),
                                      labels = c("Control", "12.5%", "25%", "50%", "100%"))) %>%
  group_by(Concentration_Label) %>%
  summarise(
    NMDS1_mean = mean(NMDS1),
    NMDS2_mean = mean(NMDS2),
    NMDS1_se = sd(NMDS1)/sqrt(n()),
    NMDS2_se = sd(NMDS2)/sqrt(n()),
    .groups = "drop"
  )

# Create plot
p_nmds_glyph <- p3_glyph + 
  geom_point(aes(color = factor(Concentration, 
                                levels = c(0, 0.125, 0.25, 0.5, 1),
                                labels = c("Control", "12.5%", "25%", "50%", "100%"))),
             size = 3, alpha = 0.6) +  # Individual points
  geom_point(data = centroid_data, 
             aes(x = NMDS1_mean, y = NMDS2_mean, fill = Concentration_Label),
             size = 8, shape = 21, color = "black", stroke = 2, inherit.aes = FALSE) +  # Centroids
  labs(
    title = "Glyphosate",
    x = "nMDS1",
    y = "nMDS2"
  ) + 
  scale_fill_brewer(type = "qual", palette = "Set1", name = "% Recommended Dose") +
  scale_color_brewer(type = "qual", palette = "Set1", guide = "none")+
  theme(
    text = element_text(size = 18),            # General text size
    axis.title = element_text(size = 20),      # Axis titles
    axis.text = element_text(size = 16),       # Axis tick labels
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5), # Plot title
    legend.title = element_text(size = 18),    # Legend title
    legend.text = element_text(size = 16)      # Legend labels
  )

print(p_nmds_glyph)
ggsave("Figures/NMDS_Glyphosate_centroids.png", plot = p_nmds_glyph, width = 8, height = 5, dpi = 300)

# Combine the 2-4d and Glyphosate nMDS centroid plots, and keep same scale each
# 1) Choose common limits (adjust if needed)
xlim_nm <- c(-0.45, 0.75)
ylim_nm <- c(-0.25, 0.35)

# 2) Apply limits (coord_cartesian keeps points outside from being dropped)
p_nmds_24d_lims <- p_nmds_24d +
  coord_equal() +
  coord_cartesian(xlim = xlim_nm, ylim = ylim_nm)

p_nmds_glyph_lims <- p_nmds_glyph +
  coord_equal() +
  coord_cartesian(xlim = xlim_nm, ylim = ylim_nm)

# 3) Print (optional, for the Rmd preview)
print(p_nmds_24d_lims)
print(p_nmds_glyph_lims)

# 4) Save with new filenames so you don't overwrite the originals
ggsave("Figures/NMDS_24D_centroids_commonAxes.png",
       plot = p_nmds_24d_lims, width = 8, height = 5, dpi = 300)

ggsave("Figures/NMDS_Glyphosate_centroids_commonAxes.png",
       plot = p_nmds_glyph_lims, width = 8, height = 5, dpi = 300)

# All treatments MDS using ellipses
# MDS on your filtered phyloseq object
ord <- ordinate(ps_filtered, method = "MDS", distance = "bray")
# Base plot (points still show herbicide)
p <- plot_ordination(ps_filtered, ord, type = "samples",
                     color = "Herbicide", shape = "Herbicide")
# Ellipses: one per Herbicide × Concentration, fill = numeric concentration
p +
  stat_ellipse(
    aes(group = interaction(Herbicide, Concentration),  # separate ellipses
        fill  = factor(Concentration, 
                       levels = c(0, 0.125, 0.25, 0.5, 1),
                       labels = c("Control", "12.5%", "25%", "50%", "100%")),               # convert to factor
        color = Herbicide),
    type = "t", level = 0.90, geom = "polygon",
    alpha = 0.70, linewidth = 0.6
  ) +
  geom_point(size = 3) +
  # discrete palette
  scale_fill_brewer(type = "qual", palette = "Set1", name = "% Recommended Dose") +
  labs(
    title = "MDS by Treatment & Concentration",
    x = "MDS1",
    y = "MDS2"
  ) + 
  coord_fixed(ratio = 1.0)  # make Y taller

# 2-4D Treatment MDS and ellipses
# Filter phyloseq object for 2-4 D and Control samples only
ps_24d <- subset_samples(ps_filtered, Herbicide %in% c("Control", "2-4 D"))
# MDS on your filtered phyloseq object
ord_24d <- ordinate(ps_24d, method = "MDS", distance = "bray")
# Base plot (points still show herbicide)
p_24d <- plot_ordination(ps_24d, ord_24d, type = "samples",
                         color = "Herbicide", shape = "Herbicide")
# Ellipses: one per Herbicide × Concentration, fill = numeric concentration
p_24d +
  stat_ellipse(
    aes(group = interaction(Herbicide, Concentration),  # separate ellipses
        fill  = factor(Concentration, 
                       levels = c(0, 0.125, 0.25, 0.5, 1),
                       labels = c("Control", "12.5%", "25%", "50%", "100%")),               # convert to factor
        color = Herbicide),
    type = "t", level = 0.90, geom = "polygon",
    alpha = 0.70, linewidth = 0.6
  ) +
  geom_point(size = 3) +
  # discrete palette
  scale_fill_brewer(type = "qual", palette = "Set1", name = "% Recommended Dose") +
  labs(
    title = "MDS by 2-4 D Treatment & Concentration",
    x = "MDS1",
    y = "MDS2"
  ) + 
  coord_fixed(ratio = 1.0)

# Glyphosate Treatment MDS Plot with ellipses
# Filter phyloseq object for Glyphosate and Control samples only
ps_glyph <- subset_samples(ps_filtered, Herbicide %in% c("Control", "Glyphosate"))
# MDS on your filtered phyloseq object
ord_glyph <- ordinate(ps_glyph, method = "MDS", distance = "bray")
# Base plot (points still show herbicide)
p_glyph <- plot_ordination(ps_glyph, ord_glyph, type = "samples",
                           color = "Herbicide", shape = "Herbicide")
# Ellipses: one per Herbicide × Concentration, fill = numeric concentration
p_glyph +
  stat_ellipse(
    aes(group = interaction(Herbicide, Concentration),  # separate ellipses
        fill  = factor(Concentration, 
                       levels = c(0, 0.125, 0.25, 0.5, 1),
                       labels = c("Control", "12.5%", "25%", "50%", "100%")),               # convert to factor
        color = Herbicide),
    type = "t", level = 0.90, geom = "polygon",
    alpha = 0.70, linewidth = 0.6
  ) +
  geom_point(size = 3) +
  # discrete palette
  scale_fill_brewer(type = "qual", palette = "Set1", name = "% Recommended Dose") +
  labs(
    title = "MDS by Glyphosate Treatment & Concentration",
    x = "MDS1",
    y = "MDS2"
  ) + 
  coord_fixed(ratio = 1.8)

# 2-4D MDS with Polygons
# Filter phyloseq object for 2-4 D and Control samples only
ps_24d <- subset_samples(ps_filtered, Herbicide %in% c("Control", "2-4 D"))
# Ordinate the filtered dataset using MDS
GP.ord.24d <- ordinate(ps_24d, "MDS", "bray")
# Base plot – color & shape by Herbicide
p2_24d <- plot_ordination(ps_24d, GP.ord.24d, type = "samples",
                          color = "Herbicide", shape = "Herbicide")
# Add polygons by Concentration
p2_24d + geom_polygon(aes(fill = factor(Concentration, 
                                        levels = c(0, 0.125, 0.25, 0.5, 1),
                                        labels = c("Control", "12.5%", "25%", "50%", "100%")), 
                          group = interaction(Herbicide, Concentration)),
                      alpha = 0.7) +
  geom_point(size = 3) +
  labs(
    title = "MDS by 2-4 D Treatment & Concentration",
    x = "MDS1",
    y = "MDS2"
  ) + 
  scale_fill_brewer(type = "qual", palette = "Set1", name = "% Recommended Dose")

# Glyphosate MDS Plot with polygons
# Filter phyloseq object for Glyphosate and Control samples only
ps_glyph <- subset_samples(ps_filtered, Herbicide %in% c("Control", "Glyphosate"))
# Ordinate the filtered dataset using MDS
GP.ord.glyph <- ordinate(ps_glyph, "MDS", "bray")
# Base plot – color & shape by Herbicide
p2_glyph <- plot_ordination(ps_glyph, GP.ord.glyph, type = "samples",
                            color = "Herbicide", shape = "Herbicide")
# Add polygons by Concentration
p2_glyph + geom_polygon(aes(fill = factor(Concentration, 
                                          levels = c(0, 0.125, 0.25, 0.5, 1),
                                          labels = c("Control", "12.5%", "25%", "50%", "100%")), 
                            group = interaction(Herbicide, Concentration)),
                        alpha = 0.7) +
  geom_point(size = 3) +
  labs(
    title = "MDS by Glyphosate Treatment & Concentration",
    x = "MDS1",
    y = "MDS2"
  ) + 
  scale_fill_brewer(type = "qual", palette = "Set1", name = "% Recommended Dose")

# Creating shorter names for the ASVs rather than the full DNA sequence to simplify things for when we do the visualizations
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

plot_richness(ps, x = "Herbicide", measures = c("Shannon", "Simpson"))

# ============================================================================
# EXPORT ALPHA DIVERSITY DATA
# ============================================================================

# Create data frame for alpha diversity results:
alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
head(alpha_div)

# Combine the map data frame and alpha_div data frames
map_alpha <- cbind(map, alpha_div[rownames(map), ])

# Export the map_alpha data frame as a CSV
write.csv(map_alpha, file="Alpha_Diversity_Metrics_Herb_Soil_NATHAN.csv", row.names=TRUE)

# ============================================================================
# CREATE BOXPLOT FIGURES FOR ASV RICHNESS
# ============================================================================

# Setting the order of the treatment groups for the figure
map_alpha$Herbicide <- factor(map_alpha$Herbicide, 
                              levels = c("Control", "2-4 D", "Glyphosate"))

# Boxplot for "Observed" Metric - All treatments
Observed_ASV_boxplot <- ggplot(map_alpha, aes(x = Herbicide, y = Observed, fill = Herbicide)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75),  # Align with boxes
             size = 2, color = "black", shape = 21, stroke = 0.4, fill = "black") +
  labs(
    title = "Observed ASV Richness: All Treatments",
    x = "Herbicide Treatment",
    y = "Observed ASV Richness"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_classic() +
  theme(legend.position = "right")

Observed_ASV_boxplot

# Observed 2-4D
# Filter for 2-4 D and Control samples only
map_alpha_24d <- map_alpha[map_alpha$Herbicide %in% c("Control", "2-4 D"), ]

# Create a combined group variable with percentages
map_alpha_24d$Treatment_Group <- ifelse(map_alpha_24d$Herbicide == "Control", 
                                        "Control", 
                                        case_when(
                                          map_alpha_24d$Concentration == 0.125 ~ "12.5%",
                                          map_alpha_24d$Concentration == 0.25 ~ "25%",
                                          map_alpha_24d$Concentration == 0.5 ~ "50%",
                                          map_alpha_24d$Concentration == 1 ~ "100%"
                                        ))

# Set the order of treatment groups
map_alpha_24d$Treatment_Group <- factor(map_alpha_24d$Treatment_Group, 
                                        levels = c("Control", "12.5%", "25%", "50%", "100%"))

Observed_ASV_boxplot_24d <- ggplot(map_alpha_24d, aes(x = Treatment_Group, y = Observed, fill = Herbicide)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75),
             size = 2, color = "black", shape = 21, stroke = 0.4, fill = "black") +
  labs(
    title = "Observed ASV Richness: 2-4 D",
    x = "% of Recommended Dose",
    y = "Observed ASV Richness"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_classic() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
Observed_ASV_boxplot_24d

# Observed Glyphosate
# Filter for Glyphosate and Control samples only
map_alpha_glyph <- map_alpha[map_alpha$Herbicide %in% c("Control", "Glyphosate"), ]

# Create a combined group variable with percentages
map_alpha_glyph$Treatment_Group <- ifelse(map_alpha_glyph$Herbicide == "Control", 
                                          "Control", 
                                          case_when(
                                            map_alpha_glyph$Concentration == 0.125 ~ "12.5%",
                                            map_alpha_glyph$Concentration == 0.25 ~ "25%",
                                            map_alpha_glyph$Concentration == 0.5 ~ "50%",
                                            map_alpha_glyph$Concentration == 1 ~ "100%"
                                          ))

# Set the order of treatment groups
map_alpha_glyph$Treatment_Group <- factor(map_alpha_glyph$Treatment_Group, 
                                          levels = c("Control", "12.5%", "25%", "50%", "100%"))

Observed_ASV_boxplot_glyph <- ggplot(map_alpha_glyph, aes(x = Treatment_Group, y = Observed, fill = Herbicide)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75),
             size = 2, color = "black", shape = 21, stroke = 0.4, fill = "black") +
  labs(
    title = "Observed ASV Richness: Glyphosate",
    x = "% of Recommended Dose",
    y = "Observed ASV Richness"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_classic() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
Observed_ASV_boxplot_glyph

# Save the plot
ggsave("Observed_ASV_Boxplot.png", plot = Observed_ASV_boxplot, width = 6, height = 4)

# Boxplot for "Shannon" Metric - All Treatments
Shannon_ASV_boxplot <- ggplot(map_alpha, aes(x = Herbicide, y = Shannon, fill = Herbicide)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75),  # Align with boxes
             size = 2, color = "black", shape = 21, stroke = 0.4, fill = "black") +
  labs(
    title = "Shannon Diversity: All Treatments",
    x = "Herbicide Treatment",
    y = "Shannon ASV Richness"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_classic() +
  theme(legend.position = "right")

Shannon_ASV_boxplot

# Shannon 2-4D
# Filter for 2-4 D and Control samples only
map_alpha_24d <- map_alpha[map_alpha$Herbicide %in% c("Control", "2-4 D"), ]

# Create a combined group variable with percentages
map_alpha_24d$Treatment_Group <- ifelse(map_alpha_24d$Herbicide == "Control", 
                                        "Control", 
                                        case_when(
                                          map_alpha_24d$Concentration == 0.125 ~ "12.5%",
                                          map_alpha_24d$Concentration == 0.25 ~ "25%",
                                          map_alpha_24d$Concentration == 0.5 ~ "50%",
                                          map_alpha_24d$Concentration == 1 ~ "100%"
                                        ))

# Set the order of treatment groups
map_alpha_24d$Treatment_Group <- factor(map_alpha_24d$Treatment_Group, 
                                        levels = c("Control", "12.5%", "25%", "50%", "100%"))

Shannon_ASV_boxplot_24d <- ggplot(map_alpha_24d, aes(x = Treatment_Group, y = Shannon, fill = Herbicide)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75),
             size = 2, color = "black", shape = 21, stroke = 0.4, fill = "black") +
  labs(
    title = "Shannon Diversity: 2-4 D",
    x = "% of Recommended Dose",
    y = "Shannon Diversity Index"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_classic() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
Shannon_ASV_boxplot_24d

# Shannon Glyphosate
# Filter for Glyphosate and Control samples only
map_alpha_glyph <- map_alpha[map_alpha$Herbicide %in% c("Control", "Glyphosate"), ]

# Create a combined group variable with percentages
map_alpha_glyph$Treatment_Group <- ifelse(map_alpha_glyph$Herbicide == "Control", 
                                          "Control", 
                                          case_when(
                                            map_alpha_glyph$Concentration == 0.125 ~ "12.5%",
                                            map_alpha_glyph$Concentration == 0.25 ~ "25%",
                                            map_alpha_glyph$Concentration == 0.5 ~ "50%",
                                            map_alpha_glyph$Concentration == 1 ~ "100%"
                                          ))

# Set the order of treatment groups
map_alpha_glyph$Treatment_Group <- factor(map_alpha_glyph$Treatment_Group, 
                                          levels = c("Control", "12.5%", "25%", "50%", "100%"))

Shannon_ASV_boxplot_glyph <- ggplot(map_alpha_glyph, aes(x = Treatment_Group, y = Shannon, fill = Herbicide)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75),
             size = 2, color = "black", shape = 21, stroke = 0.4, fill = "black") +
  labs(
    title = "Shannon Diversity: Glyphosate",
    x = "% of Recommended Dose",
    y = "Shannon Diversity Index"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_classic() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
Shannon_ASV_boxplot_glyph

# Save the plot
ggsave("Shannon_ASV_Boxplot.png", plot = Shannon_ASV_boxplot, width = 6, height = 4)

# Boxplot for "Simpson" Metric - All Treatments
Simpson_ASV_boxplot <- ggplot(map_alpha, aes(x = Herbicide, y = Simpson, fill = Herbicide)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75),  # Align with boxes
             size = 2, color = "black", shape = 21, stroke = 0.4, fill = "black") +
  labs(
    title = "Simpson Diversity: All Treatments",
    x = "Herbicide Treatment",
    y = "Simpson ASV Richness"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_classic() +
  theme(legend.position = "right")

Simpson_ASV_boxplot

# Simpson 2-4D
# Filter for 2-4 D and Control samples only
map_alpha_24d <- map_alpha[map_alpha$Herbicide %in% c("Control", "2-4 D"), ]

# Create a combined group variable with percentages
map_alpha_24d$Treatment_Group <- ifelse(map_alpha_24d$Herbicide == "Control", 
                                        "Control", 
                                        case_when(
                                          map_alpha_24d$Concentration == 0.125 ~ "12.5%",
                                          map_alpha_24d$Concentration == 0.25 ~ "25%",
                                          map_alpha_24d$Concentration == 0.5 ~ "50%",
                                          map_alpha_24d$Concentration == 1 ~ "100%"
                                        ))

# Set the order of treatment groups
map_alpha_24d$Treatment_Group <- factor(map_alpha_24d$Treatment_Group, 
                                        levels = c("Control", "12.5%", "25%", "50%", "100%"))

Simpson_ASV_boxplot_24d <- ggplot(map_alpha_24d, aes(x = Treatment_Group, y = Simpson, fill = Herbicide)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75),
             size = 2, color = "black", shape = 21, stroke = 0.4, fill = "black") +
  labs(
    title = "Simpson Diversity: 2-4 D",
    x = "% of Recommended Dose",
    y = "Simpson Diversity Index"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_classic() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
Simpson_ASV_boxplot_24d

# Simpson Glyphosate
# Filter for Glyphosate and Control samples only
map_alpha_glyph <- map_alpha[map_alpha$Herbicide %in% c("Control", "Glyphosate"), ]

# Create a combined group variable with percentages
map_alpha_glyph$Treatment_Group <- ifelse(map_alpha_glyph$Herbicide == "Control", 
                                          "Control", 
                                          case_when(
                                            map_alpha_glyph$Concentration == 0.125 ~ "12.5%",
                                            map_alpha_glyph$Concentration == 0.25 ~ "25%",
                                            map_alpha_glyph$Concentration == 0.5 ~ "50%",
                                            map_alpha_glyph$Concentration == 1 ~ "100%"
                                          ))

# Set the order of treatment groups
map_alpha_glyph$Treatment_Group <- factor(map_alpha_glyph$Treatment_Group, 
                                          levels = c("Control", "12.5%", "25%", "50%", "100%"))

Simpson_ASV_boxplot_glyph <- ggplot(map_alpha_glyph, aes(x = Treatment_Group, y = Simpson, fill = Herbicide)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75),
             size = 2, color = "black", shape = 21, stroke = 0.4, fill = "black") +
  labs(
    title = "Simpson Diversity: Glyphosate",
    x = "% of Recommended Dose",
    y = "Simpson Diversity Index"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_classic() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
Simpson_ASV_boxplot_glyph

# Save the plot
ggsave("Simpson_ASV_Boxplot.png", plot = Simpson_ASV_boxplot, width = 6, height = 4)

# ============================================================================
# RELATIVE ABUNDANCE STACKED BAR CHART FIGURES
# ============================================================================

# Abundance Stacked Bar All Treatments
top20    <- names(sort(taxa_sums(ps_filtered), decreasing = TRUE))[1:20]
ps.top20 <- prune_taxa(top20, ps_filtered)
ps.top20 <- transform_sample_counts(ps.top20, function(x) x / sum(x))

# 2) Melt to a data frame
df <- psmelt(ps.top20)

# Removing NA values from the data
df <- subset(df, !is.na(Genus))

# (optional) order facets/ticks
df$Herbicide     <- factor(df$Herbicide, levels = c("Control","2-4 D","Glyphosate"))
df$Concentration <- factor(df$Concentration, 
                           levels = c("0","0.125","0.25","0.5","1"),
                           labels = c("Control", "12.5%", "25%", "50%", "100%"))

# 3) Plot with position="fill" so each stacked bar sums to 1
stacked_bar_plot <- ggplot(df, aes(x = Concentration, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Herbicide, scales = "free_x") +
  ylab("Relative Abundance") + xlab("% of Recommended Dose")+
  ggtitle("Relative Abundance of Top 20 Genera") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(df$Genus))))

stacked_bar_plot

# Save the stacked bar plot
ggsave("top20_genera_relative_abundance.png", 
       plot = stacked_bar_plot, 
       width = 8, height = 6, dpi = 300)

# Abundance Stacked Bar 2-4D
top20 <- names(sort(taxa_sums(ps_filtered), decreasing = TRUE))[1:20]
ps.top20 <- prune_taxa(top20, ps_filtered)
ps.top20 <- transform_sample_counts(ps.top20, function(x) x / sum(x))

# Filter for 2-4 D and Control samples only
ps.top20_24d <- subset_samples(ps.top20, Herbicide %in% c("Control", "2-4 D"))

# Melt to a data frame
df_24d <- psmelt(ps.top20_24d)
# Removing NA values from the data
df_24d <- subset(df_24d, !is.na(Genus))
# Order factors
df_24d$Herbicide <- factor(df_24d$Herbicide, levels = c("Control", "2-4 D"),
                           labels = c("Control", "2,4-D"))  # Add labels parameter
df_24d$Concentration <- factor(df_24d$Concentration, 
                               levels = c("0","0.125","0.25","0.5","1"),
                               labels = c("Control", "12.5%", "25%", "50%", "100%"))

# Plot
stacked_bar_plot_24d <- ggplot(df_24d, aes(x = Concentration, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Herbicide, scales = "free_x") +
  ylab("Relative Abundance") + xlab("% of Recommended Dose") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(df_24d$Genus))))+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    strip.text   = element_text(size = 15),   # facet labels
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )

stacked_bar_plot_24d

# Save the plot
ggsave("Figures/top20_genera_24d_effects.png", 
       plot = stacked_bar_plot_24d, 
       width = 8, height = 6, dpi = 300)

# Abundance Stacked Bar Glyphosate
# Filter for Glyphosate and Control samples only
ps.top20_glyph <- subset_samples(ps.top20, Herbicide %in% c("Control", "Glyphosate"))

# Melt to a data frame
df_glyph <- psmelt(ps.top20_glyph)
# Removing NA values from the data
df_glyph <- subset(df_glyph, !is.na(Genus))
# Order factors
df_glyph$Herbicide <- factor(df_glyph$Herbicide, levels = c("Control", "Glyphosate"))
df_glyph$Concentration <- factor(df_glyph$Concentration, 
                                 levels = c("0","0.125","0.25","0.5","1"),
                                 labels = c("Control", "12.5%", "25%", "50%", "100%"))

# Plot
stacked_bar_plot_glyph <- ggplot(df_glyph, aes(x = Concentration, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~Herbicide, scales = "free_x") +
  ylab("Relative Abundance") + xlab("% of Recommended Dose") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(df_glyph$Genus))))+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    strip.text   = element_text(size = 15),   # facet labels
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )

stacked_bar_plot_glyph

# Save the plot
ggsave("Figures/top20_genera_glyphosate_effects.png", 
       plot = stacked_bar_plot_glyph, 
       width = 8, height = 6, dpi = 300)

# Save the plot
ggsave("Alpha_Diversity_Stacked_Bar_chart.png", plot = stacked_bar_plot, width = 6, height = 4)

# ============================================================================
# DADA2 TO PICRUST2 OUTPUT FORMATTING
# ============================================================================

# Creates outputs formatted to be used as inputs for PICRUsT2

# Change seq headers to asv number (ASV_1, ASV_2...)
# Citation: https://github.com/picrust/picrust2/issues/136#issuecomment-743696114
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Create and write out a fasta of our final ASV seqs with the new ID's
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "output/ASVs.fa")

# ASV count table
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "output/ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# Tax table
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "output/ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

# Merge ASV abundance and taxonomy into one file
OTU_TAX_table <- merge(asv_tab, asv_tax, by=0)
write.table(OTU_TAX_table, "~/dada2/OTU_TAX_table.tsv", sep = "\t", quote=F, col.names=NA)
