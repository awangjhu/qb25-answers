# RNA-seq DE (GTEx whole blood)
# Alan Wang
# 11/5/2025

# Packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(broom)
  library(ggrepel)
})


# Path
counts_fp <- "/Users/cmdb/qb25-answers/week8/gtex_whole_blood_counts_downsample.txt"
meta_fp   <- "/Users/cmdb/qb25-answers/week8/gtex_metadata_downsample.txt"
genes_fp  <- "/Users/cmdb/qb25-answers/week8/gene_locations.txt"
out_dir   <- "/Users/cmdb/qb25-answers/week8"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# Exercise 1: Data preprocessing
# Step 1.1: Loading data and importing libraries
# Load count data 
counts_df <- read_delim(counts_fp, delim = ",", show_col_types = FALSE) %>%
  column_to_rownames(var = "GENE_NAME")
# Load metadata 
metadata_df <- read_delim(meta_fp, delim = ",", show_col_types = FALSE)
# Convert to appropriate types for DESeq2
metadata_df <- metadata_df %>%
  mutate(
    SEX = factor(SEX, levels = c("female", "male")),  # female as reference
    DTHHRDY = factor(DTHHRDY),
    AGE = as.integer(AGE)
  )
# Peek at the data
# This is suggested by Chatgpt to visually look whether the data is correct before continuing, similar to using print to check
cat("\n=== Counts Data (first 5 genes, first 5 samples) ===\n")
print(counts_df[1:5, 1:5])

cat("\n=== Metadata (first 5 samples) ===\n")
print(head(metadata_df, 5))

cat("\nDimensions: ", nrow(counts_df), " genes x ", ncol(counts_df), " samples\n")

# Step 1.2: Create a DESeq2 object
# Checks whether the column names and metadata table match in the same order before proceeding, mentioned by Chatgpt
stopifnot(all(colnames(counts_df) == metadata_df$SUBJECT_ID))
cat("Verified: Sample IDs in count matrix columns match metadata rows in correct order.\n")

# Create DESeq2 
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts_df),
  colData = metadata_df,
  design = ~ DTHHRDY + AGE + SEX
)

cat("DESeq2 object created successfully.\n")


# Step 1.3: Normalization and PCA
# Apply variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# Generate PCA plots
p_sex <- plotPCA(vsd, intgroup = "SEX") + 
  ggtitle("PCA of GTEx Whole Blood Gene Expression (colored by SEX)") +
  theme_bw()

p_age <- plotPCA(vsd, intgroup = "AGE") + 
  ggtitle("PCA of GTEx Whole Blood Gene Expression (colored by AGE)") +
  theme_bw()

p_dthhrdy <- plotPCA(vsd, intgroup = "DTHHRDY") + 
  ggtitle("PCA of GTEx Whole Blood Gene Expression (colored by DTHHRDY)") +
  theme_bw()

# Save plots
ggsave(file.path(out_dir, "pca_by_sex.png"), p_sex, width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, "pca_by_age.png"), p_age, width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, "pca_by_dthhrdy.png"), p_dthhrdy, width = 7, height = 5, dpi = 300)

cat("PCA plots saved to outputs/ directory.\n")

# Interpretation of PCA Results
# PC1 explains around 48% of the variance while PC2 explains around 7% in gene expression
# When coloring by DTHHRDY, or the cause of death, it seems like PC1 shows a very clear separation
# as ventilator samples cluster near the right side of the plot and the fast death of natural causes
# tend to cluster towards the left side. This means that PC1 is associated with death. For sex,
# there does not seem to be a strong separation so its relative negligible on global gene expression.
# For age, there is weak DTHHRDY in explaining variance. Thus, PC1 is primarily focused and captures
# death than age or sex while PC2 is not associated with any metadata.

# Exercise 2: Perform differential expression analysis
# Step 2.1: Perform a "homemade" test for differential expression between sexes
# VST-normalized expression matrix
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble(.name_repair = "minimal")

vsd_df <- bind_cols(metadata_df, vsd_df)

#WASH7P for sex-differential expression using linear regression
cat("\n=== Homemade Linear Model: WASH7P ===\n")
m1 <- lm(formula = WASH7P ~ DTHHRDY + AGE + SEX, data = vsd_df)
wash7p_results <- summary(m1) %>% tidy()
print(wash7p_results)

# Does WASH7P show significant evidence of sex-differential expression?
# There is no statistical significance in WASH7P between males and females because
# WASH7P shows a small log2fold change, at 0.09, and a high adjusted p-value. Thus,
# WASH7P is expressed at similar levels in both sexes

# Test SLC25A47 for sex-differential expression
cat("\n=== Homemade Linear Model: SLC25A47 ===\n")
m2 <- lm(formula = SLC25A47 ~ DTHHRDY + AGE + SEX, data = vsd_df)
slc25a47_results <- summary(m2) %>% tidy()
print(slc25a47_results)

# Does SLC25A47 show evidence of sex-differential expression?
# Yes, SLC25A47 is significantly upregulate in males because it has a large log2foldchange
# at 3.06 and low adjusted p-value at 8.3e-7, so it is very much upregaulated in males.


# Step 2.2: Perform differential expression analysis with DESeq2
# Run DESeq2 analysis on raw counts 
dds <- DESeq(dds)
cat("DESeq2 analysis is complete \n")


# Step 2.3: Extract and interpret the results for sex differential expression
# Check available result names
res_names <- resultsNames(dds)
cat("\nAvailable results:\n")
print(res_names)

# Results for SEX effect (male vs female)
res_sex <- results(dds, name = "SEX_male_vs_female") %>%
  as_tibble(rownames = "GENE_NAME")

# Count significant genes at 10% FDR
sig_count_sex <- res_sex %>%
  filter(!is.na(padj), padj < 0.10) %>%
  nrow()

cat("\n=== SEX Differential Expression Results ===\n")
cat("Number of genes significantly DE at 10% FDR:", sig_count_sex, "\n")

# How many genes exhibit significant differential expression between males and females at a 10% FDR?
# 262

# Load gene location data
genes_df <- read_delim(genes_fp, delim = "\t", show_col_types = FALSE)

# Merge with chromosome information
res_sex_annot <- res_sex %>%
  left_join(genes_df, by = "GENE_NAME") %>%
  arrange(padj)

# Top hits by chromosome
top_male_up <- res_sex_annot %>%
  filter(!is.na(padj), padj < 0.10, log2FoldChange > 0) %>%
  arrange(padj) %>%
  slice_head(n = 50)

top_female_up <- res_sex_annot %>%
  filter(!is.na(padj), padj < 0.10, log2FoldChange < 0) %>%
  arrange(padj) %>%
  slice_head(n = 50)

cat("\n=== Top Male-Upregulated Genes by Chromosome ===\n")
print(top_male_up %>% count(CHROM, sort = TRUE))

cat("\n=== Top Female-Upregulated Genes by Chromosome ===\n")
print(top_female_up %>% count(CHROM, sort = TRUE))

# Top hits in males vs females?
# For the top male upregulated genes, the majority ar located on the Y chromosome
# and smaller number seam to appear on autosomal genes. For females, the most upregulated
# genes are on the X chromosome.

# Interpretation:
# This pattern should be almost expected because males have a Y chromosome and only one X 
# chromosome and females have two X chromosomes, so any Y-linked genes should show higher
# expression in males while X-linked genes should be more highly expressed in females.

# Check homemade linear model
cat("\n=== Comparing DESeq2 vs Homemade LM for WASH7P and SLC25A47 ===\n")
check_genes <- res_sex_annot %>%
  filter(GENE_NAME %in% c("WASH7P", "SLC25A47"))
print(check_genes)

# Comparison of linear model tests
# Yes, I would say the DESeq2 results are consistent with linear regression analysis.
# Both approaches show a small log2 fold change for WASH7P and large positive log2 fold change for SLC25A47
# so it seems like both the "homemade" regression and DESeq2 model agree in statistical significance.


# Trade-off between false positives and false negatives
# A stringent FDR threshold, at 1%, reduces the chance of false positives although it increases the risk of 
# false negatives by filtering out genes with modest effect size. However, a more lenient FDR threshold like at 20% detects more differentially
# expressed genes, but allows for more false positives to pass. The power to detect differential expression seems to 
# increase when you have a larger sample size or genes have stronger effect size.


# Step 2.4: Extract and interpret the results for differential expression by death classification
# Extract results for death classification effect
res_dthhrdy <- results(dds, name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes") %>%
  as_tibble(rownames = "GENE_NAME")

# Count significant genes at 10% FDR
sig_count_dthhrdy <- res_dthhrdy %>%
  filter(!is.na(padj), padj < 0.10) %>%
  nrow()

# How many genes are differentially expressed according to death classification at a 10% FDR?
# 16069

cat("\n=== DTHHRDY Differential Expression Results ===\n")
cat("Number of genes significantly DE by death classification at 10% FDR:", 
    sig_count_dthhrdy, "\n")

# Step 2.5: Estimating a false positive rate under the null hypothesis
# Permute SEX labels to break true associations (null hypothesis)
set.seed(123)
metadata_permuted <- metadata_df %>%
  mutate(SEX = sample(SEX))

# Create new DESeq2 object with permuted sex
dds_perm <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts_df),
  colData = metadata_permuted,
  design = ~ DTHHRDY + AGE + SEX
)

# Run DESeq2 on permuted data
dds_perm <- DESeq(dds_perm)

# Extract permuted results
res_perm <- results(dds_perm, name = "SEX_male_vs_female") %>%
  as_tibble(rownames = "GENE_NAME")

# Count "significant" genes in permuted analysis
sig_count_perm <- res_perm %>%
  filter(!is.na(padj), padj < 0.10) %>%
  nrow()

cat("\n=== Permutation Null Analysis ===\n")
cat("Significant genes in REAL data (10% FDR):", sig_count_sex, "\n")
cat("Significant genes in PERMUTED data (10% FDR):", sig_count_perm, "\n")

# How many genes appear "significant" in the permuted analysis at a 10% FDR?
# In the actual real data, 262 genes were identified as significantly differentially expressed
# at 10%; however, in the permuted analysis, only 15 genes passed the same threshold. Since the permuted
# data contains no real biological signal, these 15 genes are likely and in fact false positives.
# This tells me that the FDR threshold controls the rate of false discoveries, especially in these 
# large scale sequencing experiments.

# Exercise 3: Visualization - Volcano Plot
# Data for volcano plot
volcano_df <- res_sex %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significance = case_when(
      is.na(padj) ~ "Not tested",
      padj < 0.10 & abs(log2FoldChange) > 1 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Significant" = "red", "Not significant" = "gray60", "Not tested" = "gray90"),
    name = NULL
  ) +
  geom_hline(yintercept = -log10(0.10), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(
    title = "Volcano Plot: Sex-Differential Gene Expression in GTEx Whole Blood",
    subtitle = "Significant = FDR < 10% AND |log2FC| > 1",
    x = "log2(Fold Change) [male vs female]",
    y = "-log10(adjusted p-value)"
  ) +
  theme_bw() +
  theme(legend.position = "top")

top_genes <- volcano_df %>%
  filter(significance == "Significant") %>%
  arrange(padj) %>%
  slice_head(n = 10)

if (nrow(top_genes) > 0) {
  p_volcano <- p_volcano +
    geom_text_repel(
      data = top_genes,
      aes(label = GENE_NAME),
      max.overlaps = 20,
      size = 3
    )
}

# Save volcano plot
ggsave(file.path(out_dir, "volcano_plot_sex.png"), p_volcano, 
       width = 8, height = 6, dpi = 300)

cat("\n Volcano plot saved to outputs/volcano_plot_sex.png\n")