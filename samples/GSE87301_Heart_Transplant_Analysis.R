# =============================================================================
# SECTION 1: Load Required Libraries
# =============================================================================

library(affy)           # Reading CEL files
library(gcrma)          # GCRMA normalization
library(hgu133acdf)     # Array CDF annotation
library(hgu133a.db)     # Probe to gene mapping
library(AnnotationDbi)  # Annotation tools
library(limma)          # Differential expression
library(pheatmap)       # Heatmap visualization


# =============================================================================
# SECTION 2: Load CEL Files
# =============================================================================

# Set path to CEL files
cel_path <- "/home/keerthana/Documents/project/1_project/sample/sample_2"

# List all CEL.gz files
celfiles <- list.files(path = cel_path,
                       pattern = "\\.CEL\\.gz$",
                       full.names = TRUE,
                       ignore.case = TRUE)

# Verify file count (should be 14)
cat("Number of CEL files found:", length(celfiles), "\n")
print(basename(celfiles))

# Read CEL files into AffyBatch object
raw_data <- read.affybatch(filenames = celfiles)
cat("Raw data loaded successfully\n")
cat("Number of samples:", ncol(raw_data), "\n")

# =============================================================================
# SECTION 3: GCRMA Normalization
# =============================================================================

# Set single thread to avoid pthread errors
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")

# Run GCRMA normalization
cat("Running GCRMA normalization...\n")
eset <- gcrma(raw_data)

# Extract expression matrix
exprs_matrix <- exprs(eset)
cat("Normalization complete\n")
cat("Expression matrix dimensions:", dim(exprs_matrix), "\n")

# =============================================================================
# SECTION 4: Assign Group Labels (AR vs NAR)
# =============================================================================

# Get sample names from eset
samples <- basename(sampleNames(eset))

# Assign group based on GSM number in filename
AR_gsm <- c("GSM2327587", "GSM2327591", "GSM2327594",
             "GSM2327596", "GSM2327600", "GSM2327602",
             "GSM2327605", "GSM2327608", "GSM2327609",
             "GSM2327611")

NAR_gsm <- c("GSM2327586", "GSM2327588", "GSM2327589",
              "GSM2327590", "GSM2327592", "GSM2327593",
              "GSM2327595", "GSM2327597", "GSM2327598",
              "GSM2327599", "GSM2327601", "GSM2327603",
              "GSM2327604", "GSM2327606", "GSM2327607",
              "GSM2327610")

# Assign group based on GSM number in filename
group <- ifelse(grepl(paste(AR_gsm, collapse="|"),
                       samples), "AR", "NAR")

# Convert to factor
group <- factor(group,
                levels = c("NAR", "AR"))

# Assign to eset
pData(eset)$group <- group

# Verify
cat("\nGroup assignments:\n")
print(pData(eset))

# Verify counts
cat("\nAR samples:", sum(group == "AR"), "\n")    # should be 10
cat("NAR samples:", sum(group == "NAR"), "\n")    # should be 16

# =============================================================================
# SECTION 5: Quality Control
# =============================================================================

# --- 5.1: Boxplot before and after normalization ---
par(mfrow = c(1, 2))

boxplot(raw_data,
        main = "Before Normalization",
        col = "lightblue",
        las = 2,
        cex.axis = 0.5)

boxplot(exprs_matrix,
        main = "After Normalization",
        col = "lightgreen",
        las = 2,
        cex.axis = 0.5)

par(mfrow = c(1, 1))

# --- 5.2: PCA Plot ---
pca <- prcomp(t(exprs_matrix))

# Calculate variance explained
var_explained <- summary(pca)$importance[2, 1:2] * 100

plot(pca$x[, 1], pca$x[, 2],
     col = c(rep("red", 7), rep("blue", 7)),
     pch = 19,
     main = "PCA Plot - GSE87301",
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"))

legend("topright",
       legend = c("AR", "NAR"),
       col = c("red", "blue"),
       pch = 19)

# --- 5.3: Sample Correlation Heatmap ---
cor_matrix <- cor(exprs_matrix)

annotation_col <- data.frame(group = group,
                             row.names = sampleNames(eset))

pheatmap(cor_matrix,
         main = "Sample Correlation Heatmap",
         annotation_col = annotation_col,
         annotation_colors = list(group = c(AR = "red", NAR = "blue")),
         show_rownames = FALSE)

# =============================================================================
# SECTION 6: Probe to Gene Symbol Mapping
# =============================================================================

# GPL570 uses hgu133plus2.db (NOT hgu133a.db)
library(hgu133plus2.db)

probe_ids <- rownames(exprs_matrix)

gene_symbols <- mapIds(hgu133plus2.db,    # ← different from GSE5967
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

cat("Total probes:", length(gene_symbols), "\n")
cat("Unmapped probes:", sum(is.na(gene_symbols)), "\n")

# Remove NA
keep <- !is.na(gene_symbols)
exprs_clean <- exprs_matrix[keep, ]
genes_clean <- gene_symbols[keep]

# Keep highest expressed probe per gene
mean_expr <- rowMeans(exprs_clean)
ord <- order(mean_expr, decreasing = TRUE)
exprs_clean <- exprs_clean[ord, ]
genes_clean <- genes_clean[ord]

unique_idx <- !duplicated(genes_clean)
exprs_final <- exprs_clean[unique_idx, ]
rownames(exprs_final) <- genes_clean[unique_idx]

cat("Final dimensions:", dim(exprs_final), "\n")
cat("Unique genes:", nrow(exprs_final), "\n")

# Verify signature genes are present
signature_genes <- c("ALAS2", "HBD", "EPB42", "FECH")
cat("\nVerifying signature genes present:\n")
for (gene in signature_genes) {
    found <- grep(paste0("^", gene, "$"),
                  rownames(exprs_final),
                  value = TRUE)
    cat(gene, ":", ifelse(length(found) > 0, "FOUND ✓", "NOT FOUND ✗"), "\n")
}

# =============================================================================
# SECTION 7: Differential Expression Analysis (limma)
# =============================================================================

# Design matrix
design_matrix <- model.matrix(~group)
colnames(design_matrix) <- c("Intercept", "AR_vs_NAR")

# Verify dimensions match
cat("\nexprs_final columns:", ncol(exprs_final), "\n")  # 26
cat("design_matrix rows:", nrow(design_matrix), "\n")   # 26

# Fit model
fit  <- lmFit(exprs_final, design_matrix)
fit2 <- eBayes(fit)

# Get all results
DEG_results <- topTable(fit2,
                        coef = "AR_vs_NAR",
                        number = Inf,
                        adjust.method = "BH",
                        sort.by = "P")

cat("Total genes tested:", nrow(DEG_results), "\n")

# Filter DEGs
DEGs <- DEG_results[DEG_results$P.Value < 0.05 &
                    abs(DEG_results$logFC) > 1, ]

at("\n--- DEG Summary ---\n")
cat("Total DEGs:", nrow(DEGs), "\n")
cat("Upregulated:", nrow(DEGs[DEGs$logFC > 0,]), "\n")
cat("Downregulated:", nrow(DEGs[DEGs$logFC < 0,]), "\n")

# Relaxed criteria for small dataset
DEGs_relaxed <- DEG_results[DEG_results$P.Value < 0.05, ]
cat("\nDEGs with p<0.05 only:", nrow(DEGs_relaxed), "\n")

# =============================================================================
# SECTION 8: Validate Signature Genes
# =============================================================================

# Check signature genes
cat("\n--- Signature Gene Results ---\n")
sig_results <- DEG_results[rownames(DEG_results) %in%
                            signature_genes, ]
print(sig_results)

cat("\nInterpretation:\n")
for (gene in signature_genes) {
    if (gene %in% rownames(sig_results)) {
        lfc <- sig_results[gene, "logFC"]
        pval <- sig_results[gene, "P.Value"]
        direction <- ifelse(lfc > 0, "UPREGULATED", "DOWNREGULATED")
        cat(gene, ":", direction,
            "| logFC =", round(lfc, 3),
            "| P.Value =", round(pval, 4), "\n")
    }
}

# =============================================================================
# SECTION 9: Save Results
# =============================================================================

output_path <- "/home/keerthana/Documents/project/1_project/sample/sample_2"

# Save all DEG results
write.table(DEG_results,
            file = file.path(output_path, "DEG_results_GSE5967.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# Save filtered DEGs
write.table(DEGs_relaxed,
            file = file.path(output_path, "DEGs_filtered_GSE5967.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# Save signature gene results
write.table(sig_results,
            file = file.path(output_path, "signature_genes_GSE5967.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

cat("\nAll results saved to:", output_path, "\n")
cat("\nAnalysis complete!\n")

# =============================================================================
# SECTION 10: Volcano Plot
# =============================================================================

library(ggplot2)

# Add significance labels
DEG_results$significance <- "Not Significant"
DEG_results$significance[DEG_results$P.Value < 0.05 &
                         DEG_results$logFC > 1]  <- "Up"
DEG_results$significance[DEG_results$P.Value < 0.05 &
                         DEG_results$logFC < -1] <- "Down"

# Check counts
cat("\n--- Significance Summary ---\n")
print(table(DEG_results$significance))

# Add labels for signature genes only
DEG_results$label <- ""
DEG_results$label[rownames(DEG_results) %in%
                  signature_genes] <- rownames(DEG_results)[
                  rownames(DEG_results) %in% signature_genes]

# Create volcano plot
p_volcano <- ggplot(DEG_results,
                    aes(x = logFC,
                        y = -log10(P.Value),
                        color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up"   = "red",
                                  "Down" = "blue",
                                  "Not Significant" = "grey")) +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed",
               color = "black") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed",
               color = "black") +
    geom_text(aes(label = label),
              size = 3,
              vjust = -0.5,
              color = "black") +
    labs(title = "Volcano Plot - GSE5967",
         x = "log2 Fold Change",
         y = "-log10(P.Value)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

print(p_volcano)

# Save volcano plot
ggsave(file.path(output_path, "volcano_plot_GSE5967.png"),
       p_volcano,
       width = 8,
       height = 6,
       dpi = 300)

cat("Volcano plot saved!\n")


# =============================================================================
# SECTION 11: Heatmap of Top 50 DEGs
# =============================================================================

library(pheatmap)

# Get top 46 by p-value + 4 signature genes = 50 total
top46_genes    <- rownames(head(DEG_results, 46))
combined_genes <- unique(c(top46_genes, signature_genes))

# Extract expression data
heatmap_data <- exprs_final[combined_genes, ]
cat("\nHeatmap data dimensions:", dim(heatmap_data), "\n")

# Column annotation (AR vs NAR)
annotation_col <- data.frame(
    Group = group,
    row.names = colnames(exprs_final)
)

# Row annotation (highlight signature genes)
annotation_row <- data.frame(
    Type = ifelse(rownames(heatmap_data) %in% signature_genes,
                  "Signature", "DEG"),
    row.names = rownames(heatmap_data)
)

# Plot heatmap
p_heatmap <- pheatmap(heatmap_data,
         main = "Top 50 DEGs - GSE5967",
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = list(
             Group = c(AR = "red", NAR = "blue"),
             Type  = c(Signature = "orange", DEG = "grey")),
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 7,
         color = colorRampPalette(
             c("blue", "white", "red"))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete")

# Save heatmap
png(file.path(output_path, "heatmap_top50_GSE5967.png"),
    width = 900,
    height = 1100)

pheatmap(heatmap_data,
         main = "Top 50 DEGs - GSE5967",
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = list(
             Group = c(AR = "red", NAR = "blue"),
             Type  = c(Signature = "orange", DEG = "grey")),
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 7,
         color = colorRampPalette(
             c("blue", "white", "red"))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete")

dev.off()
cat("Heatmap saved!\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n", rep("=", 50), "\n", sep = "")
cat("GSE5967 ANALYSIS COMPLETE\n")
cat(rep("=", 50), "\n", sep = "")
cat("Results saved to:", output_path, "\n\n")
cat("Files generated:\n")
cat("  - DEG_results_GSE5967.txt\n")
cat("  - DEGs_filtered_GSE5967.txt\n")
cat("  - signature_genes_GSE5967.txt\n")
cat("  - volcano_plot_GSE5967.png\n")
cat("  - heatmap_top50_GSE5967.png\n")
