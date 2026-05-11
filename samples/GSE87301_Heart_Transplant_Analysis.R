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



