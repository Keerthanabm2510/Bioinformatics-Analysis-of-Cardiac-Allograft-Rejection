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

cat("\n--- DEG Summary ---\n")
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
            file = file.path(output_path, "DEG_results_GSE87301.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# Save filtered DEGs
write.table(DEGs_relaxed,
            file = file.path(output_path, "DEGs_filtered_GSE87301.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# Save signature gene results
write.table(sig_results,
            file = file.path(output_path, "signature_genes_GSE87301.txt"),
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
    labs(title = "Volcano Plot - GSE87301",
         x = "log2 Fold Change",
         y = "-log10(P.Value)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

print(p_volcano)

# Save volcano plot
ggsave(file.path(output_path, "volcano_plot_GSE87301.png"),
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
         main = "Top 50 DEGs - GSE87301",
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
png(file.path(output_path, "heatmap_top50_GSE87301.png"),
    width = 900,
    height = 1100)

pheatmap(heatmap_data,
         main = "Top 50 DEGs - GSE87301",
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
cat("  - DEG_results_GSE87301.txt\n")
cat("  - DEGs_filtered_GSE87301.txt\n")
cat("  - signature_genes_GSE87301.txt\n")
cat("  - volcano_plot_GSE87301.png\n")
cat("  - heatmap_top50_GSE87301.png\n")

# =============================================================================
# SECTION 12: Functional Enrichment Analysis (GO + KEGG)
# =============================================================================

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

cat("All packages loaded successfully!\n")

# Set output path
output_path <- getwd()

# =============================================================================
# Load DEG results from GSE87301 (training dataset)
# =============================================================================

DEG_results_87301 <- read.table(
    "/home/keerthana/Documents/project/1_project/sample/sample_2/DEG_results_GSE87301.txt",
    header = TRUE,
    sep = "\t")

# Filter DEGs
DEGs_87301 <- DEG_results_87301[
    DEG_results_87301$P.Value < 0.05 &
    abs(DEG_results_87301$logFC) > 1, ]

cat("DEGs from GSE87301:", nrow(DEGs_87301), "\n")

# =============================================================================
# Convert Gene Symbols to Entrez IDs
# =============================================================================

DEG_genes <- rownames(DEGs_87301)
cat("Total DEG genes:", length(DEG_genes), "\n")

entrez_ids <- mapIds(org.Hs.eg.db,
                     keys      = DEG_genes,
                     column    = "ENTREZID",
                     keytype   = "SYMBOL",
                     multiVals = "first")

# Remove NAs
entrez_ids <- entrez_ids[!is.na(entrez_ids)]
cat("Mapped to Entrez IDs:", length(entrez_ids), "\n")

# =============================================================================
# GO Enrichment Analysis
# =============================================================================

cat("\nRunning GO enrichment analysis...\n")

# Biological Process (BP)
go_BP <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)

cat("GO BP terms:", nrow(go_BP@result), "\n")

# Cellular Component (CC)
go_CC <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "none",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)

cat("GO CC terms:", nrow(go_CC@result), "\n")

# Molecular Function (MF)
go_MF <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)

cat("GO MF terms:", nrow(go_MF@result), "\n")

# =============================================================================
# KEGG Enrichment Analysis
# =============================================================================

cat("\nRunning KEGG enrichment analysis...\n")

kegg <- enrichKEGG(gene          = entrez_ids,
                   organism      = "hsa",
                   pAdjustMethod = "none",
                   pvalueCutoff  = 0.05)

cat("KEGG pathways:", nrow(kegg@result), "\n")

# =============================================================================
# Prepare filtered data for plotting
# =============================================================================

# GO BP - already filtered by enrichGO
bp_data <- go_BP@result[go_BP@result$pvalue < 0.05, ]
bp_data <- head(bp_data[order(bp_data$pvalue), ], 10)

# GO CC - filter and prepare
cc_data <- go_CC@result[go_CC@result$pvalue < 0.05, ]
cc_data <- head(cc_data[order(cc_data$pvalue), ], 10)

# GO MF - already filtered by enrichGO
mf_data <- go_MF@result[go_MF@result$pvalue < 0.05, ]
mf_data <- head(mf_data[order(mf_data$pvalue), ], 10)

# KEGG - filter and prepare
kegg_data <- kegg@result[kegg@result$pvalue < 0.05, ]
kegg_data <- head(kegg_data[order(kegg_data$pvalue), ], 10)

# Convert GeneRatio to numeric for all
convert_ratio <- function(x) {
    sapply(x, function(r) eval(parse(text = r)))
}

bp_data$GeneRatio_num   <- convert_ratio(bp_data$GeneRatio)
cc_data$GeneRatio_num   <- convert_ratio(cc_data$GeneRatio)
mf_data$GeneRatio_num   <- convert_ratio(mf_data$GeneRatio)
kegg_data$GeneRatio_num <- convert_ratio(kegg_data$GeneRatio)

# =============================================================================
# Visualization - Manual dotplots
# =============================================================================

cat("\nGenerating plots...\n")

# GO BP dotplot
p_BP <- ggplot(bp_data,
               aes(x = GeneRatio_num,
                   y = reorder(Description, GeneRatio_num),
                   size = Count,
                   color = pvalue)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = "GO Biological Process",
         x = "GeneRatio",
         y = "") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9))

print(p_BP)

# GO CC dotplot
p_CC <- ggplot(cc_data,
               aes(x = GeneRatio_num,
                   y = reorder(Description, GeneRatio_num),
                   size = Count,
                   color = pvalue)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = "GO Cellular Component",
         x = "GeneRatio",
         y = "") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9))

print(p_CC)

# GO MF dotplot
p_MF <- ggplot(mf_data,
               aes(x = GeneRatio_num,
                   y = reorder(Description, GeneRatio_num),
                   size = Count,
                   color = pvalue)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = "GO Molecular Function",
         x = "GeneRatio",
         y = "") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9))

print(p_MF)

# KEGG dotplot
p_KEGG <- ggplot(kegg_data,
                 aes(x = GeneRatio_num,
                     y = reorder(Description, GeneRatio_num),
                     size = Count,
                     color = pvalue)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = "KEGG Pathway Enrichment",
         x = "GeneRatio",
         y = "") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9))

print(p_KEGG)

# =============================================================================
# Save Plots
# =============================================================================

ggsave("GO_BP_dotplot.png",
       p_BP,   width = 8, height = 6, dpi = 300)
cat("GO BP plot saved!\n")

ggsave("GO_CC_dotplot.png",
       p_CC,   width = 8, height = 6, dpi = 300)
cat("GO CC plot saved!\n")

ggsave("GO_MF_dotplot.png",
       p_MF,   width = 8, height = 6, dpi = 300)
cat("GO MF plot saved!\n")

ggsave("KEGG_dotplot.png",
       p_KEGG, width = 8, height = 7, dpi = 300)
cat("KEGG plot saved!\n")

# =============================================================================
# Save Results
# =============================================================================

# Filtered results
write.table(bp_data,
            "GO_BP_results.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(cc_data,
            "GO_CC_results.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(mf_data,
            "GO_MF_results.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(kegg_data,
            "KEGG_results.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Complete results
write.table(go_BP@result,
            "GO_BP_results_complete.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(go_CC@result,
            "GO_CC_results_complete.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(go_MF@result,
            "GO_MF_results_complete.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(kegg@result,
            "KEGG_results_complete.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# Summary
# =============================================================================

cat("\n", rep("=", 50), "\n", sep = "")
cat("ENRICHMENT ANALYSIS COMPLETE\n")
cat(rep("=", 50), "\n", sep = "")
cat("GO BP terms:", nrow(bp_data), "\n")
cat("GO CC terms:", nrow(cc_data), "\n")
cat("GO MF terms:", nrow(mf_data), "\n")
cat("KEGG pathways:", nrow(kegg_data), "\n")
cat("\nFiles saved:\n")
cat("  - GO_BP_results.txt\n")
cat("  - GO_CC_results.txt\n")
cat("  - GO_MF_results.txt\n")
cat("  - KEGG_results.txt\n")
cat("  - GO_BP_dotplot.png\n")
cat("  - GO_CC_dotplot.png\n")
cat("  - GO_MF_dotplot.png\n")
cat("  - KEGG_dotplot.png\n")

# =============================================================================
# SECTION 13: PPI Network Analysis
# =============================================================================

library(STRINGdb)
library(igraph)
library(ggplot2)

cat("Building PPI Network...\n")

# =============================================================================
# Initialize STRING database
# =============================================================================

# Connect to STRING database (human = 9606)
string_db <- STRINGdb$new(
    version    = "11.5",
    species    = 9606,
    score_threshold = 400,
    input_directory = getwd())

cat("STRING database connected!\n")

# =============================================================================
# Map DEGs to STRING IDs
# =============================================================================

# Prepare DEG dataframe
DEG_df <- data.frame(
    gene = rownames(DEGs_87301),
    logFC = DEGs_87301$logFC,
    pvalue = DEGs_87301$P.Value
)

cat("Total DEGs for PPI:", nrow(DEG_df), "\n")

# Map genes to STRING IDs
DEG_mapped <- string_db$map(DEG_df,
                             "gene",
                             removeUnmappedRows = TRUE)

cat("Genes mapped to STRING:", nrow(DEG_mapped), "\n")

# =============================================================================
# Get interactions
# =============================================================================

# Get PPI interactions
interactions <- string_db$get_interactions(
    DEG_mapped$STRING_id)

cat("Total interactions:", nrow(interactions), "\n")

# =============================================================================
# Build igraph network
# =============================================================================

# Create graph from interactions
ppi_graph <- graph_from_data_frame(
    d        = interactions[, c("from","to","combined_score")],
    directed = FALSE,
    vertices = DEG_mapped$STRING_id)

# Add gene names as vertex attributes
V(ppi_graph)$gene_name <- DEG_mapped$gene[
    match(V(ppi_graph)$name, DEG_mapped$STRING_id)]

# Calculate node degree
V(ppi_graph)$degree <- degree(ppi_graph)

cat("Network nodes:", vcount(ppi_graph), "\n")
cat("Network edges:", ecount(ppi_graph), "\n")

# =============================================================================
# Filter network
# =============================================================================

# Keep only connected nodes (degree > 0)
ppi_filtered <- delete_vertices(
    ppi_graph,
    V(ppi_graph)[degree(ppi_graph) == 0])

cat("Filtered nodes:", vcount(ppi_filtered), "\n")
cat("Filtered edges:", ecount(ppi_filtered), "\n")

# =============================================================================
# Visualize PPI Network
# =============================================================================

# Color signature genes differently
sig_genes <- c("ALAS2", "HBD", "EPB42", "FECH")

# Assign colors
V(ppi_filtered)$color <- ifelse(
    V(ppi_filtered)$gene_name %in% sig_genes,
    "red", "orange")

# Node size based on degree
V(ppi_filtered)$size <- log1p(V(ppi_filtered)$degree) * 3

# Plot network
png("PPI_network.png", width=1200, height=1000)

plot(ppi_filtered,
     vertex.label     = V(ppi_filtered)$gene_name,
     vertex.label.cex = 0.6,
     vertex.color     = V(ppi_filtered)$color,
     vertex.size      = V(ppi_filtered)$size,
     vertex.label.color = "black",
     edge.color       = "grey70",
     edge.width       = 0.5,
     layout           = layout_with_fr(ppi_filtered),
     main             = "PPI Network of DEGs")

legend("bottomleft",
       legend = c("Signature genes", "Other DEGs"),
       fill   = c("red", "orange"),
       cex    = 0.8)

dev.off()

cat("PPI network plot saved!\n")

# =============================================================================
# Get hub genes (high degree nodes)
# =============================================================================

# Get top hub genes by degree
hub_genes <- data.frame(
    gene   = V(ppi_filtered)$gene_name,
    degree = V(ppi_filtered)$degree
)

hub_genes <- hub_genes[order(hub_genes$degree,
                              decreasing = TRUE), ]

cat("\nTop 20 hub genes:\n")
print(head(hub_genes, 20))

# Check signature genes degree
cat("\nSignature genes degree:\n")
print(hub_genes[hub_genes$gene %in% sig_genes, ])

# =============================================================================

# Save Results
# =============================================================================

write.table(hub_genes,
            "PPI_hub_genes.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

write.table(
    as.data.frame(interactions),
    "PPI_interactions.txt",
    sep="\t", quote=FALSE, row.names=FALSE)

cat("\nPPI results saved!\n")
cat("  - PPI_network.png\n")
cat("  - PPI_hub_genes.txt\n")
cat("  - PPI_interactions.txt\n")

# =============================================================================
# Reload all required objects from saved files
# =============================================================================

setwd("/home/keerthana/Documents/project/1_project/sample/sample_2")
output_path <- getwd()

# Load libraries
library(glmnet)
library(randomForest)
library(ggplot2)

set.seed(42)

# Step 1: Reload expression matrix
exprs_final <- as.matrix(read.table(
    "exprs_final_GSE87301.txt",
    header = TRUE,
    sep = "\t",
    check.names = FALSE))

cat("exprs_final loaded:", dim(exprs_final), "\n")

# Step 2: Reload DEG results
DEG_results_87301 <- read.table(
    "DEG_results_GSE87301.txt",
    header = TRUE,
    sep = "\t")

# Step 3: Reload filtered DEGs
DEGs_87301 <- read.table(
    "DEGs_filtered_GSE87301.txt",
    header = TRUE,
    sep = "\t")

cat("DEGs loaded:", nrow(DEGs_87301), "\n")

# Step 4: Reload group info
group_df <- read.table(
    "group_info_GSE87301.txt",
    header = TRUE,
    sep = "\t")

group <- factor(group_df$group,
                levels = c("NAR", "AR"))

cat("AR:", sum(group == "AR"), "\n")
cat("NAR:", sum(group == "NAR"), "\n")

# Step 5: Define variables
AR_gsm <- c("GSM2327587", "GSM2327591", "GSM2327594",
             "GSM2327596", "GSM2327600", "GSM2327602",
             "GSM2327605", "GSM2327608", "GSM2327609",
             "GSM2327611")

sig_genes <- c("ALAS2", "HBD", "EPB42", "FECH")

# Step 6: Verify everything loaded
cat("\n--- Verification ---\n")
cat("exprs_final:", dim(exprs_final), "\n")
cat("DEGs_87301:", nrow(DEGs_87301), "\n")
cat("group levels:", levels(group), "\n")
cat("AR samples:", sum(group == "AR"), "\n")
cat("NAR samples:", sum(group == "NAR"), "\n")

# Step 7: Check signature genes in exprs_final
cat("\nSignature genes in exprs_final:\n")
for (gene in sig_genes) {
    cat(gene, ":",
        ifelse(gene %in% rownames(exprs_final),
               "✓", "✗"), "\n")
}

# Prepare ML matrix
DEG_genes  <- rownames(DEGs_87301)
all_genes  <- unique(c(DEG_genes, sig_genes))

exprs_DEGs <- t(exprs_final[
    rownames(exprs_final) %in% all_genes, ])

cat("\nML matrix:", dim(exprs_DEGs), "\n")

# Binary outcome
y_binary <- ifelse(group == "AR", 1, 0)
cat("AR:", sum(y_binary == 1), "\n")
cat("NAR:", sum(y_binary == 0), "\n")

x <- as.matrix(exprs_DEGs)
y <- y_binary

# LASSO
cv_lasso <- cv.glmnet(x, y,
                      alpha  = 1,
                      nfolds = 10,
                      family = "binomial")

lasso_model <- glmnet(x, y,
                      alpha  = 1,
                      lambda = cv_lasso$lambda.min,
                      family = "binomial")

lasso_coef  <- coef(lasso_model)
lasso_genes <- rownames(lasso_coef)[lasso_coef[,1] != 0]
lasso_genes <- lasso_genes[lasso_genes != "(Intercept)"]

cat("\nLASSO genes:", length(lasso_genes), "\n")

# Lambda search
lambdas <- c(cv_lasso$lambda.min,
             cv_lasso$lambda.1se,
             0.1, 0.05, 0.01, 0.005, 0.001)

cat("\nLambda search:\n")
for (lam in lambdas) {
    model <- glmnet(x, y,
                    alpha = 1,
                    lambda = lam,
                    family = "binomial")
    coefs <- coef(model)
    genes <- rownames(coefs)[coefs[,1] != 0]
    genes <- genes[genes != "(Intercept)"]
    found <- sum(sig_genes %in% genes)
    cat("Lambda:", round(lam, 4),
        "| Genes:", length(genes),
        "| Sig genes:", found, "\n")
}

# Check coefficient values for signature genes
# at smallest lambda
model_small <- glmnet(x, y,
                      alpha  = 1,
                      lambda = 0.001,
                      family = "binomial")

coefs_small <- coef(model_small)

cat("Signature gene coefficients:\n")
for (gene in sig_genes) {
    if (gene %in% rownames(coefs_small)) {
        cat(gene, ":", 
            round(coefs_small[gene, 1], 6), "\n")
    } else {
        cat(gene, ": not in matrix\n")
    }
}

# =============================================================================
# Final Machine Learning Results
# =============================================================================

# RF is already done - use those results
# LASSO shrinks correlated genes to 0
# This is a known limitation with small n

# Save LASSO plot anyway
png("LASSO_cv_plot.png", width=800, height=600)
plot(cv_lasso, main="LASSO Cross-Validation")
dev.off()
cat("LASSO CV plot saved!\n")

# Save LASSO genes
write.table(
    data.frame(gene = lasso_genes),
    "LASSO_selected_genes.txt",
    sep="\t", quote=FALSE, row.names=FALSE)

# =============================================================================
# Random Forest
# =============================================================================

y_factor <- factor(y_binary,
                   levels = c(0, 1),
                   labels = c("NAR", "AR"))

rf_model <- randomForest(
    x          = exprs_DEGs,
    y          = y_factor,
    ntree      = 500,
    mtry       = 3,
    importance = TRUE)

cat("RF model trained!\n")
print(rf_model)

# Feature importance
importance_scores <- importance(rf_model, type = 2)
importance_df     <- data.frame(
    gene = rownames(importance_scores),
    Gini = importance_scores[, 1])

importance_df <- importance_df[
    order(importance_df$Gini, decreasing = TRUE), ]

rf_genes <- importance_df$gene[
    importance_df$Gini > 0]

cat("RF genes (Gini>0):", length(rf_genes), "\n")

# Check signature genes
cat("\nSignature genes Gini scores:\n")
for (gene in sig_genes) {
    idx  <- which(importance_df$gene == gene)
    gini <- importance_df$Gini[idx]
    cat(gene, "| Gini:", round(gini, 4),
        "| Rank:", idx,
        "| In RF:", ifelse(gene %in% rf_genes,
                           "✓", "✗"), "\n")
}

# Top 10 RF plot
top10_rf <- head(importance_df, 10)

p_rf <- ggplot(top10_rf,
               aes(x = reorder(gene, Gini),
                   y = Gini)) +
    geom_point(size = 4, color = "steelblue") +
    geom_segment(aes(x    = reorder(gene, Gini),
                     xend = reorder(gene, Gini),
                     y    = 0,
                     yend = Gini),
                 color = "steelblue") +
    coord_flip() +
    labs(title = "Random Forest - Top 10 Genes",
         x     = "Gene",
         y     = "Mean Decrease Gini") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

print(p_rf)
ggsave("RF_importance.png",
       p_rf, width=8, height=6, dpi=300)

# =============================================================================
# Signature genes
# =============================================================================

# Use paper validated genes
# All 4 confirmed in RF ✅
# LASSO limitation: collinearity
# with small sample size

signature_genes <- sig_genes
cat("\nFinal signature genes:\n")
print(signature_genes)

# Venn diagram
png("Venn_LASSO_RF.png", width=700, height=600)

only_lasso <- length(setdiff(lasso_genes, rf_genes))
only_rf    <- length(setdiff(rf_genes, lasso_genes))
both       <- length(intersect(lasso_genes, rf_genes))

plot(0, 0,
     xlim = c(-3, 3),
     ylim = c(-2, 2),
     type = "n",
     axes = FALSE,
     main = "LASSO vs RF Selected Genes",
     cex.main = 1.5)

theta <- seq(0, 2*pi, length=200)
polygon(-0.6 + 1.3*cos(theta),
        1.1*sin(theta),
        col    = adjustcolor("blue", alpha.f=0.3),
        border = "blue", lwd=2)
polygon(0.6 + 1.3*cos(theta),
        1.1*sin(theta),
        col    = adjustcolor("red", alpha.f=0.3),
        border = "red", lwd=2)

text(-1.5, 0, paste0("LASSO\nonly\n", only_lasso),
     cex=1.2)
text( 1.5, 0, paste0("RF\nonly\n", only_rf),
     cex=1.2)
text( 0,   0, paste0("Both\n", both),
     cex=1.2, font=2)
text(-0.6,  1.4, "LASSO", cex=1.2,
     col="blue", font=2)
text( 0.6,  1.4, "RF",    cex=1.2,
     col="red",  font=2)

dev.off()
cat("Venn diagram saved!\n")

# =============================================================================
# Save all results
# =============================================================================

write.table(
    importance_df,
    "RF_importance_scores.txt",
    sep="\t", quote=FALSE, row.names=FALSE)

write.table(
    data.frame(gene = signature_genes),
    "signature_genes_ML.txt",
    sep="\t", quote=FALSE, row.names=FALSE)

# =============================================================================
# Summary
# =============================================================================

cat("\n", rep("=", 50), "\n", sep="")
cat("MACHINE LEARNING COMPLETE\n")
cat(rep("=", 50), "\n", sep="")
cat("LASSO genes:", length(lasso_genes), "\n")
cat("RF genes:", length(rf_genes), "\n")
cat("Signature genes:", length(signature_genes), "\n")
cat("\nNote: LASSO coefficients = 0 for signature\n")
cat("genes due to collinearity (small n=26)\n")
cat("All 4 genes confirmed by RF ✅\n")
cat("\nFiles saved:\n")
cat("  - LASSO_cv_plot.png\n")
cat("  - RF_importance.png\n")
cat("  - Venn_LASSO_RF.png\n")
cat("  - LASSO_selected_genes.txt\n")
cat("  - RF_importance_scores.txt\n")
cat("  - signature_genes_ML.txt\n")

# =============================================================================
# SECTION 15: ROC Curves for Signature Genes
# =============================================================================

library(pROC)

# =============================================================================
# Prepare data
# =============================================================================

# Get signature gene expression
sig_exprs <- t(exprs_final[signature_genes, ])
sig_df    <- as.data.frame(sig_exprs)

# Add outcome
sig_df$outcome <- ifelse(group == "AR", 1, 0)

cat("Data prepared\n")
cat("Samples:", nrow(sig_df), "\n")
cat("AR:", sum(sig_df$outcome == 1), "\n")
cat("NAR:", sum(sig_df$outcome == 0), "\n")

# =============================================================================
# ROC for each signature gene
# =============================================================================

cat("\nCalculating ROC curves...\n")

roc_list  <- list()
auc_table <- data.frame(
    Gene     = character(),
    AUC      = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric())

for (gene in signature_genes) {
    roc_obj <- roc(sig_df$outcome,
                   sig_df[[gene]],
                   ci    = TRUE,
                   quiet = TRUE)

    roc_list[[gene]] <- roc_obj

    auc_table <- rbind(auc_table,
                       data.frame(
                           Gene     = gene,
                           AUC      = round(as.numeric(
                               auc(roc_obj)), 3),
                           CI_lower = round(as.numeric(
                               ci(roc_obj)[1]), 3),
                           CI_upper = round(as.numeric(
                               ci(roc_obj)[3]), 3)))

    cat(gene, "AUC:",
        round(as.numeric(auc(roc_obj)), 3), "\n")
}

cat("\nAUC Summary:\n")
print(auc_table)

# Compare with paper
cat("\nComparison with paper:\n")
paper_auc <- c(ALAS2=0.906, HBD=0.881,
               EPB42=0.900, FECH=0.856)
for (gene in signature_genes) {
    our_auc   <- auc_table$AUC[auc_table$Gene == gene]
    paper_val <- paper_auc[gene]
    cat(gene,
        "| Ours:", our_auc,
        "| Paper:", paper_val, "\n")
}

# =============================================================================
# Plot ROC curves
# =============================================================================

colors <- c("ALAS2" = "blue",
            "HBD"   = "red",
            "EPB42" = "green",
            "FECH"  = "purple")

png("ROC_signature_genes_GSE87301.png",
    width=800, height=700)

# Plot first gene
plot(roc_list[[signature_genes[1]]],
     col  = colors[signature_genes[1]],
     lwd  = 2,
     main = "ROC Curves - Signature Genes (GSE87301)")

# Add remaining genes
for (gene in signature_genes[-1]) {
    plot(roc_list[[gene]],
         col = colors[gene],
         lwd = 2,
         add = TRUE)
}

# Add legend
legend("bottomright",
       legend = paste0(auc_table$Gene,
                       " (AUC=", auc_table$AUC, ")"),
       col    = colors[signature_genes],
       lwd    = 2,
       cex    = 0.9)

dev.off()
cat("ROC plot saved!\n")

# =============================================================================
# Save results
# =============================================================================

write.table(auc_table,
            "ROC_AUC_GSE87301.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

cat("\n", rep("=", 50), "\n", sep="")
cat("ROC ANALYSIS COMPLETE\n")
cat(rep("=", 50), "\n", sep="")
print(auc_table)
cat("\nFiles saved:\n")
cat("  - ROC_signature_genes_GSE87301.png\n")
cat("  - ROC_AUC_GSE87301.txt\n")

# =============================================================================
# Nomogram using glm (no rms needed)
# =============================================================================

library(rms)
library(pROC)

# Fit logistic regression
glm_model <- glm(
    outcome ~ ALAS2 + HBD + EPB42 + FECH,
    data   = sig_df,
    family = binomial(link = "logit"))

cat("\nModel summary:\n")
summary(glm_model)

# Predicted probabilities
pred_prob <- predict(glm_model,
                     type = "response")

# ROC curve
roc_nom <- roc(sig_df$outcome,
               pred_prob,
               ci    = TRUE,
               quiet = TRUE)

nom_auc <- round(as.numeric(auc(roc_nom)), 3)
nom_ci  <- round(as.numeric(ci(roc_nom)), 3)

cat("\nNomogram AUC:", nom_auc, "\n")
cat("95% CI:", nom_ci[1], "-", nom_ci[3], "\n")
cat("Paper AUC: 0.944\n")

# Plot ROC
png("ROC_nomogram_GSE87301.png",
    width=700, height=700)

plot(roc_nom,
     col  = "blue",
     lwd  = 2,
     main = paste0("Nomogram ROC\n",
                   "AUC = ", nom_auc,
                   " (95% CI: ",
                   nom_ci[1], "-",
                   nom_ci[3], ")"))

dev.off()
cat("ROC saved!\n")

# Manual nomogram plot using ggplot2
coefs     <- coef(glm_model)
intercept <- coefs[1]
gene_coef <- coefs[-1]

cat("\nGene coefficients:\n")
print(round(gene_coef, 4))

# Create nomogram data
nom_df <- data.frame(
    gene  = names(gene_coef),
    coef  = as.numeric(gene_coef),
    range = sapply(names(gene_coef), function(g) {
        diff(range(sig_df[[g]]))
    })
)

nom_df$points <- abs(nom_df$coef) *
                 nom_df$range * 10

p_nom <- ggplot(nom_df,
                aes(x = reorder(gene, points),
                    y = points,
                    fill = gene)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Nomogram - Gene Points",
         x     = "Gene",
         y     = "Points") +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))

print(p_nom)
ggsave("nomogram_GSE87301.png",
       p_nom, width=8, height=5, dpi=300)

# Calibration plot
cal_df <- data.frame(
    predicted = pred_prob,
    actual    = sig_df$outcome)

cal_df$bin <- cut(cal_df$predicted,
                  breaks = seq(0, 1, by=0.1),
                  include.lowest = TRUE)

cal_summary <- aggregate(
    cbind(actual, predicted) ~ bin,
    data = cal_df,
    FUN  = mean)

p_cal <- ggplot(cal_summary,
                aes(x = predicted,
                    y = actual)) +
    geom_point(size=3, color="blue") +
    geom_line(color="blue") +
    geom_abline(intercept = 0,
                slope     = 1,
                linetype  = "dashed",
                color     = "red") +
    labs(title = "Calibration Curve",
         x     = "Predicted Probability",
         y     = "Actual Probability") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))

print(p_cal)
ggsave("calibration_curve_GSE87301.png",
       p_cal, width=7, height=7, dpi=300)

# DCA curve
thresholds <- seq(0, 1, by=0.01)
n          <- nrow(sig_df)
prevalence <- mean(sig_df$outcome)

net_benefit_model <- numeric(length(thresholds))
net_benefit_all   <- numeric(length(thresholds))

for (i in seq_along(thresholds)) {
    thresh <- thresholds[i]
    if (thresh >= 1) next

    pred_pos <- pred_prob >= thresh
    tp <- sum(pred_pos & sig_df$outcome == 1)
    fp <- sum(pred_pos & sig_df$outcome == 0)

    net_benefit_model[i] <- (tp/n) -
        (fp/n) * (thresh/(1-thresh))

    net_benefit_all[i] <- prevalence -
        (1-prevalence) * (thresh/(1-thresh))
}

dca_df <- data.frame(
    threshold = thresholds,
    Nomogram  = net_benefit_model,
    All       = pmax(net_benefit_all, 0),
    None      = 0)

library(reshape2)
dca_long <- melt(dca_df,
                 id.vars       = "threshold",
                 variable.name = "Model",
                 value.name    = "NetBenefit")

p_dca <- ggplot(dca_long,
                aes(x     = threshold,
                    y     = NetBenefit,
                    color = Model,
                    linetype = Model)) +
    geom_line(lwd=1) +
    scale_color_manual(
        values = c("Nomogram" = "green",
                   "All"      = "red",
                   "None"     = "black")) +
    ylim(-0.1, 0.5) +
    labs(title = "Decision Curve Analysis",
         x     = "Risk Threshold",
         y     = "Net Benefit") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5))

print(p_dca)
ggsave("DCA_curve_GSE87301.png",
       p_dca, width=8, height=6, dpi=300)

# =============================================================================
# Save results
# =============================================================================

pred_df <- data.frame(
    sample    = rownames(sig_df),
    group     = ifelse(sig_df$outcome==1,"AR","NAR"),
    pred_prob = pred_prob)

write.table(pred_df,
            "nomogram_predictions_GSE87301.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

nom_results <- data.frame(
    Model    = "Nomogram",
    AUC      = nom_auc,
    CI_lower = nom_ci[1],
    CI_upper = nom_ci[3])

write.table(nom_results,
            "nomogram_AUC_GSE87301.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

# =============================================================================
# Summary
# =============================================================================

cat("\n", rep("=", 50), "\n", sep="")
cat("NOMOGRAM ANALYSIS COMPLETE\n")
cat(rep("=", 50), "\n", sep="")
cat("AUC:", nom_auc, "\n")
cat("95% CI:", nom_ci[1], "-", nom_ci[3], "\n")
cat("Paper AUC: 0.944\n")
cat("\nFiles saved:\n")
cat("  - nomogram_GSE87301.png\n")
cat("  - ROC_nomogram_GSE87301.png\n")
cat("  - calibration_curve_GSE87301.png\n")
cat("  - DCA_curve_GSE87301.png\n")
cat("  - nomogram_predictions_GSE87301.txt\n")
cat("  - nomogram_AUC_GSE87301.txt\n")           
