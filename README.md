# Bioinformatics-Analysis-of-Cardiac-Allograft-Rejection

This repository documents my attempt to reproduce the bioinformatics analysis from the paper “Bioinformatics analysis characterizes immune infiltration landscape and identifies potential blood biomarkers for heart transplantation.”

The project includes processing of GEO datasets (GSE87301, GSE33970, GSE5967), differential gene expression analysis, functional enrichment (GO/KEGG), machine learning approaches (LASSO and Random Forest), immune cell infiltration analysis (ssGSEA), and validation of key biomarkers.

## Pipeline
```
CEL files (.cel.gz)
      ↓
read.affybatch()        # Load raw microarray data
      ↓
gcrma()                 # Background correction + normalization
      ↓
Quality Control         # Boxplot + PCA + Correlation heatmap
      ↓
mapIds(hgu133a.db)      # Probe ID → Gene symbol mapping
      ↓
lmFit() + eBayes()      # Differential expression (limma)
      ↓
topTable()              # Extract DEG results
      ↓
Signature gene          # Validate ALAS2, HBD, EPB42, FECH
validation
```
### R packages
```r
# Bioconductor
BiocManager::install(c("affy",
                       "gcrma",
                       "hgu133acdf",
                       "hgu133a.db",
                       "AnnotationDbi",
                       "limma"))

# CRAN
install.packages(c("pheatmap"))
```
### System dependencies (Linux)
```bash
sudo apt-get install -y \
    libcairo2-dev \
    libxml2-dev \
    libfontconfig1-dev
```

---

## Output Files
| File | Description |
|---|---|
| `DEG_results_GSE5967.txt` | All gene results from limma |
| `DEGs_filtered_GSE5967.txt` | Filtered DEGs (p < 0.05) |
| `signature_genes_GSE5967.txt` | Results for 4 signature genes |

## Reference
Wang, Y., & Peng, X. (2024). Bioinformatics analysis characterizes immune infiltration landscape and identifies potential blood biomarkers for heart transplantation. Transplant Immunology, 84, 102036.

## Author
Keerthana
