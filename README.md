# Bioinformatics-Analysis-of-Cardiac-Allograft-Rejection

Replication of microarray analysis from: "Bioinformatics analysis characterizes immune infiltration landscape and identifies potential blood biomarkers for heart transplantation" Wang &amp; Peng, Transplant Immunology 84 (2024) 102036

This project validates four signature genes (**ALAS2, HBD, EPB42, FECH**) as potential peripheral blood biomarkers for cardiac allograft rejection (AR) diagnosis using the GSE5967 validation dataset.

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
Wang Y, Peng X. Bioinformatics analysis characterizes immune infiltration
landscape and identifies potential blood biomarkers for heart transplantation.
*Transplant Immunology*. 2024;84:102036.
doi: 10.1016/j.trim.2024.102036

## Author
Keerthana
