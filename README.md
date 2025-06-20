# Biomarker Discovery & Survival Analysis in Lung Adenocarcinoma (LUAD)

This project identifies and validates prognostic biomarkers in LUAD by integrating **transcriptomic (RNA-seq)** and **clinical** data from **TCGA**. It uses statistical and machine learning approaches including **DEG analysis**, **Cox regression**, **LASSO**, and **external validation with GEO datasets**. The goal is to pinpoint genes that significantly affect patient survival and may serve as prognostic biomarkers.

---

## Data Sources

- **Primary Dataset**: TCGA-LUAD (RNA-Seq STAR Counts + Clinical data)
  - Retrieved using the `TCGAbiolinks` R package
- **External Validation Datasets**:
  - GSE31210
  - GSE68465

---

## Methodology Overview

### 🔹 1. **Data Retrieval and Preprocessing**
- Downloaded LUAD RNA-seq counts and clinical data using `TCGAbiolinks`
- Cleaned and merged expression and clinical metadata
- Normalized count data using `edgeR` (TMM) and logCPM transformation

### 🔹 2. **Exploratory & Survival Analysis**
- Kaplan-Meier curves by clinical variables (e.g., stage, smoking)
- Univariate Cox regression on clinical covariates

### 🔹 3. **Differential Gene Expression Analysis**
- Using both `edgeR` (LRT) and `limma-voom`
- Volcano plots highlighting top 10 upregulated and downregulated genes
- Gene annotation using `biomaRt` (HGNC and Entrez IDs)

### 🔹 4. **Functional Enrichment**
- GO Biological Process and KEGG Pathway analysis using `clusterProfiler`
- Dotplots and barplots saved as images

### 🔹 5. **Survival Analysis**
- Univariate Cox regression on all genes
- Multivariate Cox regression adjusting for age, gender, and stage
- Forest plot for top genes
- Kaplan-Meier plots for top prognostic genes

### 🔹 6. **LASSO-Cox Regression**
- LASSO model to identify predictive genes
- 10-fold cross-validation with `glmnet`
- ROC analysis for each selected gene
- AUC scores saved and visualized

### 🔹 7. **External Validation (GEO Datasets)**
- Validated top LASSO and Cox genes using:
  - GSE31210
  - GSE68465
- Genes with **adj. P < 0.05** in external datasets flagged as “Validated”

### 🔹 8. **Consensus Biomarker Identification**
- Overlap analysis between:
  - LASSO-selected genes
  - Univariate significant genes
  - Multivariate significant genes
- Final list exported for downstream enrichment

---

## Outputs & Visualizations

- PCA plots
- Volcano plots (`edgeR`, `limma-voom`)
- GO/KEGG enrichment dotplots & barplots
- Kaplan-Meier survival plots (top 10 genes)
- Forest plot (Cox regression)
- ROC plots (LASSO-selected genes)
- CSVs:
  - DEGs (annotated)
  - Cox regression results (uni/multi)
  - LASSO-selected genes (with AUC)
  - GEO validation results
  - Overlapping gene tables

---

## Tools & Packages

- **Data Access**: `TCGAbiolinks`
- **DEG Analysis**: `edgeR`, `limma`, `DESeq2`
- **Survival**: `survival`, `survminer`, `forestplot`, `glmnet`, `survivalROC`
- **Annotation**: `biomaRt`
- **Enrichment**: `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`
- **Visualization**: `ggplot2`, `ggrepel`

---

## Available Data Files

The following curated result files are included in the `data/` folder:

| File | Description |
|------|-------------|
| `LASSO_Genes_ROC_AUC_Annotated.csv` | LASSO-selected genes with Ensembl IDs, HGNC symbols, and time-dependent ROC AUC values |
| `Validated_Unique_Genes.csv` | Unique HGNC symbols from externally validated univariate + multivariate gene sets, formatted for DAVID enrichment analysis |

> **Note:** Raw data (e.g., TCGA, GEO) and intermediate outputs are excluded. You can reproduce all results by running the provided R script.

---

## Visualizations

The `photo/` folder contains all key result visualizations:

| Folder | Description |
|--------|-------------|
| `KM_Plots/` | Kaplan-Meier survival curves for top prognostic genes |
| `LASSO_ROC_Plots/` | Time-dependent ROC curves for LASSO-selected genes |
| `Photo/` | Contains enrichment plots, volcano plots, forest plots, etc. |

These plots help validate and interpret the statistical significance of selected biomarkers.

---

## Citation & Credit

- Data retrieved from [TCGA](https://portal.gdc.cancer.gov/) and [GEO](https://www.ncbi.nlm.nih.gov/geo/)
- This project is part of an MSc Bioinformatics dissertation focused on cancer genomics and survival biomarker discovery.
