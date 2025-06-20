#Installing Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks", force = TRUE)
BiocManager::install("SummarizedExperiment", force = TRUE)


library(TCGAbiolinks)
library(SummarizedExperiment) 
library(tidyverse)
library(survival)
library(survminer)
library(DESeq2)
library(edgeR)
library(tibble) 
library(readr)   
library(dplyr)
library(ggplot2)
library(limma)
library(biomaRt)
library(ggrepel)  
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(forestplot)
library(glmnet) 
library(survivalROC)


#Data Retrieving:TCGA-LUAD Dataset
#RNA-Seq Data Retrieving
# Step 1: Query LUAD RNA-seq (STAR-Counts)
query_exp <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Step 2: Download using GDC Data Transfer Tool
GDCdownload(query = query_exp, method = "client") 
# Step 3: Prepare the data as SummarizedExperiment object
data_exp <- GDCprepare(query = query_exp)
# Expression matrix
expression_matrix <- assay(data_exp)  # rows = genes, cols = samples
# Metadata (barcodes, sample types, etc.)
metadata <- colData(data_exp)

saveRDS(data_exp, file = "TCGA_LUAD_expression_data.rds")

# Save expression matrix as CSV
write.csv(as.data.frame(expression_matrix), "TCGA_LUAD_STAR_Counts.csv", row.names = TRUE)

# Convert to data frame, then flatten list columns
metadata_df <- as.data.frame(metadata)

# Keep only atomic (non-list) columns to write cleanly
metadata_flat <- metadata_df[, sapply(metadata_df, function(x) !is.list(x))]

# Save to CSV
write.csv(metadata_flat, "TCGA_LUAD_metadata.csv", row.names = TRUE)


#Clinical Data Retrieving 
# Get preprocessed clinical data directly from GDC
clinical_df <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

# View structure
str(clinical_df)

sapply(clinical_df, class) 

head(clinical_df$sites_of_involvement, 5)

# Remove the 'sites_of_involvement' column
clinical_df_clean <- clinical_df[, !sapply(clinical_df, is.list)]

# Now write to CSV
write.csv(clinical_df_clean, "TCGA_LUAD_clinical_data.csv", row.names = FALSE)



#SURVIVAL ANALYSIS - Using Clinical Data 
# Install if not already installed
install.packages(c("tidyverse", "survival", "survminer")) 

# Set your file path
clinical <- read_csv("TCGA_LUAD_clinical_data.csv")
colnames(clinical)
clinical <- clinical %>%
  mutate(
    days_to_death = as.numeric(days_to_death),
    days_to_last_follow_up = as.numeric(days_to_last_follow_up),
    survival_time = if_else(
      vital_status == "Dead" & !is.na(days_to_death), days_to_death, days_to_last_follow_up
    ),
    event = if_else(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(survival_time))
km_fit <- survfit(Surv(survival_time, event) ~ 1, data = clinical)
ggsurvplot(km_fit, data = clinical, risk.table = TRUE, pval = TRUE)

# Example: Survival by stage
fit <- survfit(Surv(survival_time, event) ~ ajcc_pathologic_stage, data = clinical)
class(fit)
# Should return "survfit"
ggsurvplot(
  fit,
  data = clinical,
  risk.table = TRUE,
  pval = TRUE,
  risk.table.height = 0.5
)

km_smoke <- survfit(Surv(survival_time, event) ~ tobacco_smoking_status, data = clinical)

ggsurvplot(
  km_smoke,
  data = clinical,
  risk.table = TRUE,
  pval = TRUE,
  risk.table.height = 0.35,        
  risk.table.y.text = FALSE,        
  legend.title = "Smoking Status",  
  legend.labs = levels(as.factor(clinical$tobacco_smoking_status)), 
  font.legend = 12                  
)


cox_model <- coxph(Surv(survival_time, event) ~ ajcc_pathologic_stage + age_at_diagnosis + gender + tobacco_smoking_status, data = clinical)
summary(cox_model)


# Extract coefficients table
cox_summary <- summary(cox_model)
coef_table <- as.data.frame(cox_summary$coefficients)

# Save as CSV
write.csv(coef_table, "cox_model_coefficients.csv", row.names = TRUE)



#NORMALISATION OF RAW COUNT DATA - RNA-Seq DATA

# Step 1: Read STAR count data
counts <- read.csv("TCGA_LUAD_STAR_Counts.csv", row.names = 1, check.names = FALSE)

# Step 2: Read metadata (barcode + sample_type: TP/NT)
metadata <- read.csv("TCGA_LUAD_metadata.csv")

# Step 3: Ensure barcodes match - get intersection
common_samples <- intersect(colnames(counts), metadata$barcode)

# Step 4: Filter both datasets to common samples
counts_filtered <- counts[, common_samples]
metadata_filtered <- metadata %>% filter(barcode %in% common_samples)

colnames(metadata)


# Step 5: Ensure same order of samples
metadata_filtered <- metadata_filtered[match(colnames(counts_filtered), metadata_filtered$barcode), ]

# Step 6: Create DGEList
group <- metadata_filtered$shortLetterCode  # TP or NT
dge <- DGEList(counts = counts_filtered, group = group)

# Step 7: Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Step 8: Normalize with TMM
dge <- calcNormFactors(dge)

dim(dge)

dge$samples

head(dge$counts)
dge

# Now `dge` is ready for PCA, DEG, or limma/voom or edgeR pipelines

# Calculate CPM (counts per million) normalized counts with TMM factors and log-transform
logCPM <- cpm(dge, log=TRUE, normalized.lib.sizes=TRUE)

# PCA on logCPM matrix (genes x samples)
pca <- prcomp(t(logCPM))  # transpose to samples x genes

# Create a data frame for plotting
pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     Group = dge$samples$group,
                     Sample = rownames(pca$x))

# Plot PCA
ggplot(pca_df, aes(x=PC1, y=PC2, color=Group)) +
  geom_point(size=3) +
  theme_minimal() +
  labs(title="PCA of samples", x="PC1", y="PC2")

#Differential Expression Analysis with edgeR

# Create design matrix for two groups (TP vs NT)
design <- model.matrix(~ group, data = dge$samples)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the negative binomial generalized log-linear model
fit <- glmFit(dge, design)

# Perform likelihood ratio test for the group coefficient (difference TP vs NT)
lrt <- glmLRT(fit, coef=2)  

# View top differentially expressed genes
topTags(lrt)

# Extract table of results
deg_results <- topTags(lrt, n = nrow(dge))$table

# You can filter significant DEGs by adjusted p-value, e.g.
sig_degs <- deg_results[deg_results$FDR < 0.05, ]

# See number of significant DEGs
cat("Number of DEGs with FDR < 0.05:", nrow(sig_degs), "\n") 

write.csv(deg_results, "edgeR_DEG_results_all.csv", row.names = TRUE)
write.csv(sig_degs, "edgeR_DEG_results_FDR_0.05.csv", row.names = TRUE)

#Volcano Plot for edgeR DEG Results

# Add a significance label
deg_results$significant <- "Not Sig"
deg_results$significant[deg_results$FDR < 0.05 & deg_results$logFC > 1] <- "Up"
deg_results$significant[deg_results$FDR < 0.05 & deg_results$logFC < -1] <- "Down"

# Plot
ggplot(deg_results, aes(x = logFC, y = -log10(FDR), color = significant)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot (edgeR)",
       x = "Log2 Fold Change",
       y = "-Log10 FDR") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")



# Differential Expression Analysis using limma-voom

# Calculate voom weights
v <- voom(dge, design)

# Fit linear model
fit <- lmFit(v, design)

# Empirical Bayes moderation
fit <- eBayes(fit)

# Get top table of DEGs
topTable(fit, coef=2)  # coef=2 corresponds to the group effect (TP vs NT)
voom_deg_results <- topTable(fit, coef=2, number=Inf)

write.csv(voom_deg_results, "limma_voom_DEG_results_all.csv", row.names = TRUE)

#Save only significant DEGs (adjusted p-value < 0.05)
voom_sig_degs <- voom_deg_results[voom_deg_results$adj.P.Val < 0.05, ]
write.csv(voom_sig_degs, "limma_voom_DEG_results_adjP_0.05.csv", row.names = TRUE)

#Volcano Plot for limma-voom Results
voom_deg_results$significant <- "Not Sig"
voom_deg_results$significant[voom_deg_results$adj.P.Val < 0.05 & voom_deg_results$logFC > 1] <- "Up"
voom_deg_results$significant[voom_deg_results$adj.P.Val < 0.05 & voom_deg_results$logFC < -1] <- "Down"

ggplot(voom_deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot (limma-voom)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")

#convert Ensembl IDs to gene symbols 
#Annotate Entrez IDs in addition to HGNC symbols

# Connect to Ensembl
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Extract clean Ensembl gene IDs from your DEGs (remove version number)
deg_results$ensembl_gene_id <- gsub("\\..*", "", rownames(deg_results))

# Query Ensembl to get HGNC symbols + Entrez Gene IDs
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = deg_results$ensembl_gene_id,
  mart = mart
)

# Merge annotation with DEG results
deg_results_annotated <- merge(deg_results, gene_map, by = "ensembl_gene_id")

# Reorder columns
deg_results_annotated <- deg_results_annotated %>%
  dplyr::select(hgnc_symbol, entrezgene_id, everything())

# Save full annotated DEGs with Entrez ID and HGNC symbol
write.csv(deg_results_annotated, "DEG_results_annotated_with_entrez.csv", row.names = FALSE)

#Save only significant DEGs (FDR < 0.05)
sig_degs <- deg_results_annotated %>% filter(FDR < 0.05)
write.csv(sig_degs, "Significant_DEGs_annotated_with_entrez.csv", row.names = FALSE)



# For example, using edgeR results:
deg_results$ensembl_gene_id <- gsub("\\..*", "", rownames(deg_results))

# Map HGNC gene symbols
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = deg_results$ensembl_gene_id,
  mart = mart
)
deg_annotated <- merge(deg_results, gene_map, by = "ensembl_gene_id")


deg_annotated$significant <- "Not Sig"
deg_annotated$significant[deg_annotated$FDR < 0.05 & deg_annotated$logFC > 1] <- "Up"
deg_annotated$significant[deg_annotated$FDR < 0.05 & deg_annotated$logFC < -1] <- "Down"


top_genes <- deg_annotated %>%
  filter(FDR < 0.05, hgnc_symbol != "") %>%
  arrange(FDR) %>%
  head(10)


ggplot(deg_annotated, aes(x = logFC, y = -log10(FDR), color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes,
                  aes(label = hgnc_symbol),
                  size = 3,
                  box.padding = 0.4,
                  max.overlaps = Inf) +
  theme_minimal() +
  labs(title = "Volcano Plot with Top 10 DEGs Labeled",
       x = "Log2 Fold Change",
       y = "-Log10 FDR")

install.packages("ggrepel")

# top 10 most upregulated and top 10 most downregulated DEGs (based on log2 fold change and significance) separately on a volcano plot

# Filter only DEGs that are significant (FDR < 0.05) and have HGNC symbol
sig_degs <- deg_annotated %>%
  filter(FDR < 0.05, hgnc_symbol != "")

top_up <- sig_degs %>%
  arrange(desc(logFC)) %>%
  head(10)

top_down <- sig_degs %>%
  arrange(logFC) %>%
  head(10)

# Combine for labeling
top_genes <- bind_rows(top_up, top_down)


ggplot(deg_annotated, aes(x = logFC, y = -log10(FDR), color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes,
                  aes(label = hgnc_symbol),
                  size = 3,
                  box.padding = 0.4,
                  max.overlaps = Inf) +
  theme_minimal() +
  labs(title = "Volcano Plot: Top 10 Up & Downregulated DEGs",
       x = "Log2 Fold Change",
       y = "-Log10 FDR")

ggsave("volcano_plot_top10_up_down_labeled.png", width = 7, height = 5, dpi = 300)

# Add a column indicating direction
top_genes <- top_genes %>%
  mutate(direction = ifelse(logFC > 0, "Up", "Down"))

# Save to CSV
write.csv(top_genes, "Top10_Up_Down_DEGs.csv", row.names = FALSE)


#Functional Enrichment Analysis

if (!requireNamespace("clusterProfiler", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "pathview"))


# Use previously annotated significant DEGs (FDR < 0.05)
sig_degs <- deg_results_annotated %>%
  filter(FDR < 0.05, !is.na(entrezgene_id))  # Ensure Entrez IDs exist

# Get unique Entrez IDs
entrez_ids <- unique(sig_degs$entrezgene_id)

#GO Enrichment Analysis (Biological Process)
ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE  # converts Entrez to gene symbols
)

#KEGG pathway enrichment
ekegg <- enrichKEGG(
  gene         = entrez_ids,
  organism     = 'hsa',     # 'hsa' for human
  pvalueCutoff = 0.05
)

# Convert Entrez IDs to symbols in KEGG results
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# View top GO terms
head(ego)

# View top KEGG pathways
head(ekegg)

# Save to CSV
write.csv(as.data.frame(ego), "GO_enrichment_results.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "KEGG_enrichment_results.csv", row.names = FALSE)

# Dotplot for GO
dotplot(ego, showCategory = 15, title = "GO Enrichment - Biological Process")

# Barplot for KEGG
barplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment")

# Save GO dotplot
png("GO_enrichment_dotplot.png", width = 2000, height = 1600, res = 300)
dotplot(ego, showCategory = 15, title = "GO Enrichment - Biological Process")
dev.off()

# Save KEGG barplot
png("KEGG_enrichment_barplot.png", width = 2000, height = 1600, res = 300)
barplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment")
dev.off()


#Survival Analysis - using merged RNA-Seq and Clinical Data

#Read the clinical data
clinical_df <- read.csv("TCGA_LUAD_clinical_data.csv", stringsAsFactors = FALSE)

#normalization
logCPM <- cpm(dge, log=TRUE, normalized.lib.sizes=TRUE)
summary(logCPM)

# Convert logCPM matrix to data frame
logCPM_df <- as.data.frame(logCPM)

# Transpose logCPM: now samples are rows
logCPM_t <- as.data.frame(t(logCPM))  
logCPM_t$patient_id <- substr(rownames(logCPM_t), 1, 12)

# Remove duplicated patient IDs (keep first occurrence)
logCPM_t_unique <- logCPM_t[!duplicated(logCPM_t$patient_id), ]

# Store patient_id as rownames
rownames(logCPM_t_unique) <- logCPM_t_unique$patient_id
logCPM_t_unique$patient_id <- NULL

# Transpose back to get genes × samples format
logCPM_unique <- t(logCPM_t_unique)

# Sanity check
dim(logCPM_unique)  


# Extract patient IDs (first 12 chars of barcode)
patient_ids <- substr(colnames(logCPM_unique), 1, 12)

# Update column names to patient IDs
colnames(logCPM_unique) <- patient_ids

# Ensure the patient ID column exists and is in same format
clinical_df$patient_id <- substr(clinical_df$bcr_patient_barcode, 1, 12)
# Keep only clinical data for which we have expression data
clinical_matched <- clinical_df %>% filter(patient_id %in% colnames(logCPM_unique))
length(unique(clinical_matched$patient_id)) 
ncol(logCPM_unique)

# Ensure same order as logCPM_unique
clinical_matched <- clinical_matched[match(colnames(logCPM_unique), clinical_matched$patient_id), ]

# Final check — should be TRUE
all(clinical_matched$patient_id == colnames(logCPM_unique))

colnames(clinical_matched)
# Create survival time: death if available, otherwise last follow-up
clinical_matched$time <- as.numeric(ifelse(
  !is.na(clinical_matched$days_to_death),
  clinical_matched$days_to_death,
  clinical_matched$days_to_last_follow_up
))

# Create event variable: 1 = died, 0 = censored (alive)
clinical_matched$event <- ifelse(
  clinical_matched$vital_status == "Dead", 1,
  ifelse(clinical_matched$vital_status == "Alive", 0, NA)
)
table(is.na(clinical_matched$time))
table(is.na(clinical_matched$event))

# Keep only patients with valid survival time
clinical_matched <- clinical_matched %>% filter(!is.na(time))

nrow(clinical_matched)  
table(clinical_matched$event)
str(clinical_matched[, c("patient_id", "time", "event")])
summary(clinical_matched$time)

expr_transposed <- as.data.frame(t(logCPM_unique))  # Now rows = samples
expr_transposed$patient_id <- rownames(expr_transposed)
full_merged <- merge(clinical_matched, expr_transposed, by = "patient_id")


# Check dimensions
dim(full_merged)  # Should be: [number of patients] x [clinical + gene columns]

# Check for missing values
sum(is.na(full_merged))  

# Check for missing values in survival-related columns only
colSums(is.na(full_merged[, c("patient_id", "time", "event")]))

# Check for any missing values
sum(is.na(expr_data))  

# Assuming expression data starts from column 101
expr_matrix <- full_merged[, 101:ncol(full_merged)]

# Total number of missing values
sum(is.na(expr_matrix))  
all(sapply(expr_matrix, is.numeric))
constant_genes <- sapply(expr_matrix, function(x) length(unique(x)) <= 1)
sum(constant_genes)  # Number of genes with no variation
summary(expr_matrix[, 1:5])  # Preview summary for first 5 genes

ncol(full_merged)

# Extract expression matrix and clinical data
expr_matrix <- full_merged[, 101:ncol(full_merged)]
clinical_data <- full_merged[, 1:100]  # Assuming clinical vars are before 101

# Initialize a data frame to store results
cox_results <- data.frame(
  gene = colnames(expr_matrix),
  HR = NA,
  HR_conf_low = NA,
  HR_conf_high = NA,
  p_value = NA
)

# Univariate Cox Regression Analysis
# Loop through each gene and run Cox regression
for (i in seq_along(expr_matrix)) {
  gene_expr <- expr_matrix[[i]]
  model <- try(coxph(Surv(clinical_data$time, clinical_data$event) ~ gene_expr), silent = TRUE)
  
  if (!inherits(model, "try-error")) {
    summary_model <- summary(model)
    cox_results$HR[i] <- summary_model$coefficients[1, "exp(coef)"]
    cox_results$HR_conf_low[i] <- summary_model$conf.int[1, "lower .95"]
    cox_results$HR_conf_high[i] <- summary_model$conf.int[1, "upper .95"]
    cox_results$p_value[i] <- summary_model$coefficients[1, "Pr(>|z|)"]
  }
}



# Filter significant results (e.g., p < 0.05)
significant_genes <- cox_results[cox_results$p_value < 0.05, ]
# Example: HR must be either > 1.5 or < 0.67 and p < 0.05
significant_genes <- cox_results[
  cox_results$p_value < 0.05 &
    (cox_results$HR > 1.5 | cox_results$HR < 0.67), ]


# Save results
write.csv(cox_results, "cox_regression_all_genes.csv", row.names = FALSE)
write.csv(significant_genes, "significant_prognostic_genes.csv", row.names = FALSE)

# View top results
head(significant_genes[order(significant_genes$p_value), ])

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract Ensembl IDs from your cox results
ensembl_ids <- gsub("\\..*", "", significant_genes$gene)  # Remove version numbers

# Map to gene symbols
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

# Merge gene symbols into your results
significant_genes$ensembl_gene_id <- gsub("\\..*", "", significant_genes$gene)
significant_genes <- merge(significant_genes, gene_map, by = "ensembl_gene_id", all.x = TRUE)

# Reorder columns if needed
significant_genes <- significant_genes[, c("gene", "hgnc_symbol", "HR", "HR_conf_low", "HR_conf_high", "p_value")]

write.csv(significant_genes, "significant_prognostic_genes.csv", row.names = FALSE)

# Install if not present
if (!require("forestplot")) install.packages("forestplot")

# Prepare top genes for plot
top_genes <- head(significant_genes[order(significant_genes$p_value), ], 10)

# Format table
tabletext <- cbind(
  c("Gene", top_genes$hgnc_symbol),
  c("HR", formatC(top_genes$HR, digits = 2, format = "f")),
  c("95% CI", paste0(
    formatC(top_genes$HR_conf_low, digits = 2),
    " - ",
    formatC(top_genes$HR_conf_high, digits = 2)))
)

# Forest plot
forestplot(
  labeltext = tabletext,
  mean = c(NA, top_genes$HR),
  lower = c(NA, top_genes$HR_conf_low),
  upper = c(NA, top_genes$HR_conf_high),
  zero = 1, boxsize = 0.2, lineheight = unit(0.5, "cm"),
  col = forestplot::fpColors(box = "royalblue", lines = "darkblue", zero = "gray50"),
  title = "Top Prognostic Genes (Cox HR)"
)


# Create a directory to save the plots
dir.create("KM_Plots", showWarnings = FALSE)

# Loop through top N genes 
top_genes <- head(significant_genes[order(significant_genes$p_value), ], 10)

for (i in 1:nrow(top_genes)) {
  
  gene_id <- top_genes$gene[i]
  gene_name <- ifelse(top_genes$hgnc_symbol[i] != "", top_genes$hgnc_symbol[i], gene_id)
  
  # Create a binary group: high vs low expression (median split)
  full_merged$group <- ifelse(full_merged[[gene_id]] >= median(full_merged[[gene_id]], na.rm = TRUE), "High", "Low")
  
  # Survival fit
  fit <- survfit(Surv(time, event) ~ group, data = full_merged)
  
  # Plot and save
  plot <- ggsurvplot(
    fit,
    data = full_merged,
    pval = TRUE,
    risk.table = TRUE,
    title = paste("KM Curve for", gene_name),
    legend.title = "Expression",
    legend.labs = c("Low", "High"),
    palette = c("blue", "red")
  )
  
  # Save to file
  ggsave(
    filename = paste0("KM_Plots/", i, "_", gene_name, "_KM_plot.png"),
    plot = plot$plot,
    width = 6, height = 6, dpi = 300
  )
}

str(full_merged[, 1:100])  

#Multivariate Cox Regression Analysis
# Store results
multi_cox_results <- data.frame()

# Loop through each gene column
for (gene in colnames(full_merged)[101:ncol(full_merged)]) {
  formula <- as.formula(paste0("Surv(time, event) ~ `", gene, "` + age_at_diagnosis + gender + ajcc_pathologic_stage"))
  try({
    fit <- coxph(formula, data = full_merged)
    s <- summary(fit)
    multi_cox_results <- rbind(multi_cox_results, data.frame(
      gene = gene,
      HR = s$coefficients[1, "exp(coef)"],
      pvalue = s$coefficients[1, "Pr(>|z|)"]
    ))
  }, silent = TRUE)
}

# Load biomaRt
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")

# Assume your dataframe is multi_cox_results and it has a column "gene" with Ensembl IDs + version
# Step 1: Remove version numbers
multi_cox_results$ensembl_id <- sub("\\.\\d+$", "", multi_cox_results$gene)

# Step 2: Connect to Ensembl BioMart
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Step 3: Fetch gene symbols
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique(multi_cox_results$ensembl_id),
  mart = mart
)

# Step 4: Merge annotations back to your data
annotated_results <- merge(multi_cox_results, gene_annotations,
                           by.x = "ensembl_id", by.y = "ensembl_gene_id",
                           all.x = TRUE)

# Reorder columns (e.g., gene symbol first)
annotated_results <- annotated_results[, c("hgnc_symbol", setdiff(names(annotated_results), "hgnc_symbol"))]

# Step 5: Save to CSV
write.csv(annotated_results, "multi_cox_results_with_gene_symbols.csv", row.names = FALSE)


# Step 1: Filter significant genes by p-value and hazard ratio
# multi_cox_results file contains columns: 'gene' (Ensembl IDs), 'p.value', and 'HR'

sig_genes <- subset(multi_cox_results, pvalue < 0.05 & (HR > 1.5 | HR < 0.67))

# Clean Ensembl IDs (remove version numbers like ".15")
sig_genes$ensembl_clean <- sub("\\..*", "", sig_genes$gene)

# Step 2: Annotate gene symbols using biomaRt
# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch annotations
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = sig_genes$ensembl_clean,
  mart = mart
)

# Step 3: Merge annotations with Cox results
annotated_sig_genes <- merge(
  sig_genes,
  gene_annotations,
  by.x = "ensembl_clean",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)

# Reorder columns (place hgnc_symbol next to gene)
annotated_sig_genes <- annotated_sig_genes[, c("gene", "hgnc_symbol", setdiff(names(annotated_sig_genes), c("gene", "hgnc_symbol")))]

# Step 4: Save to CSV
write.csv(annotated_sig_genes, "significant_multivariate_cox_genes_annotated.csv", row.names = FALSE)

cat("Annotated significant multivariate Cox results saved as 'significant_multivariate_cox_genes_annotated.csv'\n")


# LASSO Cox Regression - to predict top predictive genes
# Install & load required package
install.packages("glmnet")

# Using processed data
# expr_matrix: samples × genes 
# clinical_data: includes `time` and `event` 

# Combining expression and clinical data
combined_df <- cbind(clinical_data[, c("time", "event")], expr_matrix_clean)

# Filter out patients with time ≤ 0 or missing values
filtered_df <- combined_df %>%
  filter(time > 0 & !is.na(time) & !is.na(event))

# Separate clinical and expression data again
y <- Surv(filtered_df$time, filtered_df$event)
x <- as.matrix(filtered_df[, -(1:2)])


# Run LASSO with 10-fold cross-validation
set.seed(42)  # for reproducibility
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)

# Plot CV curve
png("LASSO_CV_Plot.png", width = 800, height = 600, res = 150)
plot(cvfit)
dev.off()

# Get selected genes at optimal lambda
lasso_coefs <- coef(cvfit, s = "lambda.min")
selected_lasso_genes <- rownames(lasso_coefs)[as.numeric(lasso_coefs) != 0]

# View and save results
selected_lasso_genes <- selected_lasso_genes[selected_lasso_genes != "(Intercept)"]
write.csv(selected_lasso_genes, "LASSO_Selected_Genes.csv", row.names = FALSE)

# Print how many genes selected
cat("Number of LASSO-selected prognostic genes:", length(selected_lasso_genes), "\n")

#Lasso-selected prognostic genes (Gene Annotation - Gene Symbols to Ensemble ids)
# Clean version numbers
lasso_genes <- read.csv("LASSO_Selected_Genes.csv", stringsAsFactors = FALSE)
lasso_genes_clean <- gsub("\\..*", "", lasso_genes$x)  # remove version suffix

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get HGNC symbols
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = lasso_genes_clean,
  mart = mart
)
# Merge
annotated_lasso_genes <- data.frame(ensembl_id = lasso_genes_clean)
annotated_lasso_genes <- merge(annotated_lasso_genes, gene_annotations,
                               by.x = "ensembl_id", by.y = "ensembl_gene_id", all.x = TRUE)
# Save
write.csv(annotated_lasso_genes, "LASSO_Selected_Genes_Annotated.csv", row.names = FALSE)


#ROC Curve for top LASSO genes

# Create folder to save ROC plots
dir.create("LASSO_ROC_Plots", showWarnings = FALSE)

# Clean Ensembl IDs (remove version numbers)
lasso_genes_clean <- gsub("\\..*", "", selected_lasso_genes)

# Map HGNC gene symbols
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = lasso_genes_clean,
  mart = mart
)

# Prepare AUC result dataframe
auc_results <- data.frame(
  Gene_Ensembl_ID_Version = selected_lasso_genes,
  Ensembl_ID = lasso_genes_clean,
  HGNC_Symbol = NA,
  AUC = NA,
  stringsAsFactors = FALSE
)

# Loop through genes
for (i in seq_along(selected_lasso_genes)) {
  gene <- selected_lasso_genes[i]
  clean_id <- lasso_genes_clean[i]
  
  # Get expression marker
  if (!gene %in% colnames(expr_matrix_clean)) next  # Skip if not found
  marker <- expr_matrix_clean[, gene]
  
  # Skip if bad marker
  if (any(is.na(marker)) || length(unique(marker)) < 2) next
  
  # Compute ROC
  roc <- survivalROC(
    Stime = clinical_data$time,
    status = clinical_data$event,
    marker = marker,
    predict.time = 1095,
    method = "KM"
  )
  
  # Save PNG plot
  filename <- paste0("LASSO_ROC_Plots/", gsub("\\.", "_", gene), "_ROC.png")
  png(filename, width = 800, height = 800, res = 150)
  plot(roc$FP, roc$TP, type = "l", col = "blue", lwd = 2,
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = paste("Time-dependent ROC for", gene))
  abline(0, 1, lty = 2)
  legend("bottomright", paste("AUC =", round(roc$AUC, 3)))
  dev.off()
  
  # Store AUC
  auc_results$AUC[i] <- roc$AUC
}

# Merge HGNC symbols into AUC table
auc_results <- merge(auc_results, gene_map, by.x = "Ensembl_ID", by.y = "ensembl_gene_id", all.x = TRUE)

# Reorder columns
auc_results <- auc_results[, c("Gene_Ensembl_ID_Version", "Ensembl_ID", "hgnc_symbol", "AUC")]
colnames(auc_results)[3] <- "HGNC_Symbol"

# Save to CSV
write.csv(auc_results, "LASSO_ROC_AUC_Annotated.csv", row.names = FALSE)

cat("All ROC plots saved in 'LASSO_ROC_Plots/' and AUC results saved in 'LASSO_ROC_AUC_Annotated.csv'\n")



# External Validation - Using GEO Datasets
#lasso top prognostic genes validation with GSE31210 GEO dataset
# Step 1: Read both input files
geo_df <- read.csv("GSE31210.top.table.csv", stringsAsFactors = FALSE)
lasso_df <- read.csv("LASSO_ROC_AUC_Annotated.csv", stringsAsFactors = FALSE)

# Step 2: Ensure proper column names
colnames(lasso_df) <- trimws(colnames(lasso_df))  # Clean column headers

# Check column names 
print(colnames(lasso_df))

# Now clean the HGNC symbols
if ("HGNC_Symbol" %in% colnames(lasso_df)) {
  lasso_df$HGNC_Symbol <- trimws(as.character(lasso_df$HGNC_Symbol))
  lasso_df <- lasso_df[!is.na(lasso_df$HGNC_Symbol) & lasso_df$HGNC_Symbol != "", ]
} else {
  stop("'HGNC_Symbol' column not found in your file. Please check the header name.")
}

# Step 3: Merge based on HGNC Symbol
merged_df <- merge(lasso_df, geo_df, by.x = "HGNC_Symbol", by.y = "Gene.symbol", all.x = TRUE)

# Step 4: Create validation column
merged_df$Validation <- ifelse(
  is.na(merged_df$adj.P.Val), "Not Found",
  ifelse(merged_df$adj.P.Val < 0.05, "Validated", "Not Significant")
)

# Step 5: Save full merged results
write.csv(merged_df, "LASSO_GSE31210_Validation_All_Matched.csv", row.names = FALSE)

# Step 6: Filter and save only validated genes
validated_df <- merged_df[merged_df$Validation == "Validated", ]
write.csv(validated_df, "LASSO_GSE31210_Validated_Genes.csv", row.names = FALSE)

cat("LASSO validation with GSE31210 completed. Files saved.\n")

##lasso top prognostic genes validation with GSE68465 GEO Dataset
# Step 1: Read both input files
geo_df <- read.csv("GSE68465.top.table.csv", stringsAsFactors = FALSE)
lasso_df <- read.csv("LASSO_ROC_AUC_Annotated.csv", stringsAsFactors = FALSE)

# Step 2: Ensure proper column names
colnames(lasso_df) <- trimws(colnames(lasso_df))  # Clean column headers

# Now clean the HGNC symbols
if ("HGNC_Symbol" %in% colnames(lasso_df)) {
  lasso_df$HGNC_Symbol <- trimws(as.character(lasso_df$HGNC_Symbol))
  lasso_df <- lasso_df[!is.na(lasso_df$HGNC_Symbol) & lasso_df$HGNC_Symbol != "", ]
} else {
  stop("'HGNC_Symbol' column not found in your LASSO file.")
}

# Step 3: Merge based on HGNC Symbol
merged_df <- merge(lasso_df, geo_df, by.x = "HGNC_Symbol", by.y = "Gene.symbol", all.x = TRUE)

# Step 4: Create validation column
merged_df$Validation <- ifelse(
  is.na(merged_df$adj.P.Val), "Not Found",
  ifelse(merged_df$adj.P.Val < 0.05, "Validated", "Not Significant")
)

# Step 5: Save full merged results
write.csv(merged_df, "LASSO_GSE68465_Validation_All_Matched.csv", row.names = FALSE)

# Step 6: Save only validated genes
validated_df <- merged_df[merged_df$Validation == "Validated", ]
write.csv(validated_df, "LASSO_GSE68465_Validated_Genes.csv", row.names = FALSE)

cat("LASSO validation with GSE68465 completed. Files saved.\n")


#finding overlapping genes between lasso_top_prognostic_genes,univariate,multivariate significant genes
# Step 1: Load files
lasso_df <- read.csv("LASSO_ROC_AUC_Annotated.csv", stringsAsFactors = FALSE)
uni_df <- read.csv("significant_prognostic_genes.csv", stringsAsFactors = FALSE)
multi_df <- read.csv("significant_multivariate_cox_genes_annotated.csv", stringsAsFactors = FALSE)

# Step 2: Select relevant columns
lasso_genes <- lasso_df[, c("HGNC_Symbol", "AUC")]
uni_genes <- uni_df[, c("hgnc_symbol", "HR", "HR_conf_low", "HR_conf_high", "p_value")]
multi_genes <- multi_df[, c("hgnc_symbol", "HR", "pvalue")]

# Step 3: Remove empty or NA HGNC symbols
lasso_genes <- lasso_genes[!is.na(lasso_genes$HGNC_Symbol) & lasso_genes$HGNC_Symbol != "", ]
uni_genes <- uni_genes[!is.na(uni_genes$hgnc_symbol) & uni_genes$hgnc_symbol != "", ]
multi_genes <- multi_genes[!is.na(multi_genes$hgnc_symbol) & multi_genes$hgnc_symbol != "", ]

# Step 4: Capitalize all gene symbols for consistent matching
lasso_genes$HGNC_Symbol <- toupper(lasso_genes$HGNC_Symbol)
uni_genes$hgnc_symbol <- toupper(uni_genes$hgnc_symbol)
multi_genes$hgnc_symbol <- toupper(multi_genes$hgnc_symbol)

# Step 5: Find common genes across all 3 datasets
common_genes <- Reduce(intersect, list(
  lasso_genes$HGNC_Symbol,
  uni_genes$hgnc_symbol,
  multi_genes$hgnc_symbol
))

# Step 6: Filter to common genes
lasso_common <- lasso_genes[lasso_genes$HGNC_Symbol %in% common_genes, ]
uni_common <- uni_genes[uni_genes$hgnc_symbol %in% common_genes, ]
multi_common <- multi_genes[multi_genes$hgnc_symbol %in% common_genes, ]

# Rename for merging
colnames(uni_common)[1] <- "HGNC_Symbol"
colnames(multi_common)[1] <- "HGNC_Symbol"

# Step 7: Merge all into one final table
merged <- merge(lasso_common, uni_common, by = "HGNC_Symbol", suffixes = c("_LASSO", "_UNI"))
merged <- merge(merged, multi_common, by = "HGNC_Symbol", suffixes = c("", "_MULTI"))

# Step 8: Save to CSV
write.csv(merged, "Overlapping_Prognostic_LASSO_Cox_Genes.csv", row.names = FALSE)

cat("Done! Saved to Overlapping_Prognostic_LASSO_Cox_Genes.csv\n")


#Finding Overalapping genes between Univariate and Multivariate Significant Genes
# Step 1: Load files
uni_df <- read.csv("significant_prognostic_genes.csv", stringsAsFactors = FALSE)
multi_df <- read.csv("significant_multivariate_cox_genes_annotated.csv", stringsAsFactors = FALSE)

# Step 2: Select relevant columns
uni_genes <- uni_df[, c("hgnc_symbol", "HR", "HR_conf_low", "HR_conf_high", "p_value")]
multi_genes <- multi_df[, c("hgnc_symbol", "HR", "pvalue")]

# Step 3: Remove empty or NA HGNC symbols
uni_genes <- uni_genes[!is.na(uni_genes$hgnc_symbol) & uni_genes$hgnc_symbol != "", ]
multi_genes <- multi_genes[!is.na(multi_genes$hgnc_symbol) & multi_genes$hgnc_symbol != "", ]

# Step 4: Capitalize all gene symbols for consistent matching
uni_genes$hgnc_symbol <- toupper(uni_genes$hgnc_symbol)
multi_genes$hgnc_symbol <- toupper(multi_genes$hgnc_symbol)

# Step 5: Find overlapping genes
common_genes <- intersect(uni_genes$hgnc_symbol, multi_genes$hgnc_symbol)

# Step 6: Filter to common genes
uni_common <- uni_genes[uni_genes$hgnc_symbol %in% common_genes, ]
multi_common <- multi_genes[multi_genes$hgnc_symbol %in% common_genes, ]

# Step 7: Rename for merging
colnames(uni_common)[1] <- "HGNC_Symbol"
colnames(multi_common)[1] <- "HGNC_Symbol"

# Step 8: Merge and save
merged <- merge(uni_common, multi_common, by = "HGNC_Symbol", suffixes = c("_UNI", "_MULTI"))
write.csv(merged, "Overlapping_Uni_Multi_Cox_Genes.csv", row.names = FALSE)

cat("Done! Saved to Overlapping_Uni_Multi_Cox_Genes.csv\n")


#Exyernal Validation using GEO dataset with "Univariate and Multivariate overlapped gene set"
# GEO Dataset 1: GSE31210.top.table.csv
# Univariate + Multivariate Validation with GSE31210

# Step 1: Load data
geo_df <- read.csv("GSE31210.top.table.csv", stringsAsFactors = FALSE)
uni_multi_df <- read.csv("Overlapping_Uni_Multi_Cox_Genes.csv", stringsAsFactors = FALSE)

# Step 2: Clean column names
colnames(uni_multi_df) <- trimws(colnames(uni_multi_df))

# Step 3: Normalize HGNC column name
if ("hgnc_symbol" %in% colnames(uni_multi_df)) {
  uni_multi_df$hgnc_symbol <- trimws(as.character(uni_multi_df$hgnc_symbol))
} else if ("HGNC_Symbol" %in% colnames(uni_multi_df)) {
  uni_multi_df$hgnc_symbol <- trimws(as.character(uni_multi_df$HGNC_Symbol))
} else {
  stop("No hgnc_symbol or HGNC_Symbol column found in Uni+Multi file.")
}

# Step 4: Remove missing symbols
uni_multi_df <- uni_multi_df[!is.na(uni_multi_df$hgnc_symbol) & uni_multi_df$hgnc_symbol != "", ]

# Step 5: Merge
merged_df <- merge(uni_multi_df, geo_df, by.x = "hgnc_symbol", by.y = "Gene.symbol", all.x = TRUE)

# Step 6: Create validation flag
merged_df$Validation <- ifelse(
  is.na(merged_df$adj.P.Val), "Not Found",
  ifelse(merged_df$adj.P.Val < 0.05, "Validated", "Not Significant")
)

# Step 7: Save outputs
write.csv(merged_df, "UniMulti_GSE31210_Validation_All_Matched.csv", row.names = FALSE)
validated_df <- merged_df[merged_df$Validation == "Validated", ]
write.csv(validated_df, "UniMulti_GSE31210_Validated_Genes.csv", row.names = FALSE)


#Univariate + Multivariate Validation with GSE68465 GEO dataset

# Step 1: Load data
geo_df <- read.csv("GSE68465.top.table.csv", stringsAsFactors = FALSE)
uni_multi_df <- read.csv("Overlapping_Uni_Multi_Cox_Genes.csv", stringsAsFactors = FALSE)

# Step 2: Clean column names
colnames(uni_multi_df) <- trimws(colnames(uni_multi_df))

# Step 3: Normalize HGNC column name
if ("hgnc_symbol" %in% colnames(uni_multi_df)) {
  uni_multi_df$hgnc_symbol <- trimws(as.character(uni_multi_df$hgnc_symbol))
} else if ("HGNC_Symbol" %in% colnames(uni_multi_df)) {
  uni_multi_df$hgnc_symbol <- trimws(as.character(uni_multi_df$HGNC_Symbol))
} else {
  stop("No hgnc_symbol or HGNC_Symbol column found in Uni+Multi file.")
}

# Step 4: Remove missing symbols
uni_multi_df <- uni_multi_df[!is.na(uni_multi_df$hgnc_symbol) & uni_multi_df$hgnc_symbol != "", ]

# Step 5: Merge
merged_df <- merge(uni_multi_df, geo_df, by.x = "hgnc_symbol", by.y = "Gene.symbol", all.x = TRUE)

# Step 6: Create validation flag
merged_df$Validation <- ifelse(
  is.na(merged_df$adj.P.Val), "Not Found",
  ifelse(merged_df$adj.P.Val < 0.05, "Validated", "Not Significant")
)

# Step 7: Save outputs
write.csv(merged_df, "UniMulti_GSE68465_Validation_All_Matched.csv", row.names = FALSE)
validated_df <- merged_df[merged_df$Validation == "Validated", ]
write.csv(validated_df, "UniMulti_GSE68465_Validated_Genes.csv", row.names = FALSE)


#Functional Enrichment Analysis
#Combining both externally validated(with GEO dataset) gene set of (Univariate and multivariate overlapped gene set)
#Removing duplicates - finding unique gene symbol for functional enrichment analysis using (DAVID)
# Load your CSV files
df1 <- read.csv("UniMulti_GSE68465_Validated_Genes.csv")
df2 <- read.csv("UniMulti_GSE31210_Validated_Genes.csv")  

genes1 <- df1$hgnc_symbol
genes2 <- df2$hgnc_symbol

# Combine both and remove duplicates
combined_genes <- unique(c(genes1, genes2))

# Save as CSV (no header or row names, perfect for enrichment tools)
write.table(combined_genes, "Validated_Unique_Genes.csv", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

length(combined_genes)
