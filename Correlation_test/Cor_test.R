# This is a script related to Correlation analysis in R. It contains the links to main sources,
# explanations and hints. The example is based on the data extracted from Kaggle:
# https://www.kaggle.com/datasets/amandam1/colorectal-cancer-patients?resource=download&select=Colorectal+Cancer+Patient+Data.csv
# The sources for coding part are:
# https://r-statistics.co/Importing-Data-in-R.html
# https://www.datanovia.com/learn/biostatistics/correlation/correlation-test-in-r#google_vignette

# IMPORT the data and look into it ------
library(readr)
library(dplyr)
cancer_meta <- read_csv("Colorectal Cancer Patient Data.csv")
cancer_expr <- read_csv("Colorectal Cancer Gene Expression Data.csv")

# Clean the data
cancer_meta[1] <- NULL
cancer_expr[1] <- NULL

cancer_meta <- cancer_meta %>% filter(row_number() <= n()-1)

# Explanation of the metadata: 
# Dukes Stage: A to D (development/progression of disease)
# Location (surgery): Left, Right, Colon or Rectum
# DFS: Disease-free survival, months (survival without the disease returning)
# DFS event: 0 or 1 (with 1 = event)
# Adj_Radio: If the patient also received radiotherapy
# Adj_Chem: If the patient also received chemotherapy

# Variables could be used for practicing:
# numerical: gene expression, Age, DFS
# ordinal: Dukes Stage

# Explanation of the expression data:
# ID_REF is an Affymetrix Microarray Probe Set ID. Specifically this identifier belongs to the 
# Affymetrix Human Genome U133 Plus 2.0 array (often cataloged as platform GPL570 in public repositories like NCBI GEO).
# Affymetrix uses systematic suffixes at the end of their probe IDs to describe how uniquely or precisely the probes target a gene transcript.
# So we need to convert the probes into the real biology - genes

# ANNOTATION of the probes ------
# Install once
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")

library(hgu133plus2.db)
library(AnnotationDbi)

# Map probes
annotation <- AnnotationDbi::select(hgu133plus2.db, keys = cancer_expr$ID_REF, keytype = "PROBEID", columns = c("SYMBOL", "GENENAME", "ENTREZID"))

# Check the percentage mapped
mean(!is.na(annotation$SYMBOL)) * 100

# Add the annotations to expression data
expression_ann <- cancer_expr %>% left_join(annotation, by = c("ID_REF" = "PROBEID"))

# Find probes mapping to multiple genes
annotation %>% group_by(PROBEID) %>% filter(n() > 1) %>% arrange(PROBEID) %>% dplyr::select(PROBEID, SYMBOL, GENENAME, ENTREZID)
annotation %>% count(PROBEID) %>% filter(n > 1) %>% nrow()

# Find multiple probes for one gene
annotation %>% filter(!is.na(SYMBOL)) %>% distinct(PROBEID, SYMBOL) %>% count(SYMBOL, name = "n_probes") %>% arrange(desc(n_probes))

# Remove ambiguous probes
unique_probes <- annotation %>% filter(!is.na(SYMBOL)) %>% group_by(PROBEID) %>% filter(n() == 1) %>% ungroup()

# Join clean annotation back to expression_data
clean_cancer_expr <- cancer_expr %>% inner_join(unique_probes, by = c("ID_REF" = "PROBEID"))

# Deal with several probes mapping to the same gene
# Strategy: use the mean expression across probes belonging to the same gene
gene_expr <- clean_cancer_expr %>% dplyr::select(-ID_REF, -GENENAME, -ENTREZID) %>% group_by(SYMBOL) %>% summarise(
  across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>% ungroup()

# Check whether there is no sample mismatch between meta and expr.data
expr_samples <- colnames(gene_expr)[-1]
meta_samples <- cancer_meta$ID_REF

length(expr_samples)
length(meta_samples)

# Check if they are exactly the same samples in the same order
setequal(expr_samples, meta_samples)
identical(expr_samples, meta_samples)

# FILTERING ------
# Let's simplify the workflow decreasing the number of genes based on variance
gene_variance <- apply(expr_mat, 1, var, na.rm = TRUE)

# Sort by variance
gene_variance <- sort(gene_variance, decreasing = TRUE)

# Inspect the variance distribution
hist(gene_variance, breaks = 50, main = "Distribution of gene-expression variance", xlab = "Variance")

# Extract top 10 genes
top_10_var_genes <- names(gene_variance)[1:10]

# Combine the genes with metadata to create one dataframe
expr_top10 <- expr_mat[top_10_var_genes, ]
expr_top10_t <- as.data.frame(t(expr_top10))
expr_top10_t$ID_REF <- rownames(expr_top10_t)

cor_data <- cancer_meta %>% left_join(expr_top10_t, by = "ID_REF")

colnames(cor_data)[which(names(cor_data) == "Age (in years)")] <- "Age"
colnames(cor_data)[which(names(cor_data) == "DFS (in months)")] <- "DFS"

# ONE GENE EXAMPLE ------
# Look at the data first, lets choose for now Age and gene REG4
library(ggpubr)

ggscatter(
  cor_data, x = "Age", y = "REG4",
  add = "reg.line", conf.int = TRUE,
  color = "#3a86d4",
  add.params = list(color = "#1f4e79"),
  xlab = "Age", ylab = "Gene REG4"
) +
  stat_cor(method = "pearson")          # prints R and the p-value on the panel

# Check if the variables normal (in case we choose Pearson)
library(rstatix)
cor_data %>% shapiro_test(Age, REG4)

ggarrange(ggqqplot(cor_data, "Age", title = "Age"), ggqqplot(cor_data, "REG4", title = "REG4"), ncol = 2)

# Conclusion: 
# Age did not significantly deviate from normality (Shapiro–Wilk, p = 0.054), whereas REG4 
# expression showed a significant deviation from normality (p < 0.001). Therefore, 
# the assumption of normality was not satisfied for both variables.

# Run the correlation test
# with rstatix 
cor_data %>% cor_test(Age, REG4, method = "spearman")
# with base R
res <- cor.test(cor_data$Age, cor_data$REG4, method = "spearman")

res$estimate
res$p.value
# Results: 
# S = 38029, p-value = 0.7437
# alternative hypothesis: true rho is not equal to 0
# sample estimates: rho 0.04236553 (weak)

# Conclusion: there is no evidence of an Age-REG4 monotonic association

# CORRELATION MATRIX ------
# Convert Duke stages to a meaningful numeric progression
cor_data$Dukes_numeric <- as.numeric(factor(cor_data$`Dukes Stage`, levels = c("A", "B", "C", "D"), ordered = TRUE))

# Prepare the final table
mydata <- cor_data %>% select(Age, Dukes_numeric, DFS, all_of(top_10_var_genes))

# matrix with rstatix
corr_mat <- mydata %>% cor_mat(method = "spearman")
corr_mat

# with base R
cor_res <- cor(mydata, method = "spearman")
round(cor_res, 2)

# get the p-values
p_val <- mydata %>% cor_mat(method = "spearman") %>% cor_get_pval()
# each cell is the p-value for that pair, and anything below 0.05 flags a correlation that is statistically significant at the 5% level. 
library(Hmisc)
res2 <- rcorr(as.matrix(mydata), type = "spearman")
res2$P
res2$r

# Reshape to a tidy table
cor_table <- mydata %>% cor_mat(method = "spearman") %>% pull_lower_triangle() %>% cor_gather()

# Mark the significant correlations
mydata %>% cor_mat(method = "spearman") %>% cor_mark_significant()
# Stars follow the usual convention: * p<0.05, ** p<0.01, *** p<0.0001.

# Visualize the correlation matrix
library(corrplot)

corrplot(cor_res, type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("#b2182b", "white", "#3a86d4"))(200))

# Interpretation: Positive correlations are blue, negative red; bigger, darker circles 
# indicate stronger relationships.

# Show only nominally significant correlations
p.mat <- res2$P
diag(p.mat) <- 0

corrplot(cor_res, p.mat = p.mat, sig.level = 0.05, insig = "blank",
         diag = FALSE, tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("#b2182b", "white", "#3a86d4"))(200))
# OTHER views on results ------
# Scatterplot matrix
library(GGally)
ggpairs(mydata) + theme_minimal()

# Heatmap
col <- colorRampPalette(c("#b2182b", "white", "#3a86d4"))(20)
par(mar = c(10, 10, 2, 2))
heatmap(cor_res, col = col, symm = TRUE, cexCol = 0.7, cexRow = 0.7)

# Replaces coefficients with symbols for a compact text view
symnum(cor_res, abbr.colnames = FALSE)

# FDR correction (if it is needed for another data) ------
cor_table <- cor_table %>%
  mutate(padj = p.adjust(p, method = "BH"))
cor_table

# Filter
cor_table %>% filter(padj < 0.05)

# SEPARATION of genes and clinical variables ------
# Gene expression matrix
gene_data <- mydata %>% select(all_of(top_10_var_genes))

# Clinical variables
clinical_data <- mydata %>% select(Age, Dukes_numeric, DFS)

# Obtain the rectangular matrix
cor_gene_clinical <- cor(gene_data, clinical_data, method = "spearman")  # use = "pairwise.complete.obs" for missing values
round(cor_gene_clinical, 2)

# p-values again
results<- rcorr(as.matrix(cbind(gene_data, clinical_data)), type = "spearman")
results$r   
results$P 

# Extract necessary portion
p_gene_clinical <- results$P[top_10_var_genes, c("Age", "Dukes_numeric", "DFS")]
round(p_gene_clinical, 4)

# Heatmap
library(pheatmap)

pheatmap(cor_gene_clinical, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
         number_format = "%.2f", breaks = seq(-1, 1, length.out = 101), color = colorRampPalette(
           c("#b2182b", "white", "#3a86d4"))(100), border_color = "white", angle_col = 0)





