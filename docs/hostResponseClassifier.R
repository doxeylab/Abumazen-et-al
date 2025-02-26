## The following code continues from the analysis described in Rmarkdown-V2.html

# Get a reduced df of genes for bacterial and viral

# Contains viral upDEGs with viral-only and bacterial-only patients (n = 64)
viral_df <- assay(prld[viral_upDEGs , ])
viral_df_subsetted <- viral_df[ , onlyvirusNames]

# Contains bacterial upDEGs with viral-only and bacterial-only patients (n = 64)
bac_df <- assay(prld[bac_upDEGs , ])
bac_df_subsetted <- bac_df[ , onlyvirusNames]

# Create new metadata object that will be used
new_metadata <- metadata
new_metadata <- metadata[rownames(metadata) %in% onlyvirusNames, ]

# Ensure that the rowname order matches
matched_indices <- match(colnames(viral_df_subsetted), rownames(new_metadata))
new_metadata <- new_metadata[matched_indices, ]

# Subset the metadata into the two groups -> viral and bacterial only
only_viral <- subset(new_metadata, bv_both_none == 2)$Filename
only_bacterial <- subset(new_metadata, bv_both_none == 1)$Filename

# Subset Viral upDEGs df, with viral-only and bacterial only groups
only_viral_Vupdeg <- viral_df_subsetted[ , only_viral]
only_bacterial_Vupdeg <- viral_df_subsetted[ , only_bacterial]

# Subset Bacterial upDEGs df, with viral-only and bacterial only groups
only_viral_Bupdeg <- bac_df_subsetted[ , only_viral]
only_bacterial_Bupdeg <- bac_df_subsetted[ , only_bacterial]

# GENE SELECTION PART 1 -> Perform t-tests to obtain the top 25 genes (by p-value)

###### VIRAL FIRST ###### 

# Viral upDEGs: perform Welch’s two-sample t-tests + bonferroni correction
pvals_viral <- mapply(function(row_viral, row_bacterial) {
  t.test(row_viral, row_bacterial)$p.value
}, as.data.frame(t(only_viral_Vupdeg)), as.data.frame(t(only_bacterial_Vupdeg)))

adjusted_pvals_viral <- sort(p.adjust(pvals_viral, method = "bonferroni"), decreasing = FALSE)
viral_genes_filtered <- adjusted_pvals_viral[1:25]

###### BACTERIAL SECOND ###### 

# Bacterial upDEGs: perform Welch’s two-sample t-tests + bonferroni correction
pvals_bacterial <- mapply(function(row_viral, row_bacterial) {
  t.test(row_viral, row_bacterial)$p.value
}, as.data.frame(t(only_viral_Bupdeg)), as.data.frame(t(only_bacterial_Bupdeg)))

adjusted_pvals_bacterial <- sort(p.adjust(pvals_bacterial, method = "bonferroni"), decreasing = FALSE)
bacterial_genes_filtered <- adjusted_pvals_bacterial[1:25]

# Current updated infection response scores (for reference, not used again)
viralResponseScore_filtered = apply(temp[names(viral_genes_filtered),],2,mean)
bacterialResponseScore_filtered = apply(temp[names(bacterial_genes_filtered),],2,mean)

# FEATURE SELECTION PART 2 -> Remove redundancy

library(Hmisc)
library(caret)

###### VIRAL FIRST ###### 

cor_viral = cor(t(as.matrix(temp[names(viral_genes_filtered),])))
cor_p_vals_viral <- rcorr(t(as.matrix(temp[names(viral_genes_filtered),])))

# For gene pairs where r > 0.75, removes gene with the largest mean correlation
cor_matrix <- cor(t(as.matrix(temp[names(viral_genes_filtered),])))
cor_threshold <- 0.75
genes_to_remove_idx <- findCorrelation(cor_matrix, cutoff = cor_threshold, exact = TRUE)
viral_genes_to_keep <- names(viral_genes_filtered)[-genes_to_remove_idx]
print(viral_genes_to_keep)

###### BACTERIAL SECOND ###### 

# Generates the correlation matrix for these genes
cor_bac = cor(t(as.matrix(temp[names(bacterial_genes_filtered),])))
cor_p_vals_bac <- rcorr(t(as.matrix(temp[names(bacterial_genes_filtered),])))

# For gene pairs where r > 0.75, removes gene with the largest mean correlation
cor_matrix <- cor(t(as.matrix(temp[names(bacterial_genes_filtered),])))
cor_threshold <- 0.75
genes_to_remove_idx <- findCorrelation(cor_matrix, cutoff = cor_threshold, exact = TRUE)
bacterial_genes_to_keep <- names(bacterial_genes_filtered)[-genes_to_remove_idx]
print(bacterial_genes_to_keep)


# Use the caret package to generate the ML models with CV

# To reproduce the results
set.seed(1234) # Arbitrary value, 1234

# Remove psudogenes and lnc genes
remove_prefixes <- c("^CKAP", "^LINC", "^LUCA", "^MIR2", "^TMEM", "^RPL", "^GTF", "^TOMM", "^CXCR2P1", "^MAN1B1", "^USP30")
viral_genes_to_keep <- viral_genes_to_keep[!grepl(paste(remove_prefixes, collapse = "|"), viral_genes_to_keep)]
bacterial_genes_to_keep <- bacterial_genes_to_keep[!grepl(paste(remove_prefixes, collapse = "|"), bacterial_genes_to_keep)]

###### VIRAL FIRST ###### 

# response_var1 will contain RNA-seq binary classifications
response_var1 = vector(length = nrow(metadata))
response_var1[] = 0
response_var1[which(metadata$high_pathogen == 1 | metadata$high_pathogen == 3)] = 1

# response_var2 will contain clinical binary classifications
response_var2 = vector(length = nrow(metadata))
response_var2[] = 0
response_var2[which(metadata$bv_both_none == 2 | metadata$bv_both_none == 3)] = 1

# a) Train the Random Forest with cross-validation (RNA-seq classifications)

data <- t(temp[viral_genes_to_keep, ])
labels <- factor(response_var1, levels = c(0, 1), labels = c("Class0", "Class1"))
control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)

rf_cv1 <- train(
  x = data, 
  y = labels, 
  method = "rf", 
  ntree = 500, 
  trControl = control, 
  metric = "ROC"
)
print(rf_cv1)

# The best performing model, predicting viral positive as defined by our RNA-seq classifications
print(max(rf_cv1$results$ROC))

# b) Train the Random Forest with cross-validation (Clinical labels)

data <- t(temp[viral_genes_to_keep, ])
labels <- factor(response_var2, levels = c(0, 1), labels = c("Class0", "Class1"))
control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)

rf_cv2 <- train(
  x = data, 
  y = labels, 
  method = "rf", 
  ntree = 500, 
  trControl = control, 
  metric = "ROC", 
)
print(rf_cv2)

# The best performing model, predicting viral positive as defined by our clinical labels
print(max(rf_cv2$results$ROC))

###### BACTERIAL SECOND ###### 

# response_bac1 will contain RNA-seq binary classifications
response_bac1 = vector(length = nrow(metadata))
response_bac1[] = 0
response_bac1[which(metadata$high_pathogen == 2 | metadata$high_pathogen == 3)] = 1

# response_bac2 will contain clinical binary classifications
response_bac2 = vector(length = nrow(metadata))
response_bac2[] = 0
response_bac2[which(metadata$bv_both_none == 1 | metadata$bv_both_none == 3)] = 1

# a) Train the Random Forest with cross-validation (RNA-seq classifications)

data <- t(temp[bacterial_genes_to_keep, ])
labels <- factor(response_bac1, levels = c(0, 1), labels = c("Class0", "Class1"))
control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)

rf_cv3 <- train(
  x = data, 
  y = labels, 
  method = "rf", 
  ntree = 500, 
  trControl = control, 
  metric = "ROC", 
)
print(rf_cv3)

# The best performing model, predicting bacterial positive as defined by our RNA-seq classifications
print(max(rf_cv3$results$ROC))

# b) Train the Random Forest with cross-validation (clinical labels)

data <- t(temp[bacterial_genes_to_keep, ])
labels <- factor(response_bac2, levels = c(0, 1), labels = c("Class0", "Class1"))
control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)

rf_cv4 <- train(
  x = data, 
  y = labels, 
  method = "rf", 
  ntree = 500, 
  trControl = control, 
  metric = "ROC", 
)
print(rf_cv4)

# The best performing model, predicting bacterial positive as defined by our clinical labels
print(max(rf_cv4$results$ROC))
