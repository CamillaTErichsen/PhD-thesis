# RNA-Seq DESeq2 for differentially expressed genes

# Load libraries
library(tidyverse)
library(clusterProfiler)
library(DESeq2)

# Set wd
setwd("C:\\Users\\au532203\\OneDrive - Aarhus universitet\\2. PhD\\Data analysis\\RNA-Seq")


# ------------------------------------------------------------------------------
# Load count data
# ------------------------------------------------------------------------------
# Load all csv files by common prefix
files <- list.files(pattern = "EksternData_GeneCounts-.*\\.csv")
# Read an prepare each file (contains gene names and count data)
prep_files <- function(file) {
  # Read file
  df <- read.csv(file)[, 1:2]
  # Extract ID from filename and rename count column for clarity
  id <- sub("EksternData_GeneCounts-", "", sub("\\.csv$", "", file))   # Extract SID from file suffix
  colnames(df)[2] <- id                                                # Rename gene counts to SID
  return(df)
}
# Read files into a large list
data_list <- lapply(files, prep_files)

# Match files by gene name and combine to one count matrix
count_matrix <- Reduce(function(x, y) merge(x, y, by ="Gene_name", all = FALSE), data_list)
print(colnames(count_matrix))
# Set gene names as index column
rownames(count_matrix) <- count_matrix$Gene_name
# Remove the gene names column
count_matrix$Gene_name <- NULL
# Check for duplicated genes and NA values
#sum(duplicated(count_matrix[0]))
sum(duplicated(rownames(count_matrix)))
which(is.na(count_matrix))

# Remove . for IDs
colnames(count_matrix) <- gsub("\\.$", "", colnames(count_matrix))
print(colnames(count_matrix))  # View SIDs
head(count_matrix)

# Remove SID from header
colnames(count_matrix) <- sub(".*-", "", colnames(count_matrix))
head(data_combined)

# Remove outlier
count_matrix <- count_matrix %>%
  select(-'[TYPE ID HERE]')



# ------------------------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------------------------
features <- read.csv("FeatureMatrix.csv", header = T)
# Set PNR as index
rownames(features) <- features$PNR
# Remove PNR column
features$PNR <- NULL
head(features)
# Remove outlier
features <- features[!rownames(features) %in% "OUTLIER_ID", ]    # GDPR
features$Death_Cause[is.na(features$Death_Cayse)] <- "unknown"
head(features)

# Sanity checks: both must be true
all(colnames(count_matrix)) == rownames(features)     # T
all(colnames(count_matrix)) %in% rownames(features)   # T



# ------------------------------------------------------------------------------
# Construct deseq object
# ------------------------------------------------------------------------------
# If you only want to correlate to diagnosis:
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = features,
                              design = ~ diag_ICD10)

# For MDD
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = features,
                              design = ~ diag_ICD10 + age + RIN + gender + Drugs_3Months)

# For SCZ
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = features,
                              design = ~ diag_ICD10 + age + RIN + gender)

# Filter low counts
keep <- rowSums(counts(dds)) >= 10
dds_filt <- dds[keep, ]
dds_filt

# Set factor level (contrl = ref)
dds_filt$diag_IDC10 <- relevel(dds_filt$diag_IDC10, ref = "CON")
# MDD vs control and SCZ vs control



# ------------------------------------------------------------------------------
# Perform DESeq2
# ------------------------------------------------------------------------------
dds_run <- DESeq(dds_filt)
res <- results(dds_run)
res005 <- results(dds_run, alpha = 0.05)

# Contrast
resultsNames(dds_run)
res_scz <- results(dds_run, contrast = c("diag_ICD10", "SCZ", "CON"))    # SCZ compared to control
res_mdd <- results(dds_run, contrast = c("diag_ICD10", "MDD", "CON"))

# Save
write.csv(res_scz, "DESeq2_SCZ.csv")
write.csv(res_mdd, "DESeq2_MDD.csv")



# ------------------------------------------------------------------------------
# Visualise
# ------------------------------------------------------------------------------
# Clusters
vsdata <- vst(dds_run, blind = F)
plotPCA(vsdata, intgroup = "gender")

# Dispersion
plotDispEsts(dds_run)

# Sign. genes
plotMA(res_scz)
plotMA(res_mdd)

DEGs_scz <- na.omit(res_scz)
DEGs_scz <- DEGs_scz[DEGs_scz$padj < 0.05, ]

DEGs_mdd <- na.omit(res_mdd)
DEGs_mdd <- DEGs_mdd[DEGs_mdd$padj < 0.05, ]

