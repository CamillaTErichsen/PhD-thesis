# Proteomics limma analysis test
# https://support.bioconductor.org/p/64484/#64554
# https://omicsplayground.readthedocs.io/en/latest/faq/ 

# Load libraries
library(tidyverse)   # also loads dplyr and ggplot2
library(limma)
library(ComplexHeatmap)
library(circlize)
library(grid)  
library(ggrepel)

# Set wd
setwd("C:\\Users\\au532203\\OneDrive - Aarhus universitet\\2. PhD\\Data analysis\\Proteomics")


# ------------------------------------------------------------------------------
# Load peak intensity data
# ------------------------------------------------------------------------------
# Load peak intensity data (not log2FC)
peaks_data <- read.csv("Data_PeakIntensity_cleaned_v1.csv", header = T, check.names = FALSE)

id_col <- "ProteinGroups"
stopifnot(id_col %in% names(peaks_data))

# Define columns
annot_cols <- c("ProteinGroups","ProteinAccessions","GeneName","ProteinDescriptions",
                "UniProtIds","ProteinName","CellularComponent","BiologicalProcess","MolecularFunction")
sample_cols <- setdiff(names(peaks_data), annot_cols)

# Replace NaN with NA
peaks_data[sample_cols] <- lapply(peaks_data[sample_cols], function(x) {
  x <- as.numeric(x)
  x[is.nan(x)] <- NA
  x
})
stopifnot(all(sapply(peaks_data[sample_cols], is.numeric)))

# Build matrix
expr <- as.matrix(peaks_data[, sample_cols])
rownames(expr) <- peaks_data$ProteinGroups

# Clean sample names
colnames(expr) <- sub("^X", "", colnames(expr))
head(expr)
stopifnot(!any(is.na(peaks_data$ProteinGroups)))

# Log-transform data (data normalisation)
logexpr <- log2(expr + 1)



# ------------------------------------------------------------------------------
# Load meta data
# ------------------------------------------------------------------------------
meta_data <- read.csv("SURVIVE_Variables.csv", header = T, check.names = FALSE)
meta_data <- meta_data[, c("SID", "diag", "age", "gender", "PMI")]


# Set SID as index
meta_data$SID <- as.character(meta_data$SID)
rownames(meta_data) <- meta_data$SID
meta_data$SID <- NULL
head(meta_data)

# Align meta data with bulk columns
common <- colnames(logexpr)[colnames(logexpr) %in% rownames(meta_data)]
logexpr <- logexpr[, common, drop = FALSE]
meta <- meta_data[common, , drop = FALSE]
stopifnot(identical(colnames(logexpr), rownames(meta)))



# ------------------------------------------------------------------------------
# Adjustments
# ------------------------------------------------------------------------------
# Factor diag
meta$diag <- factor(meta$diag, levels = c("C", "SCZ", "MDD"))
meta$gender <- factor(meta$gender)   # two levels: 0 and 1

# Scale/center age and PMI to improve model interpretability and reduce
# non-essential collinearity. Not strictly necessary but it's a good idea
# to remove some of the heterogeneity often observed in post-mortem studies
meta$age_c <- scale(meta$age, center = TRUE, scale = FALSE)
meta$PMI_c <- scale(meta$PMI, center = TRUE, scale = FALSE)


# Data filtering (already done in Perseus so doesn't change anything)
# Keep protein if present in at least 50 % of samples in one group
keep <- apply(logexpr, 1, function(x) {
  max(tapply(!is.na(x), meta$diag, mean), na.rm = TRUE) >= 0.5
})
logexpr_final <- logexpr[keep, , drop = FALSE]



# ------------------------------------------------------------------------------
# Design matrix
# ------------------------------------------------------------------------------
# Include what the model should adjust for
design <- model.matrix(~ 0 + diag + age_c + gender + PMI_c, data = meta)
colnames(design) <- gsub("^diag", "", colnames(design))



# ------------------------------------------------------------------------------
# Differential expression (DE) analysis: Limma
# ------------------------------------------------------------------------------
# Run linear model (limma)
fit <- lmFit(logexpr_final, design)
contr <- makeContrasts(
  SCZ_vs_C = SCZ - C,
  MDD_vs_C = MDD - C,
  levels = design
)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)


# Results
res_SCZ <- topTable(fit2, coef = "SCZ_vs_C", number = Inf, sort.by = "P")
res_MDD <- topTable(fit2, coef = "MDD_vs_C", number = Inf, sort.by = "P")
# Check
head(res_SCZ)
head(res_MDD)

# Save
write.csv(res_SCZ, "limma_SCZ_ProteinGroups_diag-age-gender-PMI.csv", row.names = TRUE)
write.csv(res_MDD, "limma_MDD_ProteinGroups_diag-age-gender-PMI.csv", row.names = TRUE)

# QC plots
plotMDS(logexpr_final, labels = meta$diag)
plotSA(fit2)



# ------------------------------------------------------------------------------
# Visualise: volcano plot
# ------------------------------------------------------------------------------
# Build protein group --> gene name annotation lookup
annot_lookup <- peaks_data %>%
  select(ProteinGroups, GeneName) %>%
  distinct()

# Merge results with gene names from the annotation lookup using proteingroups
# as common column
res_SCZ2 <- res_SCZ %>%
  tibble::rownames_to_column("ProteinGroups") %>%
  left_join(annot_lookup, by = "ProteinGroups")
# 'GeneName' can sometimes contains multiple names. Pick the first
res_SCZ2 <- res_SCZ2 %>%
  mutate(GeneLabel = sub(";.*$", "", GeneName))
write.csv(res_SCZ2, "DEPs_SCZ_ProteinGroups_diag-age-gender-PMI.csv", row.names = TRUE)

res_MDD2 <- res_MDD %>%
  tibble::rownames_to_column("ProteinGroups") %>%
  left_join(annot_lookup, by = "ProteinGroups")
res_MDD2 <- res_MDD2 %>%
  mutate(GeneLabel = sub(";.*$", "", GeneName))
write.csv(res_MDD2, "DEPs_MDD_ProteinGroups_diag-age-gender-PMI.csv", row.names = TRUE)


# Define thresholds and label based on gene names
# No DEPs padj < 0.05 --> explorative approach with pvalue + log
logFC_cutoff <- 0.3
p_cutoff <- 0.05

res_SCZ2 <- res_SCZ2 %>%
  mutate(status = case_when(
    P.Value < p_cutoff & logFC > logFC_cutoff ~ "Up",
    P.Value < p_cutoff & logFC < -logFC_cutoff ~ "Down",
    TRUE ~ "NS"
  ),
  neglog10p = -log10(P.Value),
  # Use gene where possible, otherwise use proteingroups
  label = ifelse(is.na(GeneLabel) | GeneLabel == "", ProteinGroups, GeneLabel)
  )

res_MDD2 <- res_MDD2 %>%
  mutate(status = case_when(
    P.Value < p_cutoff & logFC > logFC_cutoff ~ "Up",
    P.Value < p_cutoff & logFC < -logFC_cutoff ~ "Down",
    TRUE ~ "NS"
  ),
  neglog10p = -log10(P.Value),
  # Use gene where possible, otherwise use proteingroups
  label = ifelse(is.na(GeneLabel) | GeneLabel == "", ProteinGroups, GeneLabel)
  )


# Set labelling thresholds by p-values
top_n <- 12
to_label_scz <- res_SCZ2 %>%
  filter(status != "NS") %>%
  arrange(P.Value) %>%
  head(top_n)

to_label_mdd <- res_MDD2 %>%
  filter(status != "NS") %>%
  arrange(P.Value) %>%
  head(top_n)


# Volcano plot with thresholds and gene names
ggplot(res_SCZ2, aes(x = logFC, y = neglog10p)) + 
  geom_point(aes(colour = status), alpha = 0.8, size = 1.6) + 
  geom_hline(yintercept = -log10(p_cutoff), linetype = 2) + 
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = 3) + 
  geom_text_repel(
    data = to_label_scz,
    aes(label = label),
    size = 3,
    max.overlaps = Inf, 
    box.padding = 0.35,
    point.padding = 0.2
  ) + 
  scale_colour_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) + 
  labs(
    title = "Schizophrenia vs control",
    x = "Log2FC",
    y = "-log10(p-value)",
    colour = ""
  ) + 
  theme_classic()
ggsave("VolcanoPlot_SCZ.png", width = 8, height = 5)

ggplot(res_MDD2, aes(x = logFC, y = neglog10p)) + 
  geom_point(aes(colour = status), alpha = 0.8, size = 1.6) + 
  geom_hline(yintercept = -log10(p_cutoff), linetype = 2) + 
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = 3) + 
  geom_text_repel(
    data = to_label_mdd,
    aes(label = label),
    size = 3,
    max.overlaps = Inf, 
    box.padding = 0.35,
    point.padding = 0.2
  ) + 
  scale_colour_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) + 
  labs(
    title = "Depression vs control",
    x = "Log2FC",
    y = "-log10(p-value)",
    colour = ""
  ) + 
  theme_classic()
ggsave("VolcanoPlot_MDD.png", width = 8, height = 5)



# ------------------------------------------------------------------------------
# Visualise: heatmap. Distance metric = pearson, linking method = average
# ------------------------------------------------------------------------------
# Build function for complex heatmap
plot_dep_heatmap <- function(res2, logmat, meta, title,
                             keep_diags = c("C", "SCZ"),
                             p_cutoff = 0.05,
                             logFC_cutoff = 0.3,
                             cluster_columns = TRUE,
                             show_row_names = TRUE) {
  
  # Subset samples to desired diagnoses
  meta_sub <- meta %>% dplyr::filter(diag %in% keep_diags)
  
  logmat_sub <- logmat[, rownames(meta_sub), drop = FALSE]
  
  # Select DEPs
  dep_tbl <- res2 %>%
    dplyr::filter(P.Value < p_cutoff,
                  abs(logFC) > logFC_cutoff) %>%
    dplyr::arrange(P.Value)
  
  message(title, ": DEPs selected = ", nrow(dep_tbl))
  if (nrow(dep_tbl) < 2) stop("Too few DEPs to plot with current thresholds.")
  
  # Subset expression matrix
  dep_ids <- dep_tbl$ProteinGroups
  dep_ids <- dep_ids[dep_ids %in% rownames(logmat_sub)]
  mat <- logmat_sub[dep_ids, , drop = FALSE]
  
  # Row z-score
  mat_z <- t(scale(t(mat)))
  # Replace remaining NA values after scaling
  mat_z[is.na(mat_z)] <- 0
  
  # Row labels: GeneLabel (fallback to ProteinGroups) where possible
  dep_tbl <- dep_tbl %>%
    dplyr::mutate(
      GeneLabel = ifelse(
        is.na(GeneLabel) | GeneLabel == "",
        ProteinGroups,
        GeneLabel))
  
  label_map <- dep_tbl$GeneLabel
  names(label_map) <- dep_tbl$ProteinGroups
  rownames(mat_z) <- make.unique(label_map[rownames(mat_z)])
  
  # Column annotations (aligned)
  meta_ann <- meta_sub[colnames(mat_z), , drop = FALSE]
  
  ha <- HeatmapAnnotation(
    diag = meta_ann$diag,
    gender = meta_ann$gender,
    age = meta_ann$age,
    PMI = meta_ann$PMI,
    annotation_name_side = "left"
  )
  
  # For column split in heatmap
  direction <- dep_tbl$logFC
  names(direction) <- dep_tbl$ProteinGroups
  row_split <- ifelse(direction[dep_ids] > 0, "Upregulated", "Downregulated")
  
  # Define the complex heatmap
  Heatmap(
    mat_z,
    name = "z-score",
    top_annotation = ha,
    # Clustering
    cluster_columns = cluster_columns,
    cluster_rows = TRUE,
    # Distance metric
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    # Linkage method
    clustering_method_rows = "average",
    clustering_method_columns = "average",
    
    # Make split between clusters and proteins
    row_split = row_split,
    
    # Display
    show_column_names = FALSE,
    show_row_names = show_row_names,
    row_names_gp = grid::gpar(fontsize = 7),  # good for ~60–70 genes
    column_title = title
  )
}


# Heatmaps clustered by samples
ht_scz_cluster <- plot_dep_heatmap(
  res2 = res_SCZ2,
  logmat = logexpr_final,
  meta = meta,
  keep_diags = c("C", "SCZ"),
  title = "SCZ vs Control DEPs (Pearson correlation + average linkage)",
  cluster_columns = TRUE,
  show_row_names = TRUE
)
draw(ht_scz_cluster)  # 69 DEPs
# Save manually: 650 x 800

ht_mdd_cluster <- plot_dep_heatmap(
  res2 = res_MDD2,
  logmat = logexpr_final,
  meta = meta,
  keep_diags = c("C", "MDD"),
  title = "MDD vs Control DEPs (Pearson correlation + average linkage)",
  cluster_columns = TRUE,
  show_row_names = TRUE
)
draw(ht_mdd_cluster)  # 60 DEPs

# Save the heatmaps
pdf("Heatmap_SCZ_SampleClustered.pdf", width = 10, height = 8)
draw(ht_scz_cluster)
dev.off()
pdf("Heatmap_MDD_SampleClustered.pdf", width = 10, height = 8)
draw(ht_scz_ordered)
dev.off()



# ------------------------------------------------------------------------------
# Secondary (DE) analysis: ROTS
# ------------------------------------------------------------------------------
# A secondary analysis with ROTS only to check the sensitivity/robustness of limma
# results; it is only for ranking concordance, and it does not handle covariates.
library(ROTS)

# Make annotation lookup
annot_lookup <- peaks_data %>%
  dplyr::select(ProteinGroups, GeneName) %>%
  dplyr::distinct()

# ROTS needs features in rows and samples in columns
# Helper function
run_rots_contrast <- function(logmat, meta, group1, group2,
                              annot_lookup,
                              B = 1000, K = 500, seed = 123) {
  
  # Subset to two groups only
  keep_samples <- rownames(meta)[meta$diag %in% c(group1, group2)]
  mat_sub <- logmat[, keep_samples, drop = FALSE]
  meta_sub <- meta[keep_samples, , drop = FALSE]
  
  # Ensure group order: group1 = reference/control, group2 = case
  meta_sub$group <- factor(meta_sub$diag, levels = c(group1, group2))
  
  keep_complete <- rowSums(is.na(mat_sub)) == 0
  mat_sub_complete <- mat_sub[keep_complete, , drop = FALSE]
  
  message(group2, " vs ", group1, ": proteins retained after complete-case filter = ",
          nrow(mat_sub_complete))
  
  # Groups must be numeric or factor-like vector
  groups <- as.numeric(meta_sub$group) - 1  # group1 = 0, group2 = 1
  
  set.seed(seed)
  rots_out <- ROTS(
    data = mat_sub_complete,
    groups = groups,
    B = B,
    K = K,
    seed = seed,
    log = FALSE,
    progress = TRUE
  )
  
  rots_res <- data.frame(
    ProteinGroups = rownames(mat_sub_complete),
    ROTS_stat = rots_out$d,
    ROTS_pvalue = rots_out$pvalue,
    ROTS_FDR = rots_out$FDR,
    stringsAsFactors = FALSE
  ) %>%
    # Add annotation
    dplyr::left_join(annot_lookup, by = "ProteinGroups") %>%
    dplyr::mutate(
      GeneLabel = sub(";.*$", "", GeneName)
    ) %>%
    dplyr::arrange(ROTS_pvalue)
  
  return(list(
    rots_object = rots_out,
    results = rots_res
  ))
}


# Run contrasts
rots_SCZ <- run_rots_contrast(
  logmat = logexpr_final,
  meta = meta,
  group1 = "C",
  group2 = "SCZ",
  annot_lookup = annot_lookup
)

rots_MDD <- run_rots_contrast(
  logmat = logexpr_final,
  meta = meta,
  group1 = "C",
  group2 = "MDD",
  annot_lookup = annot_lookup
)

# Save
write.csv(rots_SCZ$results, "ROTS_SCZ_vs_C.csv", row.names = FALSE)
write.csv(rots_MDD$results, "ROTS_MDD_vs_C.csv", row.names = FALSE)



# ------------------------------------------------------------------------------
# Compare ROTS with limma
# ------------------------------------------------------------------------------
# Add ProteinGroups column to limma instead of being rownames
limma_SCZ <- res_SCZ %>%
  tibble::rownames_to_column("ProteinGroups") %>%
  dplyr::select(ProteinGroups, limma_logFC = logFC, limma_t = t,
                limma_p = P.Value, limma_FDR = adj.P.Val)

limma_MDD <- res_MDD %>%
  tibble::rownames_to_column("ProteinGroups") %>%
  dplyr::select(ProteinGroups, limma_logFC = logFC, limma_t = t,
                limma_p = P.Value, limma_FDR = adj.P.Val)

# Merge
comp_SCZ <- limma_SCZ %>%
  inner_join(rots_SCZ$results, by = "ProteinGroups")

comp_MDD <- limma_MDD %>%
  inner_join(rots_MDD$results, by = "ProteinGroups")

# Save
write.csv(comp_SCZ, "LimmavsROTS_SCZ_vs_C.csv", row.names = FALSE)
write.csv(comp_MDD, "LimmavsROTS_MDD_vs_C.csv", row.names = FALSE)


# Rank concordance
cor(comp_SCZ$limma_t, comp_SCZ$ROTS_stat, method = "spearman", use = "complete.obs")   # Spearman correlation = 0.8993136
cor(comp_MDD$limma_t, comp_MDD$ROTS_stat, method = "spearman", use = "complete.obs")   # Spearman correlation = 0.8978318

# Overlap between top 100 proteins
top_n <- 100
# Schizophrenia
top_limma_SCZ <- comp_SCZ %>%
  arrange(limma_p) %>%
  slice_head(n = top_n) %>%
  pull(ProteinGroups)

top_rots_SCZ <- comp_SCZ %>%
  arrange(ROTS_pvalue) %>%
  slice_head(n = top_n) %>%
  pull(ProteinGroups)
length(intersect(top_limma_SCZ, top_rots_SCZ))    # 51

# Depression
top_limma_MDD <- comp_MDD %>%
  arrange(limma_p) %>%
  slice_head(n = top_n) %>%
  pull(ProteinGroups)
 
top_rots_MDD <- comp_MDD %>%
  arrange(ROTS_pvalue) %>%
  slice_head(n = top_n) %>%
  pull(ProteinGroups)
length(intersect(top_limma_MDD, top_rots_MDD))   # 49


# Plot concordance
library(ggplot2)
ggplot(comp_SCZ, aes(x = limma_t, y = ROTS_stat)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "SCZ vs C: limma vs ROTS ranking concordance",
    x = "limma moderated t",
    y = "ROTS statistic"
  )
ggsave("Limma-ROTS-concordance_SCZ.png", width = 11, height = 8)

ggplot(comp_MDD, aes(x = limma_t, y = ROTS_stat)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "MDD vs C: limma vs ROTS ranking concordance",
    x = "limma moderated t",
    y = "ROTS statistic"
  )
ggsave("Limma-ROTS-concordance_MDD.png", width = 11, height = 8)

