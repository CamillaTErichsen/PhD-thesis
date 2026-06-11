# limma

# Load libraries
library(tidyverse)
library(pheatmap)
library(limma)

# Set working directory
setwd("C:\\Users\\au532203\\OneDrive - Aarhus universitet\\2. PhD\\Data analysis\\Nanostring")


# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
# Load data
counts <- read.csv("CountData_SCZ.csv", check.names = FALSE)
names(counts)[1] <- "GeneName"
count_mat <- counts %>%
  column_to_rownames("GeneName") %>% 
  as.matrix()
colnames(count_mat) <- as.character(colnames(count_mat))

# Load meta data
meta <- read.csv("SURVIVE_Variables.csv", header = TRUE) %>%
  select(SID, diag) %>%
  mutate(SID = as.character(SID)) %>%
  filter(diag %in% c("C", "SCZ")) %>%
  mutate(diag = factor(diag, levels = c("C", "SCZ")))

# Keep only C and SCZ in count data: Subset using meta data
keep_ids <- intersect(colnames(count_mat), meta$SID)
count_mat_scz <- count_mat[, keep_ids, drop = FALSE]
meta_scz <- meta %>%
  filter(SID %in% keep_ids) %>%
  arrange(match(SID, colnames(count_mat_scz)))
stopifnot(all(meta_scz$SID == colnames(count_mat_scz)))
stopifnot(!any(is.na(meta_scz$diag)))



# ------------------------------------------------------------------------------
# Limma
# ------------------------------------------------------------------------------
mat <- log2 (count_mat_scz)
meta_scz$diag <- factor(meta_scz$diag, levels = c("C", "SCZ"))

# Design matrix
design <- model.matrix(~ diag, data = meta_scz)

# Fit model
fit <- lmFit(mat, design)
fit <- eBayes(fit)

results <- topTable(
  fit,
  coef = "diagSCZ",
  number = Inf,
  sort.by = "none"
)

results <- results %>%
  rownames_to_column("gene")

# Degrees of freedom
df <- fit$df.total[1]
# Critical t value for 95% CI
t_crit <- qt(0.975, df)

# Calculate standard error
results <- results %>%
  mutate(
    SE = logFC / t,
    CI_lower = logFC - t_crit * SE,
    CI_upper = logFC + t_crit * SE,
    FoldChange = 2^logFC,
    FC_lower = 2^CI_lower,
    FC_upper = 2^CI_upper
  )

write.csv(results, "LimmaResults_SCZvsC.csv", row.names = FALSE)



# ------------------------------------------------------------------------------
# Visualise data
# ------------------------------------------------------------------------------
mat <- log2 (count_mat_scz)
z_mat <- t(scale(t(mat)))

ann_col <- data.frame(Diagnosis = meta_scz$diag)
rownames(ann_col) <- meta_scz$SID
stopifnot(all(rownames(ann_col) == colnames(z_mat)))
ann_colours <- list(
  Diagnosis = c(C = "grey70", SCZ = "firebrick3")
)

# Heatmap: all samples
pheatmap(
  z_mat,
  annotation_col = ann_col,
  annotation_colors = ann_colours,
  show_colnames = FALSE,
  fontsize_row = 7,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "SCZ vs Control – log2 normalised counts (row z-score)",
  filename = "Heapmap_AllSamples_SCZ.png",
  width = 8,
  height = 10
)


# Heatmap: all samples - ordered by group
# Order samples: all C first, then all SCZ
ord <- order(meta_scz$diag)
mat_ord <- mat[, ord]
meta_ord <- meta_scz[ord, ]
# Row z-score for colour scaling
z_mat <- t(scale(t(mat_ord)))

ann_col <- data.frame(Diagnosis = meta_ord$diag)
rownames(ann_col) <- meta_ord$SID
ann_colors <- list(
  Diagnosis = c(C = "grey70", SCZ = "firebrick3")
)

pheatmap(
  z_mat,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  cluster_cols = FALSE,                 # <-- key change
  clustering_distance_rows = "correlation",
  show_colnames = FALSE,
  fontsize_row = 7,
  main = "SCZ vs Control (samples ordered by diagnosis)",
  filename = "Heapmap_AllSamplesReordered_SCZ.png",
  width = 8,
  height = 10
)


# Heatmap: Grouped
# Mean gene expression/group
group_means <- sapply(levels(meta_scz$diag), function(g) {
  cols <- meta_scz$SID[meta_scz$diag == g]
  rowMeans(mat[, cols, drop = FALSE])
})
# z-score across the two group means (per gene)
group_means_z <- t(scale(t(group_means)))

pheatmap(
  group_means_z,
  cluster_cols = FALSE,
  fontsize_row = 7,
  main = "Group means (C vs SCZ) – log2 normalised counts",
  filename = "Heapmap_Grouped_SCZ.png",
  width = 8,
  height = 10
)


# Gene summary
df_long <- as.data.frame(mat) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "SID", values_to = "value") %>%
  left_join(meta_scz, by = "SID")

gene_diff <- df_long %>%
  group_by(gene) %>%
  summarise(
    log2FC = mean(value[diag == "SCZ"]) - mean(value[diag == "C"]),
    .groups = "drop"
  ) %>%
  arrange(log2FC)

ggplot(gene_diff, aes(x = reorder(gene, log2FC), y = log2FC)) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "log2 fold change (SCZ – C)",
    title = "All genes – group difference"
  ) +
  theme_minimal()

