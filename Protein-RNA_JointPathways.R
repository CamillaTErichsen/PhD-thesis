# Joint pathway analyses between proteomics and sequencing data

# Load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggrepel)

setwd("C:\\Users\\au532203\\OneDrive - Aarhus universitet\\2. PhD\\Data analysis\\Joint-pathway analyses")


# ------------------------------------------------------------------------------
# Load pathway data
# ------------------------------------------------------------------------------
# Load pathway data for sequencing and proteomics.
prot <- read.csv("GSEA-proteomics_scz.csv", header = T, check.names = FALSE)
# Use the non-simplified sequencing pathways
seq <- read.csv("GSEA-sequencing_all_scz.csv", header = T, check.names = FALSE)

# Filter columns and rename
prot_gsea <- prot %>%
  dplyr::select(
    ID,
    Description,
    NES_protein = NES,
    pvalue_protein = pvalue,
    padj_protein = p.adjust,
    qvalue_protein = qvalue
  )

seq_gsea <- seq %>%
  dplyr::select(
    ID,
    Description,
    NES_RNA = NES,
    pvalue_RNA = pvalue,
    padj_RNA = p.adjust,
    qvalue_RNA = qvalue
  )



# ------------------------------------------------------------------------------
# Join pathways
# ------------------------------------------------------------------------------
# Join by GO ID - the common ID
join_gsea <- seq_gsea %>%
  dplyr::full_join(
    prot_gsea,
    by = "ID",
    suffix = c("_RNA_desc", "_protein_desc")
  ) %>%
  dplyr::mutate(
    Description = dplyr::coalesce(
      Description_RNA_desc, Description_protein_desc
    ),
    
    # Define significance
    sig_RNA = !is.na(padj_RNA) & padj_RNA < 0.05,
    sig_protein = !is.na(padj_protein) & padj_protein < 0.05,
    
    # Define directions
    direction_RNA = dplyr::case_when(
      NES_RNA > 0 ~ "Positive",
      NES_RNA < 0 ~ "Negative",
      TRUE ~ NA_character_
    ),
    direction_protein = dplyr::case_when(
      NES_protein > 0 ~ "Positive",
      NES_protein < 0 ~ "Negative",
      TRUE ~ NA_character_
    ),
    
    # Concordance
    concordance = dplyr::case_when(
      sig_RNA & sig_protein & sign(NES_RNA) == sign(NES_protein) ~ "Significant in both, same direction",
      sig_RNA & sig_protein & sign(NES_RNA) != sign(NES_protein) ~ "Significant in both, opposite direction",
      sig_RNA & !sig_protein ~ "RNA only",
      !sig_RNA & sig_protein ~ "Protein only",
      TRUE ~ "Not significant"
    )
  )

# Save joined table
write.csv(join_gsea, "Joint_GSEA_SCZ.csv", row.names = FALSE)



# ------------------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------------------
table(join_gsea$concordance)

join_sig_summary <- join_gsea %>%
  dplyr::filter(sig_RNA | sig_protein) %>%
  dplyr::count(concordance)
join_sig_summary

write.csv(join_sig_summary, "Joint_GSEA-summary_SCZ.csv", row.names = FALSE)



# ------------------------------------------------------------------------------
# Visualisation: NES scatterplot - plot only sign
# ------------------------------------------------------------------------------
# RNA vs protein: only sign pathways in either omics layer
joint_sig <- join_gsea %>%
  dplyr::filter(sig_RNA | sig_protein) %>%
  dplyr::mutate(
    is_shared_same = concordance == "Significant in both, same direction",
    point_size = ifelse(is_shared_same, 4.5, 1.8),
    point_alpha = ifelse(is_shared_same, 1, 0.45),
    label = ifelse(is_shared_same, Description, NA)
  )

# Plot
ggplot(joint_sig, aes(x = NES_RNA, y = NES_protein)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  
  # Background: RNA-only / protein-only
  geom_point(
    data = joint_sig %>% filter(!is_shared_same),
    aes(colour = concordance),
    size = 1.5,
    alpha = 0.45
  ) +
  
  # Foreground: shared same-direction pathways
  geom_point(
    data = joint_sig %>% filter(is_shared_same),
    aes(colour = concordance),
    size = 3,
    alpha = 1
  ) +
  
  # Labels only for shared same-direction pathways
  ggrepel::geom_text_repel(
    data = joint_sig %>% filter(is_shared_same),
    aes(label = Description),
    colour = "#D81B60",
    size = 3.2,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.25,
    max.overlaps = Inf,
    segment.colour = "#D81B60"
  ) +
  
  scale_colour_manual(
    values = c(
      "RNA only" = "#009E73",
      "Protein only" = "orange",
      "Significant in both, same direction" = "#D81B60",
      "Significant in both, opposite direction" = "purple"
    )
  ) +
  
  labs(
    title = "SCZ vs Control: RNA–protein pathway concordance",
    subtitle = "Pathways significant in at least one omics layer",
    x = "RNA-seq NES",
    y = "Proteomics NES",
    colour = ""
  ) +
  theme_classic()

ggsave("Joint_GSEA_NES_scatter-SigOnly_scz.pdf", width = 8, height = 6)



# ------------------------------------------------------------------------------
# Correlate NES values in shared pathways
# ------------------------------------------------------------------------------
shared <- join_gsea %>%
  filter(!is.na(NES_RNA),
         !is.na(NES_protein))

cor(shared$NES_RNA,
    shared$NES_protein,
    method = "spearman")
# Correlation = 0.1078398.



# ------------------------------------------------------------------------------
# Table
# ------------------------------------------------------------------------------
# Pathways significant in either omics layer
joint_sig <- join_gsea %>%
  dplyr::filter(sig_RNA | sig_protein) %>%
  dplyr::arrange(concordance, padj_RNA, padj_protein)

write.csv(joint_sig, "Joint_GSEA_EitherSig_scz.csv", row.names = FALSE)

