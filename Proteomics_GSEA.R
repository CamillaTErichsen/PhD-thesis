# GSEA through ClusterProfiler

# Load libraries
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(data.table)
library(org.Hs.eg.db)
library(tidyr)
library(stringr)

# Set wd
setwd("C:\\Users\\au532203\\OneDrive - Aarhus universitet\\2. PhD\\Data analysis\\Proteomics\\Normalisation_log2FCOnly - use")


#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
# Load limma data in csv format. Use the 't' column
scz_raw <- read.csv("DEPs_SCZ_ProteinGroups_diag-age-gender-PMI.csv", header = TRUE, check.names = FALSE)

# Remove weak data, keep genenames and t values from limma results
scz <- scz_raw %>%
  transmute(
    GeneName = GeneLabel,
    t = as.numeric(t)
  ) %>%
  filter(!is.na(GeneName), GeneName != "", !is.na(t)) %>%
  group_by(GeneName) %>%
  slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  ungroup()

# Map gene names to entrez IDs
entrez_map <- bitr(
  scz$GeneName, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db) %>%
  distinct(SYMBOL, .keep_all = TRUE)  # remove any duplicates

# Join with t column
entrez_scz <- scz %>%
  left_join(entrez_map, by = c("GeneName" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  ungroup()


# Convert to vector for GSEA and rank t column (highest to lowest)
vect_scz_entrez <- entrez_scz$t
names(vect_scz_entrez) <- entrez_scz$ENTREZID
vect_scz_entrez <- sort(vect_scz_entrez, decreasing = TRUE)
head(vect_scz_entrez)
tail(vect_scz_entrez)
length(vect_scz_entrez)



#-------------------------------------------------------------------------------
# GSEA (ClusterProfiler)
#-------------------------------------------------------------------------------
# Run GSEA from ClusterProfiler
gsea_scz <- gseGO(
  geneList = vect_scz_entrez,      # ranked t values
  OrgDb = org.Hs.eg.db,            # human gene database
  keyType = "ENTREZID",            # type of gene identifier
  ont = "BP",                      # biological process
  pvalueCutoff = 1,             
  verbose = TRUE
)
gsea_results_scz <- as.data.frame(gsea_scz@result)
sig_gsea_scz <- gsea_results_scz %>%    # sign pathways only
  filter(p.adjust < 0.05) %>%
  arrange(desc(NES))

write.csv(as.data.frame(gsea_scz), "GSEA_results_scz-diag-age-gender-PMI.csv", row.names = F)
write.csv(sig_gsea_scz, "GSEA_sig_results_scz-diag-age-gender-PMI.csv", row.names = F)

# Plot top enriched pathways
setwd("C:\\Users\\au532203\\OneDrive - Aarhus universitet\\2. PhD\\Data analysis\\Proteomics\\Normalisation_log2FCOnly - use\\SCZ figures")

# Dotplot
dotplot(gsea_scz, showCategory = 17)
ggsave("DotPlot_20EnrichedPathways_scz-diag-age-gender-PMI.png", width = 12, height = 8)



#-------------------------------------------------------------------------------
# Plot: NES bar plots
#-------------------------------------------------------------------------------
gsea_df <- as.data.frame(gsea_scz@result) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::mutate(
    Description = stringr::str_trunc(Description, 45),
    direction = ifelse(NES > 0, "Up", "Down")
  )

ggplot(gsea_df, aes(
  x = reorder(Description, NES),
  y = NES,
  fill = direction
)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("Up" = "#CD3333", "Down" = "#4F94CD")
  ) +
  labs(
    x = NULL,
    y = "Normalized Enrichment Score (NES)",
    fill = "",
    title = "Top GO:BP pathways (SCZ vs Control)"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  theme_classic()
ggsave("GSEAplot_SCZ-diag-age-gender-PMI.png", width = 11, height = 5)



#-------------------------------------------------------------------------------
# Plot: Individual enrichment plots
#-------------------------------------------------------------------------------
library(enrichplot)

# The some of the top terms (by adjusted p or by |NES|)
top_ids <- as.data.frame(gsea_scz@result) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::pull(ID)

gseaplot2(gsea_scz, geneSetID = top_ids[1], title = top_ids[1])
gseaplot2(gsea_scz, geneSetID = top_ids[2], title = top_ids[2])
gseaplot2(gsea_scz, geneSetID = top_ids[3], title = top_ids[3])



#-------------------------------------------------------------------------------
# Plot: Ridgeplot
#-------------------------------------------------------------------------------
ridgeplot(gsea_scz)
ggsave("GSEARidgePlot_SCZ-diag-age-gender-PMI.png", width = 14, height = 15)



#-------------------------------------------------------------------------------
# Plot: Enrichment map
#-------------------------------------------------------------------------------
library(enrichplot)

# Term similarity
gsea_scz_sim <- pairwise_termsim(gsea_scz)

# Enrichment plot
emapplot(
  gsea_scz_sim,
  showCategory = 26,    # only have nine terms
  layout = "kk"         # can also try fr layout
)
ggsave("GSEAEnrichmentMap_MostEnriched_SCZ-diag-age-gender-PMI.png", width = 15, height = 10)



# ------------------------------------------------------------------------------
# Plot results: Upset plot
# https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
# ------------------------------------------------------------------------------
gsea_df <- as.data.frame(gsea_scz@result)

# Plot (part of the enrichplot library)
upsetplot(gsea_scz)
upsetplot(gsea_scz, showCategory = 17)    # Choose this option if too many gene sets to display
# Save manually with size 800x550



#-------------------------------------------------------------------------------
# Plot: Treemap
#-------------------------------------------------------------------------------
treemap_df <- gsea_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 30) %>%
  mutate(
    neglog10padj = -log10(p.adjust),
    Description = stringr::str_trunc(Description, 40)
  )

ggplot(treemap_df, aes(
  area = setSize,
  fill = NES,
  label = Description,
  subgroup = ID
)) +
  treemapify::geom_treemap() +
  treemapify::geom_treemap_text(reflow = TRUE, place = "centre", min.size = 8) +
  labs(title = "GSEA (GO:BP) Treemap", fill = "NES") +
  theme_minimal()
ggsave("GSEATreemap_MostEnriched_SCZ-diag-age-gender-PMI.png", width = 15, height = 10)



#-------------------------------------------------------------------------------
# Plot: Gene-concept network
#-------------------------------------------------------------------------------
# Plot
cnetplot(
  gsea_scz,
  showCategory = 17,
  circular = FALSE,
  color.params = list(edge = TRUE)
)
ggsave("GSEAGeneConceptNetwork_MostEnriched_SCZ-diag-age-gender-PMI.png", width = 15, height = 10)



#-------------------------------------------------------------------------------
# Extract leading-edge
#-------------------------------------------------------------------------------
# Convert significant pathways to dataframe
df_sig <- as.data.frame(sig_gsea_scz)
head(df_sig)

# Extract leading edge
LEGs_long <- df_sig %>%
  dplyr::select(ID, Description, NES, p.adjust, core_enrichment) %>%
  dplyr::mutate(ENTREZID = str_split(core_enrichment, "/")) %>%
  unnest(ENTREZID) %>%
  dplyr::select(-core_enrichment)
head(LEGs_long)


# Convert entrez ID to gene name
gene_map <- bitr(
  unique(LEGs_long$ENTREZID),
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

LEGs_long_join <- LEGs_long %>%
  left_join(gene_map, by = "ENTREZID")
head(LEGs_long_join)
write.csv(LEGs_long_join, "LEGs_SCZ_AllCovariates.csv", row.names = FALSE)


# Unique LEGs across significant pathways
LEGs_uniques <- LEGs_long_join %>%
  distinct(ENTREZID, SYMBOL) %>%
  arrange(SYMBOL)
head(LEGs_uniques)
write.csv(LEGs_uniques, "LEGsUniques_SCZ_AllCovariates.csv", row.names = FALSE)


# Recurrent LEGs across pathways
LEGs_recurrent <- LEGs_long_join %>%
  count(ENTREZID, Description, SYMBOL, sort = TRUE)
head(LEGs_recurrent)
write.csv(LEGs_recurrent, "LEGsRecurrent_SCZ_AllCovariates.csv", row.names = FALSE)



#-------------------------------------------------------------------------------
# Recurrent LEGs in clusters
#-------------------------------------------------------------------------------
gsea_df <- as.data.frame(sig_gsea_scz)
head(df_sig)

# Define clusters
cluster_map <- tibble::tribble(
  ~Description, ~Cluster,
  # RNA processing and splicing
  "mRNA processing", "RNA processing and splicing",
  "regulation of RNA splicing", "RNA processing and splicing",
  "RNA splicing", "RNA processing and splicing",
  "RNA processing", "RNA processing and splicing",
  "RNA splicing, via transesterification reactions", "RNA processing and splicing",
  "mRNA metabolic process", "RNA processing and splicing",
  
  # Ion transport and calcium/cGMP signalling
  "inorganic ion transmembrane transport", "Ion transport and calcium/cGMP signalling",
  "inorganic cation transmembrane transport", "Ion transport and calcium/cGMP signalling",
  "monoatomic cation transmembrane transport", "Ion transport and calcium/cGMP signalling",
  "monoatomic ion transmembrane transport", "Ion transport and calcium/cGMP signalling",
  "regulation of cytosolic calcium ion concentration", "Ion transport and calcium/cGMP signalling",
  "cGMP-mediated signaling", "Ion transport and calcium/cGMP signalling",
  
  # Intermediate filament / cytoskeletal organisation
  "keratinocyte differentiation", "Intermediate filament and cytoskeletal organisation",
  "epidermal cell differentiation", "Intermediate filament and cytoskeletal organisation",
  "intermediate filament-based process", "Intermediate filament and cytoskeletal organisation",
  "intermediate filament cytoskeleton organisation", "Intermediate filament and cytoskeletal organisation"
)


# Add cluster labels to GSEA results
gsea_clustered <- gsea_df %>%
  left_join(cluster_map, by = "Description") %>%
  filter(!is.na(Cluster))


# Extract core proteins
core_long <- gsea_clustered %>%
  dplyr::select(ID, Description, Cluster, NES, p.adjust, core_enrichment) %>%
  dplyr::mutate(ENTREZID = stringr::str_split(core_enrichment, "/")) %>%
  tidyr::unnest(ENTREZID) %>%
  dplyr::select(-core_enrichment)


# Convert IDs to symbols
gene_map <- bitr(
  unique(core_long$ENTREZID),
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

core_long <- core_long %>%
  left_join(gene_map, by = "ENTREZID")


# Count recurrent core genes within each cluster
core_recurrent_by_cluster <- core_long %>%
  distinct(Cluster, Description, ENTREZID, SYMBOL) %>%
  count(Cluster, ENTREZID, SYMBOL, name = "n_pathways") %>%
  arrange(Cluster, desc(n_pathways), SYMBOL)
head(core_recurrent_by_cluster)


# Save recurrent LEGs
write.csv(core_long, "LEGs_GenesPathway_SCZ.csv", row.names = FALSE)
write.csv(core_recurrent_by_cluster, "LEGs_Cluster_SCZ.csv", row.names = FALSE)

