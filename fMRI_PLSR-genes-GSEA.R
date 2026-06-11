# GSEA: PLSR genes for left and right dlPFC

# Load libraries
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(data.table)
library(org.Hs.eg.db)
library(ggplot2)
library(simplifyEnrichment)

# Set working directory
setwd("C:/Users/au532203/OneDrive - Aarhus universitet/Exchange/Codes/FC-genes/GSEA")


#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
# OBS! Data should be ranked highest to lowest zscore
data_L <- read.table("Genes_pass_zscore_all_L.csv", sep=",", header=T)   # 15633
data_R <- read.table("Genes_pass_zscore_all_R.csv", sep=",", header=T)
is.data.frame(data_L)     # Data is stored in a df
head(data_L)              # Data is divided into two columns: GeneName, zscore

# Convert gene names to Entrez IDs
entrez_L <- bitr(data_L$GeneName, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%   # 14973
  distinct(SYMBOL, .keep_all = TRUE)  # remove any duplicates

entrez_R <- bitr(data_R$GeneName, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
  distinct(SYMBOL, .keep_all = TRUE)  # remove any duplicates

# Merge Entrez IDs with stat column
L_entrez <- merge(data_L, entrez_L, by.x = "GeneName", by.y = "SYMBOL")    # 14973
R_entrez <- merge(data_R, entrez_R, by.x = "GeneName", by.y = "SYMBOL")
# Rank based on zscores
vect_L_entrez <- setNames(L_entrez$zscore, L_entrez$ENTREZID) %>%          # 14973
  na.omit() %>%
  sort(decreasing = TRUE)
vect_R_entrez <- setNames(R_entrez$zscore, R_entrez$ENTREZID) %>%
  na.omit() %>%
  sort(decreasing = TRUE)



#-------------------------------------------------------------------------------
# GSEA (ClusterProfiler): Left dlPFC
#-------------------------------------------------------------------------------
# Run GSEA from ClusterProfiler
gsea_L <- gseGO(
  geneList = vect_L_entrez,      # ranked stat values
  OrgDb = org.Hs.eg.db,            # human gene database
  keyType = "ENTREZID",            # type of gene identifier
  ont = "BP",                      # biological processes
  pvalueCutoff = 0.05,             # 
  verbose = TRUE
)

gsea_results_L <- as.data.frame(gsea_L@result)
sig_gsea_L <- gsea_results_L %>%    # sign pathways only
  filter(p.adjust < 0.05) %>%
  arrange(desc(NES))

# Save
write.csv(as.data.frame(gsea_L), "GSEA_results_L.csv", row.names = F)
write.csv(sig_gsea_L, "GSEA_sig_results_L.csv", row.names = F)



#-------------------------------------------------------------------------------
# GSEA (ClusterProfiler): Right dlPFC
#-------------------------------------------------------------------------------
# Run GSEA from ClusterProfiler
gsea_R <- gseGO(
  geneList = vect_R_entrez,      # ranked stat values
  OrgDb = org.Hs.eg.db,            # human gene database
  keyType = "ENTREZID",            # type of gene identifier
  ont = "BP",                      # biological processes
  pvalueCutoff = 0.05,             # 
  verbose = TRUE
)

gsea_results_R <- as.data.frame(gsea_R@result)
sig_gsea_R <- gsea_results_R %>%    # sign pathways only
  filter(p.adjust < 0.05) %>%
  arrange(desc(NES))

# Save
write.csv(as.data.frame(gsea_R), "GSEA_results_R.csv", row.names = F)
write.csv(sig_gsea_R, "GSEA_sig_results_R.csv", row.names = F)



#-------------------------------------------------------------------------------
# Simplify results: left dlPFC
#-------------------------------------------------------------------------------
# Too many significant pathways
# Extract GO terms with p.adjust < 0.05
sig_terms_L <- gsea_L@result %>%
  filter(p.adjust < 0.05) %>%
  pull(ID)  # Extract GO term IDs

# Compute GO similarity matrix
simMatrix_L <- GO_similarity(sig_terms_L, ont = "BP")
print(dim(simMatrix_L))  # Should be square (120x120)
pdf("SimilarityMatrix_L.pdf", width = 9, height = 5)   # 13 terms
simplifyGO(simMatrix_L, plot = TRUE)
dev.off()

# Cluster GO terms (for further simplification plots): binary cut
simMatrix_clusters_L <- simplifyGO(simMatrix_L, plot = FALSE)

# Extract clusters
terms_cluster_L <- simMatrix_clusters_L %>%
  dplyr::rename(GO_ID = id) %>%
  dplyr::select(GO_ID, cluster)

# Function to find general term descriptors
terms_gen <- function(go_ids) {
  terms_info <- Term(go_ids)
  return(names(which.min(nchar(terms_info))))
}

# Representative term for each cluster
terms_repres_L <- terms_cluster_L %>%
  group_by(cluster) %>%
  summarise(representative_term = terms_gen(GO_ID))

# Add GO term lvls to the similarity matrices
terms_cluster_L <- terms_cluster_L %>%
  left_join(terms_repres_L, by = "cluster") %>%
  mutate(cluster = as.character(cluster))



#-------------------------------------------------------------------------------
# Plot sign. simplified pathways + NES: left dlPFC
#-------------------------------------------------------------------------------
# Map GO term IDs to names
GO_term_map <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, gs_exact_source) %>%
  distinct()  # Keep unique mappings of gene set names to GO IDs

# Merge the cluster mapping with the GSEA results (display average NES/cluster)
merged_gsea_terms_L <- gsea_L@result %>%
  dplyr::select(ID, NES, Description, p.adjust) %>%  # Select necessary columns, including p.adjust
  filter(p.adjust < 0.05) %>%  # Filter for significant results (adjusted p-value < 0.05)
  left_join(terms_cluster_L, by = c("ID" = "GO_ID"))  # Merge based on GO term ID

# Calculate the mean NES for each cluster
cluster_gsea_L <- merged_gsea_terms_L %>%
  group_by(cluster, representative_term) %>%
  summarise(mean_NES = mean(NES, na.rm = TRUE)) %>%   # Calculate mean NES for each cluster
  arrange(desc(mean_NES))  # Sort by mean NES value

# Merge GO map with GO IDs to add pathway names
gsea_simple_L <- inner_join(cluster_gsea_L, GO_term_map, by = c("representative_term" = "gs_exact_source"))
write.csv(gsea_simple_L, "GSEA_simplified_L.csv", row.names = TRUE)

# Plot clusters and mean NES
gsea_simple_L %>%
  mutate(gs_name = gsub("GOBP_", "", gs_name),
         gs_name = gsub("_", " ", gs_name)) %>%
  ggplot(aes(x = reorder(gs_name, mean_NES), y = mean_NES, fill = mean_NES > 0)) + 
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#CD3333", "FALSE" = "#4F94CD")) + 
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black")) +
  coord_flip() +
  labs(y = "Mean normalised enrichment score",
       x = "Simplified GO terms",
       title = "Simplified GSEA") +
  guides(fill = guide_legend(title = "NES Direction"))
ggsave("GSEAplot_L.png", width = 11, height = 5)



#-------------------------------------------------------------------------------
# Simplify results: right dlPFC
#-------------------------------------------------------------------------------
# Too many significant pathways
# Extract GO terms with p.adjust < 0.05
sig_terms_R <- gsea_R@result %>%
  filter(p.adjust < 0.05) %>%
  pull(ID)  # Extract GO term IDs

# Compute GO similarity matrix
simMatrix_R <- GO_similarity(sig_terms_R, ont = "BP")
print(dim(simMatrix_R))  # Should be square (120x120)
pdf("SimilarityMatrix_R.pdf", width = 9, height = 5)   # 9 terms
simplifyGO(simMatrix_R, plot = TRUE)
dev.off()

# Cluster GO terms (for further simplification plots): binary cut
simMatrix_clusters_R <- simplifyGO(simMatrix_R, plot = FALSE)

# Extract clusters
terms_cluster_R <- simMatrix_clusters_R %>%
  dplyr::rename(GO_ID = id) %>%
  dplyr::select(GO_ID, cluster)

# Representative term for each cluster
terms_repres_R <- terms_cluster_R %>%
  group_by(cluster) %>%
  summarise(representative_term = terms_gen(GO_ID))

# Add GO term lvls to the similarity matrices
terms_cluster_R <- terms_cluster_R %>%
  left_join(terms_repres_R, by = "cluster") %>%
  mutate(cluster = as.character(cluster))



#-------------------------------------------------------------------------------
# Plot sign. simplified pathways + NES: right dlPFC
#-------------------------------------------------------------------------------
# Merge the cluster mapping with the GSEA results (display average NES/cluster)
merged_gsea_terms_R <- gsea_R@result %>%
  dplyr::select(ID, NES, Description, p.adjust) %>%  # Select necessary columns, including p.adjust
  filter(p.adjust < 0.05) %>%  # Filter for significant results (adjusted p-value < 0.05)
  left_join(terms_cluster_R, by = c("ID" = "GO_ID"))  # Merge based on GO term ID

# Calculate the mean NES for each cluster
cluster_gsea_R <- merged_gsea_terms_R %>%
  group_by(cluster, representative_term) %>%
  summarise(mean_NES = mean(NES, na.rm = TRUE)) %>%   # Calculate mean NES for each cluster
  arrange(desc(mean_NES))  # Sort by mean NES value

# Merge GO map with GO IDs to add pathway names
gsea_simple_R <- inner_join(cluster_gsea_R, GO_term_map, by = c("representative_term" = "gs_exact_source"))
write.csv(gsea_simple_R, "GSEA_simplified_R.csv", row.names = TRUE)

# Plot clusters and mean NES
gsea_simple_R %>%
  mutate(gs_name = gsub("GOBP_", "", gs_name),
         gs_name = gsub("_", " ", gs_name)) %>%
  ggplot(aes(x = reorder(gs_name, mean_NES), y = mean_NES, fill = mean_NES > 0)) + 
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#CD3333", "FALSE" = "#4F94CD")) + 
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black")) +
  coord_flip() +
  labs(y = "Mean normalised enrichment score",
       x = "Simplified GO terms",
       title = "Simplified GSEA") +
  guides(fill = guide_legend(title = "NES Direction"))
ggsave("GSEAplot_R.png", width = 11, height = 5)

