# GSEA through ClusterProfiler

# Load libraries
library(tidyverse)
library(clusterProfiler)     
library(data.table)
library(simplifyEnrichment)
library(enrichplot)
library(GO.db)
library(org.Hs.eg.db)
library(msigdbr)            

# set working directory
setwd("C:\\Users\\au532203\\OneDrive - Aarhus universitet\\2. PhD\\Data analysis\\GSEA_simplified")



#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
# Load deseq data in csv format. Use the 'stat' column, ranked from lowest to
# highest. Include all genes
scz <- read.table("DESeq2_SCZ_age-RIN-gender.csv", sep = ",", header = T) %>%
  dplyr::select(GeneName = X, stat)

# Convert gene names to Entrez IDs
entrez_scz <- bitr(scz$GeneName, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
  distinct(SYMBOL, .keep_all = TRUE)  # remove any duplicates

# Merge Entrez IDs with stat column
scz_entrez <- merge(scz, entrez_scz, by.x = "GeneName", by.y = "SYMBOL")
vect_scz_entrez <- setNames(scz_entrez$stat, scz_entrez$ENTREZID) %>%
  na.omit() %>%
  sort(decreasing = TRUE)
vect_scz_entrez



#-------------------------------------------------------------------------------
# GSEA (ClusterProfiler)
#-------------------------------------------------------------------------------
# Run GSEA from ClusterProfiler
gsea_scz <- gseGO(
  geneList = vect_scz_entrez,      # ranked stat values
  OrgDb = org.Hs.eg.db,            # human gene database
  keyType = "ENTREZID",            # type of gene identifier
  ont = "BP",                      # biological process
  pvalueCutoff = 0.05,         
  verbose = TRUE
)

gsea_results_scz <- as.data.frame(gsea_scz@result)
sig_gsea_scz <- gsea_results_scz %>%    # sign pathways only
  filter(p.adjust < 0.05) %>%
  arrange(desc(NES))

write.csv(as.data.frame(gsea_scz), "GSEA_results_scz.csv", row.names = F)
write.csv(sig_gsea_scz, "GSEA_sig_results_scz.csv", row.names = F)


# Plot top 10 enriched pathways
dotplot(gsea_scz, showCategory = 20)
ggsave("DotPlot_20EnrichedPathways_scz.png", width = 12, height = 8)



#-------------------------------------------------------------------------------
# Extract leading-edge genes (LEGs)
#-------------------------------------------------------------------------------
# Convert results to df
df_res <- as.data.frame(gsea_results_scz)
head(df_res)

# Cluster 1 (cognition and plasticity): GO:0048167 – regulation of synaptic plasticity, GO:0050808 – synapse organization, GO:0050890 – cognition
# Cluster 2 (neurotransmitter release and vesicles): GO:0016079 – synaptic vesicle exocytosis, GO:0007269 – neurotransmitter secretion
# Cluster 3 (excitability and ion transport): GO:0071805 – potassium ion transmembrane transport, GO:0042391 – regulation of membrane potential
# Are all these GO terms present in the significant data? --> yes
df_res[df_res$ID %in% c("GO:0050890",
                        "GO:0048167",
                        "GO:0050808",
                        "GO:0016079",
                        "GO:0007269",
                        "GO:0042391",
                        "GO:0071805"), ]
# Gets LEGs
get_LEGs <- function(df, go_id){
  x <- df$core_enrichment[df$ID == go_id]
  unlist(strsplit(x, "/"))
}
# Convert IDs to gene names
library(org.Hs.eg.db)
library(AnnotationDbi)
convert_symbols <- function(entrez){
  syms <- mapIds(org.Hs.eg.db,
                 keys = entrez,
                 keytype = "ENTREZID",
                 column = "SYMBOL",
                 multiVals = "first")
  na.omit(unique(as.character(syms)))
}


# CLUSTER 1
legs_cluster1 <- list(
  plasticity = convert_symbols(get_LEGs(df_res, "GO:0048167")),
  organisation  = convert_symbols(get_LEGs(df_res, "GO:0050808")),
  cognition = convert_symbols((get_LEGs(df_res, "GO:0050890")))
)
# Overlap
core_cluster1 <- Reduce(intersect, legs_cluster1)
# Frequency ranking
freq_cluster1 <- sort(table(unlist(legs_cluster1)), decreasing = TRUE)
head(freq_cluster1, 20)


# CLUSTER 2
legs_cluster2 <- list(
  exocytosis = convert_symbols(get_LEGs(df_res, "GO:0016079")),
  secretion  = convert_symbols(get_LEGs(df_res, "GO:0007269"))
)
core_cluster2 <- Reduce(intersect, legs_cluster2)
freq_cluster2 <- sort(table(unlist(legs_cluster2)), decreasing = TRUE)
head(freq_cluster2, 20)


# CLUSTER 3
legs_cluster3 <- list(
  transport = convert_symbols(get_LEGs(df_res, "GO:0071805")),
  potential  = convert_symbols(get_LEGs(df_res, "GO:0042391"))
)
core_cluster3 <- Reduce(intersect, legs_cluster3)
freq_cluster3 <- sort(table(unlist(legs_cluster3)), decreasing = TRUE)
head(freq_cluster3, 20)

#----------------
# Make LEGs tables
library(tibble)
library(stringr)
# --- Helper functions ---
get_LE_entrez <- function(df, go_id){
  x <- df$core_enrichment[df$ID == go_id]
  if (length(x) == 0) stop(paste("GO ID not found:", go_id))
  unlist(strsplit(x, "/"))
}

entrez_to_symbol <- function(entrez_ids){
  syms <- mapIds(
    org.Hs.eg.db,
    keys = entrez_ids,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  syms <- na.omit(unique(as.character(syms)))
  syms
}

# Representative GO terms
rep_go <- c(
  "GO:0050890", # cognition
  "GO:0048167", # regulation of synaptic plasticity
  "GO:0050808", # synapse organization
  "GO:0016079", # synaptic vesicle exocytosis
  "GO:0007269", # neurotransmitter secretion
  "GO:0042391", # regulation of membrane potential
  "GO:0071805"  # potassium ion transmembrane transport
)


# Table for each representative GO term
RepresLEGs <- df_res %>%
  dplyr::filter(ID %in% rep_go) %>%
  dplyr::select(GO_ID = ID,
         GO_Description = Description,
         NES,
         Adjusted_p_value = p.adjust,
         core_enrichment) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    LEGs_entrez = list(get_LE_entrez(df_res, GO_ID)),
    LEGs_symbol = list(entrez_to_symbol(LEGs_entrez)),
    LEGs = paste(LEGs_symbol, collapse = ", ")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(GO_ID, GO_Description, NES, Adjusted_p_value, LEGs) %>%
  dplyr::arrange(Adjusted_p_value, desc(abs(NES)))
# View in console
print(RepresLEGs, n = Inf)
# Save to file
write.csv(RepresLEGs, "SupplementaryTable__RepresLEGs_SCZvsC.csv", row.names = FALSE)


# Core recurrent LEGs per functional cluster
cluster_map <- list(
  Cluster1_Cognition_Plasticity_SynapseOrg = c("GO:0050890", "GO:0048167", "GO:0050808"),
  Cluster2_Vesicle_Release                 = c("GO:0016079", "GO:0007269"),
  Cluster3_Excitability_IonTransport       = c("GO:0042391", "GO:0071805")
)

# Build a long table of (cluster, GO term, gene)
cluster_gene_long <- bind_rows(lapply(names(cluster_map), function(cl){
  terms <- cluster_map[[cl]]
  bind_rows(lapply(terms, function(go){
    genes <- entrez_to_symbol(get_LE_entrez(df_res, go))
    tibble(
      Cluster = cl,
      GO_ID = go,
      Gene = genes
    )
  }))
}))

# Summarise recurrence within each cluster
sup_table_Y <- cluster_gene_long %>%
  distinct(Cluster, GO_ID, Gene) %>%   # avoid duplicates within a term
  group_by(Cluster, Gene) %>%
  summarise(
    n_terms_present = n_distinct(GO_ID),
    Terms_present = paste(sort(unique(GO_ID)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Cluster, desc(n_terms_present), Gene)

# Keep only "core" genes (present in all selected terms in that cluster)
cluster_sizes <- tibble(
  Cluster = names(cluster_map),
  n_terms_total = sapply(cluster_map, length)
)

sup_table_Y_core_only <- sup_table_Y %>%
  left_join(cluster_sizes, by = "Cluster") %>%
  dplyr::filter(n_terms_present == n_terms_total) %>%
  dplyr::select(Cluster, Gene, n_terms_present, Terms_present)
# Save both (full recurrence + core-only)
write.csv(sup_table_Y, "SupplementaryTable_RecurrentLEGsCluster_SCZvsC.csv", row.names = FALSE)
write.csv(sup_table_Y_core_only, "SupplementaryTable_CoreLEGsCluster_SCZvsC.csv", row.names = FALSE)

# Quick sanity check: top recurrent genes per cluster
top_per_cluster <- sup_table_Y %>%
  group_by(Cluster) %>%
  slice_max(order_by = n_terms_present, n = 20, with_ties = TRUE) %>%
  ungroup()

print(top_per_cluster, n = 60)



#-------------------------------------------------------------------------------
# Simplify GSEA results
#-------------------------------------------------------------------------------
# Too many significant pathways
# Extract GO terms with p.adjust < 0.05
sig_terms_scz <- gsea_scz@result %>%
  filter(p.adjust < 0.05) %>%
  pull(ID)  # Extract GO term IDs


# Compute GO similarity matrix
simMatrix_scz <- GO_similarity(sig_terms_scz, ont = "BP")
print(dim(simMatrix_scz))  # Should be square (952x952)

pdf("SimilarityMatrix_scz-NoDrugs.pdf", width = 9, height = 5)
simplifyGO(simMatrix_scz, plot = TRUE)
dev.off()


# Cluster GO terms (for further simplification plots): binary cut
simMatrix_clusters_scz <- simplifyGO(simMatrix_scz, plot = FALSE)

# Extract clusters
terms_cluster_scz <- simMatrix_clusters_scz %>%
  dplyr::rename(GO_ID = id) %>%
  dplyr::select(GO_ID, cluster)


# Function to find general term descriptors
terms_gen <- function(go_ids) {
  terms_info <- Term(go_ids)
  return(names(which.min(nchar(terms_info))))
}

# Representative term for each cluster
terms_repres_scz <- terms_cluster_scz %>%
  group_by(cluster) %>%
  summarise(representative_term = terms_gen(GO_ID))

# Add GO term lvls to the similarity matrices
terms_cluster_scz <- terms_cluster_scz %>%
  left_join(terms_repres_scz, by = "cluster") %>%
  mutate(cluster = as.character(cluster))



#-------------------------------------------------------------------------------
# Plot simplified terms + NES
#-------------------------------------------------------------------------------
# Map GO term IDs to names
GO_term_map <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, gs_exact_source) %>%
  distinct()  # Keep unique mappings of gene set names to GO IDs

# Merge the cluster mapping with the GSEA results (display average NES/cluster)
merged_gsea_terms <- gsea_scz@result %>%
  dplyr::select(ID, NES, Description, p.adjust) %>%  # Select necessary columns, including p.adjust
  filter(p.adjust < 0.05) %>%  # Filter for significant results (adjusted p-value < 0.05)
  left_join(terms_cluster_scz, by = c("ID" = "GO_ID"))  # Merge based on GO term ID

# Calculate the mean NES for each cluster
cluster_gsea_scz <- merged_gsea_terms %>%
  group_by(cluster, representative_term) %>%
  summarise(mean_NES = mean(NES, na.rm = TRUE)) %>%   # Calculate mean NES for each cluster
  arrange(desc(mean_NES))  # Sort by mean NES value

# Merge GO map with GO IDs to add pathway names
gsea_simple_scz <- inner_join(cluster_gsea_scz, GO_term_map, by = c("representative_term" = "gs_exact_source"))
write.csv(gsea_simple_scz, "GSEA_simplified_scz.csv", row.names = TRUE)

# Plot clusters and mean NES
gsea_simple_scz %>%
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
       title = "Simplified GSEA: Schizophrenia vs control") +
  guides(fill = guide_legend(title = "NES Direction"))
ggsave("GSEAplot_SCZ.png", width = 11, height = 5)



#-------------------------------------------------------------------------------
# Plot: Enrichment map
#-------------------------------------------------------------------------------
library(enrichplot)

# Too many data points (go terms): show only strong enrichment: NES > 2 or < -2
strong_terms <- gsea_scz@result %>%
  filter(abs(NES) > 2 & p.adjust < 0.01)
filtered_gsea <- gsea_scz
filtered_gsea@result <- strong_terms

emapplot(pairwise_termsim(filtered_gsea), showCategory = 25, layout = "kk")
ggsave("GSEAEnrichmentMap_MostEnriched_SCZ.png", width = 15, height = 10)


# Plot only representative terms
rep_terms <- terms_repres_scz$representative_term

# Filter GSEA result to representative terms
rep_results <- gsea_scz@result %>%
  filter(ID %in% rep_terms)
gsea_simple <- gsea_scz
gsea_simple@result <- rep_results

emapplot(pairwise_termsim(gsea_simple), showCategory = nrow(rep_results), layout = "kk")
ggsave("GSEAEnrichmentMap_RepresentativeTerms_SCZ.png", width = 15, height = 10)



# ------------------------------------------------------------------------------
# Plot results: Upset plot
# https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
# ------------------------------------------------------------------------------
library(ComplexUpset)

# Plot the most enriched terms
upsetplot(filtered_gsea)
# Save manually with size 800x550



# ------------------------------------------------------------------------------
# Plot results: Treemap
# ------------------------------------------------------------------------------
library(treemapify)

gsea_results_df <- as.data.frame(gsea_scz@result)

treemap_data <- merged_gsea_terms %>%
  group_by(cluster, representative_term) %>%
  summarise(mean_NES = mean(NES), .groups = "drop") %>%
  left_join(
    gsea_results_df %>% dplyr::select(ID, Description),
    by = c("representative_term" = "ID")
  )


ggplot(treemap_data, aes(area = abs(mean_NES), fill = mean_NES, label = Description)) +
  geom_treemap() +
  geom_treemap_text(colour = "white", place = "centre", grow = TRUE) +
  scale_fill_gradient2(low = "#4F94CD", mid = "white", high = "#CD3333", midpoint = 0) +
  theme_void() +
  labs(title = "Simplified GSEA Treemap by Cluster")
ggsave("Treemap_scz.png", width = 10, height = 8)

