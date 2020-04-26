#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Find markers of each tissue relative to all other tissues and perform enrichment analysis
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data")
  place <- "wolfpack"
}

# Libraries
library('Seurat')
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('clusterProfiler')
library('ReactomePA')
library('org.Hs.eg.db')
library('msigdbr')
library('gtools')


# Load prefiltered and clustered Seurat Object
if(place == "local") {
  filtered_exp <- readRDS("Prefiltered_experiment_practice_seurat_integrated_markers_by_cluster_tissueVStissue.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("Prefiltered_experiment_all_seurat_integrated_markers_by_cluster_tissueVStissue.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

Idents(object = filtered_exp) <- "Tissue"
DefaultAssay(object = filtered_exp) <- "SCT"

# Find markers for all clusters and perform GO enrichment
# --------------------------------------------------------------------------
all_markers <- FindAllMarkers(filtered_exp, assay = "SCT", test.use = "LR", slot = "data", latent.vars = c("Replicate", "Lib_size", "S.Score", "G2M.Score", "digest_stress1"))
filtered_exp@misc$tissueVSall_markers <- all_markers

# Make gene universe(s)
universe <- as.character(unique(filtered_exp@assays$RNA@meta.features[filtered_exp@assays$RNA@meta.features$Any_Active, "Ensembl"]))
universe_entrez <- as.character(unique(filtered_exp@assays$RNA@meta.features[filtered_exp@assays$RNA@meta.features$Any_Active, "EntrezID"]))
universe_genesymbol <- as.character(unique(filtered_exp@assays$RNA@meta.features[filtered_exp@assays$RNA@meta.features$Any_Active, "GeneSymbol"]))

# Add ENSEMBL gene ids to Marker dataframe
idx <- match(all_markers$gene, filtered_exp@assays$RNA@meta.features$Row.names)
all_markers$Ensembl <- filtered_exp@assays$RNA@meta.features$Ensembl [idx]
all_markers$EntrezID <- filtered_exp@assays$RNA@meta.features$EntrezID [idx]
all_markers$GeneSymbol <- filtered_exp@assays$RNA@meta.features$GeneSymbol [idx]

# Split into clusters and find GO annotations all, up or down regulated
all_markers$Upregulated <- all_markers$avg_logFC > 0
all_markers$Downregulated <- all_markers$avg_logFC < 0
markers.filtered_exp <- split(all_markers, all_markers$Tissue)

# Makes hallmark geneset for enrichment testing
h_df <- msigdbr(species = "Homo sapiens", category = "H")
h_t2g <- h_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()

# All-regulated markers in each cluster
for(i in names(markers.filtered_exp)){
  dir.create(paste0("markers/Tissue_", i))
  interesting <- markers.filtered_exp[[i]]
  write.csv(interesting, paste0("markers/Tissue_", i,"/Tissue_", i, "_allregulated_markers_tissueVSall.csv"))
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/Tissue_", i,"/GOBP_markers_all_Tissue_", i, "_tissueVSall.csv"))
  
  MFenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("markers/Tissue_", i,"/GOMF_markers_all_Tissue_", i, "_tissueVSall.csv"))
  
  CCenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("markers/Tissue_", i,"/GOCC_markers_all_Tissue_", i, "_tissueVSall.csv"))

  # Perform HALLMARK enrichment analysis
  HallmarkEnrich <- enricher(gene = interesting[interesting$p_val_adj < 0.05, "GeneSymbol"], 
                             TERM2GENE = h_t2g,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             pAdjustMethod = "bonferroni",
                             minGSSize = 15,
                             maxGSSize = 500,
                             universe = universe_genesymbol)
  Hallmark_GOI <- as.data.frame(HallmarkEnrich@result)
  write.csv(Hallmark_GOI, paste0("markers/Tissue_", i,"/HALLMARK_markers_all_Tissue_", i, "_tissueVSall.csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/Tissue_", i,"/KEGG_markers_all_Tissue_", i, "_tissueVSall.csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/Tissue_", i,"/Reactome_markers_all_Tissue_", i, "_tissueVSall.csv"))
}

# Up-regulated markers in each cluster
for(i in names(markers.filtered_exp)){
  interesting <- markers.filtered_exp[[i]]
  interesting <- interesting[interesting$Upregulated, ]
  write.csv(interesting, paste0("markers/Tissue_", i,"/Tissue_", i, "_upregulated_markers_tissueVSall.csv"))
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/Tissue_", i,"/GOBP_markers_up_Tissue_", i, "_tissueVSall.csv"))
  
  MFenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("markers/Tissue_", i,"/GOMF_markers_up_Tissue_", i, "_tissueVSall.csv"))
  
  CCenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("markers/Tissue_", i,"/GOCC_markers_up_Tissue_", i, "_tissueVSall.csv"))

  # Perform HALLMARK enrichment analysis
  HallmarkEnrich <- enricher(gene = interesting[interesting$p_val_adj < 0.05, "GeneSymbol"], 
                             TERM2GENE = h_t2g,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             pAdjustMethod = "bonferroni",
                             minGSSize = 15,
                             maxGSSize = 500,
                             universe = universe_genesymbol)
  Hallmark_GOI <- as.data.frame(HallmarkEnrich@result)
  write.csv(Hallmark_GOI, paste0("markers/Tissue_", i,"/HALLMARK_markers_up_Tissue_", i, "_tissueVSall.csv"))
    
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/Tissue_", i,"/KEGG_markers_up_Tissue_", i, "_tissueVSall.csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/Tissue_", i,"/Reactome_markers_up_Tissue_", i, "_tissueVSall.csv"))
}

# Down-regulated markers in each cluster
for(i in names(markers.filtered_exp)){
  interesting <- markers.filtered_exp[[i]]
  interesting <- interesting[interesting$Downregulated, ]
  write.csv(interesting, paste0("markers/Tissue_", i,"/Tissue_", i, "_downregulated_markers_tissueVSall.csv"))
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/Tissue_", i,"/GOBP_markers_down_Tissue_", i, "_tissueVSall.csv"))
  
  MFenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("markers/Tissue_", i,"/GOMF_markers_down_Tissue_", i, "_tissueVSall.csv"))
  
  CCenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("markers/Tissue_", i,"/GOCC_markers_down_Tissue_", i, "_tissueVSall.csv"))
  
  # Perform HALLMARK enrichment analysis
  HallmarkEnrich <- enricher(gene = interesting[interesting$p_val_adj < 0.05, "GeneSymbol"], 
                             TERM2GENE = h_t2g,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             pAdjustMethod = "bonferroni",
                             minGSSize = 15,
                             maxGSSize = 500,
                             universe = universe_genesymbol)
  Hallmark_GOI <- as.data.frame(HallmarkEnrich@result)
  write.csv(Hallmark_GOI, paste0("markers/Tissue_", i,"/HALLMARK_markers_down_Tissue_", i, "_tissueVSall.csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/Tissue_", i,"/KEGG_markers_down_Tissue_", i, "_tissueVSall.csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/Tissue_", i,"/Reactome_markers_down_Tissue_", i, "_tissueVSall.csv"))
}

if(place == "local") {
  saveRDS(filtered_exp, "Prefiltered_experiment_practice_seurat_integrated_markers_by_cluster_tissueVStissue_VSall.rds")
} else {
  saveRDS(filtered_exp, "Prefiltered_experiment_all_seurat_integrated_markers_by_cluster_tissueVStissue_VSall.rds")
}

if(place == "local") {
  saveRDS(all_markers, "Whole_experiment_tissueVSall_markers.rds")
} else {
  saveRDS(all_markers, "Whole_experiment_tissueVSall_markers.rds")
}
