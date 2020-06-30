#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Cluster samples merge samples by replicate within each tissue. Calculate site specificity scores for each cell.
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCPDXmets_chrcha/project_results/seurat/by_organ/") # Uses full dataset as is only small
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCPDXmets_chrcha/project_results/seurat/by_organ")
  place <- "wolfpack"
}

# Libraries
library('Seurat')
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('ggplot2')
library('readr')
library('Matrix')
library('phateR')
library('cowplot')
library('factoextra')
library('gtools')
library('doParallel')
library('foreach')
library('clusterProfiler')
library('ReactomePA')
library('org.Hs.eg.db')
library('Rmagic')
library('ggridges')
library('colorblindr')


# Load prefiltered and clustered Seurat Object
if(place == "local") {
  filtered_exp_sce <- readRDS("/Users/mac/cloudstor/sarah_projects/SCPDXmets_chrcha/project_results/prefiltering/all_data/Prefiltered_QC_experiment_all.rds") # uses practice data if local
} else {
  filtered_exp_sce <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCPDXmets_chrcha/project_results/prefiltering/all_data/Prefiltered_QC_experiment_all.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Rename colData Sample slot. Add Tissue, Met, Replicate information.
filtered_exp_sce$old_Sample <- filtered_exp_sce$Sample

new_ids <- data.frame(filtered_exp_sce$Sample) %>%
  separate(1, c("TISSUE", "REP"), "_", remove = FALSE)
colnames(new_ids) <- c("Sample", "Tissue", "Replicate")

idx <- match(filtered_exp_sce$Sample, new_ids$Sample)
filtered_exp_sce$Tissue <- as.character(new_ids$Tissue) [idx]
filtered_exp_sce$Replicate <- as.character(new_ids$Replicate) [idx]
filtered_exp_sce$cellIDs <- rownames((colData(filtered_exp_sce)))

# Set up prefiltered object, add cell cycle difference info and split by sample
# --------------------------------------------------------------------------
filtered_exp_seurat <- as.Seurat(filtered_exp_sce, counts = "counts", data = NULL) # convert to Seurat

# Add Entrez IDs
rowData(filtered_exp_sce)$EntrezID <- mapIds(org.Hs.eg.db, keys=rowData(filtered_exp_sce)$Ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
rowMetaData <- data.frame(rowData(filtered_exp_sce))
filtered_exp_seurat@assays$RNA@meta.features <- merge(filtered_exp_seurat@assays$RNA@meta.features, rowMetaData, by.x = 0, by.y = 0)

# Add cell cycle
s.genes <- cc.genes.updated.2019$s.genes
s.genes <- c(rownames(rowMetaData[is.element(rowMetaData$GeneSymbol, s.genes), ]))
g2m.genes <- cc.genes.updated.2019$g2m.genes
g2m.genes <- c(rownames(rowMetaData[is.element(rowMetaData$GeneSymbol, g2m.genes), ]))
filtered_exp_seurat <- CellCycleScoring(filtered_exp_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# Add score to idenitfy dying cells and stress scores with low mito content (as per DOI 10.1186/s13059-019-1830-0)
digest_stress <- c("FOS", "CXCL2", "ZFP36", "FOSB", "DUSP1", "ATF3", "CXCL8", "NR4A1", "CXCL3", "PPP1R15A", "JUNB", "EGR1", "HSPA1A", "HSPA1B", "SOCS3", "KLF6", "JUN", "IER2", "CXCL1", "NKFBIA", "HSPA6", "DNAJB1", "IER3", "CCNL1", "MTRNR2L2", "IER5", "ID1", "CEBPD", "KRT6A", "CYR61", "DEPP1", "CLDN4", "IRF1", "DUSP2", "BTG2", "PLAUR", "MAFF", "KLF4", "PHLDA2", "TNFAIP3")
digest_stress <- list(c(rownames(rowMetaData[is.element(rowMetaData$GeneSymbol, digest_stress), ])))
filtered_exp_seurat <- AddModuleScore(object = filtered_exp_seurat, features = digest_stress, name = 'digest_stress')

dying <- c("HLA-A", "HLA-B", "HLA-C", "B2M")
dying <- list(c(rownames(rowMetaData[is.element(rowMetaData$GeneSymbol, dying), ])))
filtered_exp_seurat <- AddModuleScore(object = filtered_exp_seurat, features = dying, name = 'dying')

# Calculate tissue specificity score (correlation) for all cells. Correlation between signature and scaled normalised counts calculated across the whole experiment.
# --------------------------------------------------------------------------

# Load site specificity signatures
if(place == "local") {
  tissue_signatures <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Pseudo-bulk_tissue_specific_signature.rds") # uses practice data if local
} else {
  tissue_signatures <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Pseudo-bulk_tissue_specific_signature.rds") # uses whole dataset if wolfpack
}

# Calculate normalised, scaled counts across the whole experiment
DEG_GOIs <- as.character(unique(unlist(lapply(tissue_signatures, `[[`, "Gene_name"))))
DEG_GOIs <- DEG_GOIs[is.element(DEG_GOIs, rownames(filtered_exp_seurat@assays$RNA@counts))]
filtered_exp_SCT <- SCTransform(filtered_exp_seurat, verbose = TRUE, vars.to.regress = c("Replicate", "Lib_size", "S.Score", "G2M.Score", "digest_stress1"), variable.features.n = 5000, return.only.var.genes = FALSE, new.assay.name = "SCT_whole")
filtered_exp_magic <- Rmagic::magic(filtered_exp_SCT, assay = "SCT_whole", genes = c(DEG_GOIs), t = "auto")
filtered_exp_magic@active.assay <- "MAGIC_SCT_whole"
filtered_exp_magic <- ScaleData(filtered_exp_magic)

for(sig in names(tissue_signatures)) {
  sig_tmp <- tissue_signatures[[sig]]
  rownames(sig_tmp) <- as.character(sig_tmp$Gene_name)
  signature_GOIs <- as.character(sig_tmp$Gene_name)[is.element(sig_tmp$Gene_name, rownames(filtered_exp_magic@assays$MAGIC_SCT_whole@scale.data))]
  SCT_signature_GOIs <- Matrix::as.matrix(filtered_exp_magic@assays$MAGIC_SCT_whole@scale.data)[signature_GOIs, ]
  unique(signature_GOIs == rownames(SCT_signature_GOIs))
  sig_tmp <- sig_tmp[rownames(SCT_signature_GOIs),]
  correlation <- apply(SCT_signature_GOIs, 2, cor.test, sig_tmp$Signature, method = "spearman")
  filtered_exp_seurat[[paste0(sig, "_corr")]] <- unlist(lapply(correlation, `[[`, "estimate"))
  filtered_exp_seurat[[paste0(sig, "_p.adj")]] <- p.adjust(unlist(lapply(correlation, `[[`, "p.value")), method = "bonf")
}

# Save and remove large objects
if(place == "local") {
  saveRDS(filtered_exp_seurat, "Prefiltered_experiment_practice_seurat_integrated_w_SiteSpecSig.rds")
} else if(place == "wolfpack") {
  saveRDS(filtered_exp_seurat, "Prefiltered_experiment_all_seurat_integrated_w_SiteSpecSig.rds")
} else {
  print("Not overwritten")
}

# Save and remove large objects
if(place == "local") {
  saveRDS(filtered_exp_SCT, "Prefiltered_experiment_practice_seurat_integrated_whole_SCT.rds")
} else if(place == "wolfpack") {
  saveRDS(filtered_exp_SCT, "Prefiltered_experiment_all_seurat_integrated_whole_SCT.rds")
} else {
  print("Not overwritten")
}
rm(filtered_exp_SCT)

if(place == "local") {
  saveRDS(filtered_exp_magic, "Prefiltered_experiment_practice_seurat_integrated_MAGIC_DGEGOI.rds")
} else if(place == "wolfpack") {
  saveRDS(filtered_exp_magic, "Prefiltered_experiment_all_seurat_integrated_MAGIC_DGEGOI.rds")
} else {
  print("Not overwritten")
}
rm(filtered_exp_magic)

# Identify clusters within each tissue
# --------------------------------------------------------------------------
filtered_exp.list <- SplitObject(filtered_exp_seurat, split.by = "Sample") # split into individual samples

for (tissue in names(filtered_exp.list)) {
  tissue_exp <- filtered_exp.list[[tissue]]
  print(tissue)
  dir.create(tissue)
  
  # Normalise transform counts within each experiment. Note: will not overwrite if already exists.
  # --------------------------------------------------------------------------
  
  # Perform SCT normalisation on each dataset individually
  tissue_exp <- SCTransform(tissue_exp, verbose = TRUE, vars.to.regress = c("Lib_size", "S.Score", "G2M.Score", "digest_stress1"), variable.features.n = 5000, return.only.var.genes = TRUE)

  
  tissue_exp.ind <- tissue_exp
  # Seurat integration pipeline
  tissue_exp.ind <- RunPCA(tissue_exp.ind, npcs = 20, assay = "SCT", ndims.print = 1:5, nfeatures.print = 5)
  tissue_exp.ind <- FindNeighbors(tissue_exp.ind, reduction = "pca", dims = 1:20)
  tissue_exp.ind <- FindClusters(tissue_exp.ind, resolution = 0.03)
  tissue_exp.ind <- RunUMAP(tissue_exp.ind, reduction = "pca", dims = 1:20)
  tissue_exp.ind[[paste0(tissue, "_seurat_PCA_clusters")]] <- Idents(object = tissue_exp.ind)
  
  
  # Run PHATE on intergated dataset and determine clusters using kmeans
  # --------------------------------------------------------------------------
  
  # Run phate, find number of clusters in embeddings using silhouette
  phate.out <- phate(Matrix::t(GetAssayData(tissue_exp.ind, assay = "SCT", slot = "data")), ndim = 10) # Runs PHATE diffusion map
  tissue_exp.ind[[paste0(tissue, "_phate")]] <- CreateDimReducObject(embeddings = phate.out$embedding, key = "PHATE_", assay = "SCT")
  ProjectDim(tissue_exp.ind, reduction = paste0(tissue, "_phate"))
  tissue_exp.ind <- FindNeighbors(tissue_exp.ind, reduction = paste0(tissue, "_phate"), dims = 1:10)
  tissue_exp.ind <- FindClusters(tissue_exp.ind, resolution = 0.03)
  tissue_exp.ind[[paste0(tissue, "_PHATE_clusters")]] <- Idents(object = tissue_exp.ind)
  
  phate_embed <- data.frame(Embeddings(tissue_exp.ind, reduction = paste0(tissue, "_phate")))
  
  for(i in c("pca", "umap")) {
    p1 <- DimPlot(tissue_exp.ind, reduction = i, group.by = "Tissue")
    p2 <- DimPlot(tissue_exp.ind, reduction = i, group.by = "Replicate")
    p3 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
    p4 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("dying1"), sort.cell = TRUE)
    p5 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("Mito_percent"), sort.cell = TRUE)
    p6 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("Lib_size"), sort.cell = TRUE)
    p7 <- DimPlot(tissue_exp.ind, reduction = i, group.by = "Phase")
    p8 <- DimPlot(tissue_exp.ind, reduction = i, label = TRUE)
    p9 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("Total_Liver_sig_corr"), sort.cell = TRUE)
    p9 <- p9 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p11 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("Total_Lung_sig_corr"), sort.cell = TRUE)
    p11  <- p11 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p12 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("Total_Primary_sig_corr"), sort.cell = TRUE)
    p12  <- p12 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    gridit1 <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
    gridit2 <- plot_grid(p9, p11, p12, nrow = 2)
    
    ggsave(paste0(tissue, "/PHATE_clusters_", tissue, "_", i, ".png"), plot = gridit1, device = "png")
    ggsave(paste0(tissue, "/PHATE_specificity_", tissue, "_", i, ".png"), plot = gridit2, device = "png")
  }
  
  # No idea why phate needs coordinates set
  for(i in c(paste0(tissue, "_phate"))) {
    p1 <- DimPlot(tissue_exp.ind, reduction = i, group.by = "Tissue")
    p2 <- DimPlot(tissue_exp.ind, reduction = i, group.by = "Replicate")
    p3 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
    p3 <- p3 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p4 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("dying1"), sort.cell = TRUE)
    p4 <- p4 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p5 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("Mito_percent"), sort.cell = TRUE)
    p5 <- p5 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p6 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("Lib_size"), sort.cell = TRUE)
    p6 <- p6 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p7 <- DimPlot(tissue_exp.ind, reduction = i, group.by = "Phase")
    p8 <- DimPlot(tissue_exp.ind, reduction = i, label = TRUE)
    p9 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("Total_Liver_sig_corr"), sort.cell = TRUE)
    p9 <- p9 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p9 <- p9 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p11 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("Total_Lung_sig_corr"), sort.cell = TRUE)
    p11  <- p11 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p11 <- p11 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p12 <- FeaturePlot(tissue_exp.ind, reduction = i, features = c("Total_Primary_sig_corr"), sort.cell = TRUE)
    p12  <- p12 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p12 <- p12 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    
    gridit1 <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
    gridit2 <- plot_grid(p9, p11, p12, nrow = 2)
    
    ggsave(paste0(tissue, "/PHATE_clusters_", tissue, "_", i, ".png"), plot = gridit1, device = "png")
    ggsave(paste0(tissue, "/PHATE_specificity_", tissue, "_", i, ".png"), plot = gridit2, device = "png")
  }
  
  tissue_exp.ind$Primary_like_cells <- tissue_exp.ind$Total_Primary_sig_corr > 0 & tissue_exp.ind$Total_Primary_sig_p.adj < 0.05
  tissue_exp.ind$Liver_like_cells <- tissue_exp.ind$Total_Liver_sig_corr > 0 & tissue_exp.ind$Total_Liver_sig_p.adj < 0.05
  tissue_exp.ind$Lung_like_cells <- tissue_exp.ind$Total_Lung_sig_corr > 0 & tissue_exp.ind$Total_Lung_sig_p.adj < 0.05
  
  tissue_exp.ind$Primary_like_cells[tissue_exp.ind$Liver_like_cells | tissue_exp.ind$Lung_like_cells] <- FALSE
  tissue_exp.ind$Liver_like_cells[tissue_exp.ind$Primary_like_cells | tissue_exp.ind$Lung_like_cells] <- FALSE
  tissue_exp.ind$Lung_like_cells[tissue_exp.ind$Liver_like_cells | tissue_exp.ind$Primary_like_cells] <- FALSE
  
  table(tissue_exp.ind$Primary_like_cells, tissue_exp.ind$Tissue)
  table(tissue_exp.ind$Lung_like_cells, tissue_exp.ind$Tissue)
  table(tissue_exp.ind$Liver_like_cells, tissue_exp.ind$Tissue)
  
  tissue_exp.ind_MD <- data.frame(tissue_exp.ind@meta.data)
  
  for(i in c("Liver", "Lung", "Primary")) {
    plotData <- tissue_exp.ind_MD[, c("Tissue", paste0("Total_", i, "_sig_corr"), paste0(i, "_like_cells"))]
    colnames(plotData) <- c("Tissue", "Corr", "Likeness")
    ggplot(plotData, aes(x = Corr, y = Tissue, fill = Tissue)) +
      geom_density_ridges() +
      scale_fill_OkabeIto() +
      theme_classic() +
      ggsave(paste0(tissue, "/Cell_identity_ridge_", i, ".pdf"))
    
    sample_likeness <- as.data.frame.matrix(table(plotData$Tissue, plotData$Likeness))
    colnames(sample_likeness) <- c("Not", "Likeness")
    sample_likeness$Total_cells <- rowSums(sample_likeness)
    sample_likeness$Percent_likeness <- sample_likeness$Likeness/sample_likeness$Total_cells*100
    sample_likeness$Tissue <- rownames(sample_likeness)
    sample_likeness$Tissue <- factor(sample_likeness$Tissue, levels = c("Primary", "Lung", "BM"))
    write.csv(sample_likeness, paste0(tissue, "/Percent_likeness_",i, ".csv"))
    
    ggplot(data=sample_likeness, aes(x=Tissue, y=Percent_likeness, fill = Tissue)) +
      geom_bar(stat="identity") +
      ggtitle(paste0(i, "-adapted cells (%)")) +
      scale_fill_OkabeIto() +
      ylim(0, 100) +
      coord_flip() +
      theme_classic() +
      ggsave(paste0(tissue, "/Percent_likeness_", i, ".pdf"))
  }
  
  
  # Save seurat objects
  # --------------------------------------------------------------------------
  
  if(place == "local" & exists("phate.out")) {
    # saveRDS(phate.out, paste0(tissue, "/Prefiltered_experiment_practice_phate.out_", tissue,".rds"))
  } else if(place == "wolfpack" & exists("phate.out")) {
    saveRDS(phate.out, paste0(tissue, "/Prefiltered_experiment_all_phate.out_", tissue,".rds"))
  } else {
    print("Not overwritten")
  }
  
  if(place == "local" & exists("tissue_exp.ind.list")) {
    # saveRDS(tissue_exp.ind.list, paste0(tissue, "/Prefiltered_experiment_practice_seurat_list_", tissue,".rds"))
  } else if(place == "wolfpack" & exists("phate.out")) {
    saveRDS(tissue_exp.ind.list, paste0(tissue, "/Prefiltered_experiment_all_seurat_list_", tissue,".rds"))
  } else {
    print("Not overwritten")
  }
  
  if(place == "local" & exists("tissue_exp.ind.anchors")) {
    # saveRDS(tissue_exp.ind.anchors, paste0(tissue, "/Prefiltered_experiment_practice_seurat_anchors_", tissue,".rds"))
  } else if(place == "wolfpack" & exists("phate.out")) {
    saveRDS(tissue_exp.ind.anchors, paste0(tissue, "/Prefiltered_experiment_all_seurat_anchors_", tissue,".rds"))
  } else {
    print("Not overwritten")
  }
  
  if(place == "local" & exists("tissue_exp.ind")) {
    saveRDS(tissue_exp.ind, paste0(tissue, "/Prefiltered_experiment_practice_seurat_integrated_", tissue,".rds"))
  } else if(place == "wolfpack" & exists("phate.out")) {
    saveRDS(tissue_exp.ind, paste0(tissue, "/Prefiltered_experiment_all_seurat_integrated_", tissue,".rds"))
  } else {
    print("Not overwritten")
  }
  }
}
