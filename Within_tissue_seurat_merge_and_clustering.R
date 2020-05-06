#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Cluster samples merge samples by replicate within each tissue. Calculate site specificity scores for each cell.
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_by_organ/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/by_organ")
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


# Load prefiltered and clustered Seurat Object
if(place == "local") {
  filtered_exp <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data/Prefiltered_experiment_all_seurat_integrated.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

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
DEG_GOIs <- DEG_GOIs[is.element(DEG_GOIs, rownames(data.frame(filtered_exp@assays$RNA@counts)))]
filtered_exp <- SCTransform(filtered_exp, verbose = TRUE, vars.to.regress = c("Replicate", "Lib_size", "S.Score", "G2M.Score", "digest_stress1"), variable.features.n = 5000, return.only.var.genes = FALSE, new.assay.name = "SCT_whole")
filtered_exp_magic <- magic(filtered_exp, assay = "SCT_whole", genes = c(DEG_GOIs), t = "auto")
filtered_exp_magic@active.assay <- "MAGIC_SCT_whole"
filtered_exp_magic <- ScaleData(filtered_exp_magic)

for(sig in names(tissue_signatures)) {
  signature_GOIs <- as.character(tissue_signatures[[sig]]$Gene_name)[is.element(tissue_signatures[[sig]]$Gene_name, rownames(data.frame(filtered_exp_magic@assays$MAGIC_SCT_whole@scale.data)))]
  SCT_signature_GOIs <- data.frame(filtered_exp_magic@assays$MAGIC_SCT_whole@scale.data)[signature_GOIs, ]
  unique(signature_GOIs == rownames(SCT_signature_GOIs))
  sig_tmp <- tissue_signatures[[sig]]
  rownames(sig_tmp) <- sig_tmp$Gene_name
  sig_tmp <- sig_tmp[rownames(SCT_signature_GOIs),]
  correlation <- apply(SCT_signature_GOIs, 2, cor.test, sig_tmp$Signature)
  filtered_exp[[paste0(sig, "_corr")]] <- unlist(lapply(correlation, `[[`, "estimate"))
  filtered_exp[[paste0(sig, "_p.adj")]] <- p.adjust(unlist(lapply(correlation, `[[`, "p.value")), method = "bonf")
}

# Identify clusters within each tissue
# --------------------------------------------------------------------------
filtered_exp.list <- SplitObject(filtered_exp, split.by = "Tissue") # split into individual samples

for (tissue in names(filtered_exp.list)) {
  tissue_exp <- filtered_exp.list[[tissue]]
  print(tissue)
  dir.create(tissue)
  
  # Normalise transform counts within each experiment. Note: will not overwrite if already exists.
  # --------------------------------------------------------------------------
  
  tissue_exp.list <- SplitObject(tissue_exp, split.by = "Replicate") # split into individual samples
  
  # Perform SCT normalisation on each dataset individually
  for (i in 1:length(tissue_exp.list)) {
    tissue_exp.list[[i]] <- SCTransform(tissue_exp.list[[i]], verbose = TRUE, vars.to.regress = c("Lib_size", "S.Score", "G2M.Score", "digest_stress1"), variable.features.n = 5000, return.only.var.genes = FALSE)
  }
  
  # Integrate datasets based on highly correlated features
  tissue_exp.features <- SelectIntegrationFeatures(object.list = tissue_exp.list, nfeatures = 5000)
  tissue_exp.list <- PrepSCTIntegration(object.list = tissue_exp.list, verbose = TRUE, anchor.features = tissue_exp.features)
  # tissue_exp.list <- lapply(X = tissue_exp.list, FUN = RunPCA, verbose = TRUE, features = tissue_exp.features) # Perform PCA on each object individually (needed for rpca)
  reference_datasets <- which(names(tissue_exp.list) == "3")
  tissue_exp.anchors <- FindIntegrationAnchors(object.list = tissue_exp.list, normalization.method = "SCT", anchor.features = tissue_exp.features, verbose = TRUE, reduction = "cca", reference = reference_datasets)
  tissue_exp.integrated <- IntegrateData(anchorset = tissue_exp.anchors, normalization.method = "SCT", verbose = TRUE)
  
  # Run PCA on intergated dataset and determine clusters using Seurat
  # --------------------------------------------------------------------------
  
  # Seurat integration pipeline
  tissue_exp.integrated <- RunPCA(tissue_exp.integrated, dims = 1:50, assay = "integrated", ndims.print = 1:5, nfeatures.print = 5)
  tissue_exp.integrated <- FindNeighbors(tissue_exp.integrated, reduction = "pca", dims = 1:50)
  tissue_exp.integrated <- FindClusters(tissue_exp.integrated, resolution = 0.03)
  tissue_exp.integrated <- RunUMAP(tissue_exp.integrated, reduction = "pca", dims = 1:50)
  tissue_exp.integrated[[paste0(tissue, "_seurat_PCA_clusters")]] <- Idents(object = tissue_exp.integrated)
  
  # for(i in c("pca", "umap")) {
  #   p1 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Tissue")
  #   p2 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Replicate")
  #   p3 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
  #   p4 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("dying1"), sort.cell = TRUE)
  #   p5 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Phase")
  #   p6 <- DimPlot(tissue_exp.integrated, reduction = i, label = TRUE)
  #   gridit <- plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3)
  # ggsave(paste0(tissue, "/Seurat_clusters_", tissue, "_", i, ".png"), plot = gridit, device = "png")
  # }
  
  # Run PHATE on intergated dataset and determine clusters using kmeans
  # --------------------------------------------------------------------------
  
  # Run phate, find number of clusters in embeddings using silhouette
  phate.out <- phate(Matrix::t(GetAssayData(tissue_exp.integrated, assay = "integrated", slot = "data")), ndim = 10) # Runs PHATE diffusion map
  tissue_exp.integrated[[paste0(tissue, "_phate")]] <- CreateDimReducObject(embeddings = phate.out$embedding, key = "PHATE_", assay = "integrated")
  ProjectDim(tissue_exp.integrated, reduction = paste0(tissue, "_phate"))
  tissue_exp.integrated <- FindNeighbors(tissue_exp.integrated, reduction = paste0(tissue, "_phate"), dims = 1:10)
  tissue_exp.integrated <- FindClusters(tissue_exp.integrated, resolution = 0.03)
  tissue_exp.integrated[[paste0(tissue, "_PHATE_clusters")]] <- Idents(object = tissue_exp.integrated)
  
  phate_embed <- data.frame(Embeddings(tissue_exp.integrated, reduction = paste0(tissue, "_phate")))
  
  for(i in c("pca", "umap")) {
    p1 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Tissue")
    p2 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Replicate")
    p3 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
    p4 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("dying1"), sort.cell = TRUE)
    p5 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Mito_percent"), sort.cell = TRUE)
    p6 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Lib_size"), sort.cell = TRUE)
    p7 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Phase")
    p8 <- DimPlot(tissue_exp.integrated, reduction = i, label = TRUE)
    p9 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Total_Liver_sig_corr"), sort.cell = TRUE)
    p9 <- p9 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p10 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Total_LN_sig_corr"), sort.cell = TRUE)
    p10  <- p10 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p11 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Total_Lung_sig_corr"), sort.cell = TRUE)
    p11  <- p11 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p12 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Total_Primary_sig_corr"), sort.cell = TRUE)
    p12  <- p12 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    gridit1 <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
    gridit2 <- plot_grid(p9, p10, p11, p12, nrow = 2)
    
    ggsave(paste0(tissue, "/PHATE_clusters_", tissue, "_", i, ".png"), plot = gridit1, device = "png")
    ggsave(paste0(tissue, "/PHATE_specificity_", tissue, "_", i, ".png"), plot = gridit2, device = "png")
  }
  
  # No idea why phate needs coordinates set
  for(i in c(paste0(tissue, "_phate"))) {
    p1 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Tissue")
    p2 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Replicate")
    p3 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
    p3 <- p3 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p4 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("dying1"), sort.cell = TRUE)
    p4 <- p4 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p5 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Mito_percent"), sort.cell = TRUE)
    p5 <- p5 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p6 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Lib_size"), sort.cell = TRUE)
    p6 <- p6 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p7 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Phase")
    p8 <- DimPlot(tissue_exp.integrated, reduction = i, label = TRUE)
    p9 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Total_Liver_sig_corr"), sort.cell = TRUE)
    p9 <- p9 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p9 <- p9 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p10 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Total_LN_sig_corr"), sort.cell = TRUE)
    p10  <- p10 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p10 <- p10 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p11 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Total_Lung_sig_corr"), sort.cell = TRUE)
    p11  <- p11 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p11 <- p11 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    p12 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("Total_Primary_sig_corr"), sort.cell = TRUE)
    p12  <- p12 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    p12 <- p12 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
    
    gridit1 <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
    gridit2 <- plot_grid(p9, p10, p11, p12, nrow = 2)
    
    ggsave(paste0(tissue, "/PHATE_clusters_", tissue, "_", i, ".png"), plot = gridit1, device = "png")
    ggsave(paste0(tissue, "/PHATE_specificity_", tissue, "_", i, ".png"), plot = gridit2, device = "png")
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
  
  if(place == "local" & exists("tissue_exp.list")) {
    # saveRDS(tissue_exp.list, paste0(tissue, "/Prefiltered_experiment_practice_seurat_list_", tissue,".rds"))
  } else if(place == "wolfpack" & exists("phate.out")) {
    saveRDS(tissue_exp.list, paste0(tissue, "/Prefiltered_experiment_all_seurat_list_", tissue,".rds"))
  } else {
    print("Not overwritten")
  }
  
  if(place == "local" & exists("tissue_exp.anchors")) {
    # saveRDS(tissue_exp.anchors, paste0(tissue, "/Prefiltered_experiment_practice_seurat_anchors_", tissue,".rds"))
  } else if(place == "wolfpack" & exists("phate.out")) {
    saveRDS(tissue_exp.anchors, paste0(tissue, "/Prefiltered_experiment_all_seurat_anchors_", tissue,".rds"))
  } else {
    print("Not overwritten")
  }
  
  if(place == "local" & exists("tissue_exp.integrated")) {
    saveRDS(tissue_exp.integrated, paste0(tissue, "/Prefiltered_experiment_practice_seurat_integrated_", tissue,".rds"))
  } else if(place == "wolfpack" & exists("phate.out")) {
    saveRDS(tissue_exp.integrated, paste0(tissue, "/Prefiltered_experiment_all_seurat_integrated_", tissue,".rds"))
  } else {
    print("Not overwritten")
  }
  
}
