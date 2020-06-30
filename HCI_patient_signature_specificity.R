#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Create SC experiment from HCI pdx data. Calculate site specificity scores for each cell.
# --------------------------------------------------------------------------

# Working directory
setwd("/Users/mac/cloudstor/sarah_projects/SCHCImets_chrcha/project_results/seurat/signatures/") # Uses practice data (5% of cells from each sample) if running locally
place <- "local"

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

# Load site specificity signatures
tissue_signatures <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Pseudo-bulk_tissue_specific_signature.rds") # uses practice data if local

rowData_summed <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Pseudo-bulk_rowData_all.csv", header = TRUE, row.names = 1)

# Load HCI data and create Seurat Object
HCI001 <- read.table("/Users/mac/cloudstor/sarah_projects/SCHCImets_chrcha/data/GSE123837_HCI001_combinedmatrix.geneIDensemblID_FPKM_Submission.txt", header = TRUE, row.names = 1)
HCI002 <- read.csv("/Users/mac/cloudstor/sarah_projects/SCHCImets_chrcha/data/GSE123837_HCI002_combinedmatrix.geneIDensemblID_FPKM_Submission_curated_cleaned.csv", header = TRUE, row.names = 1)
HCI010 <- read.table("/Users/mac/cloudstor/sarah_projects/SCHCImets_chrcha/data/GSE123837_HCI010_combinedmatrix.geneIDensemblID_FPKM_Submission.txt", header = TRUE, row.names = 1)

NAMES <- Reduce(intersect, list(rownames(HCI001), rownames(HCI002), rownames(HCI010)))
ALL <- merge(HCI001[NAMES, ], HCI002[NAMES, ], by.x = 0, by.y = 0, )
rownames(ALL) <- ALL[, 1]
ALL <- ALL[, -1]
ALL <- merge(ALL, HCI010[NAMES, ], by.x = 0, by.y = 0)
rownames(ALL) <- ALL[, 1]
ALL <- ALL[, -1]

PDX_patients <- list("ALL" = Matrix::as.matrix(ALL), "HCI001" = Matrix::as.matrix(HCI001), "HCI002" = Matrix::as.matrix(HCI002), "HCI010" = Matrix::as.matrix(HCI010))

# Prepare Seurat object
# --------------------------------------------------------------------------

for(pdx in names(PDX_patients)) {
  p10 <- CreateSeuratObject(counts = PDX_patients[[pdx]], min.cells = 8, min.features = 1000, project = "HCI010")
  DefaultAssay(p10) <- "RNA"
  p10[["percent.mt"]] <- PercentageFeatureSet(p10, pattern = "^MT-")
  p10 <- subset(p10, subset = nFeature_RNA > 2500 & percent.mt < 50)
  p10[["RNA"]]@data <- Matrix::as.matrix(log(p10[["RNA"]]@counts + 1))  # Normalise FPKMs
  
  new_ids <- data.frame(colnames(PDX_patients[[pdx]])) %>%
    separate(1, c("PDX", "MOUSE", "TYPE", "TISSUE", "CELL"), "_", remove = FALSE)
  colnames(new_ids) <- c("Cell_ID", "PDX", "MOUSE", "TYPE", "TISSUE", "CELL")
  
  PDX <- new_ids$PDX
  names(PDX) <- new_ids$Cell_ID
  p10$PDX <- PDX

  TISSUE <- new_ids$TISSUE
  names(TISSUE) <- new_ids$Cell_ID
  p10$Tissue <- TISSUE
  
  MOUSE <- new_ids$MOUSE
  names(MOUSE) <- new_ids$Cell_ID
  p10$Replicate <- MOUSE

  # Regress out cell cycle
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  p10 <- CellCycleScoring(p10, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  p10$CC.Difference <- p10$S.Score - p10$G2M.Score
  if(pdx == "ALL"){
    p10 <- ScaleData(p10, vars.to.regress = c("percent.mt", "nFeature_RNA", "CC.Difference", "PDX"))
  } else if (pdx == "HCI002" ) {
    p10 <- ScaleData(p10, vars.to.regress = c("percent.mt", "nFeature_RNA", "CC.Difference"))
  } else{
    p10 <- ScaleData(p10, vars.to.regress = c("percent.mt", "nFeature_RNA"))
  }
  
  # Calculate tissue specificity score (correlation) for all cells. Correlation between signature and scaled normalised counts calculated across the whole experiment.
  # --------------------------------------------------------------------------
  
  # Calculate normalised, scaled counts across the whole experiment
  DEG_GOIs <- data.frame("Gene_name" = as.character(unique(unlist(lapply(tissue_signatures, `[[`, "Gene_name")))))
  idx <- match(DEG_GOIs$Gene_name, rownames(rowData_summed))
  DEG_GOIs$GeneSymbol <- rowData_summed$GeneSymbol [idx]
  DEG_GOIs <- as.character(DEG_GOIs[is.element(DEG_GOIs$GeneSymbol, rownames(p10@assays$RNA@counts)), "GeneSymbol"])
  p10_magic <- Rmagic::magic(p10, assay = "RNA", genes = c(DEG_GOIs), t = "auto")
  p10_magic@active.assay <- "MAGIC_RNA"
  p10_magic <- ScaleData(p10_magic)

  for(sig in names(tissue_signatures)) {
    sig_tmp <- tissue_signatures[[sig]]
    idx <- match(sig_tmp$Gene_name, rownames(rowData_summed))
    sig_tmp$GeneSymbol <- rowData_summed$GeneSymbol [idx]
    rownames(sig_tmp) <- sig_tmp$GeneSymbol
    
    signature_GOIs <- as.character(sig_tmp$GeneSymbol)
    signature_GOIs <- signature_GOIs[is.element(signature_GOIs, rownames(p10_magic@assays$MAGIC_RNA@scale.data))]
    SCT_signature_GOIs <- Matrix::as.matrix(p10_magic@assays$MAGIC_RNA@scale.data)[signature_GOIs, ]
    unique(signature_GOIs == rownames(SCT_signature_GOIs))
    
    sig_tmp <- sig_tmp[rownames(SCT_signature_GOIs),]
    correlation <- apply(SCT_signature_GOIs, 2, cor.test, sig_tmp$Signature, method = "spearman")
    p10[[paste0(sig, "_corr")]] <- unlist(lapply(correlation, `[[`, "estimate"))
    p10[[paste0(sig, "_p.adj")]] <- p.adjust(unlist(lapply(correlation, `[[`, "p.value")), method = "bonf")
  }
  
  p10 <- FindVariableFeatures(p10)
  p10 <- RunPCA(p10, dims = 1:50, assay = "RNA", ndims.print = 1:5, nfeatures.print = 5)
  p10 <- FindNeighbors(p10, reduction = "pca", dims = 1:50)
  p10 <- FindClusters(p10, resolution = 0.03)
  p10 <- RunUMAP(p10, reduction = "pca", dims = 1:50)
  
  for(i in c("pca", "umap")) {
    plot1 <- DimPlot(p10, reduction = i, group.by = "Tissue")
    plot2 <- DimPlot(p10, reduction = i, group.by = "Replicate")
    plot3 <- FeaturePlot(p10, reduction = i, features = c("nFeature_RNA"), sort.cell = TRUE)
    plot5 <- DimPlot(p10, reduction = i, group.by = "Phase")
    plot6 <- FeaturePlot(p10, reduction = i, features = c("Total_Liver_sig_corr"), sort.cell = TRUE)
    plot6 <- plot6 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    plot7 <- FeaturePlot(p10, reduction = i, features = c("Total_Lung_sig_corr"), sort.cell = TRUE)
    plot7  <- plot7 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    plot8 <- FeaturePlot(p10, reduction = i, features = c("Total_Primary_sig_corr"), sort.cell = TRUE)
    plot8  <- plot8 + scale_colour_gradientn(colours = hcl.colors(7, palette = "TealRose"), limits = c(-1, 1))
    gridit1 <- plot_grid(plot1, plot2, plot3, plot5, plot6, plot7, plot8, nrow = 4)
    ggsave(paste0("Site_specificity_QC_", pdx,"_", i, ".png"), plot = gridit1, device = "png")
  }
  
  p10$Primary_like_cells <- p10$Total_Primary_sig_corr > 0 & p10$Total_Primary_sig_p.adj < 0.05
  p10$Liver_like_cells <- p10$Total_Liver_sig_corr > 0 & p10$Total_Liver_sig_p.adj < 0.05
  p10$Lung_like_cells <- p10$Total_Lung_sig_corr > 0 & p10$Total_Lung_sig_p.adj < 0.05
  
  p10$Primary_like_cells[p10$Liver_like_cells | p10$Lung_like_cells] <- FALSE
  p10$Liver_like_cells[p10$Primary_like_cells | p10$Lung_like_cells] <- FALSE
  p10$Lung_like_cells[p10$Liver_like_cells | p10$Primary_like_cells] <- FALSE
  
  table(p10$Primary_like_cells, p10$Tissue)
  table(p10$Lung_like_cells, p10$Tissue)
  table(p10$Liver_like_cells, p10$Tissue)
  
  p10_MD <- data.frame(p10@meta.data)
  
  for(i in c("Liver", "Lung", "Primary")) {
    plotData <- p10_MD[, c("Tissue", paste0("Total_", i, "_sig_corr"), paste0(i, "_like_cells"))]
    colnames(plotData) <- c("Tissue", "Corr", "Likeness")
    ggplot(plotData, aes(x = Corr, y = Tissue, fill = Tissue)) +
      geom_density_ridges() +
      scale_fill_OkabeIto() +
      theme_classic() +
      ggsave(paste0("Cell_identity_ridge_", pdx,"_", i, ".pdf"))
    
    sample_likeness <- as.data.frame.matrix(table(plotData$Tissue, plotData$Likeness))
    colnames(sample_likeness) <- c("Not", "Likeness")
    sample_likeness$Total_cells <- rowSums(sample_likeness)
    sample_likeness$Percent_likeness <- sample_likeness$Likeness/sample_likeness$Total_cells*100
    sample_likeness$Tissue <- rownames(sample_likeness)
    sample_likeness$Tissue <- factor(sample_likeness$Tissue, levels = c("PT", "LN", "LU"))
    write.csv(sample_likeness, paste0("Percent_likeness_", pdx,"_",i, ".csv"))
    
    ggplot(data=sample_likeness, aes(x=Tissue, y=Percent_likeness, fill = Tissue)) +
      geom_bar(stat="identity") +
      ggtitle(paste0(i, "-adapted cells (%)")) +
      scale_fill_OkabeIto() +
      ylim(0, 100) +
      coord_flip() +
      theme_classic() +
      ggsave(paste0("Percent_likeness_", pdx,"_",i, ".pdf"))
  }
  
  saveRDS(p10, paste0(pdx,"_signature_seurat_object.rds"))

}

