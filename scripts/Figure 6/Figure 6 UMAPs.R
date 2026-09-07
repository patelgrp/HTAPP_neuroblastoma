library(Seurat)
library(dplyr)
library(SingleR)
library(ggplot2)
library(RColorBrewer)
library(rliger)
library(SCpubr)
library(UCell)
library(scCustomize)
library(MetBrewer)
library(viridis)

#Data for this figure were generated using the Broad Terra pipeline
#This pipeline takes Cell Ranger aligned data, and uses Cellbender to remove
#technical artifacts. Cells/nuclei underwent QC and filtering, and were integrated
#using LIGER. See the pipeline code to investigate further.

######################CODE TO GENERATE FIGURE 6 UMAP PANELS######################
#load complete dataset with all single-cell and single-nucleus data (only malignant cells/nuclei)
# gcdata <-
#   readRDS(
#     "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/malignant_combined_dataset_k12.Rds"
#   )
output.dir <- output.dir <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/Fig3"
# 
# #cluster 7 is a minor cluster contaminated with doublets expressing tumor + macrophage markers (PTPRC+, CD68+, LYZ+)
# gcdata <-
#   subset(gcdata, subset = inmfNorm.cluster == 7, invert = TRUE)
# 
# #reprocesss data
# gcdata <- normalize(gcdata)
# gcdata <- selectGenes(gcdata)
# gcdata <- scaleNotCenter(gcdata)
# gcdata
# 
# gcdata <- runINMF(gcdata, k = 8)
# gcdata <- quantileNorm(gcdata)
# gcdata
# 
# gcdata <-
#   RunUMAP(gcdata, reduction = "inmfNorm", dims = 1:8)
# 
# saveRDS(gcdata,
#         file = paste0(output.dir, "malignant_combined_dataset_k8_reprocess.Rds"))

gcdata <- readRDS(paste0(output.dir, "malignant_combined_dataset_k8_reprocess.Rds"))
###################### PANEL A ######################
pdf(paste0(output.dir, "/panel 6G 20220714 HTAPP malignant liger.pdf"))
print(
  DimPlot_scCustom(
    seurat_object = gcdata,
    figure_plot = TRUE,
    group.by = "inmfNorm.cluster",
    colors_use = met.brewer("Signac", 12)
  )
)
dev.off()

pdf(paste0(output.dir, "/panel 6G 20220714 HTAPP malignant liger with labels.pdf"))
print(
  DimPlot_scCustom(
    seurat_object = gcdata,
    figure_plot = TRUE, 
    group.by = "inmfNorm.cluster", 
    colors_use = met.brewer("Signac", 12), 
    label = TRUE
  )
)
dev.off()

#pull top 50 feature embeddings for each NMF cluster for supplemental tables
nmf.embeddings <- matrix(nrow = 50, ncol = 8)
colnames(nmf.embeddings) <- paste0("iNMF", 1:8)
for (i in 1:8)
{
  nmf.embeddings[, i] <-
    TopFeatures(object = gcdata[["inmfNorm"]],
                dim = i,
                nfeatures = 50)
}
write.csv(
  nmf.embeddings,
  file = paste0(output.dir, "/NMF embeddings.csv"),
  row.names = FALSE
)

gcdata.flat <- JoinLayers(gcdata)
#calculate DE genes for each iNMF cluster
de.markers <-
  FindAllMarkers(
    object = gcdata.flat,
    group.by = "inmfNorm.cluster",
    logfc.threshold = 2,
    only.pos = TRUE
  )
de.markers.filter <- subset(de.markers, subset = p_val_adj < 1E-100)
write.csv(de.markers.filter, file = paste0(output.dir, "/Table S5 DE genes.csv"))


###################### PANEL 6F ######################
signatures <-
  read.csv("ADRN_MES_signatures.csv")  #signatures from von Groningen 2017 Nat Genet
ADRN.genes <- signatures$ADRENERGIC
MES.genes <- signatures$MESENCHYMAL
ADRN.genes <- ADRN.genes[ADRN.genes != ""]
MES.genes <- MES.genes[MES.genes != ""]

signatures_dyer <-
  read.csv("dyer_state_sigs.csv")  #signatures from Dyer re-analysis
ADRN.genes.dyer <- signatures_dyer$ADRN
MES.genes.dyer <- signatures_dyer$MES
ADRN.genes.dyer <- ADRN.genes.dyer[ADRN.genes.dyer != ""]
MES.genes.dyer <- MES.genes.dyer[MES.genes.dyer != ""]

Idents(gcdata) <- "inmfNorm.cluster"
gcdata@active.ident <-
  factor(gcdata@active.ident, levels = ordering)

gcdata <- JoinLayers(gcdata)
gcdata <-
  AddModuleScore_UCell(obj = gcdata,
                       features = list(ADRN = ADRN.genes, MES = MES.genes))
gcdata <-
  AddModuleScore_UCell(obj = gcdata,
                       features = list(ADRN.dyer = ADRN.genes.dyer, MES.dyer = MES.genes.dyer))


pdf(paste0(output.dir, "/Mike scores.pdf"))

FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = c("ADRN_UCell"),
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024), 
  na_cutoff = NULL,
  order = TRUE, max.cutoff  = 0.3, min.cutoff = 0.2,
  colors_use = viridis(n = 10, option = "D", alpha = 0.7)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = c("MES_UCell"),
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024), 
  na_cutoff = NULL,
  order = TRUE, max.cutoff  = 0.25, min.cutoff = 0.15,
  colors_use = viridis(n = 10, option = "D", alpha = 0.7)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = c("ADRN.dyer_UCell"),
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024), 
  na_cutoff = NULL,
  order = TRUE, max.cutoff  = 0.2, min.cutoff = 0.05,
  colors_use = viridis(n = 10, option = "D", alpha = 0.7)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = c("MES.dyer_UCell"),
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024), 
  na_cutoff = NULL,
  order = TRUE, max.cutoff  = 0.05, min.cutoff = 0.01,
  colors_use = viridis(n = 10, option = "D", alpha = 0.7)
)
VlnPlot_scCustom(
  seurat_object = gcdata,
  features = c("ADRN_UCell", "MES_UCell", "ADRN.dyer_UCell", "MES.dyer_UCell"),
  pt.size = 0,
  group.by = "inmfNorm.cluster",
  num_columns = 1,
  colors_use = met.brewer("Signac", 12),
  plot_median = TRUE,
  median_size = 5
)

dev.off()

new.cluster.ids <- c(
  "ADRN",
  "ADRN",
  "ADRN",
  "ADRN",
  "ADRN",
  "ADRN",
  "ADRN",
  "ADRN",
  "MES",
  "MES",
  "MES",
  "sympathoblast"
)
names(new.cluster.ids) <- levels(gcdata@active.ident)
gcdata <- RenameIdents(gcdata, new.cluster.ids)
gcdata$annot.ids <- gcdata@active.ident

pdf(paste0(
  output.dir ,
  "/panel 3D-G 20220722 HTAPP malignant van Groningen plots.pdf"
))
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "ADRN1",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024),
  na_cutoff = NULL,
  max.cutoff = 0.4,
  order = TRUE,
  colors_use = viridis(n = 10, option = "D", alpha = 0.5)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "MES1",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024),
  na_cutoff = NULL,
  order = TRUE,
  colors_use = viridis(n = 10, option = "D", alpha = 0.5)
)
print(
  SCpubr::do_ViolinPlot(
    sample = gcdata,
    features = "ADRN1",
    plot_boxplot = T,
    group.by = "annot.ids",
    colors.use = c(
      "MES" = "#8A181A",
      "ADRN" = "#1C83BA",
      "sympathoblast" = "#009E73"
    )
  )
)
print(
  SCpubr::do_ViolinPlot(
    sample = gcdata,
    features = "MES1",
    plot_boxplot = T,
    group.by = "annot.ids",
    colors.use = c(
      "MES" = "#8A181A",
      "ADRN" = "#1C83BA",
      "sympathoblast" = "#009E73"
    )
  )
)
# print(FeatureScatter(gcdata, feature1 = "ADRN1", feature2 = "MES1", shuffle = F, raster = F,
#                      pt.size = 1, cols = c("#2E307B", "#009D73", "#89191B"), jitter = T))
dev.off()
