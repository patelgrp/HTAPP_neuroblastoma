library(Seurat)
library(dplyr)
library(tidyr)
library(scCustomize)
library(SCpubr)
library(MetBrewer)
library(ggplot2)
library(viridis)

#Data for this figure were generated using the Broad Terra pipeline
#This pipeline takes Cell Ranger aligned data, and uses Cellbender to remove
#technial artifacts. Cells/nuclei underwent QC and filtering, and were integrated
#using LIGER. See the pipeline code to investigate further.

######################CODE TO GENERATE FIGURE S2 PANELS######################
#load complete dataset with all single-cell and single-nucleus data
gcdata <-
  readRDS(
    "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/combined_dataset_k20.Rds"
  )
output.dir <- output.dir <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/FigS2"

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
write.csv(de.markers.filter, file = paste0(output.dir, "/Table S4 DE genes.csv"))

  
###################### PANEL A-C ########################
pdf(paste0(output.dir ,"/panel S2A-C HTAPP UMAP marker feature plots.pdf"))
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "PTPRC",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(2048, 2048),
  na_cutoff = NULL,
  colors_use = viridis(n = 10, option = "D", alpha = 0.7)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "PECAM1",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(2048, 2048),
  na_cutoff = NULL,
  colors_use = viridis(n = 10, option = "D", alpha = 0.7)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "NCAM1",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(2048, 2048),
  na_cutoff = NULL,
  colors_use = viridis(n = 10, option = "D", alpha = 0.7)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "VIM",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(2048, 2048),
  na_cutoff = NULL,
  colors_use = viridis(n = 10, option = "D", alpha = 0.7)
)
dev.off()

###################### PANEL D ########################
gcdata.subset <-
  subset(gcdata, subset = malignant_calling == "SNV+CNV" &
           annotated_coarse != 'erythrocyte')

pdf(paste0(output.dir, "/malignant calls.pdf"))
Meta_Highlight_Plot(
  seurat_object = gcdata,
  meta_data_column = "annotated_coarse",
  raster = T,
  figure_plot = TRUE,
  raster.dpi = c(1024, 1024),
  meta_data_highlight = "malignant",
  highlight_color = "black", 
  background_color = "lightgrey"
)
dev.off()