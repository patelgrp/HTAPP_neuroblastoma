library(Seurat)
library(dplyr)
library(tidyr)
library(scCustomize)
library(MetBrewer)
library(viridis)
library(SCpubr)

#read in combined dataset of integrated neuroblastoma data
gcdata <-
  readRDS(
    "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/combined_dataset_k20.Rds"
  )

#set an output directory
output.dir <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/Fig2"

#print a pdf with a umap colored by iNMF cluster
pdf(paste0(output.dir, "/UMAP by iNMF cluster.pdf"))
print(
  DimPlot_scCustom(
    seurat_object = gcdata,
    figure_plot = TRUE,
    group.by = "inmfNorm.cluster",
    colors_use = met.brewer("Signac", 20)
  )
)
dev.off()

#print a pdf with a umap colored by annotated cluster
pdf(paste0(output.dir, "/UMAP by coarse annotation.pdf"))
print(
  DimPlot_scCustom(
    seurat_object = gcdata,
    figure_plot = TRUE,
    group.by = "annotated_coarse",
    colors_use = met.brewer("Hokusai3", 5)
  )
)
dev.off()

#pull top 50 feature embeddings for each NMF cluster for supplemental table S1
nmf.embeddings <- matrix(nrow = 50, ncol = 20)
colnames(nmf.embeddings) <- paste0("iNMF", 1:20)
for (i in 1:20)
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

#reorder clusters for heatmap
ordering <- c(
  "2",
  "3",
  "7",
  "18",
  "12",
  "14",
  "19",
  "9",
  "10",
  "11",
  "16",
  "17",
  "8" ,
  "6",
  "5",
  "4",
  "1",
  "20",
  "15",
  "13"
)
nmf.embeddings <-
  nmf.embeddings[, paste0("iNMF", ordering)]

#limit genes for heatmap to top 25 from each NMF (and remove duplicates)
gene.list <- unique(as.vector(nmf.embeddings[25:50, ]))

#average expression for each NMF cluster
Idents(gcdata) <- "inmfNorm.cluster"
avg.neuroblastoma <- AverageExpression(gcdata,
                                       return.seurat = T)
custom.order <- paste0("g", ordering)
avg.neuroblastoma@active.ident <-
  factor(avg.neuroblastoma@active.ident, levels = custom.order)

color.scheme <-
  colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(256)

#plot heatmap of major genes
pdf(paste0(output.dir, "/GEP heatmap.pdf"))
print(
  DoHeatmap(
    avg.neuroblastoma,
    features = gene.list,
    draw.lines = F,
    group.colors = met.brewer("Signac", 20)[as.numeric(ordering)]
  ) +
    ggplot2::scale_fill_gradientn(colours = rev(color.scheme))
)
dev.off()

#pulls cell ids for every sample which had numbat calling and remove RBCs
gcdata.subset <-
  subset(gcdata, subset = malignant_calling == "SNV+CNV" &
           annotated_coarse != 'erythrocyte')

#feature plot of allel probabilities for those cells that had numbat calls.
pdf(paste0(output.dir, "/malignant state calls.pdf"))
FeaturePlot_scCustom(
  seurat_object = gcdata.subset,
  features = "allele",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(2048, 2048),
  colors_use = viridis(n = 10, option = "D", alpha = 0.7)
)
dev.off()

pdf(paste0(output.dir, "/malignant state violin plot.pdf"))
VlnPlot_scCustom(
  seurat_object = gcdata.subset,
  features = "allele",
  group.by = "annotated_coarse",
  split.by = "condition",
  plot_boxplot = FALSE,
  pt.size = 0,
  plot_median = FALSE,
  colors_use = c("#6FA7B9", "#D1B07C")
)
dev.off()

pdf(paste0(output.dir, "/barplots.pdf"))
do_BarPlot(
  sample = gcdata.subset,
  split.by = "condition",
  group.by = "annotated_coarse",
  position = "fill",
  colors.use = c(
    "endothelium" = met.brewer("Hokusai3", 5)[1],
    "immune" = met.brewer("Hokusai3", 5)[3],
    "malignant" = met.brewer("Hokusai3", 5)[4],
    "stroma" = met.brewer("Hokusai3", 5)[5]
  )
)
dev.off()

#write table for scCODA analyses
write.csv(
  table(gcdata.subset$annotated_coarse, gcdata.subset$orig.ident),
  file = paste0(output.dir, "/cell counts.csv")
)