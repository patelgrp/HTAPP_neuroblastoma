library(Seurat)
library(dplyr)
library(tidyr)
library(scCustomize)
library(MetBrewer)
library(viridis)

#a function for calculating jaccard similarity for gene sets
jaccard_index <- function(a, b) {
  unique_a <- unique(a)
  unique_b <- unique(b)
  intersection <- length(intersect(unique_a, unique_b))
  union <- length(union(unique_a, unique_b))
  if (union == 0)
    return (NA)
  return(intersection / union)
}

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

#pull top 100 feature embeddings for each NMF cluster
nmf.embeddings <- matrix(nrow = 50, ncol = 20)
colnames(nmf.embeddings) <- paste0("iNMF", 1:20)
for (i in 1:20)
{
  nmf.embeddings[, i] <-
    TopFeatures(object = gcdata[["inmfNorm"]],
                dim = i,
                nfeatures = 50)
}

# jaccard.nmf <-
#   outer(
#     1:ncol(nmf.embeddings),
#     1:ncol(nmf.embeddings),
#     Vectorize(function(i, j)
#       jaccard_index(nmf.embeddings[, i], nmf.embeddings[, j]))
#   )
#
# jaccard.dist.nmf <- 1 - jaccard.nmf
#
# hc.nmf <- hclust(as.dist(jaccard.dist.nmf), method = "average")
# plot(hc.nmf, main = "Hierarchial Clustering Dendogram", xlab = "NMF", sub = "", cex = 0.8)
#
# nmf.embeddings <- nmf.embeddings[, hc.nmf$order]
#pull the top 25 feature loadings for each NMF and merge them (removing duplicates)

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
  "13",
  "20",
  "15"
)
nmf.embeddings <-
  nmf.embeddings[, paste0("iNMF", ordering)]
gene.list <- unique(as.vector(nmf.embeddings[25:50,]))

Idents(gcdata) <- "inmfNorm.cluster"
avg.neuroblastoma <- AverageExpression(gcdata,
                                       return.seurat = T)
custom.order <- paste0("g", ordering)

avg.neuroblastoma@active.ident <-
  factor(avg.neuroblastoma@active.ident, levels = custom.order)
color.scheme <-
  colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(256)

pdf(paste0(output.dir, "\GEP heatmap.pdf"))
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

#pulls cell ids for every sample which had numbat calling
gcdata.subset <-
  subset(gcdata, subset = malignant_calling == "SNV+CNV")

FeaturePlot_scCustom(
  seurat_object = gcdata.subset,
  features = "allele",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(2048, 2048),
  colors_use = viridis(n=10, option = "D", alpha = 0.7)
)
