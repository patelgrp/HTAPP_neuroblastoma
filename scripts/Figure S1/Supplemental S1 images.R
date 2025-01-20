library(Seurat)
library(dplyr)
library(tidyr)
library(scCustomize)
library(SCpubr)
library(MetBrewer)
library(ggplot2)
library(ggeasy)

#Data for this figure were generated using the Broad Terra pipeline
#This pipeline takes Cell Ranger aligned data, and uses Cellbender to remove
#technial artifacts. Cells/nuclei underwent QC and filtering, and were integrated
#using LIGER. See the pipeline code to investigate further.

######################CODE TO GENERATE FIGURE S1 PANELS######################
#load complete dataset with all single-cell and single-nucleus data
gcdata <-
  readRDS(
    "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/combined_dataset_k20.Rds"
  )
output.dir <- output.dir <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/FigS1"

###Panel A is a schematic

################### PANEL B ########################
#UMAP plot colored by sample
color.scheme <- MetBrewer::met.brewer("Hokusai1", n = 56)
color.scheme <- sample(color.scheme)
#note sample 11 and 12 are duplicates
color.scheme[11] <- color.scheme[12]

pdf(paste0(output.dir, "/panel S1B All HTAPP plot by sample.pdf"))
print(
  DimPlot_scCustom(
    seurat_object = gcdata,
    figure_plot = TRUE,
    group.by = "orig.ident",
    colors_use = color.scheme
  )
)
print(
  DimPlot_scCustom(
    seurat_object = gcdata,
    figure_plot = TRUE,
    group.by = "orig.ident",
    colors_use = color.scheme,
  ) & NoLegend()
)
dev.off()

################### PANEL C ########################
sample.ids <- sort(unique(gcdata$orig.ident))
Idents(gcdata) <- "orig.ident"

meta.data <- Fetch_Meta(gcdata)

plot.list <- list()
for (i in c(1:11, 13:56))
  #12 is a duplicate value
{
  cell.highlight <- WhichCells(gcdata, idents = sample.ids[i])
  if (meta.data[cell.highlight[1], ]$condition == "cell")
  {
    title.color <- "#6FA8BA"
  }
  else{
    title.color <- "#D2B17C"
  }
  #cell.highlight <- list(sample = cell.highlight)
  plot.list[[i]] <-
    Meta_Highlight_Plot(
      seurat_object = gcdata,
      meta_data_column = "orig.ident",
      raster = F,
      figure_plot = TRUE,
      meta_data_highlight = sample.ids[i],
      highlight_color = color.scheme[i]
    )
  plot.list[[i]][[1]] <-
    plot.list[[i]][[1]] + NoLegend() + ggtitle(sample.ids[i]) + easy_all_text_color(color = title.color) + easy_center_title()
  plot.list[[i]] <- plot.list[[i]][[1]]
}

png(
  paste0(output.dir, "/panel S1C single sample UMAPs 1.png"),
  units = "in",
  res = 600,
  width = 24,
  height = 15
)
#patchwork::wrap_plots(plot.list[c(1:11, 13:56)])
patchwork::wrap_plots(plot.list[c(1:11, 17, 23:30, 36:43, 49:56)], ncol = 8)
dev.off()

png(
  paste0(output.dir, "/panel S1C single sample UMAPs 2.png"),
  units = "in",
  res = 600,
  width = 15,
  height = 9
)
patchwork::wrap_plots(plot.list[c(18:22, 31:35, 44:48)], ncol = 5)
dev.off()


###################### PANEL D ######################
plot.cells <-
  Meta_Highlight_Plot(
    seurat_object = gcdata,
    meta_data_column = "condition",
    raster = T,
    figure_plot = TRUE,
    raster.dpi = c(1024, 1024),
    meta_data_highlight = "cell",
    highlight_color = "#6FA8BA"
  )
plot.cells <-
  plot.cells[[1]] + NoLegend() + ggtitle("single-cell RNA-seq") + easy_all_text_color(color = "#6FA8BA") + easy_center_title()

plot.nucleus <-
  Meta_Highlight_Plot(
    seurat_object = gcdata,
    meta_data_column = "condition",
    raster = T,
    figure_plot = TRUE,
    raster.dpi = c(1024, 1024),
    meta_data_highlight = "nucleus",
    highlight_color = "#D2B17C"
  )
plot.nucleus <-
  plot.nucleus[[1]] + NoLegend() + ggtitle("single-nucleus RNA-seq") + easy_all_text_color(color = "#D2B17C") + easy_center_title()

png(
  paste0(output.dir, "/panel S1D UMAP plots by condition.png"),
  units = "in",
  res = 600,
  width = 12,
  height = 6
)
patchwork::wrap_plots(list(plot.cells, plot.nucleus), ncol = 2)
dev.off()

###################### PANELS E-H ######################
#These lines of code  apply stress signatures from Biermann et al Cell (2022) and dissociation signature from van den Brink et al Nat
#Methods (2017)
stress.sigs <-
  read.csv(file = "./stress_signatures.csv",
           fill = F,
           strip.white = T)

gcdata.flat <- JoinLayers(gcdata)
gcdata.flat <-
  AddModuleScore(
    gcdata.flat,
    features = list(stress.sigs$stress_module_Biermann[1:31]),
    name = "stress_Biermann"
  )
gcdata.flat <-
  AddModuleScore(gcdata.flat,
                 features = list(stress.sigs$stress_vanBrink),
                 name = "stress_vanBrink")

pdf(paste0(output.dir, "/panel S1E-H Violin plot All HTAPP QC metrics.pdf"))
VlnPlot_scCustom(
  seurat_object = gcdata.flat,
  features = "nFeature_RNA",
  group.by = "condition",
  pt.size = 0,
  plot_boxplot = T,
  colors_use = c("cell" = "#6FA8BA", "nucleus" = "#D2B17C")
)
VlnPlot_scCustom(
  seurat_object = gcdata.flat,
  features = "percent_mito",
  group.by = "condition",
  pt.size = 0,
  plot_boxplot = T,
  colors_use = c("cell" = "#6FA8BA", "nucleus" = "#D2B17C")
)
VlnPlot_scCustom(
  seurat_object = gcdata.flat,
  features = "stress_Biermann1",
  group.by = "condition",
  pt.size = 0,
  plot_boxplot = T,
  colors_use = c("cell" = "#6FA8BA", "nucleus" = "#D2B17C")
)
VlnPlot_scCustom(
  seurat_object = gcdata.flat,
  features = "stress_vanBrink1",
  group.by = "condition",
  pt.size = 0,
  plot_boxplot = T,
  colors_use = c("cell" = "#6FA8BA", "nucleus" = "#D2B17C")
)
dev.off()

#statistics using Wilcox rank sum to compare cell vs nucleus QC metrics
stats.frame <-
  data.frame(
    genes = gcdata.flat$nFeature_RNA,
    mito = gcdata.flat$percent_mito,
    stress = gcdata.flat$stress_Biermann1,
    dissoc = gcdata.flat$stress_vanBrink1,
    condition = gcdata.flat$condition
  )
stats.cell <- subset(stats.frame, subset = condition == "cell")
stats.nucleus <-
  subset(stats.frame, subset = condition == "nucleus")
test1 <- wilcox.test(stats.cell$genes, stats.nucleus$genes)
test2 <- wilcox.test(stats.cell$mito, stats.nucleus$mito)
test3 <- wilcox.test(stats.cell$stress, stats.nucleus$stress)
test4 <- wilcox.test(stats.cell$dissoc, stats.nucleus$dissoc)
