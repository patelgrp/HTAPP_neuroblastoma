library(Seurat)
library(dplyr)
library(tidyr)
library(scCustomize)
library(SCpubr)
library(MetBrewer)

#Data for this figure were generated using the Broad Terra pipeline
#This pipeline takes Cell Ranger aligned data, and uses Cellbender to remove
#technial artifacts. Cells/nuclei underwent QC and filtering, and were integrated
#using LIGER. See the pipeline code to investigate further.

######################CODE TO GENERATE FIGURE S1 PANELS######################
#load complete dataset with all single-cell and single-nucleus data
gcdata <- readRDS("/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/combined_dataset_k20.Rds")
output.dir <- output.dir <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/FigS1"

###Panel A is a schematic

################### PANEL B ########################
#UMAP plot colored by sample
color.scheme <- MetBrewer::met.brewer("Hokusai1", n=56)
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

################### PANEL B ########################
#note: this code will not run in Seurat v5!
sample.ids <- sort(unique(gcdata$orig.ident))
Idents(gcdata) <- "orig.ident"

plot.list <- list()
for (i in c(1:11, 13:57))#12 is a duplicate value 
{
  cell.highlight <- WhichCells(gcdata, idents = sample.ids[i])
  plot.list[[i]] <- Meta_Highlight_Plot(seurat_object = gcdata, meta_data_column = "orig.ident", 
                                        meta_data_highlight = sample.ids[i], highlight_color = color.scheme[i]) + NoLegend() + ggplot2::ggtitle(sample.ids[i])
}

pdf("panel S1G single sample UMAPs.pdf", width = 30, height = 30)
patchwork::wrap_plots(plot.list[c(1:11,13:57)])
dev.off()



###################### PANEL B ######################
Idents(gcdata) <- "condition"

pdf("panel S1B UMAP plots by condition.pdf")
print(DimPlot(gcdata, cells = WhichCells(gcdata, idents = "cell"), cols = c("#6FA789"), 
              pt.size = 1, raster = T, raster.dpi = c(1024, 1024)) + NoLegend())
print(DimPlot(gcdata, cells = WhichCells(gcdata, idents = "nucleus"), cols = c("#D1B07C"),
              pt.size = 1, raster = T, raster.dpi = c(1024, 1024)) + NoLegend())
dev.off()

###################### PANELS C-E ######################
#These lines of code  apply stress signatures from Biermann et al Cell (2022) and dissociation signature from van den Brink et al Nat 
#Methods (2017)
stress.sigs <- read.csv(file = "./stress_signatures.csv", fill = F, strip.white = T)
gcdata <- AddModuleScore(gcdata, features = list(stress.sigs$stress_module_Biermann[1:31]), name = "stress_Biermann")
gcdata <- AddModuleScore(gcdata, features = list(stress.sigs$stress_vanBrink), name = "stress_vanBrink")

pdf("panel S1C-E Violin plot All HTAPP QC metrics.pdf")
SCpubr::do_ViolinPlot(sample = gcdata, feature = c("nFeature_RNA"), 
                      group.by = "condition",
                      colors.use = c("cell" = "#70A9BA","nucleus" = "#D2B27C"))
SCpubr::do_ViolinPlot(sample = gcdata, feature = c("percent_mito"), 
                      group.by = "condition",
                      colors.use = c("cell" = "#70A9BA","nucleus" = "#D2B27C"))
SCpubr::do_ViolinPlot(sample = gcdata, feature = c("stress_Biermann1"), 
                      group.by = "condition",
                      colors.use = c("cell" = "#70A9BA","nucleus" = "#D2B27C"))
SCpubr::do_ViolinPlot(sample = gcdata, feature = c("stress_vanBrink1"), 
                      group.by = "condition",
                      colors.use = c("cell" = "#70A9BA","nucleus" = "#D2B27C"))
dev.off()

#statistics using Wilcox rank sum to compare cell vs nucleus QC metrics
stats.frame <- data.frame(genes = gcdata$nFeature_RNA, mito = gcdata$percent_mito, 
                          stress = gcdata$stress_Biermann1, dissoc = gcdata$stress_vanBrink1, condition = gcdata$condition)
stats.cell <- subset(stats.frame, subset = condition == "cell")
stats.nucleus <- subset(stats.frame, subset = condition == "nucleus")
test1 <- wilcox.test(stats.cell$genes, stats.nucleus$genes)
test2 <- wilcox.test(stats.cell$mito, stats.nucleus$mito)
test3 <- wilcox.test(stats.cell$stress, stats.nucleus$stress)
test4 <- wilcox.test(stats.cell$dissoc, stats.nucleus$dissoc)




