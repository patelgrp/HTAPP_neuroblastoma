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
htapp.124 <- readRDS("HTAPP-124-SMP-61_convert_cumulus_Seurat3.Rds")
all.metadata <- read.csv("HTAPP-124-SMP-61_metadata.csv", row.names = 1)
all.metadata <- all.metadata[match(Cells(htapp.124), rownames(all.metadata)),]

htapp.124 <- AddMetaData(htapp.124, metadata = all.metadata$annotate_round3, col.name = "annot.ids")

cols.annot3 <- c('B cell' = "#CC99FF", #lavender
                 'sympathoblast' = "#88CCE9", #cyan
                 'Endothelial' = "#D08C4A", #brown
                 'Erythrocyte' = "#FFCC00", #yellow
                 'Myeloid' = "#996699", #light purple
                 'MES' = "#8A181A", #darkred
                 'Stroma' = "#CC6600", #orange
                 'ADRN' = "#3953A4", #blue
                 'T cell' = "#660099", #deep purple
                 'Zona' = "#00FF33", #green
                 'NA' = "#999999") #grey

pdf("panel S2D HTAPP-124 UMAP.pdf")
print(DimPlot(htapp.124, group.by = "annot.ids", cols = cols.annot3, pt.size = 1))
dev.off()

###################### PANEL E ########################
#the heatmap was generated as part of the inferCNV pipeline (see png file)
#to generate the cell id bar on the observation section, see code below

dend.labels <- read.dendrogram("infercnv.observations_dendrogram.txt") #generated from the inferCNV pipeline

htapp.124 <- subset(gcdata, subset = orig.ident == "HTAPP-124-SMP-61")
htapp.124.labels <- htapp.124$annotate_refine_coarse
htapp.124.labels2 <- data.frame(annot = htapp.124.labels[labels(dend.labels)], cell.id = labels(dend.labels), col = "NA")
htapp.124.labels2$col[htapp.124.labels2$annot == "SCP"] <- "#8A181A"
htapp.124.labels2$col[htapp.124.labels2$annot == "sympathoblast"] <- "#1C83BA"
htapp.124.labels2$col[htapp.124.labels2$annot == "cycling sympathoblast"] <- "#88CCE9"
htapp.124.labels2$col[htapp.124.labels2$annot == "Stroma"] <- "#CC6000"
htapp.124.labels2$col[is.na(htapp.124.labels2$annot)] <- "lightgrey"

htapp.124.labels2$cell.id <- factor(htapp.124.labels2$cell.id, levels = htapp.124.labels2$cell.id)

p1 <- ggdendrogram(dend.labels, labels = F, leaf_labels = F, rotate = F) +
  scale_y_continuous(limits = c(0,75), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,4162), expand = c(0,0))

p2 <- ggplot(htapp.124.labels2, aes(cell.id, y=1, fill = factor(annot))) + geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#8A181A", "#88CCE9", "#1C83BA", "#CC6000", "lightgrey"))


gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)

maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])

gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)

#this produces a dendogram which can be overlaid over the png
pdf("panel S2E HTAPP-124 cell id strip.pdf")
grid.arrange(gp1,gp2,ncol = 1, heights = c(4/5, 1/5))
dev.off()

###################### PANEL F ########################
# This code is adapted from the Geeleher lab. To see the complete implementation, please
# refer to their article - Chapple RH et al BioRxiv (2023) https://doi.org/10.1101/2023.04.13.536639 

# Outer tracks are putative CNV regions (gains = red, losses = blue) for each patient<br>
# Inner track is expression data divided by malignant cell type

hg38.lengths = list("1" = 248956422, #GRCh38
                    "2" = 242193529,
                    "3" = 198295559,
                    "4" = 190214555,
                    "5" = 181538259,
                    "6" = 170805979,
                    "7" = 159345973,
                    "8" = 145138636,
                    "9" = 138394717,
                    "10" = 133797422,
                    "11" = 135086622,
                    "12" = 133275309,
                    "13" = 114364328,
                    "14" = 107043718,
                    "15" = 101991189,
                    "16" = 90338345,
                    "17" = 83257441,
                    "18" = 80373285,
                    "19" = 58617616,
                    "20" = 64444167,
                    "21" = 46709983,
                    "22" = 50818468,
                    "X" = 156040895)

#Read in inferCNV results for HTAPP-124
regions <- read.delim2(file = "HMM_CNV_predictions.HMMi6.qnorm.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat")
regions$sample_id <- "HTAPP-124-SMP-61"

# Remove losses on Chr6 in test cells
# This is a result of using immune cells as the reference due to HLA copy number differences in these cells
regions <- subset(regions, !(chr == "6" & state < 3))
regions <- subset(regions, !(chr == "Y"))

#Define colors for CNV arcs
#Losses (i.e. states 1 and 2) are blue
#Gains (i.e. states > 3) are red
#Opacities for arcs
#Greater CN gains/loss are less opaque

regions$cols <- NA
try(regions[regions$state == 1,]$cols <- "#2166AC")
try(regions[regions$state == 2,]$cols <- "#92C5DE")
try(regions[regions$state == 3,]$cols <- "#F7F7F7")
try(regions[regions$state == 4,]$cols <- "#F4A582")
try(regions[regions$state == 5,]$cols <- "#D6604D")
try(regions[regions$state == 6,]$cols <- "#B2182B")

#Opacities for arcs
#Greater CN gains/loss are less opaque
regions$ops <- NA
try(regions[regions$state == 1,]$ops <- 1)
try(regions[regions$state == 2,]$ops <- 0.5)
try(regions[regions$state == 3,]$ops <- 0.3)
try(regions[regions$state == 4,]$ops <- 0.3)
try(regions[regions$state == 5,]$ops <- 0.5)
try(regions[regions$state == 6,]$ops <- 1)

#Define ARC Tracklists
tracklist <- NULL

#Sets the position for the first track
minRadius = 1.16

maxRadius <- minRadius + 0.04

track <- BioCircosArcTrack(trackname = "HTAPP-124 arc",
                           chromosomes = as.character(regions$chr),
                           starts = regions$start, ends = regions$end,
                           minRadius = minRadius, maxRadius = maxRadius, 
                           colors =  regions$cols, opacities = regions$ops)
tracklist <- tracklist + track

widgetBioCircos <- BioCircos(tracklist, yChr = FALSE, genome = hg38.lengths,
                             genomeFillColor = viridis::magma(23, direction = -1),
                             chrPad = 0.05,
                             displayGenomeBorder = F, genomeTicksDisplay = F, genomeLabelDisplay = F, zoom = T, 
                             ARCMouseOverDisplay = T, ARCMouseOverTooltipsHtml04 = "<br/>Patient ID: ")
#Obtain list of genes (filtered by those used in inferCNV analysis), along with their chromosomal positions

infercnv_obj <- readRDS("run.final.infercnv_obj")
gene_order <- infercnv_obj@gene_order
gene_order$chr_num <- as.numeric(gene_order$chr)
chr6_gene_ind <- gene_order %>% filter(chr == "6" & start > 26055967 & stop < 56951643)
gene_order <- gene_order[!(rownames(gene_order) %in% rownames(chr6_gene_ind)), ]
gene_order <- subset(gene_order, subset = chr != "Y")

ref_ind <- unlist(infercnv_obj@reference_grouped_cell_indices)
obs_ind <- unlist(infercnv_obj@observation_grouped_cell_indices)

exprmat <- infercnv_obj@expr.data 

Idents(htapp.124) <- "annotate_refine_coarse"
cell.types <- unique(htapp.124$annotate_refine_coarse)

topProgCells <- list()
topProgCell_expr <- list()
topProgCell_expr_vals <- list()

for(j in 1:length(cell.types))
{
  topProgCells[[j]] <- WhichCells(htapp.124, idents = cell.types[j], downsample = 1000)
  topProgCell_expr[[j]] <- t(exprmat[, topProgCells[[j]]])
  topProgCell_expr_vals[[j]] <- apply(topProgCell_expr[[j]], 2, mean)
}

Idents(htapp.124) <- "annotate_refine_1"
randomCells <- WhichCells(htapp.124, 
                          idents = c("Malignant"),
                          invert = T)
randomCell_expr <- t(exprmat[, randomCells])
randomCell_expr_vals <- data.frame(random = apply(randomCell_expr, 2, mean))
Idents(htapp.124) <- "annotate_refine_coarse"

topProgCell_expr_vals <- as.data.frame(do.call(cbind, topProgCell_expr_vals))
colnames(topProgCell_expr_vals) <- cell.types

expr_vals <- as.data.frame(cbind(topProgCell_expr_vals, randomCell_expr_vals))

expr_vals <- expr_vals[, names(expr_vals) %in% c("adrenergic", "mesenchymal", "sympathoblast", "random")]

expr_vals <- expr_vals[rownames(expr_vals) %in% rownames(gene_order),]

chrs <- as.character(gene_order$chr)
starts <- gene_order$start
ends <- gene_order$stop
chrs[chrs == 23] <- "X"

tracklist <- track

for (k in c(4,3,1,2))
{
  if (colnames(expr_vals)[k] == "adrenergic") {
    snpcol = "darkblue"
    opacity = 1
  }
  if (colnames(expr_vals)[k] == "sympathoblast") {
    snpcol = "cyan"
    opacity = 1
  }
  if (colnames(expr_vals)[k] == "mesenchymal") {
    snpcol = "darkred"
    opacity = 1
  }
  if (colnames(expr_vals)[k] == "random") {
    snpcol = "grey"
    opacity = 0.3
  }
  
  tracklist <- tracklist + BioCircosSNPTrack("expr", chrs, starts, expr_vals[,k], 
                                             size = 1,
                                             maxRadius = 1.0, minRadius = 0.6, range = c(0.8,1.2),
                                             color = snpcol, opacities = opacity)
}

tracklist <- tracklist + BioCircosBackgroundTrack("myBackgroundTrack_testcell", 
                                                  maxRadius = 1.0, minRadius = 0.6, 
                                                  borderColors = "#E8FDE8", borderSize = 0.3, fillColors = "#F7F7F7")  

widgetBioCircos <- BioCircos(tracklist, yChr = FALSE, genome = hg38.lengths,
                             genomeFillColor = viridis::magma(23, direction = -1),
                             chrPad = 0.05,
                             displayGenomeBorder = F, genomeTicksDisplay = F, genomeLabelDisplay = F, zoom = T, 
                             ARCMouseOverDisplay = T, ARCMouseOverTooltipsHtml04 = "<br/>Patient ID: ")

htmlwidgets::saveWidget(widgetBioCircos, "panel S2F BioCircos.html")

