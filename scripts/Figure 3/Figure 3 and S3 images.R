library(Seurat)
library(dplyr)
library(SingleR)
library(ggplot2)
library(RColorBrewer)
library(rliger)
library(SCpubr)
library(UCell)

#Data for this figure were generated using the Broad Terra pipeline
#This pipeline takes Cell Ranger aligned data, and uses Cellbender to remove
#technical artifacts. Cells/nuclei underwent QC and filtering, and were integrated
#using LIGER. See the pipeline code to investigate further.

######################CODE TO GENERATE FIGURE 3 and S3 PANELS######################
#load complete dataset with all single-cell and single-nucleus data (only malignant cells/nuclei)
gcdata <-
  readRDS(
    "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/malignant_combined_dataset_k12.Rds"
  )
output.dir <- output.dir <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/Fig3"

#cluster 7 is a minor cluster contaminated with doublets expressing tumor + macrophage markers (PTPRC+, CD68+, LYZ+)
gcdata <-
  subset(gcdata, subset = inmfNorm.cluster == 7, invert = TRUE)

#reprocesss data
gcdata <- normalize(gcdata)
gcdata <- selectGenes(gcdata)
gcdata <- scaleNotCenter(gcdata)
gcdata

gcdata <- runINMF(gcdata, k = 12)
gcdata <- quantileNorm(gcdata)
gcdata

gcdata <-
  RunUMAP(gcdata, reduction = "inmfNorm", dims = 1:12)

saveRDS(gcdata,
        file = paste0(output.dir, "malignant_combined_dataset_k12_reprocess.Rds"))

###################### PANEL A ######################
pdf(paste0(output.dir, "/panel 3A 20220714 HTAPP malignant liger.pdf"))
print(
  DimPlot_scCustom(
    seurat_object = gcdata,
    figure_plot = TRUE,
    group.by = "inmfNorm.cluster",
    colors_use = met.brewer("Signac", 12)
  )
)
dev.off()

#pull top 50 feature embeddings for each NMF cluster for supplemental table S4
nmf.embeddings <- matrix(nrow = 50, ncol = 12)
colnames(nmf.embeddings) <- paste0("iNMF", 1:12)
for (i in 1:12)
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



###################### PANEL B ######################
#reorder clusters for heatmap
ordering <- c("3",
              "7",
              "2",
              "8" ,
              "6",
              "12",
              "4",
              "10",
              "5",
              "9",
              "11",
              "1")
nmf.embeddings <-
  nmf.embeddings[, paste0("iNMF", ordering)]

#limit genes for heatmap to top 25 from each NMF (and remove duplicates)
gene.list <- unique(as.vector(nmf.embeddings[25:50,]))

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
pdf(paste0(output.dir, "/panel 3B GEP heatmap.pdf"))
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


###################### PANEL C ######################
pdf(paste0(output.dir, "/panel 3C HTAPP malignant feature plots.pdf"))
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "VIM",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024),
  na_cutoff = NULL,
  max.cutoff = 1.5,
  order = TRUE,
  colors_use = viridis(n = 10, option = "D", alpha = 0.5)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "B2M",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024),
  na_cutoff = NULL,
  max.cutoff = 3,
  order = TRUE,
  colors_use = viridis(n = 10, option = "D", alpha = 0.5)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "TH",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024),
  na_cutoff = NULL,
  max.cutoff = 1.5,
  order = TRUE,
  colors_use = viridis(n = 10, option = "D", alpha = 0.5)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "DBH",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024),
  na_cutoff = NULL,
  max.cutoff = 2,
  order = TRUE,
  colors_use = viridis(n = 10, option = "D", alpha = 0.5)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "MKI67",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024),
  na_cutoff = NULL,
  max.cutoff = 1.5,
  order = TRUE,
  colors_use = viridis(n = 10, option = "D", alpha = 0.5)
)
FeaturePlot_scCustom(
  seurat_object = gcdata,
  features = "TOP2A",
  figure_plot = T,
  pt.size = 2,
  raster = T,
  raster.dpi = c(1024, 1024),
  na_cutoff = NULL,
  max.cutoff = 2,
  order = TRUE,
  colors_use = viridis(n = 10, option = "D", alpha = 0.5)
)
dev.off()

###################### PANEL D-G ######################
signatures <-
  read.csv("ADRN_MES_signatures.csv")  #signatures from von Groningen 2017 Nat Genet
ADRN.genes <- signatures$ADRENERGIC
MES.genes <- signatures$MESENCHYMAL
ADRN.genes <- ADRN.genes[ADRN.genes != ""]
MES.genes <- MES.genes[MES.genes != ""]

Idents(gcdata) <- "inmfNorm.cluster"
gcdata@active.ident <-
  factor(gcdata@active.ident, levels = ordering)

gcdata <- JoinLayers(gcdata)
gcdata <-
  AddModuleScore_UCell(obj = gcdata,
                       features = list(ADRN = ADRN.genes, MES = MES.genes))
# gcdata <- AddModuleScore(object = gcdata,
#                          features = list(ADRN.genes),
#                          name = "ADRN")
# gcdata <- AddModuleScore(object = gcdata,
#                          features = list(MES.genes),
#                          name = "MES")

p <-
  do_EnrichmentHeatmap(
    sample = gcdata,
    input_gene_list = list("ADRN" = ADRN.genes, "MES" = MES.genes),
    viridis.direction = -1,
    group.by = "inmfNorm.cluster",
    flip = TRUE
  )
p2 <-
  do_EnrichmentHeatmap(
    sample = gcdata,
    input_gene_list = list("ADRN" = ADRN.genes, "MES" = MES.genes),
    viridis.direction = -1, flavor = "UCell,"
    group.by = "inmfNorm.cluster",
    flip = TRUE
  )

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

###Stats
adrn.frame <-
  data.frame(adrn = gcdata$ADRN1, cluster = gcdata$annot.ids)
adrn.MES <- subset(adrn.frame, subset = cluster == "MES")
adrn.ADRN <- subset(adrn.frame, subset = cluster == "ADRN")
adrn.symp <- subset(adrn.frame, subset = cluster == "sympathoblast")

test1 <- wilcox.test(adrn.MES$adrn, adrn.ADRN$adrn, exact = FALSE)
test2 <- wilcox.test(adrn.MES$adrn, adrn.symp$adrn, exact = FALSE)
test3 <- wilcox.test(adrn.ADRN$adrn, adrn.symp$adrn, exact = FALSE)

mes.frame <-
  data.frame(mes = gcdata$MES1, cluster = gcdata$annot.ids)
mes.MES <- subset(mes.frame, subset = cluster == "MES")
mes.ADRN <- subset(mes.frame, subset = cluster == "ADRN")
mes.symp <- subset(mes.frame, subset = cluster == "sympathoblast")

test4 <- wilcox.test(mes.MES$mes, mes.ADRN$mes, exact = FALSE)
test5 <- wilcox.test(mes.MES$mes, mes.symp$mes, exact = FALSE)
test6 <- wilcox.test(mes.ADRN$mes, mes.symp$mes, exact = FALSE)

###################### PANEL H ######################
#pySCENIC was run in the command line (see Methods section)

#subset cells from malignant NB data
cells.gcdata <- subset(gcdata, subset = condition == "cell")

#read in SCENIC matrix for cells
scenic.cells <-
  t(read.csv(file = "cells_mal_auc_nomask.csv", row.names = 1))
cells.gcdata <- subset(cells.gcdata,
                       subset = annot.ids == "MES" |
                         annot.ids == "ADRN" |
                         annot.ids == "sympathoblast")

cells.gcdata[["TF"]] <- CreateAssayObject(data = scenic.cells)

DefaultAssay(cells.gcdata) <- "TF"
tf.markers <-
  FindAllMarkers(
    object = cells.gcdata,
    assay = "TF",
    logfc.threshold = 0.05,
    only.pos = T
  )

tf.pick <-
  subset(tf.markers, subset = p_val_adj < 1E-100 & avg_log2FC > 0)
tf.pick <- as.data.frame(tf.pick)
tf.pick <- tf.pick[!duplicated(tf.pick$gene), ]
tf.pick <- tf.pick[order(match(tf.pick$cluster, custom.order)), ]

avg.tf.neuroblastoma <-
  AverageExpression(
    object = cells.gcdata,
    return.seurat = T,
    group.by = "annot.ids",
    slot = "data"
  )

group.colors <- c("#2E307B", "#009D73", "#89191B")

pdf("panel 3H SCENIC TF heatmap.pdf",
    width = 10,
    height = 10)
DoHeatmap(
  avg.tf.neuroblastoma,
  features = c(
    "ELK4...",
    "EGR3...",
    "SOX10...",
    "MAFF...",
    "IRF1...",
    "IRF2...",
    "IRF3....",
    "FLI1...",
    "RUNX2...",
    "FOSL1...",
    "FOSL2...",
    "HAND2...",
    "PHOX2A...",
    "SOX11...",
    "GATA2...",
    "KLF7...",
    "TFAP2B...",
    "ISL1...",
    "PHOX2B...",
    "PBX3...",
    "EZH2...",
    "MYCN..."
  ),
  disp.min = 0,
  assay = "TF",
  draw.lines = F,
  group.colors = group.colors
) + ggplot2::scale_fill_gradient2(
  low = "white",
  mid = "yellow",
  high = "darkred",
  midpoint = 0.6
)
dev.off()


### Panel I are IHC images

###################### PANEL J and K ######################
#imports data from Jansky et al Nat Genetics (2021)
Westermann.adrenal <- readRDS("adrenal_medulla_Seurat.RDS")

Westermann.matrix <- Westermann.adrenal@assays$RNA@data
Westermann.labels <- Westermann.adrenal@active.ident
gcdata.matrix <- gcdata@assays$RNA@data

#the cell-by-cell transfer takes a soul-crushing amount of time. I've saved the output to an RDS object to speed things along.
#pred <- SingleR(test = gcdata.sample.matrix, ref = Westermann.matrix, labels = Westermann.labels, de.method = "wilcox")
pred <- readRDS("predicted cell by cell status.Rds")

#We also ran SingleR on a cluster-by-cluster basis for panel K
pred2 <-
  SingleR(
    test = gcdata.matrix,
    ref = Westermann.matrix,
    labels = Westermann.labels,
    de.method = "wilcox",
    clusters = gcdata$liger_clusters
  )

gcdata <-
  AddMetaData(gcdata,
              metadata = pred$pruned.labels,
              col.name = "SingleR.mal")
gcdata <- AddMetaData(gcdata, metadata = pred$scores)

pdf("panel 3J SingleR predictions using Jansky data.pdf")
print(DimPlot(
  gcdata,
  group.by = "SingleR.mal",
  pt.size = 0.5,
  cols = c(
    "#FAA51A",
    "#718191",
    "#010101",
    "#059E74",
    "#D4434A",
    "#424B52",
    "#2D317B",
    "#8B191B",
    "#88CCE9",
    "#ED2224",
    "#A7A9AC"
  ),
  raster.dpi = c(1024, 1024)
))
dev.off()

write.csv(prop.table(table(gcdata$SingleR.mal)), "panel 3J percentages.csv")

pdf("panel 3K SingleR heatmap.pdf")
print(plotScoreHeatmap(
  pred2,
  show.pruned = F,
  show_colnames = T,
  cluster_cols = T
))
dev.off()

###################### PANEL L ######################
HTAPP.194.geneedit <- readRDS("HTAPP-194-SMP-251.Rds")

#now, we'll import the PDX dataset
MAST97 <- readRDS("MAST97_1_1_sn_sn_nocellbender.Rds")

#remove cluster 9, which are mouse cells
MAST97 <- subset(MAST97, subset = RNA_snn_res.0.4 != "9")
MAST97$orig.ident <- "MAST97"

#use Seurat's transfer mapping functions
features <-
  SelectIntegrationFeatures(object.list = list(HTAPP.194.geneedit, MAST97),
                            nfeatures = 3000)
anchors <-
  FindTransferAnchors(
    reference = HTAPP.194.geneedit,
    query = MAST97,
    features = features,
    reference.reduction = "pca",
    reduction = "rpca"
  )

MAST97 <-
  MapQuery(
    anchorset = anchors,
    reference = HTAPP.194.geneedit,
    query = MAST97,
    refdata = "annotate_round3",
    reference.reduction = "pca",
    reduction.model = "umap",
    transferdata.args = list(prediction.assay = T)
  )

transfer.scores <-
  TransferData(anchorset = anchors, refdata = HTAPP.194.geneedit$annotate_round3)
MAST97 <- AddMetaData(MAST97, metadata = transfer.scores)

p1 <-
  DimPlot(HTAPP.194.geneedit, group.by = "annotate_round3") + xlim(c(-12, 17)) + ylim(c(-16, 6))
p2 <-
  DimPlot(
    MAST97,
    reduction = "ref.umap",
    group.by = "predicted.id",
    cols = c("orange", "orange")
  ) + xlim(c(-12, 17)) + ylim(c(-16, 6))
p3 <-
  DimPlot(HTAPP.194.geneedit,
          group.by = "orig.ident",
          cols = c("black")) + xlim(c(-12, 17)) + ylim(c(-16, 6))

pdf("panel 3L comparison of HTAPP-194 and MAST97.pdf",
    height = 10,
    width = 20)
print(p1 + p2)
print(p3 + p2)
dev.off()

###################### FIGURE S3 ######################
###################### PANEL A ########################
#UMAP plot colored by sample
color.scheme <- MetBrewer::met.brewer("Hokusai1", n = 57)
color.scheme <- sample(color.scheme)
#note sample 12 and 13 are duplicates
color.scheme[12] <- color.scheme[13]
pdf("panel S3A All HTAPP plot by sample.pdf")
print(
  DimPlot(
    gcdata,
    group.by = "orig.ident",
    cols = color.scheme,
    label = F,
    pt.size = 2,
    raster = T,
    shuffle = T,
    raster.dpi = c(2048, 2048)
  ) + NoLegend()
)
print(
  DimPlot(
    gcdata,
    group.by = "orig.ident",
    cols = color.scheme,
    label = F,
    pt.size = 2,
    raster = T,
    shuffle = T,
    raster.dpi = c(2048, 2048)
  )
)
dev.off()

#confusion matrix
all.metadata <- readRDS("combined unfiltered metadata full.Rds")
all.metadata <- all.metadata[Cells(gcdata), ]
gcdata <-
  AddMetaData(gcdata,
              metadata = all.metadata$liger_clusters,
              col.name = "original.NMF")


confusion_matrix <-
  table(gcdata$liger_clusters, gcdata$original.NMF)
confusion_matrix <- confusion_matrix[c("4", "6", "11",
                                       "12", "2", "8", "7", "10", "3", "9", "5",
                                       "1"),
                                     c(
                                       "1",
                                       "8",
                                       "16",
                                       "20",
                                       "9",
                                       "4",
                                       "18",
                                       "11",
                                       "7",
                                       "14",
                                       "12",
                                       "15",
                                       "10",
                                       "17",
                                       "13",
                                       "2",
                                       "6",
                                       "3",
                                       "5",
                                       "19"
                                     )]
confusion_prop <- prop.table(confusion_matrix, margin = 1)
write.csv(confusion_matrix, "NMF confusion matrix.csv")
write.csv(confusion_prop, "NMF confusion matrix proportions.csv")

###################### PANEL B ########################
Westermann.adrenal <- UpdateSeuratObject(Westermann.adrenal)
pdf("panel S3A Westermann adrenal medulla.pdf",
    width = 10,
    height = 10)
print(DimPlot(
  Westermann.adrenal,
  cols = c(
    "#8A181A",
    "red",
    "#D3424A",
    "orange",
    "black",
    "#708090",
    "#414A51",
    "#009E73",
    "#2A2D7C",
    "#88CCE9"
  ),
  pt.size = 2
))
dev.off()

###################### PANEL C ########################
pdf("panel S3B predicted Westermann scores.pdf",
    width = 10,
    height = 10)
print(
  FeaturePlot(
    gcdata,
    features = c("SCPs"),
    pt.size = 1,
    min.cutoff = 0,
    max.cutoff = 0.5
  ) + viridis::scale_color_viridis(alpha = 0.7)
)
print(
  FeaturePlot(
    gcdata,
    features = c("late SCPs"),
    pt.size = 1,
    min.cutoff = 0,
    max.cutoff = 0.5
  ) + viridis::scale_color_viridis(alpha = 0.7)
)
print(
  FeaturePlot(
    gcdata,
    features = c("cycling SCPs"),
    pt.size = 1,
    min.cutoff = 0,
    max.cutoff = 0.5
  ) + viridis::scale_color_viridis(alpha = 0.7)
)
print(
  FeaturePlot(
    gcdata,
    features = c("Bridge"),
    pt.size = 1,
    min.cutoff = 0,
    max.cutoff = 0.5
  ) + viridis::scale_color_viridis(alpha = 0.7)
)
print(
  FeaturePlot(
    gcdata,
    features = c("connecting Chromaffin cells"),
    pt.size = 1,
    min.cutoff = 0,
    max.cutoff = 0.5
  ) + viridis::scale_color_viridis(alpha = 0.7)
)
print(
  FeaturePlot(
    gcdata,
    features = c("Chromaffin cells"),
    pt.size = 1,
    min.cutoff = 0,
    max.cutoff = 0.5
  ) + viridis::scale_color_viridis(alpha = 0.7)
)
print(
  FeaturePlot(
    gcdata,
    features = c("late Chromaffin cells"),
    pt.size = 1,
    min.cutoff = 0,
    max.cutoff = 0.5
  ) + viridis::scale_color_viridis(alpha = 0.7)
)
print(
  FeaturePlot(
    gcdata,
    features = c("Neuroblasts"),
    pt.size = 1,
    min.cutoff = 0,
    max.cutoff = 0.5
  ) + viridis::scale_color_viridis(alpha = 0.7)
)
print(
  FeaturePlot(
    gcdata,
    features = c("late Neuroblasts"),
    pt.size = 1,
    min.cutoff = 0,
    max.cutoff = 0.5
  ) + viridis::scale_color_viridis(alpha = 0.7)
)
print(
  FeaturePlot(
    gcdata,
    features = c("cycling Neuroblasts"),
    pt.size = 1,
    min.cutoff = 0,
    max.cutoff = 0.5
  ) + viridis::scale_color_viridis(alpha = 0.7)
)
dev.off()

###################### PANEL D ########################
gcdata$late_SCPs <-
  gcdata$`late SCPs` #renames the metadata because empty spaces much up the FeatureScatter command

pdf(
  "panel S3C scatter plots of similarity score vs van Groningen scores.pdf",
  width = 10,
  height = 10
)
print(FeatureScatter(
  gcdata,
  feature1 = "late_SCPs",
  feature2 = "MES1",
  shuffle = T,
  cols = c("black", "black", "black")
))
print(FeatureScatter(
  gcdata,
  feature1 = "late_SCPs",
  feature2 = "ADRN1",
  shuffle = T,
  cols = c("black", "black", "black")
))
dev.off()


###################### PANEL E-F ########################
adr_all <- readRDS(file = "Kildisiute_adrenal.rds")

adr_all <- UpdateSeuratObject(adr_all)
DimPlot(adr_all)

#all of this code syncs the annotations with the published ones in Kildisiute et al
adr_all@meta.data$new_clust = as.character(adr_all@active.ident)
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust %in% c("25"))] =
  "SCPs"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust %in% c("13"))] =
  "Chromaffin"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust %in% c("24"))] =
  "Bridge"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust %in% c("19", "20"))] =
  "Sympathoblastic"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust %in% c("3", "15", "14", "17"))] =
  "Endothelium"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust %in% c("23", "11", "18"))] =
  "Mesenchyme"
adr_all@meta.data$new_clust[which(
  adr_all@meta.data$new_clust %in% c("22", "5", "28", "0", "1", "2", "10", "6", "8", "9", "21")
)] = "Cortex"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust %in% c("16", "26"))] =
  "Leukocytes"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust %in% c("30", "7", "4", "12"))] =
  "Erythroblasts"
adr_all@meta.data$new_clust[which(adr_all@meta.data$new_clust %in% c("27", "31", "29"))] =
  "Other"
adr_all@meta.data$new_clust = factor(
  adr_all@meta.data$new_clust,
  levels = c(
    "SCPs",
    "Bridge",
    "Sympathoblastic",
    "Chromaffin",
    "Endothelium",
    "Mesenchyme",
    "Cortex",
    "Leukocytes",
    "Erythroblasts",
    "Other"
  )
)
adr_all@meta.data$sample_name = "x"
adr_all@meta.data$sample_name[which(sapply(strsplit(names(
  adr_all@active.ident
), "_"), "[", 6) %in% c("Adr8710632", "Adr8710633"))] = "w21_1"
adr_all@meta.data$sample_name[which(sapply(strsplit(names(
  adr_all@active.ident
), "_"), "[", 6) %in% c("Adr8710634", "Adr8710635"))] = "w21_2"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident == "babyAdrenal1")] =
  "w8"
adr_all@meta.data$sample_name[which(adr_all@meta.data$orig.ident == "babyAdrenal2")] =
  "w8d6"
adr_all@meta.data$sample_name[which(
  adr_all@meta.data$orig.ident %in% c(
    "5388STDY7717452",
    "5388STDY7717453",
    "5388STDY7717454",
    "5388STDY7717455"
  )
)] = "w10d5_1"
adr_all@meta.data$sample_name[which(
  adr_all@meta.data$orig.ident %in% c("5388STDY7717456", "5388STDY7717458",
                                      "5388STDY7717459")
)] = "w10d5_2"
adr_all@meta.data$sample_name[which(
  adr_all@meta.data$orig.ident %in% c("5698STDY7839907", "5698STDY7839909",
                                      "5698STDY7839917")
)] = "w11"

#now subset the adrenal medulla data for downstream analysis
adr_all_med = subset(adr_all, idents = c(25, 24, 13, 19, 20))

pdf("panel S3D Behjati fetal adrenal medulla.pdf",
    width = 10,
    height = 10)
print(DimPlot(
  adr_all_med,
  cols = c("darkred", "orange", "darkblue", "#708090"),
  pt.size = 2,
  group.by = "new_clust"
))
dev.off()

Behjati.matrix <- adr_all_med@assays$RNA@data
Behjati.labels <- adr_all_med$new_clust

# pred <- SingleR(test = gcdata.matrix, ref = Behjati.matrix, labels = Behjati.labels, de.method = "wilcox")
# saveRDS(pred, "Behjati prediction scores.Rds")
pred.Behjati <- readRDS("Behjati prediction scores.Rds")

#We also ran SingleR on a cluster-by-cluster basis for panel K
# pred2 <- SingleR(test = gcdata.matrix, ref = Westermann.matrix, labels = Westermann.labels, de.method = "wilcox",
#                  clusters = gcdata$liger_clusters)

gcdata <-
  AddMetaData(gcdata,
              metadata = pred.Behjati$pruned.labels,
              col.name = "SingleR.mal.Behjati")

pdf("panel S3E SingleR predictions using Kildisiute data.pdf")
print(DimPlot(
  gcdata,
  group.by = "SingleR.mal.Behjati",
  pt.size = 1,
  cols = c("orange", "#708090", "darkred", "darkblue"),
  raster.dpi = c(1024, 1024)
))
dev.off()

###################### PANEL G ########################
all.data <-
  readRDS("../20220117 final HTAPP neuroblastoma LIGER.Rds")
cell.types <- data.frame(
  sample = all.data$orig.ident,
  fine = all.data$annotate_refine_coarse,
  coarse = all.data$annotate_coarse,
  malignant = all.data$Malignant
)
#cell.types$new.ids <- all.data$annotate_refine_coarse
load("../2020_08_18_all_Malignant_batchcorrect_liger_k12.Rda")
liger.ids <- gcdata$liger_clusters

cell.types$new.ids <- as.character(cell.types$fine)
cell.types[names(liger.ids), ]$new.ids <- as.character(liger.ids)

write.csv(table(cell.types$sample, cell.types$new.ids), file = "panel S3G raw cell counts using iNMF.csv")
write.csv(prop.table(table(cell.types$sample, cell.types$new.ids), margin = 1), file = "panel S3G raw cell props using iNMF.csv")

mal.cell.types <-
  subset(cell.types, subset = malignant == "Malignant")
write.csv(table(mal.cell.types$sample, mal.cell.types$new.ids), file = "panel S3H raw malignant cell counts using iNMF.csv")
write.csv(prop.table(
  table(mal.cell.types$sample, mal.cell.types$new.ids),
  margin = 1
), file = "panel S3H malignant cell props using iNMF.csv")

#write.csv(prop.table(table(all.data$orig.ident, all.data$annotate_refine_coarse), margin = 1), file = "panel s3f cell type proportions.csv")

#this data was input into GraphPad Prism to generate bar plots

###################### PANEL H ########################
write.csv(table(gcdata$orig.ident, gcdata$liger_clusters), file = "panel s3h ")



####ADDITIONAL ANALYSIS####
Idents(gcdata) <- "annot.ids"
x <-
  FindMarkers(
    gcdata,
    ident.1 = "ADRN",
    ident.2 = "sympathoblast",
    logfc.threshold = 1
  )
