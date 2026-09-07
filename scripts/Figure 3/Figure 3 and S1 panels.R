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
library(scales)
library(readxl)

options("SCpubr.ColorPaletteEnds" = FALSE) #fixes the issue where diverging palettes end up white on the extremes

output.dir <- output.dir <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/Fig3"

#load complete dataset with all single-cell and single-nucleus data
gcdata <-
  readRDS(
    "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/combined_dataset_k20.Rds"
  )

gcdata.all <- gcdata #keep combined dataset object

#subset all immune cells
gcdata <-
  subset(gcdata, subset = annotated_coarse == "immune")

#gcdata <- readRDS(paste0(output.dir, "immune_combined_dataset_k12_process.Rds"))

#reprocesss data
gcdata <- normalize(gcdata)
gcdata <- selectGenes(gcdata)
gcdata <- scaleNotCenter(gcdata)
gcdata

num_nmfs <- 12
gcdata <- runINMF(gcdata, k = num_nmfs)
gcdata <- quantileNorm(gcdata)
gcdata

gcdata <-
  RunUMAP(gcdata, reduction = "inmfNorm", dims = 1:num_nmfs)

saveRDS(
  gcdata,
  file = paste0(
    output.dir,
    "/immune_combined_dataset_k",
    num_nmfs,
    "_reprocess.Rds"
  )
)

######################CODE TO GENERATE FIGURE 3 PANELS######################
pdf(paste0(output.dir, "/panel 3A 20220714 HTAPP immune liger.pdf"))
print(
  DimPlot_scCustom(
    seurat_object = gcdata,
    figure_plot = TRUE,
    group.by = "inmfNorm.cluster",
    colors_use = met.brewer("Signac", num_nmfs),
  )
)
dev.off()

#UMAP plot colored by sample
color.scheme <- MetBrewer::met.brewer("Hokusai1", n = 56)
color.scheme <- sample(color.scheme)
#note sample 11 and 12 are duplicates
color.scheme[11] <- color.scheme[12]


pdf(paste0(
  output.dir,
  "/20220714 HTAPP immune liger with labels.pdf"
))
print(
  DimPlot_scCustom(
    seurat_object = gcdata,
    figure_plot = TRUE,
    group.by = "inmfNorm.cluster",
    colors_use = met.brewer("Signac", num_nmfs),
    label = TRUE
  )
)
dev.off()

#pull top 50 feature embeddings for each NMF cluster for supplemental tables
nmf.embeddings <- matrix(nrow = 50, ncol = num_nmfs)
colnames(nmf.embeddings) <- paste0("iNMF", 1:num_nmfs)

for (i in 1:num_nmfs)
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
write.csv(de.markers.filter, file = paste0(output.dir, "/Immune DE genes.csv"))

# write.csv(
#   GetAssayData(gcdata.flat, slot = "counts"),
#   file = paste0(output.dir, "/immune_count_matrix.csv")
# )

#CODE FOR RUNNING CELLTYPIST IN COMMAND LINE
#celltypist --indata /mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/Fig4/immune_count_matrix.csv --model Immune_All_Low.pkl --transpose-input --majority-voting --outdir /mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/Fig4/ --prefix low --xlsx
#celltypist --indata /mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/Fig4/immune_count_matrix.csv --model Immune_All_High.pkl --transpose-input --majority-voting --outdir /mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/Fig4/ --prefix high`` --xlsx

low_annotations <-
  read_xlsx(path = paste0(output.dir, "/lowannotation_result.xlsx"),
            sheet = 1)
rownames(low_annotations) <- low_annotations$...1
low_annotations$...1 <- NULL

high_annotations <-
  read_xlsx(path = paste0(output.dir, "/highannotation_result.xlsx"),
            sheet = 1)
rownames(high_annotations) <- high_annotations$...1
high_annotations$...1 <- NULL

gcdata.flat$celltypist_low <- low_annotations$majority_voting
gcdata.flat$celltypist_high <- high_annotations$majority_voting

#generates panel 3A
pdf(paste0(output.dir, "/20250623 HTAPP immune celltypist.pdf"))
print(
  DimPlot_scCustom(
    seurat_object = gcdata.flat,
    figure_plot = TRUE,
    group.by = "celltypist_low",
    colors_use = met.brewer("Signac", length(unique(
      gcdata.flat$celltypist_low
    ))),
    label = TRUE
  )
)
print(
  DimPlot_scCustom(
    seurat_object = gcdata.flat,
    figure_plot = TRUE,
    group.by = "celltypist_high",
    colors_use = met.brewer("Signac", length(unique(
      gcdata.flat$celltypist_high
    ))),
    label = TRUE
  )
)

color.list <- c(
  "B cells" = "#EEC76E",
  "B-cell lineage" = "#D8C673",
  "DC" = "#6E9869",
  "Double-positive thymocytes" = "#5C67A8",
  "Endothelial cells" = "#CF8B4A",
  "Epithelial cells" = "#545454",
  "Erythroid" = "#000000",
  "Fibroblasts" = "#BEBDBD",
  "ILC" = "#464A74",
  "Macrophages" = "#ACC67E",
  "Mast cells" = "#A9B4D4",
  "Monocytes" = "#8AA69D",
  "pDC" = "#9AAEB9",
  "Plasma cells" = "#C2C678",
  "T cells" = "#0F31BC"
)
print(
  DimPlot_scCustom(
    seurat_object = gcdata.flat,
    figure_plot = TRUE,
    group.by = "celltypist_high",
    colors_use = color.list,
    label = TRUE
  )
)
print(
  DimPlot_scCustom(
    seurat_object = gcdata.flat,
    figure_plot = TRUE,
    group.by = "celltypist_high",
    colors_use = color.list,
    label = FALSE
  )
)
print(
  DimPlot_scCustom(
    seurat_object = gcdata.flat,
    figure_plot = TRUE,
    group.by = "celltypist_high",
    colors_use = color.list,
    label = FALSE
  ) + NoLegend()
)
dev.off()

pdf(paste0(output.dir, "/20250623 HTAPP immune umap by condition.pdf"))
print(
  DimPlot_scCustom(
    seurat_object = gcdata.flat,
    figure_plot = TRUE,
    group.by = "condition",
    colors_use = c("cell" = "#6FA7B9",
                   "nucleus" = "#D1B07C"),
    label = FALSE
  )
)
dev.off()

gcdata2 <- gcdata.flat #stow data in case of a mishap
gcdata.flat$update_annot <- gcdata.flat$annotated_coarse
gcdata.all$update_annot <- gcdata.all$annotated_coarse

immune.metadata <- gcdata.flat[[]]
all.metadata <- gcdata.all[[]]
################MYELOID########################
myeloid <-
  subset(
    gcdata.flat,
    subset = celltypist_high %in% c("DC", "Macrophages", "Mast cells", "Monocytes", "pDC") &
      condition == "nucleus"
  )
myeloid.counts <- table(myeloid$orig.ident)

#subset to myeloid nucleus data (39/44 datasets had > 100 myeloid nuclei)
myeloid.subset <- names(myeloid.counts[myeloid.counts >= 100])

myeloid <-
  subset(myeloid, subset = orig.ident %in% myeloid.subset &
           condition == "nucleus")

Idents(myeloid) <- "orig.ident"
myeloid <- normalize(myeloid)
myeloid <- selectGenes(myeloid)
myeloid <- scaleNotCenter(myeloid)
myeloid

num_nmfs <- 8
myeloid <- runINMF(myeloid, k = num_nmfs)
myeloid <- quantileNorm(myeloid)
myeloid

myeloid <-
  RunUMAP(myeloid, reduction = "inmfNorm", dims = 1:num_nmfs)

saveRDS(
  myeloid,
  file = paste0(
    output.dir,
    "/myeloid_combined_dataset_k",
    num_nmfs,
    "_reprocess.Rds"
  )
)

myeloid <-
  readRDS(file = paste0(
    output.dir,
    "/myeloid_combined_dataset_k",
    num_nmfs,
    "_reprocess.Rds"
  ))

color.scheme <- MetBrewer::met.brewer("Hokusai1", n = 39)
color.scheme <- sample(color.scheme)

#generates a bunch of plots for QC and cross-validation
pdf(paste0(output.dir, "/Myeloid nuclei UMAP.pdf"),
    width = 10,
    height = 10)
print(
  DimPlot_scCustom(
    myeloid,
    group.by = "inmfNorm.cluster",
    colors_use = c(
      "#62324E",
      "#A54D7A",
      "#E19F4A",
      "#292C6A",
      "#0CB7A6",
      "#90AC3D",
      "#A21D21",
      "#0A6946"
    ),
    shuffle = T,
    figure_plot = TRUE,
    label = F,
    pt.size = 1,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    myeloid,
    group.by = "celltypist_high",
    shuffle = T,
    label = T,
    figure_plot = TRUE,
    pt.size = 1,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    myeloid,
    group.by = "celltypist_low",
    shuffle = T,
    label = T,
    figure_plot = TRUE,
    pt.size = 1,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    myeloid,
    group.by = "Phase",
    shuffle = T,
    label = T,
    figure_plot = TRUE,
    pt.size = 1,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    myeloid,
    group.by = "orig.ident",
    colors_use = color.scheme,
    label = F,
    figure_plot = TRUE,
    pt.size = 0.5,
    shuffle = T,
    raster = F
  ) + NoLegend()
)
print(
  DimPlot_scCustom(
    myeloid,
    group.by = "orig.ident",
    colors_use = color.scheme,
    label = F,
    figure_plot = TRUE,
    pt.size = 0.5,
    shuffle = T,
    raster = F
  )
)
dev.off()

Idents(myeloid) <- "inmfNorm.cluster"
annot.ids.myeloid = c(
  "cDC1",
  "cDC2",
  "CD68_macrophage",
  "CD68_macrophage",
  "monocyte",
  "CD163_macrophage",
  "CD163_macrophage",
  "pDC"
)
names(annot.ids.myeloid) <- levels(myeloid$inmfNorm.cluster)

myeloid <-
  Rename_Clusters(object = myeloid, new_idents = annot.ids.myeloid)
myeloid$annot_ids <- myeloid@active.ident

custom.order <-
  c("cDC1",
    "cDC2",
    "pDC",
    "monocyte",
    "CD68_macrophage",
    "CD163_macrophage")

myeloid$annot_ids <-
  factor(myeloid$annot_ids, levels = custom.order)

Idents(myeloid) <- "annot_ids"
color.scheme <-
  colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(256)

#generates panel 3B
pdf(
  paste0(output.dir, "/Dim plot myeloid nucleus.pdf"),
  width = 20,
  height = 20
)
print(
  DimPlot_scCustom(
    myeloid,
    colors_use = c(
      "#A11D21",
      "#292C6A",
      "#8FAB3D",
      "#E09E4A",
      "#A44D7A",
      "#62324E"
    ),
    shuffle = T,
    label = F,
    pt.size = 2,
    figure_plot = TRUE,
    raster = F
  )
)
dev.off()

#generates panel 3C
pdf(
  paste0(output.dir, "/Dot plot myeloid nucleus.pdf"),
  width = 20,
  height = 20
)
print(
  do_DotPlot(
    myeloid,
    features = c(
      "MRC1",
      "CD163",
      "CD68",
      "FCER1G",
      "VCAN",
      "FCN1",
      "S100A8",
      "S100A9",
      "CLEC4C",
      "GZMB",
      "IL3RA",
      "CX3CR1",
      "CD1C",
      "FCER1A",
      "CLEC9A",
      "XCR1"
    ),
    use_viridis = F,
    flip = T,
    sequential.palette = "RdBu",
    sequential.direction = -1,
    scale = T,
    dot.scale = 25,
    dot_border = F
  ) + scale_color_gradient2(
    low = "darkblue",
    mid = "white",
    high = "darkred",
    midpoint = 0.5, limits = c(-1, 2), oob = squish
  ) + guides(color = guide_colorbar())
)
dev.off()

#generates panel S1F
pdf(
  paste0(output.dir, "/Bar plot myeloid nucleus.pdf"),
  width = 20,
  height = 20
)
print(
  do_BarPlot(
    sample = myeloid,
    group.by = "annot_ids",
    split.by = "orig.ident",
    position = "fill",
    colors.use = c(
      "cDC1" = "#A11D21",
      "cDC2" = "#292C6A",
      "pDC" = "#8FAB3D",
      "monocyte" = "#E09E4A",
      "CD68_macrophage" = "#A44D7A",
      "CD163_macrophage" = "#62324E"
    )
  )
)
dev.off()

###################### T CELLS ######################
t.cell <-
  subset(
    gcdata.flat,
    subset = celltypist_high %in% c("Double-positive thymocytes", "ILC", "T cells") &
      condition == "nucleus"
  )
t.counts <- table(t.cell$orig.ident)


#subset to t nucleus data (29/41 datasets had > 100 t nuclei)
t.subset <- names(t.counts[t.counts >= 100])

t.cell <-
  subset(t.cell, subset = orig.ident %in% t.subset &
           condition == "nucleus")

Idents(t.cell) <- "orig.ident"
t.cell <- normalize(t.cell)
t.cell <- selectGenes(t.cell)
t.cell <- scaleNotCenter(t.cell)
t.cell

num_nmfs <- 8
t.cell <- runINMF(t.cell, k = num_nmfs)
t.cell <- quantileNorm(t.cell)
t.cell

t.cell <-
  RunUMAP(t.cell, reduction = "inmfNorm", dims = 1:num_nmfs)

saveRDS(t.cell,
        file = paste0(
          output.dir,
          "/tcell_combined_dataset_k",
          num_nmfs,
          "_reprocess.Rds"
        ))

t.cell <-
  readRDS(file = paste0(
    output.dir,
    "/tcell_combined_dataset_k",
    num_nmfs,
    "_reprocess.Rds"
  ))

color.scheme <- MetBrewer::met.brewer("Hokusai1", n = 29)
color.scheme <- sample(color.scheme)

#generates a bunch of plots for QC and cross-validation
pdf(paste0(output.dir, "/T cell nuclei UMAP.pdf"),
    width = 10,
    height = 10)
print(
  DimPlot_scCustom(
    t.cell,
    group.by = "inmfNorm.cluster",
    colors_use = c(
      "#A11D21",
      "#27377A",
      "#FFCD12",
      "#067E40",
      "#B86192",
      "#03B6A5",
      "#BEBD8D",
      "#721C3F"
    ),
    shuffle = T,
    label = F,
    pt.size = 1,
    figure_plot = TRUE,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    t.cell,
    group.by = "celltypist_high",
    shuffle = T,
    label = T,
    pt.size = 1,
    figure_plot = TRUE,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    t.cell,
    group.by = "celltypist_low",
    shuffle = T,
    label = T,
    pt.size = 1,
    figure_plot = TRUE,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    t.cell,
    group.by = "Phase",
    shuffle = T,
    label = T,
    pt.size = 1,
    figure_plot = TRUE,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    t.cell,
    group.by = "orig.ident",
    colors_use = color.scheme,
    label = F,
    pt.size = 0.5,
    shuffle = T,
    figure_plot = TRUE,
    raster = F
  ) + NoLegend()
)
print(
  DimPlot_scCustom(
    t.cell,
    group.by = "orig.ident",
    colors_use = color.scheme,
    label = F,
    pt.size = 0.5,
    shuffle = T,
    figure_plot = TRUE,
    raster = F
  )
)
dev.off()

Idents(t.cell) <- "inmfNorm.cluster"
annot.ids.t = c("ILC",
                "CD8_Teff",
                "naive_T",
                "naive_T",
                "CD4_Treg",
                "CD4_Th",
                "NK",
                "CD8_Teff")
names(annot.ids.t) <- levels(t.cell$inmfNorm.cluster)

t.cell <- Rename_Clusters(object = t.cell, new_idents = annot.ids.t)
t.cell$annot_ids <- t.cell@active.ident

custom.order <-
  c("naive_T", "CD8_Teff", "CD4_Treg", "CD4_Th", "NK", "ILC")

t.cell$annot_ids <- factor(t.cell$annot_ids, levels = custom.order)

Idents(t.cell) <- "annot_ids"
color.scheme <-
  colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(256)

#generates panel 3D
pdf(paste0(output.dir, "/Dim plot T nucleus.pdf"),
    width = 20,
    height = 20)
print(
  DimPlot_scCustom(
    t.cell,
    colors_use = c(
      "#A11D21",
      "#27377A",
      "#067E40",
      "#FFCD12",
      "#B86192",
      "#721C3F"
    ),
    shuffle = T,
    label = F,
    pt.size = 2,
    figure_plot = TRUE,
    raster = F
  )
)
dev.off()

#generates panel 3E
pdf(paste0(output.dir, "/Dot plot T nucleus.pdf"),
    width = 20,
    height = 20)
print(
  do_DotPlot(
    t.cell,
    features = c(
      "LIF",
      "RORC",
      "IL23R",
      "PRF1",
      "GNLY",
      "NKG7",
      "RORA",
      "CD226",
      "CTLA4",
      "FOXP3",
      "CD4",
      "CD8A",
      "CCL5",
      "GZMA",
      "LEF1",
      "CCR7",
      "SELL"
    ),
    use_viridis = F,
    flip = T,
    sequential.palette = "RdBu",
    sequential.direction = -1,
    scale = T,
    dot.scale = 25,
    dot_border = F
  )+ scale_color_gradient2(
    low = "darkblue",
    mid = "white",
    high = "darkred",
    midpoint = 0, limits = c(-1, 1), oob = squish
  ) + guides(color = guide_colorbar())
)
dev.off()

#geneates panel S1G
pdf(paste0(output.dir, "/Bar plot T nucleus.pdf"),
    width = 20,
    height = 20)
print(
  do_BarPlot(
    sample = t.cell,
    group.by = "annot_ids",
    split.by = "orig.ident",
    position = "fill",
    colors.use = c(
      "naive_T" = "#A11D21",
      "CD8_Teff" = "#27377A",
      "CD4_Treg" = "#067E40",
      "CD4_Th" = "#FFCD12",
      "NK" = "#B86192",
      "ILC" = "#721C3F"
    )
  )
)
dev.off()

###################### B CELLS ######################
b.cell <-
  subset(
    gcdata.flat,
    subset = celltypist_high %in% c("B cells", "B-cell lineage", "Plasma cells") &
      condition == "nucleus"
  )
b.counts <- table(b.cell$orig.ident)


#subset to b nucleus data (7/37 datasets had > 100 b nuclei)
b.subset <- names(b.counts[b.counts >= 100])

b.cell <-
  subset(b.cell, subset = orig.ident %in% b.subset &
           condition == "nucleus")

Idents(b.cell) <- "orig.ident"
b.cell <- normalize(b.cell)
b.cell <- selectGenes(b.cell)
b.cell <- scaleNotCenter(b.cell)
b.cell

num_nmfs <- 8
b.cell <- runINMF(b.cell, k = num_nmfs)
b.cell <- quantileNorm(b.cell)
b.cell

b.cell <-
  RunUMAP(b.cell, reduction = "inmfNorm", dims = 1:num_nmfs)

saveRDS(b.cell,
        file = paste0(
          output.dir,
          "/bcell_combined_dataset_k",
          num_nmfs,
          "_reprocess.Rds"
        ))

b.cell <-
  readRDS(file = paste0(
    output.dir,
    "/bcell_combined_dataset_k",
    num_nmfs,
    "_reprocess.Rds"
  ))

color.scheme <- MetBrewer::met.brewer("Hokusai1", n = 7)
color.scheme <- sample(color.scheme)

#generates a bunch of plots for QC and cross-validation
pdf(paste0(output.dir, "/B cell nuclei UMAP.pdf"),
    width = 10,
    height = 10)
print(
  DimPlot_scCustom(
    b.cell,
    group.by = "inmfNorm.cluster",
    colors_use = c(
      "#A11D21",
      "#27377A",
      "#FFCD12",
      "#067E40",
      "#B86192",
      "#03B6A5",
      "#BEBD8D",
      "#721C3F"
    ),
    shuffle = T,
    label = F,
    pt.size = 1,
    figure_plot = TRUE,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    b.cell,
    group.by = "celltypist_high",
    shuffle = T,
    label = T,
    pt.size = 1,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    b.cell,
    group.by = "celltypist_low",
    shuffle = T,
    label = T,
    pt.size = 1,
    figure_plot = TRUE,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    b.cell,
    group.by = "Phase",
    shuffle = T,
    label = T,
    pt.size = 1,
    figure_plot = TRUE,
    raster = F
  )
)
print(
  DimPlot_scCustom(
    b.cell,
    group.by = "orig.ident",
    cols = color.scheme,
    label = F,
    pt.size = 0.5,
    figure_plot = TRUE,
    shuffle = T,
    raster = F
  ) + NoLegend()
)
print(
  DimPlot_scCustom(
    b.cell,
    group.by = "orig.ident",
    cols = color.scheme,
    label = F,
    pt.size = 0.5,
    figure_plot = TRUE,
    shuffle = T,
    raster = F
  )
)
dev.off()

Idents(b.cell) <- "inmfNorm.cluster"
annot.ids.b = c(
  "naive_B",
  "memory_B",
  "plasma",
  "memory_B",
  "memory_B",
  "naive_B",
  "memory_B",
  "naive_B"
)
names(annot.ids.b) <- levels(b.cell$inmfNorm.cluster)

b.cell <- Rename_Clusters(object = b.cell, new_idents = annot.ids.b)
b.cell$annot_ids <- b.cell@active.ident

custom.order <- c("naive_B", "memory_B", "plasma")

b.cell$annot_ids <- factor(b.cell$annot_ids, levels = custom.order)

Idents(b.cell) <- "annot_ids"
color.scheme <-
  colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(256)

#generates panel 3F
pdf(paste0(output.dir, "/Dim plot B nucleus.pdf"),
    width = 20,
    height = 20)
print(
  DimPlot_scCustom(
    b.cell,
    colors_use = c("#919059",
                   "#67332F",
                   "#184E49"),
    shuffle = T,
    label = F,
    pt.size = 4,
    figure_plot = TRUE,
    raster = F
  )
)
dev.off()

#generates panel 3G
pdf(paste0(output.dir, "/Dot plot B nucleus.pdf"),
    width = 20,
    height = 20)
print(
  do_DotPlot(
    b.cell,
    features = c(
      "CD38",
      "SDC1",
      "CD27",
      "CD24",
      "CD69",
      "CD79B",
      "IGHD",
      "MS4A1",
      "CD19"
    ),
    use_viridis = F,
    flip = T,
    sequential.palette = "RdBu",
    sequential.direction = -1,
    scale = T,
    dot.scale = 25,
    dot_border = F
  ) + scale_color_gradient2(
    low = "darkblue",
    mid = "white",
    high = "darkred",
    midpoint = 0, limits = c(-1, 1), oob = squish
  ) + guides(color = guide_colorbar())
)
dev.off()

#generates panel S1H
pdf(paste0(output.dir, "/Bar plot B nucleus.pdf"),
    width = 20,
    height = 20)
print(
  do_BarPlot(
    sample = b.cell,
    group.by = "annot_ids",
    split.by = "orig.ident",
    position = "fill",
    colors.use = c(
      "naive_B" = "#919059",
      "memory_B" = "#67332F",
      "plasma" = "#184E49"
    )
  )
)
dev.off()

immune.metadata[immune.metadata$celltypist_high %in% c("DC", "Macrophages", "Mast cells", "Monocytes", "pDC"), ]$update_annot <-
  "myeloid"

immune.metadata[immune.metadata$celltypist_high %in% c("Double-positive thymocytes", "ILC", "T cells"), ]$update_annot <-
  "t cell"

immune.metadata[immune.metadata$celltypist_high %in% c("B cells", "B-cell lineage", "Plasma cells"), ]$update_annot <-
  "b cell"

immune.metadata[immune.metadata$update_annot == "immune", ]$update_annot <- NA


is_immune <- all.metadata$annotated_coarse == "immune"
donor <- immune.metadata$update_annot[match(rownames(all.metadata), rownames(immune.metadata))]

fallback <- all.metadata$annotated_coarse

all.metadata$update_annot <- ifelse(is_immune & !is.na(donor), donor, fallback)
all.metadata[all.metadata$update_annot == "immune", ]$update_annot <- NA

gcdata.all$update_annot <- all.metadata$update_annot
gcdata.flat$update_annot <- immune.metadata$update_annot

##########DATASET HANDLING#############
saveRDS(gcdata.all, paste0(output.dir, "/final annotated HTAPP object.Rds"))
saveRDS(gcdata.flat, paste0(output.dir, "/final annotated immune object.Rds"))
saveRDS(myeloid, paste0(output.dir, "/final myeloid object.Rds"))
saveRDS(t.cell, paste0(output.dir, "/final t cell object.Rds"))
saveRDS(b.cell, paste0(output.dir, "/final b cell object.Rds"))

