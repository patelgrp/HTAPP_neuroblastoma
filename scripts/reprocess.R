library(Seurat)
library(rliger)
library(SeuratData)
library(patchwork)
library(numbat)
library(dplyr)
library(scCustomize)
library(stringi)
library(ggplot2)

#the purpose of this code is to import datasets from HTAPP, exclude HTAPP-124-SMP-61 (which unfortunately
#was a sample contaminated with a portion of adrenocortical adenoma. This has been validated by comparing
#pseudobulked data to reference bulk sequencing data.
source.dir <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/data/"
output.dir <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/"

#read in original Seurat object with the offending object
#NOTE - liger throws a weird error - it does not affect the analysis (rerun the command and it'll work!)
gcdata <-
  readRDS(paste0(source.dir, "20220117 final HTAPP neuroblastoma LIGER.Rds"))

#create a new Seurat object removing HTAPP-124-SMP-61
Seurat.new <- subset(gcdata, orig.ident != "HTAPP-124-SMP-61")

#remove old corrupted object to avoid mischief
remove(gcdata)

#break up the Seurat object by channel so we can integrate
Seurat.new[["RNA"]] <-
  split(Seurat.new[["RNA"]], f = Seurat.new$Channel)
sample.names <- unique(Seurat.new$Channel)

#import numbat data so that we can annotate SNV and CNV probability scores for each sample
numbat.dir <-
  "/media/ResearchHome/dyergrp/projects/ALSF_Pediatric_Atlas/common/REFERENCE_PUBLISHED/HTAPP/numbat/processed2/"

#create empty lists to dump data into
#numbat_probs <- list()
Seurat.obj <- list()
numbat.probs <- list()

cluster.annos <-
  read.csv(file = "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/scripts/cluster_annotations.csv",
           header = TRUE,
           row.names = 1)
confidence.annos <-
  read.csv(file = "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/scripts/confidence.csv",
           header = TRUE,
           row.names = 1)

#this loop goes channel-by-channel, generates plots and incorporates annotations from the
#cluster_annotations.csv file
#after running this list, there will be a list of Seurat objects ready for merging, re-integration and subsetting
for (i in 1:length(sample.names))
{
  Seurat.obj[[i]] <-
    subset(Seurat.new, subset = Channel == sample.names[i])
  
  #reprocess each channel individually
  Seurat.obj[[i]] <- NormalizeData(Seurat.obj[[i]])
  Seurat.obj[[i]] <- FindVariableFeatures(Seurat.obj[[i]])
  Seurat.obj[[i]] <- ScaleData(Seurat.obj[[i]])
  Seurat.obj[[i]] <- CellCycleScoring(Seurat.obj[[i]],
                                      s.features = cc.genes.updated.2019$s.genes,
                                      g2m.features = cc.genes.updated.2019$g2m.genes)
  Seurat.obj[[i]] <- RunPCA(Seurat.obj[[i]])
  Seurat.obj[[i]] <- FindNeighbors(Seurat.obj[[i]], dims = 1:30)
  Seurat.obj[[i]] <-
    FindClusters(Seurat.obj[[i]],
                 resolution = 0.4,
                 cluster.name = "sample_cluster")
  Seurat.obj[[i]] <- RunUMAP(Seurat.obj[[i]], dims = 1:30)
  
  sample.dirs <- paste0(numbat.dir, sample.names[i])
  
  #apply manual annotations based off singleR + numbat probability scores
  Idents(Seurat.obj[[i]]) <- "sample_cluster"
  annos <- t(cluster.annos[i, ])
  annos <- stri_remove_empty(annos)
  try(Seurat.obj[[i]] <-
        Rename_Clusters(Seurat.obj[[i]], 
                        new_idents = annos,
                        new_ident_name = "annotated_coarse", overwrite = TRUE))
  try(Seurat.obj[[i]]$annotated_coarse <-
        Seurat.obj[[i]]@active.ident)
  Seurat.obj[[i]]$malignant_calling <-
    confidence.annos$confidence[i]
  
  if (file.exists(paste0(sample.dirs, "/joint_post_2.tsv")))
  {
    #incorporate numbat output
    numbat_output <- Numbat$new(out_dir = sample.dirs)
    numbat.probs[[i]] <-
      data.frame(
        cell = numbat_output$clone_post$cell,
        joint = numbat_output$clone_post$p_cnv,
        expression = numbat_output$clone_post$p_cnv_x,
        allele = numbat_output$clone_post$p_cnv_y
      )
    rownames(numbat.probs[[i]]) <- numbat.probs[[i]]$cell
    numbat.probs[[i]] <- numbat.probs[[i]][Cells(Seurat.obj[[i]]),]
    Seurat.obj[[i]] <-
      AddMetaData(Seurat.obj[[i]], metadata = numbat.probs[[i]])
    
    pdf(file = paste0(output.dir, sample.names[i], "_Seurat_analysis.pdf"))
    print(DimPlot(Seurat.obj[[i]], group.by = "sample_cluster"))
    print(DimPlot(Seurat.obj[[i]], group.by = "sample_cluster",  label = T) + NoLegend())
    print(DimPlot(Seurat.obj[[i]], group.by = "Phase"))
    print(DimPlot(Seurat.obj[[i]], group.by = "Phase", label = T) + NoLegend())
    print(DimPlot(Seurat.obj[[i]], group.by = "SingleR.cluster.labels"))
    print(DimPlot(Seurat.obj[[i]], group.by = "SingleR.cluster.labels", label = T) + NoLegend())
    print(numbat_output$plot_phylo_heatmap())
    print(numbat_output$plot_consensus())
    print(
      numbat_output$plot_sc_tree(
        label_size = 3,
        branch_width = 0.5,
        tip_length = 0.5,
        tip = TRUE
      )
    )
    print(numbat_output$plot_mut_history())
    print(
      FeaturePlot_scCustom(
        seurat_object = Seurat.obj[[i]],
        features = c("allele", "expression", "joint", "cnv"),
        num_columns = 2,
        colors_use = viridis_light_high
      )
    )
    try(print(DimPlot(Seurat.obj[[i]], group.by = "annotated_coarse")))
    try(print(DimPlot(Seurat.obj[[i]], group.by = "annotated_coarse", label = T) + NoLegend()))
    print(
      VlnPlot_scCustom(
        seurat_object = Seurat.obj[[i]],
        features = c("allele", "expression", "joint", "cnv"),
        num_columns = 2,
        colors_use = viridis_light_high,
        plot_median = TRUE,
        group.by = "sample_cluster",
        pt.size = 0
      )
    )
    print(
      VlnPlot_scCustom(
        seurat_object = Seurat.obj[[i]],
        features = c("allele", "expression", "joint", "cnv"),
        num_columns = 2,
        colors_use = viridis_light_high,
        plot_median = TRUE,
        pt.size = 0
      )
    )
    dev.off()
    
    remove(numbat_output)
  }
  else
  {
    numbat.probs[[i]] <-
      data.frame(
        cells = Cells(Seurat.obj[[i]]),
        joint = NA,
        expression = NA,
        allele = NA
      )
    Seurat.obj[[i]] <-
      AddMetaData(Seurat.obj[[i]], metadata = numbat.probs[[i]])
    
    pdf(file = paste0(output.dir, sample.names[i], "_Seurat_analysis.pdf"))
    print(DimPlot(Seurat.obj[[i]], group.by = "sample_cluster"))
    print(DimPlot(Seurat.obj[[i]], group.by = "sample_cluster",  label = T) + NoLegend())
    print(DimPlot(Seurat.obj[[i]], group.by = "Phase"))
    print(DimPlot(Seurat.obj[[i]], group.by = "Phase", label = T) + NoLegend())
    print(DimPlot(Seurat.obj[[i]], group.by = "SingleR.cluster.labels"))
    print(DimPlot(Seurat.obj[[i]], group.by = "SingleR.cluster.labels", label = T) + NoLegend())
    print(
      FeaturePlot_scCustom(
        seurat_object = Seurat.obj[[i]],
        features = c("cnv"),
        colors_use = viridis_light_high
      )
    )
    try(print(DimPlot(Seurat.obj[[i]], group.by = "annotated_coarse")))
    try(print(DimPlot(Seurat.obj[[i]], group.by = "annotated_coarse", label = T) + NoLegend()))
    print(
      VlnPlot_scCustom(
        seurat_object = Seurat.obj[[i]],
        features = c("cnv"),
        colors_use = viridis_light_high,
        plot_median = TRUE,
        group.by = "sample_cluster",
        pt.size = 0
      )
    )
    print(
      VlnPlot_scCustom(
        seurat_object = Seurat.obj[[i]],
        features = c("cnv"),
        num_columns = 2,
        colors_use = viridis_light_high,
        plot_median = TRUE,
        pt.size = 0
      )
    )
    dev.off()
  }
  #save each object so I don't ever have to do this again
  saveRDS(Seurat.obj[[i]],
          file = paste0(output.dir, "/", sample.names[i], "_Seurat.Rds"))
  remove(annos)
}

remove(Seurat.new)
Seurat.merge <-
  merge(Seurat.obj[[1]], Seurat.obj[2:length(sample.names)])
integrated.Seurat <- Seurat.merge

#####INTEGRATION OF ALL DATA WITH LIGER
integrated.Seurat <- normalize(integrated.Seurat)
integrated.Seurat <- selectGenes(integrated.Seurat)
integrated.Seurat <- scaleNotCenter(integrated.Seurat)
integrated.Seurat

integrated.Seurat <- runINMF(integrated.Seurat, k = 20)
integrated.Seurat <- quantileNorm(integrated.Seurat)
integrated.Seurat

integrated.Seurat <-
  RunUMAP(integrated.Seurat, reduction = "inmfNorm", dims = 1:20)
gg.byDataset <-
  DimPlot(integrated.Seurat, group.by = "orig.ident", label = T) + NoLegend()
gg.byCluster <-
  DimPlot(integrated.Seurat, group.by = "inmfNorm.cluster", label = T) + NoLegend()
gg.bySingleR <-
  DimPlot(integrated.Seurat, group.by = "SingleR.cluster.labels", label = T) + NoLegend()
gg.byAnnot <-
  DimPlot(integrated.Seurat, group.by = "annotated_coarse", label = T) + NoLegend()

# Seurat.new <-
#   FindNeighbors(Seurat.new, reduction = "inmfNorm", dims = 1:20)
# Seurat.new <- FindClusters(Seurat.new, resolution = 0.4)

saveRDS(integrated.Seurat,
        file = paste0(output.dir, "combined_dataset_k20.Rds"))

# integrated.Seurat <-
#   readRDS(file = paste0(output.dir, "combined_dataset_k20.Rds"))

Idents(integrated.Seurat) <- "Channel"
pdf(
  file = paste0(output.dir, "integrated data split by Seurat_analysis.pdf"),
  width = 20,
  height = 10
)
print(DimPlot(
  integrated.Seurat,
  group.by = c("inmfNorm.cluster", "annotated_coarse"),
  label = T
) + ggtitle("All samples plot"))

for (i in 1:length(sample.names))
{
  print(
    DimPlot(
      integrated.Seurat,
      group.by = c("inmfNorm.cluster", "annotated_coarse"),
      label = T,
      cells = WhichCells(integrated.Seurat, idents = sample.names[i])
    ) + ggtitle(paste0(sample.names[i], " plot"))
  )
  print(
    DimPlot(
      integrated.Seurat,
      group.by = c("seurat_clusters", "annotated_coarse"),
      label = T,
      cells = WhichCells(integrated.Seurat, idents = sample.names[i])
    ) + ggtitle(paste0(sample.names[i], " plot"))
  )
}
dev.off()
#####INTEGRATION OF EACH SUBSET WITH LIGER
Seurat.mal <-
  subset(integrated.Seurat,
         subset = malignant_calling == "SNV+CNV" &
           annotated_coarse == "malignant")
Seurat.mal <- normalize(Seurat.mal)
Seurat.mal <- selectGenes(Seurat.mal)
Seurat.mal <- scaleNotCenter(Seurat.mal)
Seurat.mal

Seurat.mal <- runINMF(Seurat.mal, k = 12)
Seurat.mal <- quantileNorm(Seurat.mal)
Seurat.mal

Seurat.mal <-
  RunUMAP(Seurat.mal, reduction = "inmfNorm", dims = 1:12)

saveRDS(Seurat.mal,
        file = paste0(output.dir, "malignant_combined_dataset_k12.Rds"))

# Seurat.mal <-
#   readRDS(file = paste0(output.dir, "malignant_combined_dataset_k20.Rds"))