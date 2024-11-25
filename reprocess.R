library(Seurat)
library(rliger)
library(SeuratData)
library(patchwork)
library(numbat)
library(dplyr)
library(scCustomize)
library(stringi)

#the purpose of this code is to import datasets from HTAPP, exclude HTAPP-124-SMP-61 (which unfortunately)
#mapped to adrenocortical adenoma instead of neuroblastoma when comparing
#pseudobulked data to reference bulk sequencing data
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
  "/mnt/storage1/Anand_temp/HTAPP/numbat/processed2/"

#create empty lists to dump data into
#numbat_probs <- list()
Seurat.obj <- list()
numbat.probs <- list()

cluster.annos <-
  read.csv(file = "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/cluster_annotations.csv",
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
                        meta_col_name = "annotated_coarse"))
  try(Seurat.obj[[i]]$annotated_coarse <-
        Seurat.obj[[i]]@active.ident)
  
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
  merge(Seurat.obj[[1]], Seurat.ob[2:length(sample.names)])

#####INTEGRATION OF ALL DATA WITH LIGER
Seurat.new <- normalize(Seurat.new)
Seurat.new <- selectGenes(Seurat.new)
Seurat.new <- scaleNotCenter(Seurat.new)
Seurat.new

Seurat.new <- runINMF(Seurat.new, k = 20)
Seurat.new <- quantileNorm(Seurat.new)
Seurat.new

Seurat.new <-
  RunUMAP(Seurat.new, reduction = "inmfNorm", dims = 1:20)
gg.byDataset <-
  DimPlot(Seurat.new, group.by = "orig.ident", label = T) + NoLegend()
gg.byCluster <-
  DimPlot(Seurat.new, group.by = "inmfNorm.cluster", label = T) + NoLegend()
gg.bySingleR <-
  DimPlot(Seurat.new, group.by = "SingleR.cluster.labels", label = T) + NoLegend()

Seurat.new <-
  FindNeighbors(Seurat.new, reduction = "inmfNorm", dims = 1:20)
Seurat.new <- FindClusters(Seurat.new, resolution = 0.4)

saveRDS(Seurat.new, file = paste0(output.dir, "combined_dataset_k20.Rds"))


#####INTEGRATION OF EACH SUBSET WITH LIGER