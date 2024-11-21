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
  "/media/ResearchHome/dyergrp/projects/ALSF_Pediatric_Atlas/common/REFERENCE_PUBLISHED/HTAPP_core/scrna_level4/"
output.dir <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output/"

#read in original Seurat object with the offending object
load(
  paste0(
    source.dir,
    "2020_08_18_all_cluster_final_all_cluster_final_batchcorrect_liger_k20.Rda"
  )
)

#create a new Seurat object removing HTAPP-124-SMP-61
Seurat.new <- subset(gcdata, orig.ident != "HTAPP-124-SMP-61")

#remove old corrupted object to avoid mischief
remove(gcdata)

#break up the Seurat object by channel so we can integrate
Seurat.new[["RNA"]] <-
  split(Seurat.new[["RNA"]], f = Seurat.new$Channel)

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

numbat.dir <-
  "/media/ResearchHome/dyergrp/projects/ALSF_Pediatric_Atlas/common/REFERENCE_PUBLISHED/HTAPP_core/numbat/processed"
sample.dirs <- list.dirs(path = numbat.dir, recursive = F)
sample.names <-
  list.dirs(path = numbat.dir,
            recursive = F,
            full.names = F)

numbat_probs <- list()

for (i in 1:length(sample.dirs))
{
  if (file.exists(paste0(sample.dirs[i], "/joint_post_2.tsv")))
  {
    Seurat.temp <-
      subset(Seurat.new, subset = Channel == sample.names[i])
    Seurat.temp <- NormalizeData(Seurat.temp)
    Seurat.temp <- FindVariableFeatures(Seurat.temp)
    Seurat.temp <- ScaleData(Seurat.temp)
    Seurat.temp <- CellCycleScoring(
      Seurat.temp,
      s.features = cc.genes.updated.2019$s.genes,
      g2m.features = cc.genes.updated.2019$g2m.genes
    )
    Seurat.temp <- RunPCA(Seurat.temp)
    Seurat.temp <- FindNeighbors(Seurat.temp, dims = 1:30)
    Seurat.temp <-
      FindClusters(Seurat.temp,
                   resolution = 0.4,
                   cluster.name = "sample_cluster")
    Seurat.temp <- RunUMAP(Seurat.temp, dims = 1:30)
    
    numbat_output <- Numbat$new(out_dir = sample.dirs[i])
    numbat.probs[[i]] <-
      data.frame(
        cell = numbat_output$clone_post$cell,
        joint = numbat_output$clone_post$p_cnv,
        expression = numbat_output$clone_post$p_cnv_x,
        allele = numbat_output$clone_post$p_cnv_y
      )
    rownames(numbat.probs[[i]]) <- numbat.probs[[i]]$cell
    numbat.probs[[i]] <- numbat.probs[[i]][Cells(Seurat.temp), ]
    Seurat.temp <-
      AddMetaData(Seurat.temp, metadata = numbat.probs[[i]])
    
    pdf(file = paste0(output.dir, sample.names[i], "_Seurat_analysis.pdf"))
    print(numbat_output$plot_phylo_heatmap())
    print(numbat_output$plot_consensus())
    print(DimPlot(Seurat.temp))
    print(DimPlot(Seurat.temp, label = T) + NoLegend())
    print(DimPlot(Seurat.temp, group.by = "SingleR.cluster.labels"))
    print(DimPlot(Seurat.temp, group.by = "SingleR.cluster.labels", label = T) + NoLegend())
    print(
      FeaturePlot_scCustom(
        seurat_object = Seurat.temp,
        features = c("allele", "expression", "joint"),
        num_columns = 2,
        colors_use = viridis_light_high
      )
    )
    print(
      numbat_output$plot_sc_tree(
        label_size = 3,
        branch_width = 0.5,
        tip_length = 0.5,
        tip = TRUE
      )
    )
    print(numbat_output$plot_mut_history())
    dev.off()
    saveRDS(Seurat.temp,
            file = paste0(output.dir, "/", sample.names[i], "_Seurat.Rds"))
    remove(numbat_output, Seurat.temp)
  }
  else
  {
    Seurat.temp <-
      subset(Seurat.new, subset = Channel == sample.names[i])
    Seurat.temp <- NormalizeData(Seurat.temp)
    Seurat.temp <- FindVariableFeatures(Seurat.temp)
    Seurat.temp <- ScaleData(Seurat.temp)
    Seurat.temp <- CellCycleScoring(
      Seurat.temp,
      s.features = cc.genes.updated.2019$s.genes,
      g2m.features = cc.genes.updated.2019$g2m.genes
    )
    Seurat.temp <- RunPCA(Seurat.temp)
    Seurat.temp <- FindNeighbors(Seurat.temp, dims = 1:30)
    Seurat.temp <-
      FindClusters(Seurat.temp,
                   resolution = 0.4,
                   cluster.name = "sample_cluster")
    Seurat.temp <- RunUMAP(Seurat.temp, dims = 1:30)
    
    numbat.probs[[i]] <-
      data.frame(
        cell = Cells(Seurat.temp),
        joint = NA,
        expression = NA,
        allele = NA
      )
    
    pdf(file = paste0(output.dir, sample.names[i], "_Seurat_analysis.pdf"))
    print(DimPlot(Seurat.temp))
    print(DimPlot(Seurat.temp, group.by = "SingleR.cluster.labels"))
    print(DimPlot(Seurat.temp, group.by = "SingleR.cluster.labels", label = T) + NoLegend())
    dev.off()
    
    saveRDS(Seurat.temp,
            file = paste0(output.dir, "/", sample.names[i], "_Seurat.Rds"))
    remove(numbat_output, Seurat.temp)
  }
}

#annotate each sample, subset maligannt cells and generate integrated malignant, immune and stromal objects
cluster.annos <-
  read.csv(file = "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/cluster_annotations.csv",
           header = TRUE,
           row.names = 1)
output.dir2 <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/output_round2/"
Seurat.obj <- list()

for (i in 1:length(sample.names))
{
  Seurat.obj[[i]] <-
    readRDS(paste0(output.dir, "/", sample.names[i], "_Seurat.Rds"))
  Idents(Seurat.obj[[i]]) <- "sample_cluster"
  
  annos <- t(cluster.annos[i,])
  annos <- stri_remove_empty(annos)
  
  Seurat.obj[[i]] <-
    Rename_Clusters(Seurat.obj[[i]],
                    new_idents = annos,
                    meta_col_name = "annotated_coarse")
  Seurat.obj[[i]]$annotated_coarse <- Seurat.obj[[i]]@active.ident
  
  pdf(file = paste0(output.dir2, sample.names[i], "_Seurat_analysis.pdf"))
  print(DimPlot(Seurat.obj[[i]], group.by = "sample_cluster"))
  print(DimPlot(Seurat.obj[[i]], group.by = "sample_cluster",  label = T) + NoLegend())
  print(DimPlot(Seurat.obj[[i]], group.by = "SingleR.cluster.labels"))
  print(DimPlot(Seurat.obj[[i]], group.by = "SingleR.cluster.labels", label = T) + NoLegend())
  print(
    FeaturePlot_scCustom(
      seurat_object = Seurat.obj[[i]],
      features = c("allele", "expression", "joint"),
      num_columns = 2,
      colors_use = viridis_light_high
    )
  )
  print(DimPlot(Seurat.obj[[i]], group.by = "annotated_coarse"))
  print(DimPlot(Seurat.obj[[i]], group.by = "annotated_coarse", label = T) + NoLegend())
  print(
    VlnPlot_scCustom(
      seurat_object = Seurat.obj[[i]],
      features = c("allele", "expression", "joint", "cnv"),
      num_columns = 2,
      colors_use = viridis_light_high, 
      plot_median = TRUE, pt.size = 0
    )
  )
  dev.off()
  
  remove(annos)
}