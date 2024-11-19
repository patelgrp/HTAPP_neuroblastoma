library(Seurat)
library(rliger)
library(SeuratData)
library(patchwork)

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
Seurat.new[["RNA"]] <- split(Seurat.new[["RNA"]], f = Seurat.new$Channel)

Seurat.new <- normalize(Seurat.new)
Seurat.new <- selectGenes(Seurat.new)
Seurat.new <- scaleNotCenter(Seurat.new)
Seurat.new

Seurat.new <- runINMF(Seurat.new, k=12)
Seurat.new <- quantileNorm(Seurat.new)
Seurat.new

Seurat.new <- RunUMAP(Seurat.new, reduction = "inmfNorm", dims = 1:20)
gg.byDataset <- DimPlot(Seurat.new, group.by = "orig.ident", label = T) + NoLegend()

Seurat.new <- FindNeighbors(Seurat.new, reduction = "inmfNorm", dims = 1:20)
Seurat.new <- FindClusters(Seurat.new)