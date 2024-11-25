library(numbat)
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)

data.dir <-
  "/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/data/"
numbat.dir <-
  "/media/ResearchHome/dyergrp/projects/ALSF_Pediatric_Atlas/common/REFERENCE_PUBLISHED/HTAPP_core/numbat/"
#load(file = "../2020_08_18_all_cluster_final_all_cluster_final_batchcorrect_liger_k20.Rda")
gcdata <-
  readRDS(paste0(data.dir, "20220117 final HTAPP neuroblastoma LIGER.Rds"))

#remove HTAPP-124 which was an adrenal adenoma
gcdata <- subset(gcdata, orig.ident != "HTAPP-124-SMP-61")

sample.ids <- unique(gcdata$Channel)

#import fetal adrenal medulla data from Jansky et al Nat Genetics (2023) and prep it
#as a reference dataset for numbat
ref_Jansky <-
  readRDS(paste0(data.dir, "adrenal_medulla_Seurat.RDS"))
ref_Jansky <- UpdateSeuratObject(ref_Jansky)
Jansky_annots <-
  data.frame(cell = Cells(ref_Jansky),
             group = as.character(ref_Jansky@active.ident))
Jansky_counts <- ref_Jansky@assays$RNA@counts
ref_Jansky_pseudobulk <-
  aggregate_counts(count_mat = Jansky_counts, annot = Jansky_annots)

all_sample_props <- data.frame()

for (i in 3:length(sample.ids))
{
  print(sample.ids[i])
  Seurat.obj <- subset(gcdata, subset = Channel == sample.ids[i])
  counts <- Seurat.obj@assays$RNA@counts
  # annots <-
  #   data.frame(cell = Cells(Seurat.obj),
  #              annot = Seurat.obj$annotate_refine_coarse)
  
  df_allele <- fread(paste0(
    numbat.dir,
    sample.ids[i],
    "/",
    sample.ids[i],
    "_allele_counts.tsv.gz"
  ))
  df_allele$cell <- gsub(pattern = "-1",
                         replacement = "",
                         x = df_allele$cell)
  df_allele$cell <- paste(sample.ids[i], df_allele$cell, sep = "-")
  
  try(out = run_numbat(
    count_mat = counts,
    #lambdas_ref = ref_hca,
    lambdas_ref = ref_Jansky_pseudobulk,
    df_allele = df_allele,
    genome = "hg38",
    t = 1e-5,
    ncores = 48,
    min_LLR = 5,
    call_clonal_loh = T,
    plot = T,
    out_dir = paste0(
      "/mnt/storage1/Anand_temp/HTAPP/numbat/processed2/",
      sample.ids[i],
      "/"
    )
  ))
  
  
  try(nb <-
        Numbat$new(out_dir = paste0(
          "/mnt/storage1/Anand_temp/HTAPP/numbat/processed2/",
          sample.ids[i],
          "/"
        )))
  
  try(write.csv(
    nb$joint_post %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y),
    file = paste0(
      "/mnt/storage1/Anand_temp/HTAPP/numbat/processed2/",
      sample.ids[i],
      "/",
      sample.ids[i],
      " cnv props individual chromosomes.csv"
    )
  ))
  
  try(joint_props <- nb$clone_post %>%
        {
          setNames(.$p_cnv, .$cell)
        })
  try(Seurat.obj$joint_props <- joint_props)
  try(write.csv(
    joint_props,
    file = paste0(
      "/mnt/storage1/Anand_temp/HTAPP/numbat/processed2/",
      sample.ids[i],
      "/",
      sample.ids[i],
      " cnv props joint probabilities.csv"
    )
  ))
  
  try(exp_props <- nb$clone_post %>%
        {
          setNames(.$p_cnv_x, .$cell)
        })
  try(Seurat.obj$exp_props <- exp_props)
  try(write.csv(
    exp_props,
    file = paste0(
      "/mnt/storage1/Anand_temp/HTAPP/numbat/processed2/",
      sample.ids[i],
      "/",
      sample.ids[i],
      " cnv props expression probabilities.csv"
    )
  ))
  
  try(allele_props <- nb$clone_post %>%
        {
          setNames(.$p_cnv_y, .$cell)
        })
  try(Seurat.obj$allele_props <- allele_props)
  try(write.csv(
    allele_props,
    file = paste0(
      "/mnt/storage1/Anand_temp/HTAPP/numbat/processed2/",
      sample.ids[i],
      "/",
      sample.ids[i],
      " cnv props allele probabilities.csv"
    )
  ))
  
  try(combined_props <- cbind(joint_props, exp_props, allele_props))
  
  try(all_sample_props <- rbind(all_sample_props, combined_props))
  pdf(
    paste0(
      "/mnt/storage1/Anand_temp/HTAPP/numbat/processed2/",
      sample.ids[i],
      "/",
      sample.ids[i],
      ".pdf"
    )
  )
  try(print(nb$plot_consensus()))
  
  try(print(nb$plot_mut_history()))
  
  try(print(nb$bulk_clones %>%
              filter(n_cells > 50) %>%
              plot_bulks(min_LLR = 20, # filtering CNVs by evidence
                         legend = TRUE)))
  
  dev.off()
  
  try(remove(Seurat.obj, nb, df_allele, counts)
  )
}

saveRDS(all_sample_props, "numbat cell scores.Rds")
