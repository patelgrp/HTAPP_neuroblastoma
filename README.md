# HTAPP_neuroblastoma
Last updated: September 7, 2026

This repository contains code for processing and re-generating HTAPP figures related to Patel et al Cancer Cell "Spatial Omics Resolves Adrenergic and Mesenchymal Cell States in Neuroblastoma". A total of 54 patient samples, collected from St. Jude Children's Research Hospital and Dana-Farber Cancer Institute, were processed for multiple molecular (bulk DNA-seq, bulk RNA-seq, DNA methylation), single-cell/nucleus RNA-seq, and for multiple spatial omic methods (Slide-seqV2, CODEX, H&E, MIBI, and multiplexed immunofluorescence).
<img width="2297" height="607" alt="Patel et al  Fig 1 REV V5 overview" src="https://github.com/user-attachments/assets/ec78a740-a771-404c-aad7-eb7fbf02e59e" />

We aimed to create a cohort that was comprehensive, and that incorporated samples from every stage and spanned major molecular subtypes of disease, including MYCN-amplified and ALK mutant samples. Multiple samples were obtained after exposure to chemotherapy or combined chemotherapy + anti-GD2 (chemoimmunotherapy).
<img width="2273" height="1179" alt="HTAPP cohort" src="https://github.com/user-attachments/assets/c303ae8f-93e5-453b-8d69-894b180297ee" />

The code in this repository makes use of preprocessed data available through the HTAN portal (https://humantumoratlas.org/publications/hta1_2024_biorxiv_anand-g-patel). We recommend downloading the level 4 single-cell/nucleus RNA-seq Rds object and re-processing the data using the script in scripts/reprocess.R, which removes one failed sample and creates updated objects. The scripts folder houses code for each figure.  
