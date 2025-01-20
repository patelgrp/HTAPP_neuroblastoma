library(Seurat)
library(dplyr)
library(tidyr)
library(scCustomize)
library(SCpubr)
library(MetBrewer)
library(phylogram)
library(dplyr)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(htmlwidgets)
library(BioCircos)

#This addendum generates the color bars for the infercnv plots
setwd("/Volumes/groups/dyergrp/projects/ALSF_Pediatric_Atlas/common/REFERENCE_PUBLISHED/code for Github upload/")
gcdata <- readRDS("20220117 final HTAPP neuroblastoma LIGER.Rds")
###################### PANEL E ########################
#the heatmap was generated as part of the inferCNV pipeline (see png file)
#to generate the cell id bar on the observation section, see code below

#HTAPP-130
infercnv.dir <- "../HTAPP_core/infercnv/2020_09_15_neuroblastoma_inferCNV/HTAPP-130-SMP-91/"
dend.labels <- read.dendrogram(paste0(infercnv.dir,"infercnv.observations_dendrogram.txt")) #generated from the inferCNV pipeline

htapp.130 <- subset(gcdata, subset = orig.ident == "HTAPP-130-SMP-91")
htapp.130.labels <- htapp.130$annotate_refine_coarse
htapp.130.labels2 <- data.frame(annot = htapp.130.labels[labels(dend.labels)], cell.id = labels(dend.labels), col = "NA")
htapp.130.labels2$col[htapp.130.labels2$annot == "mesenchymal"] <- "#8A181A"
htapp.130.labels2$col[htapp.130.labels2$annot == "sympathoblast"] <- "#009F74"
htapp.130.labels2$col[htapp.130.labels2$annot == "adrenergic"] <- "#3953A4"
htapp.130.labels2$col[htapp.130.labels2$annot == "Myeloid"] <- "lightgrey"
htapp.130.labels2$col[htapp.130.labels2$annot == "T cell"] <- "lightgrey"
htapp.130.labels2$col[htapp.130.labels2$annot == "B cell"] <- "lightgrey"
htapp.130.labels2$col[htapp.130.labels2$annot == "Zona"] <- "lightgrey"
htapp.130.labels2$col[htapp.130.labels2$annot == "Endothelial"] <- "lightgrey"
htapp.130.labels2$col[htapp.130.labels2$annot == "Stroma"] <- "#CC6000"
htapp.130.labels2$col[htapp.130.labels2$annot == "Erythrocyte"] <- "lightgrey"
htapp.130.labels2$col[is.na(htapp.130.labels2$annot)] <- "lightgrey"

htapp.130.labels2$cell.id <- factor(htapp.130.labels2$cell.id, levels = htapp.130.labels2$cell.id)

p1 <- ggdendrogram(dend.labels, labels = F, leaf_labels = F, rotate = F) +
  scale_y_continuous(limits = c(0,200), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,length(htapp.130.labels2$annot)), expand = c(0,0))

p2 <- ggplot(htapp.130.labels2, aes(cell.id, y=1, fill = factor(annot))) + geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#8A181A", "#88CCE9", "#1C83BA", "lightgrey"))

p3 <- ggplot(htapp.130.labels2, aes(cell.id, y=1, fill = factor(annot))) + geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("dodgerblue", "dodgerblue", "dodgerblue", "goldenrod"))


gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)
gp3 <- ggplotGrob(p3)

maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)

#this produces a dendogram which can be overlaid over the png
pdf("~/Downloads/panel S2E-add HTAPP-130 cell id strip.pdf")
grid.arrange(gp1,gp2,ncol = 1, heights = c(4/5, 1/5))

maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp3$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp3$widths[2:5] <- as.list(maxWidth)

grid.arrange(gp1,gp3,ncol = 1, heights = c(4/5, 1/5))
dev.off()


remove(infercnv.dir, dend.labels, p1, p2, p3, gp1, gp2, gp3, maxWidth)


#HTAPP-171
infercnv.dir <- "../HTAPP_core/infercnv/2020_09_15_neuroblastoma_inferCNV/HTAPP-171-SMP-231/"
dend.labels <- read.dendrogram(paste0(infercnv.dir,"infercnv.observations_dendrogram.txt")) #generated from the inferCNV pipeline

htapp.171 <- subset(gcdata, subset = orig.ident == "HTAPP-171-SMP-231")
htapp.171.labels <- htapp.171$annotate_refine_coarse
htapp.171.labels2 <- data.frame(annot = htapp.171.labels[labels(dend.labels)], cell.id = labels(dend.labels), col = "NA")
htapp.171.labels2$col[htapp.171.labels2$annot == "mesenchymal"] <- "#8A181A"
htapp.171.labels2$col[htapp.171.labels2$annot == "sympathoblast"] <- "#009F74"
htapp.171.labels2$col[htapp.171.labels2$annot == "adrenergic"] <- "#3953A4"
htapp.171.labels2$col[htapp.171.labels2$annot == "Myeloid"] <- "lightgrey"
htapp.171.labels2$col[htapp.171.labels2$annot == "T cell"] <- "lightgrey"
htapp.171.labels2$col[htapp.171.labels2$annot == "B cell"] <- "lightgrey"
htapp.171.labels2$col[htapp.171.labels2$annot == "Zona"] <- "lightgrey"
htapp.171.labels2$col[htapp.171.labels2$annot == "Endothelial"] <- "lightgrey"
htapp.171.labels2$col[htapp.171.labels2$annot == "Stroma"] <- "#CC6000"
htapp.171.labels2$col[htapp.171.labels2$annot == "Erythrocyte"] <- "lightgrey"
htapp.171.labels2$col[is.na(htapp.171.labels2$annot)] <- "lightgrey"

htapp.171.labels2$cell.id <- factor(htapp.171.labels2$cell.id, levels = htapp.171.labels2$cell.id)

p1 <- ggdendrogram(dend.labels, labels = F, leaf_labels = F, rotate = F) +
  scale_y_continuous(limits = c(0,200), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,length(htapp.171.labels2$annot)), expand = c(0,0))

p2 <- ggplot(htapp.171.labels2, aes(cell.id, y=1, fill = factor(annot))) + geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#8A181A", "#88CCE9", "#1C83BA", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "#CC6000"))

p3 <- ggplot(htapp.171.labels2, aes(cell.id, y=1, fill = factor(annot))) + geom_tile() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("dodgerblue", "dodgerblue", "dodgerblue", "goldenrod", "goldenrod", "goldenrod", "goldenrod", "goldenrod"))


gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)
gp3 <- ggplotGrob(p3)

maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)

#this produces a dendogram which can be overlaid over the png
pdf("~/Downloads/panel S2E-add HTAPP-171 cell id strip.pdf")
grid.arrange(gp1,gp2,ncol = 1, heights = c(4/5, 1/5))

maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp3$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp3$widths[2:5] <- as.list(maxWidth)

grid.arrange(gp1,gp3,ncol = 1, heights = c(4/5, 1/5))
dev.off()


remove(infercnv.dir, dend.labels, p1, p2, p3, gp1, gp2, gp3, maxWidth)

