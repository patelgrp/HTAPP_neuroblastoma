library(ComplexHeatmap)
library(dplyr)

#set working directory
setwd("/mnt/storage1/Anand_temp/github_repos/HTAPP_neuroblastoma/scripts/Figure 1/")

#read oncoprint info
data <- read.csv(file = "20241227 oncoprint.csv", row.names = 1)

#rearrange data matrix to plot by stage -> MYCN status -> ALK status
data <- arrange(data, INSS.Stage, MYCN, ALK)

#convert into a matrix to play nice with the oncoprint function
mat <- as.matrix(data[, 1:29])
#reorder based on INSS stage (this is redundant)
mat <- mat[order(mat[, 6], decreasing = TRUE), ]

#pull ages to generate a dot plot
ages <- as.data.frame(x = mat[, 1], row.names = rownames(mat))
colnames(ages) <- "Age"
ages$Age <- as.numeric(ages$Age)
mat <- mat[, -1]
col.order <- colnames(mat[, c(5, 8, 9, 10, 19, 21:28)])


pdf("Oncoprint.pdf", width = 5, height = 10)
oncoPrint(
  mat = mat[, c(5, 8, 9, 10, 19, 21:28)],
  alter_fun = list(
    background = alter_graphic("rect", fill = "white"),
    '1' = alter_graphic("rect", fill = "#3954A1"),
    '2A' = alter_graphic("rect", fill = "#85CAE6"),
    '2B' = alter_graphic("rect", fill = "#5F9C9C"),
    '3' = alter_graphic("rect", fill = "#6CBD45"),
    '4' = alter_graphic("rect", fill = "#F5A11E"),
    '4S' = alter_graphic("rect", fill = "#E82829"),
    'Recurrent' = alter_graphic("rect", fill = "black"),
    'unknown' = alter_graphic("rect", fill = "#BDBABA"),
    'NB' = alter_graphic("rect", fill = "#3954A1"),
    'GNB' = alter_graphic("rect", fill = "#E9298E"),
    'amplified' = alter_graphic("rect", fill = "#80509C"),
    'non-amplified' = alter_graphic("rect", fill = "#F5A11E"),
    'wild type' = alter_graphic("rect", fill = "#F5A11E"),
    'mutant' = alter_graphic("rect", fill = "#80509C"),
    'Yes' = alter_graphic("rect", fill = "#6ABC45"),
    'No' = alter_graphic("rect", fill = "#3953A2"),
    'cell' = alter_graphic("rect", fill = "#6FA7B9"),
    'nucleus' = alter_graphic("rect", fill = "#D0B07B"),
    'yes' = alter_graphic("rect", fill = "black"),
    'both' = alter_graphic("rect", fill = "#CE7C8D")
  ),
  show_column_names = T,
  column_labels = c(
    "INSS stage",
    "histology",
    "MYCN",
    "ALK",
    "treatment",
    "sc/snRNA-seq",
    "WES",
    "RNA-seq",
    "DNA-me", 
    "multiplex IF",
    "MIBI", 
    "CODEX",
    "Slide-SeqV2"
  ),
  show_row_names = F,
  show_pct = F,
  row_order = rownames(mat),
  column_order = col.order,
  left_annotation = HeatmapAnnotation(age = anno_points(as.matrix(ages$Age), axis_param = list(direction = "reverse")), which = "row")
)
dev.off()