library(Seurat)
library(UCell)

load(
  "/Users/apatel2/Downloads/Slideseq/drive-download-20251124T141718Z-1-001/2025_11_22_neuroblastoma_slideseq_HTAPP-102-SMP-11_6_analysis_updated.Rda")
)

van_Groningen <- read.csv("../durbin_nb_states/ADRN_MES_signatures.csv")

signatures <- list(ADRN = van_Groningen$ADRENERGIC, MES = van_Groningen$MESENCHYMAL)
gcdata <- AddModuleScore_UCell(gcdata, 
                                     features=signatures, name=NULL)

pdf(file = "~/Downloads/Slideseq/spatial_data_HTAPP-102_6.pdf")
print(SpatialDimPlot(gcdata))
print(SpatialFeaturePlot(gcdata, features = c("ADRN"), min.cutoff = "q10", max.cutoff = "q90") + ggplot2::scale_fill_continuous(low = "white", high = "darkred"))
print(SpatialFeaturePlot(gcdata, features = c("MES"), min.cutoff = "q10", max.cutoff = "q90") + ggplot2::scale_fill_continuous(low = "white", high = "darkred"))

print(VlnPlot(gcdata, features = c("ADRN", "MES"), pt.size = 0))
dev.off()
