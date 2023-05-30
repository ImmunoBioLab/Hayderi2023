#Clustering and celltype allocation
#Find top variable features for PCA
scDown %<>% FindVariableFeatures(., selection.method = "vst", nfeatures = 2000)

scDown %>%
  ggpubr::ggarrange(plotlist = list(VariableFeaturePlot(.),
                                    LabelPoints(plot = VariableFeaturePlot(.), points = head(VariableFeatures(.), 10), repel = TRUE))
  )

#Scale
scDown %<>% ScaleData(., features = rownames(.))

#PCA
scDown %<>% RunPCA(., features = VariableFeatures(.))

DimPlot(scDown, reduction = "pca")
DimHeatmap(scDown, dims = 1:15, cells = 500, balanced = TRUE)

#Find useful PCs
scDown %<>% JackStraw(., num.replicate = 100)
scDown %<>% ScoreJackStraw(., dims = 1:20)

JackStrawPlot(scDown, dims = 1:20)

ElbowPlot(scDown)

#Cluster with Louvain algorithm
scDown %<>% FindNeighbors(., dims = 1:12)
scDown %<>% FindClusters(., resolution = 0.5, algorithm = 1)

scDown %<>% RunUMAP(., dims = 1:13)
DimPlot(scDown, reduction = "umap")

tiff(file.path(figDir, "UMAP_Clusters_&_Area.tiff"), width = 8, height = 4, res = 600, units = "in")
ggpubr::ggarrange(plotlist = list(DimPlot(scDown, reduction = "umap"),
                                  DimPlot(scDown, reduction = "umap", group.by = "Area", cols = c("red", "blue")))
)
dev.off()

#Find markers
markers <- FindAllMarkers(scDown, min.pct = 0.25)
markers %<>% split(., .$cluster)

save(scDown, file = file.path(projDir, "GSE159677.RData"))


#Find cell types based on common marker genes
#T cells: Clusters 0, 2, 15
VlnPlot(scDown, features = c("CD2", "CD3E", "CD3D", "CD4", "CD8A"))

#T cell subclustering
scTc <- subset(scTc, idents = "T cells")

#Cluster with Louvain algorithm
scTc %<>% FindNeighbors(., dims = 1:12)
scTc %<>% FindClusters(., resolution = 1, algorithm = 1)

scTc %<>% RunUMAP(., dims = 1:12)
DimPlot(scTc, reduction = "umap")

#T cell subsets
markerGenes <- c("CD4", "CD8A",
                 "CXCR3", "IFNG", "IL2", "TBX21", #Th1
                 "CCR4", "PTGDR2", "IL4", "IL5", "GATA3", #Th2
                 "IL9", "SPI1", "IRF4", #Th9
                 "KLRB1", "CCR6", "IL17A", "RORC", #Th17
                 "IL2RA", "IL10", "FOXP3", #Treg
                 "CCR10", "IL22", "AHR", "FOXO4" #Th22
                 ) 

tiff(file.path(figDir, "T_subset_markers.tiff"), width = 10, height = 4, res = 600, units = "in")
DotPlot(scTc, features = markerGenes) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

#NKT cells: Clusters 0 (2, 15) - based on CD3, with this clustering NKTs cluster together with T cells and do not form a separate cluster.
VlnPlot(scDown, features = c("NKG7", "CTSW", "XCL1", "KLRB1", "KLRD1"))

#Macrophages: Clusters 5, 6, 7, 12
#FCGR3A - CD16, FCGR1A - CD64, TFRC - CD71, ITGAM - CD11b, ITGAX - CD11c
VlnPlot(scDown, features = c("CD14", "CD68", "FCGR3A", "FCGR1A", "TFRC", "ITGAM", "ITGAX"))

#Macrophage subclustering: M2 - 0, 1, 3, 4, 5, 9, 11, 13
scMc <- subset(scDown, subset = seurat_clusters %in% c("5", "6", "7", "12"))

#Cluster with Louvain algorithm
scMc %<>% FindNeighbors(., dims = 1:12)
scMc %<>% FindClusters(., resolution = 1, algorithm = 1)

scMc %<>% RunUMAP(., dims = 1:12)
DimPlot(scMc, reduction = "umap")

VlnPlot(scMc, features = c("IL10", "IL12A", "FCGR3A", "FCGR1A", "MRC1", "CD14"))

m2cells <- colData(scMc) %>% .[.$seurat_clusters %in% c(0, 1, 3, 4, 5, 9, 11, 13), ] %>% rownames()
m1cells <- colData(scMc) %>% .[!.$seurat_clusters %in% c(0, 1, 3, 4, 5, 9, 11, 13), ] %>% rownames()

#VSMCs: Clusters 3, 8, 9. Modulated: Clusters 8, 9
VlnPlot(scDown, features = c("ACTA2", "MYH11", "CNN1", "TAGLN"))

tiff(file.path(figDir, "VSMC_markers.tiff"), width = 10, height = 4, res = 600, units = "in")
subset(scDown, subset = seurat_clusters %in% c("3", "8", "9")) %>% 
  VlnPlot(., features = c("TCF21", "SOST", "COL1A1", "CNN1"))
dev.off()

#Endothelial cells: Clusters 1, 10, 13. Activated: Clusters 1, 10
VlnPlot(scDown, features = c("PECAM1", "VWF", "CDH5", "ICAM1", "VCAM1", "SELE"))
subset(scDown, subset = seurat_clusters %in% c("1", "10", "13")) %>% 
  DotPlot(., features = c("PECAM1", "ICAM1", "VCAM1", "SELE"))

FeaturePlot(scDown, features = c("ICAM1", "PECAM1", "VWF", "CDH5"))

#B cells: Cluster 11
VlnPlot(scDown, features = c("CD19", "CD79A", "MS4A1", "IGKC"))

#Fibroblasts: Cluster 4
VlnPlot(scDown, features = c("DCN", "MFAP5", "APOD", "PDGFRA"))

#pDCs:
VlnPlot(scDown, features = c("CLEC4C", "KIT", "MS4A2", "CPA3"))


#Assign identities
celltypes <- c(0:16) %>% as.character()
names(celltypes) <- levels(scDown)
subtypes <- c(0:16) %>% as.character()
names(subtypes) <- levels(scDown)

celltypes[celltypes %in% c(0, 2, 15)] <- "T cells"
celltypes[celltypes %in% c(5, 6, 7, 12)] <- "Macrophages"
celltypes[celltypes %in% c(3, 8, 9)] <- "VSMCs"
celltypes[celltypes %in% c(1, 10, 13)] <- "Endotheliocytes"
celltypes[celltypes %in% c(11)] <- "B cells"
celltypes[celltypes %in% c(4)] <- "Fibroblasts"
celltypes[celltypes %in% c(14)] <- "Mast cells"
celltypes[celltypes %in% c(16)] <- "pDCs"

scDown %<>% RenameIdents(., celltypes)
scDown$CellType <- celltypes %>% .[match(colData(scDown)$seurat_clusters, as.double(names(.)))]
Idents(scDown) <- colData(scDown)$CellType

#Assign subtype identities
subtypes[subtypes %in% c(0, 2, 15)] <- "Th1"
subtypes[subtypes %in% c(5, 6, 7, 12)] <- "Macrophages"
subtypes[subtypes %in% c(3)] <- "Quiscent VSMCs"
subtypes[subtypes %in% c(8, 9)] <- "Modulated VSMCs"
subtypes[subtypes %in% c(1, 10)] <- "Activated Endotheliocytes"
subtypes[subtypes %in% c(13)] <- "Quiscent Endotheliocytes"
subtypes[subtypes %in% c(11)] <- "B cells"
subtypes[subtypes %in% c(4)] <- "Fibroblasts"
subtypes[subtypes %in% c(14)] <- "Mast cells"
subtypes[subtypes %in% c(16)] <- "pDCs"

scDown$SubType <- subtypes %>% .[match(colData(scDown)$seurat_clusters, as.double(names(.)))]
scDown$SubType[rownames(colData(scDown)) %in% m1cells] <- "M1"
scDown$SubType[rownames(colData(scDown)) %in% m2cells] <- "M2"
rm(scMc, m1cells, m2cells, subtypes, celltypes)


#Celltypes by area
tiff(file.path(figDir, "UMAP_CellTypes.tiff"), width = 6, height = 4, res = 600, units = "in")
DimPlot(scDown, reduction = "umap", label = TRUE, pt.size = 0.5) +
  theme(plot.title = element_blank())
dev.off()

tiff(file.path(figDir, "UMAP_Area.tiff"), width = 5, height = 4, res = 600, units = "in")
DimPlot(scDown, reduction = "umap", group.by = "Area", cols = c(rgb(248, 118, 109, maxColorValue = 255), rgb(7, 137, 146, maxColorValue = 255))) 
dev.off()

ggpubr::ggarrange(plotlist = list(DimPlot(subset(scDown, subset = Area == "PA"),
                                          reduction = "umap", label = TRUE, pt.size = 0.5),
                                  DimPlot(subset(scDown, subset = Area == "AC"),
                                          reduction = "umap", label = TRUE, pt.size = 0.5)
)
)

#Subtypes by area
Idents(scDown) <- scDown$SubType

tiff(file.path(figDir, "UMAP_Subypes_&_Area.tiff"), width = 15, height = 4, res = 600, units = "in")
ggpubr::ggarrange(plotlist = list(DimPlot(scDown, reduction = "umap", label = TRUE, pt.size = 0.5),
                                  DimPlot(scDown, reduction = "umap", group.by = "Area", cols = c(rgb(248, 118, 109, maxColorValue = 255), rgb(7, 137, 146, maxColorValue = 255))))
)
dev.off()


#Markers in all cell types
tiff(file.path(figDir, "UMAP_Celltype_Markers.tiff"), width = 12, height = 4, res = 600, units = "in")
DotPlot(scDown, features = c("NKG7", "CTSW", "XCL1", "KLRB1", "KLRD1", 
                             "CD2", "CD3E", "CD3D", 
                             "CD14", "CD68", "FCGR3A", "FCGR1A", "TFRC", "ITGAM", "ITGAX",
                             "ACTA2", "MYH11", "CNN1", "TAGLN",
                             "ICAM1", "PECAM1", "VWF", "CDH5",
                             "CD19", "CD79A", "MS4A1", "IGKC",
                             "DCN", "MFAP5", "APOD", "PDGFRA",
                             "CLEC4C", "KIT", "MS4A2", "CPA3")
) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
dev.off()


#What are the proportions of different cells in PA and AC?
cellData <- colData(scDown) %>% .[, c("Area", "SubType")] %>% table() %>%
  as.data.frame()
cellData$Area %<>% factor(., levels = c("PA", "AC"))
cellData$SubType %<>% factor(., levels = c("Quiscent Endotheliocytes", "Activated Endotheliocytes", "Quiscent VSMCs", "Modulated VSMCs", "Th1", 
                                           "M1", "M2", "Fibroblasts", "B cells", "Mast cells", "pDCs"))

tiff(file.path(figDir, "Subtype_frequencies.tiff"), width = 7, height = 6, res = 600, units = "in")
ggplot() +
  geom_bar(data = cellData, mapping = aes(x = Area, y = Freq, fill = SubType),
           stat = "identity") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Cell count") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 18),
        axis.line = element_line(size = 0.5), legend.title = element_blank(),
        legend.text = element_text(size = 10)
  )
dev.off()

#Save
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Subtype_frequencies")
openxlsx::writeData(wb, 1, cellData, startRow = 1, startCol = 1, colNames = TRUE, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file.path(dataDir, "Subtype_frequencies.xlsx"), TRUE)
rm(wb, cellData)