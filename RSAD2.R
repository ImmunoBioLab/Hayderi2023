#RSAD2+ cells
#Mark cells
scDown@meta.data$RSADness <- as.vector(scDown@assays$RNA["RSAD2", ]) > 0
colData(scDown)$RSADness <- ifelse(colData(scDown)$RSADness == TRUE, "Positive", "Negative")

tiff(file.path(figDir, "RSAD2_UMAP.tiff"), width = 10, height = 4, res = 600, units = "in")
ggpubr::ggarrange(plotlist = list(DimPlot(scDown, reduction = "umap", group.by = "RSADness", cols = c("gray90", "red")) +
                                    theme(plot.title = element_text(color = "white")),
                                  FeaturePlot(scDown, reduction = "umap", features = c("RSAD2"))
                                  )
                  )
dev.off()

tiff(file.path(figDir, "RSAD2_expression.tiff"), width = 5, height = 4, res = 600, units = "in")
subset(scDown, idents = c("T cells", "Endotheliocytes", "VSMCs", "Macrophages", "Mast cells", "pDCs")) %>% 
  subset(., subset = RSADness == "Positive") %>%
  VlnPlot(., features = c("RSAD2"), split.by = "Area")
dev.off()


#Which proportion of cells per celltype per area contribute to RSAD2 positivity?
cellData <- subset(scDown, subset = RSADness == "Positive") %>%
  colData(.) %>% .[, c("Area", "CellType")] %>% 
  table() %>% as.data.frame()
cellData %<>% split(., .$Area) %>% .[c(2, 1)] 

#PieChart using: https://r-graph-gallery.com/piechart-ggplot2.html
tiff(file.path(figDir, "RSAD2+_Celltype_PieChart.tiff"), width = 12, height = 6, res = 600, units = "in")
ggpubr::ggarrange(plotlist = lapply(cellData, function(Area) {
  ggplot(data = Area, mapping = aes(x = "", y = Freq, fill = CellType)) +
    geom_bar(stat = "identity", width=1, color="white") +
    coord_polar("y", start=0) +
    labs(title = unique(Area$Area)) +
    theme_void() +
    theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) 
}), common.legend = TRUE, legend = "right"
)
dev.off()


#Which proportion of cells per subtype per area contribute to RSAD2 positivity?
cellData <- subset(scDown, subset = RSADness == "Positive") %>%
  colData(.) %>% .[, c("Area", "SubType")] %>% 
  table() %>% as.data.frame()
cellData$SubType %<>% factor(., levels = c("Quiscent Endotheliocytes", "Activated Endotheliocytes", "Quiscent VSMCs", "Modulated VSMCs", "Th1", 
                                           "M1", "M2", "Fibroblasts", "B cells", "Mast cells", "pDCs"))
cellData %<>% split(., .$Area) %>% .[c(2, 1)] 

#PieChart using: https://r-graph-gallery.com/piechart-ggplot2.html
tiff(file.path(figDir, "RSAD2+_PieChart.tiff"), width = 12, height = 6, res = 600, units = "in")
ggpubr::ggarrange(plotlist = lapply(cellData, function(Area) {
  ggplot(data = Area, mapping = aes(x = "", y = Freq, fill = SubType)) +
    geom_bar(stat = "identity", width=1, color="white") +
    coord_polar("y", start=0) +
    labs(title = unique(Area$Area)) +
    theme_void() +
    theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) 
  }), common.legend = TRUE, legend = "right"
)
dev.off()

#Save
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "PieChart")
openxlsx::writeData(wb, 1, do.call("rbind", cellData), startRow = 1, startCol = 1, colNames = TRUE, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file.path(dataDir, "RSAD2+_PieChart.xlsx"), TRUE)
rm(wb, cellData)


#Count positive cells per celltype
sumMeta <- dplyr::group_by_at(colData(scDown), c("Area", "CellType", "RSADness"), drop = FALSE) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  dplyr::mutate(freq = n / sum(n) * 100)

#No RSAD+ cells were found in AC for fibroblasts or PA for B cells and pDCs - adding 0s
sumMeta %<>% rbind(., data.frame(Area = c("PA", "PA", "AC"), CellType = c("B cells", "pDCs", "Fibroblasts"), RSADness = rep("Positive", 3), n = rep(0, 3), freq = rep(0, 3)))

sumMeta$Area %<>% factor(., levels = c("PA", "AC"))

tiff(file.path(figDir, "RSAD2+_frequency.tiff"), width = 7, height = 6, res = 600, units = "in")
sumMeta %>% .[.$RSADness == "Positive",] %>%
  ggplot(data = .) +
  geom_bar(mapping = aes(x = CellType, y = freq, fill = Area),
           stat = "identity", position = position_dodge(0.7), width = 0.7) +
  scale_fill_manual(values = c(rgb(7, 137, 146, maxColorValue = 255), rgb(248, 118, 109, maxColorValue = 255))) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Percentage of RSAD2+ cells") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 14), axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 16), legend.title = element_blank(),
        legend.text = element_text(size = 12))
dev.off()

rm(metadata, sumMeta)


#Characterization of RSAD2+ cells
Idents(scDown) <- colData(scDown)$RSADness

rsad <- subset(scDown, subset = RSADness == "Positive")
rsad <- data.frame(RSADness = "Positive",
                   Area = factor(colData(rsad)$Area, levels = c("PA", "AC")),
                   RSAD2 = as.vector(rna(rsad) %>% .[rownames(.) == "RSAD2",])
                   )

tiff(file.path(figDir, "RSAD2_Area.tiff"), width = 4, height = 4, res = 600, units = "in")
ggplot(data = rsad, mapping = aes(x = Area, y = RSAD2, fill = Area)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width = 0.1, size = 1.3, show.legend = FALSE, outlier.colour = NA) +
  geom_point(size = 1, color = "gray50", alpha = 0.5, show.legend = FALSE, position = position_jitter(height = 0, width = 0.25)) +
  scale_fill_manual(values = c(rgb(7, 137, 146, maxColorValue = 255), rgb(248, 118, 109, maxColorValue = 255))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) +
  ylab("RSAD2 expression level") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16), axis.title.x = element_blank(), axis.line = element_line(size = 0.75),
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 16), legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  geom_pwc(method = "wilcox.test", label = "p = {p}")
dev.off()


#Find markers
markers <- FindMarkers(scDown, min.pct = 0.25, ident.1 = "Positive", ident.2 = "Negative")
markers %<>% dplyr::arrange(., desc(avg_log2FC))

DoHeatmap(scDown, group.by = "RSADness", features = rownames(head(markers, 20)))
VlnPlot(scDown, features = rownames(head(markers, 20)), group.by = "CellTypes")

tiff(file.path(figDir, "RSAD2+_markers.tiff"), width = 20, height = 20, res = 600, units = "in")
VlnPlot(scDown, features = rownames(head(markers, 20)), group.by = "RSADness")
dev.off()