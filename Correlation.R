#Correlation
#All cells
rsadCor <- corToGene(scDown, "RSAD2", "RSADness", c("CXCL9", "CXCL10", "CXCL11"))

#Plot correlation with RSAD2 for all top genes in all cells
corPlots <- corPlot(rsadCor, scDown, "RSAD2", c("CXCL9", "CXCL10", "CXCL11"))

pdf(file.path(figDir, "CorGenes.pdf"), width = 15, height = 5)
ggpubr::ggarrange(plotlist = corPlots, nrow = 1, ncol = 3)
dev.off()

#Save
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Correlation")
openxlsx::writeData(wb, 1, rsadCor, startRow = 1, startCol = 1, colNames = TRUE, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file.path(dataDir, "CorData.xlsx"), TRUE)
rm(rsadCor, wb)



#VSMCs
rsadCor <- subset(scDown, subset = CellType %in% "VSMCs") %>%
  corToGene(., "RSAD2", "RSADness", c("CXCL9", "CXCL10", "CXCL11"))

#Plot correlation with RSAD2 for all top genes in all cells
corPlots <- corPlot(rsadCor, subset(scDown, subset = CellType %in% "VSMCs"), "RSAD2", c("CXCL9", "CXCL10", "CXCL11"))

pdf(file.path(figDir, "CorGenesVSMC.pdf"), width = 15, height = 5)
ggpubr::ggarrange(plotlist = corPlots, nrow = 1, ncol = 3)
dev.off()

#Save
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Correlation")
openxlsx::writeData(wb, 1, rsadCor, startRow = 1, startCol = 1, colNames = TRUE, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file.path(dataDir, "CorDataVSMC.xlsx"), TRUE)
rm(wb)