#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677
projDir <- file.path(getwd(), "GSE159677")

dir.create(file.path(projDir, "Data"))
dataDir <- file.path(projDir, "Data")
dir.create(file.path(projDir, "Figures"))
figDir <- file.path(projDir, "Figures")

#Import the matrix
scData <- ReadMtx(mtx = file.path(dataDir, "matrix.mtx"),
                  cells = file.path(dataDir, "barcodes.tsv"),
                  features = file.path(dataDir, "features.tsv")
                  ) %>% CreateSeuratObject(., project = "GSE159677")

#Import metadata keys
metaDf <- readr::read_delim(file.path(dataDir, "GSM4837528_Aggregated.Sample.Meta.txt"), skip = 4, col_names = FALSE, delim = " ")
colnames(metaDf) <- c("Patient", "Area", "Assay", "File", "Index", "IdType")
metaDf$Area %<>% factor(., levels = c("PA", "AC"))
metaDf$Index %<>% as.character(.) 

scData@meta.data$Index <- str_extract(rownames(scData@meta.data), "\\-\\d")
scData@meta.data$Cell <- rownames(scData@meta.data)
scData@meta.data %<>% dplyr::left_join(., metaDf[, c("Patient", "Area", "Index")], "Index")
rownames(scData@meta.data) <- scData@meta.data$Cell 

#QC
scData[["percent.mt"]] <- PercentageFeatureSet(scData, pattern = "^MT-")

VlnPlot(scData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

qcPlots <- ggpubr::ggarrange(plotlist = list(FeatureScatter(scData, feature1 = "nCount_RNA", feature2 = "percent.mt"),
                                  FeatureScatter(scData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
                             )
qcPlots <- FeatureScatter(scData, feature1 = "nCount_RNA", feature2 = "percent.mt")
save(qcPlots, file = file.path(figDir, "qcPlots.RData"))

#Filter out cells with more then 5% mitochondrial genes
scData %<>% subset(., subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

VlnPlot(scData, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Downsample
scData@meta.data$toDown <- str_c(scData@meta.data$Patient, scData@meta.data$Area, sep = " ")
scData@meta.data$toDown %>% table()

set.seed(1000)
scDown <- DownsampleSeurat(scData, column = "toDown", cells = min(table(scData@meta.data$toDown)))
scDown@meta.data$toDown %>% table()

#Normalize
scDown %<>% NormalizeData(., normalization.method = "LogNormalize", scale.factor = 10000)