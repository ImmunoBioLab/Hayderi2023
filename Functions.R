library(Seurat)
library(patchwork)
library(magrittr)
library(stringr)
library(ggplot2)
library(pheatmap)

#Are all instances of vector X present in vector Y?
setGeneric("%v%", function(x, y) standardGeneric("%v%"))
setMethod("%v%", signature(x = "vector", y = "vector"), function(x, y) {
  all(x %in% y)
})
setMethod("%v%", signature(x = "vector", y = "list"), function(x, y) {
  all(sapply(y, function(z) x %v% z))
})

#Accessors
setGeneric("rna", function(object, features = NULL) standardGenercic("rna"))
setMethod("rna", signature(object = "Seurat"), function(object, features = NULL) {
  output <- object@assays$RNA
  
  if(is.null(features) == FALSE) {
    if(!features %v% rownames(object)) {
      features %<>% .[features %in% rownames(object)]
      
      warning("Some requested features are missing from the object!")
    }
    
    output %<>% .[rownames(.) %in% features,]
  }

  return(output)
})
setGeneric("rna<-", function(x, value) standardGenercic("rna<-"))
setMethod("rna<-", signature(x = "Seurat", value = "matrix"), function(x, value) {
  if(nrow(value) != nrow(x)) {
    warning("Inequality of genes detected!")
    
    value %<>% .[rownames(value) %in% rownames(x),]
  }
  
  x@assays$RNA <- value
  
  return(x)
})
setGeneric("colData", function(object, cells = NULL) standardGenercic("colData"))
setMethod("colData", signature(object = "Seurat"), function(object, cells = NULL) {
  output <- object@meta.data
  
  if(is.null(cells) == FALSE) {
    if(cells %v% colnames(object)) {
      cells %<>% .[cells %in% colnames(object)]
      
      warning("Some requested cells are missing from the object!")
    }
    
    output %<>% .[, colnames(.) %in% cells]
  }
  
  return(output)
})
setGeneric("colData<-", function(x, value) standardGenercic("colData<-"))
setMethod("colData<-", signature(x = "Seurat", value = "data.frame"), function(x, value) {
  if(nrow(value) != ncol(x)) {
    warning("Inequality of cells detected!")
    
    value %<>% .[rownames(value) %in% colnames(x),]
  }
  
  x@meta.data <- value
  
  return(x)
})

#Downsample a Seurat object per defined Metadata column
setGeneric("DownsampleSeurat", function(object, column, cells) standardGeneric("DownsampleSeurat"))
setMethod("DownsampleSeurat", signature(object = "Seurat"), function(object, column, cells) {
  metadata <- object@meta.data %>% split(., .[[column]])
  
  metadata %<>% lapply(., function(cellGroup) {
    #Generate a vector with cells number of TRUE, the rest - FALSE
    toKeep <- data.frame(toKeep = c(rep(TRUE, cells), rep(FALSE, nrow(cellGroup) - cells)),
                         Index = sample(c(1:nrow(cellGroup)), nrow(cellGroup))
                         ) %>%
      dplyr::arrange(., Index)
    
    cellGroup$toKeep <- toKeep$toKeep
    
    cellGroup %<>% .[.$toKeep == TRUE,]
    
    return(cellGroup)
  }) %>% do.call("rbind", .)
  rownames(metadata) <- metadata$Cell
  
  object <- object@assays$RNA %>% .[, colnames(.) %in% rownames(metadata)] %>%
    CreateSeuratObject(., meta.data = metadata, project = "GSE159677")
  
  return(object)
})

#Correlation
setGeneric("corToGene", function(object, geneName, metaGeneCol, features = NULL) standardGeneric("corToGene"))
setMethod("corToGene", signature(object = "Seurat", geneName = "character", metaGeneCol = "character"), function(object, geneName, metaGeneCol, features = NULL) {
  if(is.null(features) == TRUE) {
    features <- rownames(object)
  } else if(is.character(features) == FALSE) {
    stop("features must be a character vector!")
  } else if(!c(geneName, features) %v% rownames(object)) {
    stop("features are missing from the object!")
  }
  
  features %<>% c(geneName, .) %>% unique()
  
  geneCor <- rna(object, features = features)
  
  geneCor <- lapply(rownames(geneCor), function(gene) {
    corTmp <- cor.test(geneCor[geneName,], geneCor[gene,], method = "pearson")
    
    data.frame(Gene = gene,
               rVal = corTmp$estimate,
               pVal = corTmp$p.value
    )
  }) %>% do.call("rbind", .)
  
  geneCor %<>% dplyr::arrange(., desc(rVal))
  geneCor$fdr <- p.adjust(geneCor$pVal, method = "fdr", n = nrow(geneCor))
  geneCor$logFDR <- log(geneCor$fdr, 10) * -1
  
  #How many cells are top genes expressed in?
  geneCor %<>% dplyr::mutate(., PosCells = rna(object)[.$Gene,] %>% apply(., 1, function(gene) sum(gene > 0)))
  geneCor %<>% dplyr::mutate(., PosPercent = round(.$PosCells/ncol(rna(object))*100, 0))
  
  #How many gene+  and gene- cells are top genes expressed in?
  geneCor %<>% dplyr::mutate(., genePosCells = rna(object)[.$Gene, colData(object)[[metaGeneCol]] == "Positive"] %>% apply(., 1, function(gene) sum(gene > 0)))
  geneCor %<>% dplyr::mutate(., genePosPercent = round(.$genePosCells/ncol(rna(object)[, colData(object)[[metaGeneCol]] == "Positive"])*100, 0))
  geneCor %<>% dplyr::mutate(., geneNegCells = rna(object)[.$Gene, colData(object)[[metaGeneCol]] == "Negative"] %>% apply(., 1, function(gene) sum(gene > 0)))
  geneCor %<>% dplyr::mutate(., geneNegPercent = round(.$geneNegCells/ncol(rna(object)[, colData(object)[[metaGeneCol]] == "Negative"])*100, 0))
  geneCor %<>% dplyr::mutate(., ChDiff = .$genePosPercent - .$geneNegPercent)
  
  return(geneCor)
})
corPlot <- function(corDf, seurat, geneName, features = NULL) {
  if(is.null(features) == TRUE) {
    features <- corDf$Gene
  } else if(is.character(features) == FALSE) {
    stop("features must be a character vector!")
  } else if(!c(geneName, features) %v% corDf$Gene) {
    stop("features are missing from the corDf!")
  }
  
  corPlots <- lapply(seq_along(features), function(index) {
    gene <- data.frame(geneExp = as.vector(rna(seurat)[geneName,]),
                       GeneExp = as.vector(rna(seurat)[features[index],])
                       )
    
    pVal <- rsadCor %>% .[.$Gene == features[index], "fdr"] %>% formatC(., format = "e", digits = 2)
    rVal <- rsadCor %>% .[.$Gene == features[index], "rVal"] %>% round(., 2)
    
    ggplot(data = gene, aes(y = geneExp, x = GeneExp)) +
      geom_point(color = "black", alpha = 0.5) +
      geom_text(data = gene, mapping = aes(x = Inf, y = 2.5, vjust = 1, hjust = 1,
                                           label = paste0("r == ", rVal)), inherit.aes = FALSE, parse = TRUE) +
      geom_text(data = gene, mapping = aes(x = Inf, y = 2.5, vjust = 2.65, hjust = 1,
                                           label = paste0("p == ", pVal)), inherit.aes = FALSE, parse = TRUE) +
      xlab(features[index]) +
      ylab(geneName) +
      theme_classic() +
      theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), 
            panel.grid.major = element_line(color = "gray80", size = 0.2), panel.grid.minor = element_line(color = "gray90", size = 0.1))
  })
  
  return(corPlots)
}


