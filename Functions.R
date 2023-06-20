library(Seurat)
library(patchwork)
library(magrittr)
library(stringr)
library(ggplot2)
library(pheatmap)
library(ggpubr)

#Are all instances of vector X present in vector Y?
setGeneric("%v%", function(x, y) standardGeneric("%v%"))
setMethod("%v%", signature(x = "vector", y = "vector"), function(x, y) {
  all(x %in% y)
})
setMethod("%v%", signature(x = "vector", y = "list"), function(x, y) {
  all(sapply(y, function(z) x %v% z))
})

#Accessors
setGeneric("rna", function(object, features = NULL) standardGeneric("rna"))
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
setGeneric("rna<-", function(x, value) standardGeneric("rna<-"))
setMethod("rna<-", signature(x = "Seurat", value = "matrix"), function(x, value) {
  if(nrow(value) != nrow(x)) {
    warning("Inequality of genes detected!")
    
    value %<>% .[rownames(value) %in% rownames(x),]
  }
  
  x@assays$RNA <- value
  
  return(x)
})
setGeneric("colData", function(object, cells = NULL) standardGeneric("colData"))
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
setGeneric("colData<-", function(x, value) standardGeneric("colData<-"))
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


