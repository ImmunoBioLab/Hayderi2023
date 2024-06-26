## Overview

This repository contains code for the analysis of single cell RNA seq
data for the paper “RSAD2 is abundant in atherosclerotic plaques and promotes interferon-induced CXCR3-chemokines in human smooth muscle cells” (<a href="https://www.nature.com/articles/s41598-024-58592-9">Hayderi et al., 2024</a>). The dataset analyzed here
was published and described in
<a href="https://www.nature.com/articles/s42003-022-04056-7">Alsaigh et
al., 2020</a>
<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677">(GSE159677)</a>.
Processed data was downloaded from GEO and analysed using the code
described here. Details on the analysis are described below.

## Contents

<ol>
<li>
<a href="https://github.com/ImmunoBioLab/Hayderi2023/blob/main/Functions.R">Functions</a>:
describes packages and functions necessary for analysis below.
</li>
<li>
<a href="https://github.com/ImmunoBioLab/Hayderi2023/blob/main/Analysis.R">Analysis</a>:
describes data import, metadata construction, quality control and
filtering, downsampling and normalization.
</li>
<li>
<a href="https://github.com/ImmunoBioLab/Hayderi2023/blob/main/Celltypes.R">Celltypes</a>:
describes PCA, clustering, cell type and subtype identification, and
plotting of celltype-related data
</li>
<li>
<a href="https://github.com/ImmunoBioLab/Hayderi2023/blob/main/RSAD2.R">RSAD2</a>:
describes labeling of RSAD2+ and RSAD2- cells, their analysis and marker
gene identification, as well as plotting of data
</li>
</ol>

## Analysis

Here, Seurat pipeline was used for the downstream analysis. Cells with
less than 200 and more than 6000 features, and cells where more than 5%
of genes were of mitochondrial origin were excluded. To equalize the
contribution of every patient and area of the aorta to the dataset,
dataset was downsampled by randomly selecting 2348 cells per sample. Log
normalization was performed with scale factor 10000. Seurat function
FindVariableFeatures was used to identify 2000 most significantly
regulated genes, and PCA was performed. Data was scaled using Seurat
function ScaleData and Louvain clustering was performed based on 12 top
PCs with resolution 0.5. Uniform Manifold Approximation and Projection
(UMAP) was used to visualize clusters. Cell types were assigned based on
expression of marker genes per cluster (Hayderi et al., 202X, Figure
S2). As most cell types were assigned to more than one cluster per cell
type, cell types were further divided into subtypes based on common
marker genes. For T cells and macrophages, to achieve this cells
belonging to these cell types were isolated and separate Louvain
clustering was performed, with resolution 1. Cells with more than 0
RSAD2 UMI were labelled as RSAD2-positive. Seurat function FindMarkers
was used to identify markers of RSAD2+ and RSAD2- cells.
