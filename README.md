## Overview

This repository contains code for the analysis of single cell RNA seq
data for the paper “” (Hayderi et al., 202X). Original data was
published in
<a href="https://www.nature.com/articles/s42003-022-04056-7">Alsaigh et
al., 2020</a>
<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677">(GSE159677)</a>.
Processed data was download from GEO and anaysed using the code
described here. Details on the analysis are described below.

## Contents

<ol>
<li>
<span class="subTitle">Functions</span>: described packages and
functions necessary for analysis below.
</li>
<li>
<span class="subTitle">Analysis</span>: describes data import, metadata
construction, quality control and filtering, downsampling and
normalization.
</li>
<li>
<span class="subTitle">Celltypes</span>: describes PCA, clustering, cell
type and subtype identification, and plotting of celltype-related data
</li>
<li>
<span class="subTitle">RSAD2</span>: describes labeling of RSAD2+ and
RSAD2- cells, their analysis and marker gene identification, as well as
plotting of data
</li>
<li>
<span class="subTitle">Correlation</span>: describes correlation of
chemokine expression to RSAD2
</li>
</ol>

## Analysis

Seurat pipeline was used for the analysis. Cells with less then 200 and
more then 6000 features were excluded. Further, cells where more then 5%
of genes were of mitochondrial origin were excluded. To equalize the
contribution of every patient and area of the aorta to the dataset,
dataset was downsampled to 2348 cells per sample. Log normalization was
performed with scale factor 10000. Seurat function FindVariableFeatures
was used to identify 2000 most significantly regulated genes. These
genes were used for PCA. Data was scaled using Seurat function ScaleData
and Louvain clustering was performed based on 12 top PCs with resolution
0.5. Uniform Manifold Approximation and Projection (UMAP) was used to
visualize clusters. Cell types were assigned based on expression of
common marker genes per cluster. As most cell types were assigned to
more then one cluster per cell type, cell types were further divided
into subtypes based on common marker genes. For T cells and macrophages,
to achieve this cells belonging to these subtypes were isolated and
separate Louvain clustering was performed, resolution 1. Cells with more
then 0 RSAD2 UMI were labelled as RSAD2-positive. Seurat function
FindMarkers was used to identify markers of RSAD2+ and RSAD2- cells.
Pearson correlation was used to assess correlation between expression of
RSAD2 and genes of interest.
