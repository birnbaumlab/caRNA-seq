# caRNA-seq

Author: Caleb R. Perez

This repository contains scripts used in the processing of caRNA-seq datasets generated in Perez et al., 2024 (pre-printed as of May 2024).

caRNA-seq is a scRNA-seq-based platform that allows for the simultaneous measurement of genome-wide transcriptional responses induced by a library of diverse CAR molecules, using the detection of CAR-specific barcodes to assign CAR identity to single cells. The functions in these scripts provide basic functionality for caRNA-seq analysis, including cell calling of CAR variants, as well as wrappers for conventional scRNA-seq analysis pipelines (e.g. normalization, clustering, DEG analysis, etc.). These are written to interface with the [Seurat](https://satijalab.org/seurat/) scRNA-seq analysis package (written for v4, should be compatible with v5, although this was not tested).

For a detailed protocol on the sample preparation required to generate paired transcriptomic and CAR barcode datasets using a 10X Genomics Chromium platform, see [here](/docs/caRNA-seq_sample_prep_protocol.pdf).

## Broad overview of scripts included in this repository

### [caRNA-seq_preprocessing.R](src/caRNA-seq_preprocessing.R)
Provides basic dataset preprocessing functionality, including filtering on conventional QC metrics, calling cells by HTO, and calling cells by CAR BC. CAR identity is stored in a metadata column, labeled `CARID`.

### [caRNA-seq_analysis.R](src/caRNA-seq_analysis.R)

Provides basic analysis functionality, including normalization, clustering, DEG analysis, and geneset scoring, as well as functionality for CAR-level analysis (e.g. CAR distribution analyses).

See source code for detailed documentation for each function.


## Guided caRNA-seq analysis tutorial

In this tutorial, we will walk through a simple analysis of an example caRNA-seq dataset, in which ~29,000 T cells express one of 35 different CD19-targeted CAR variants, the identity of which is represented by 35 unique CAR-specific barcodes (CAR BCs) encoded in the CAR transcript. This library of CAR T cells was stimulated via coculture with CD19+ leukemia cells, then encapsulated and sequenced via a Chromium Next GEM Single Cell 5’ Kit v2, generating a gene expression library and a CAR BC library that were sequenced on an Illumina NovaSeq6000. The CellRanger count pipeline (10X Genomics) was used to generate cell-feature matrices from both libraries, in which the CAR BC library was treated as a “Custom” feature library type (for more information, see relevant CellRanger documentation [here](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-gex-count) and [here](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-libraries-csv#library-types), and an example of the [feature reference CSV file](data/CARBC_CellRanger_feature_reference.csv) for our library of CAR BCs that can be used for input into CellRanger). The resulting count matrices can be found uploaded to this repository, [here](/data/filtered_feature_bc_matrix). For simplicity, this data consists of cells from a single donor, making up a subset of the data shown in Figure 2-3 of our manuscript.

We begin with data import and pre-processing. First, we import our count matrices into a Seurat object in which transcript counts and CAR BC counts are stored in two different assays, `RNA` and `CARBC`, respectively. We provide support for CITE-seq data as well, which would be stored in an `ADT` assay; if cells were hashed with barcoded antibodies to multiplex samples, we also allow hash counts to be stored in an `HTO` assay for easy sample demultiplexing.


```r
# Setup working directory and source caRNA-seq scripts
repo_directory <- 'path/to/directory'
setwd(repo_directory)
source('./src/caRNA-seq_preprocessing.R')
source('./src/caRNA-seq_analysis.R')

# Import data and create Seurat object under default settings
seurat_obj <- create_seurat(mat_path = './data/filtered_feature_bc_matrix')
```


Next, we can filter out low-quality cells based on basic QC metrics drawn from the gene expression library. Here, we will use a minimum threshold of 1250 unique genes detected in each cell, and a maximum threshold of 5% of transcript counts coming from mitochondrial genes; each cell batch and sequencing dataset will differ on these thresholds, so we recommend optimization for your own applications.


```r
# Filter out low-quality cells on the basis of unique genes and mitochondrial reads
seurat_obj <- filter_seurat(seurat_obj, unique_gene_thresh = 1250, mt_thresh = 5)
```


As a final pre-processing step, we identify CARs expressed by each cell on the basis of CAR BC counts, utilizing the [MULTIseqDemux](https://satijalab.org/seurat/reference/multiseqdemux) algorithm. For each cell, this identifies which CAR BCs are expressed over background levels, allowing discrimination of cells expressing only a single CAR from those expressing multiple, or those expressing none. By default, we keep only cells expressing a single CAR for downstream analysis.


```r
# Assign each cell to individual CARs on the basis of CAR BC expression
seurat_obj <- CAR_demux(seurat_obj)
```


Following this pre-processing pipeline, we are left only with high quality cells that were confidently called to a single CAR in the library. The Seurat object we generated can be easily input into in any Seurat-compatible pipeline as desired for downstream analysis; we also provide a series of wrappers for commonly used Seurat functions to simplify this process. Here we will outline one analysis pipeline to explore CAR-intrinsic transcriptional profiles, by comparing average transcriptional profiles across CARs using pseudobulking. 

We begin by normalizing our transcriptional data; by default, we use the Seurat SCTransform algorithm, but our wrapper also provides support for log-normalization. At this point, we also identify variable genes. By default, we filter out TCR genes from these variable features to avoid clonality driving clustering or other downstream analysis.


```r
# Normalize gene expression counts and identify variable genes
seurat_obj <- normalize(seurat_obj, workflow = 'SCT', filter_variable_genes = 'TCR')
```


At this point, you could then proceed to dimensionality reduction and unsupervised clustering of cells. By default, our wrappers use the Leiden algorithm on the basis of the first 30 PCs, but we encourage optimization for each dataset. Cluster markers, alongside CAR membership within individual clusters, can then be analyzed to reveal CAR-specific signals.


```r
# PCA and cluster
seurat_obj <- seurat_obj %>% pca() %>% cluster(resolution = 0.4)
```


For the sake of this tutorial, we will take a simpler approach, in which we calculate the pseudobulked average gene expression profiles for each CAR in the library. We can do this using the `compute_pb_avg()` function. By default, we pseudobulk normalized transcript counts across all cells expressing each CAR (`group.by = ‘CARID’`), utilizing only variable genes, and we scale the results to allow for easy visualization of Z-scores in a heatmap. For visualization, we typically use the [ComplexHeatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html) package.


```r
# Pseudobulk by CAR to calculate average gene expression profiles across all cells expressing each CAR
pb_data <- compute_pb_avg(seurat_obj, group.by = 'CARID', features = VariableFeatures(seurat_obj))

# Visualize average expression profiles by CAR
ComplexHeatmap::Heatmap(pb_data, name = 'Gene Z Scores', show_row_names = F)
```


Alternatively, within a single function call, we can compute pseudobulked data and visualize the results:


```r
DoAverageHeatmap(seurat_obj, group.by = 'CARID', features = VariableFeatures(seurat_obj),
                 show_row_names = F, name = 'Gene Z Scores')
```


This analysis seems to suggest that CAR1 drives a unique transcriptional profile. We can look specifically at the genes driving this signature using differential gene expression analysis, computing genes significantly over or underexpressed in cells expressing CAR1 relative to those expressing all other CARs:



```r
# Perform DEG analysis, CAR1 vs. all other CARs
CAR1_markers <- CAR_markers(seurat_obj, CAR.1 = 'CAR1')

# Look at top 10 overexpressed genes
head(CAR1_markers, n = 10)

# Look at top 10 underexpressed genes
tail(CAR1_markers, n = 10)
```


And we find that these marker genes make sense biologically, as CAR1 is in fact the signaling-deficient negative control that we built into our library. Without any signaling domains to drive activation upon antigen binding, it follows that cells expressing this particular construct would show reduced expression of activation markers and increased expression of genes associated with naïve, unactivated T cells, when compared to the other functional CARs in the library. We can visualize these differences explicitly by comparing CAR1-expressing cells to those that express the two clinically approved CAR architectures, which are encoded in our library as CAR4 and CAR5 (1928z and 19BBz CARs, respectively).


```r
# Visualize expression of these markers compared to CAR4 and 5 (1928z and 19BBz, clinically-approved pos. ctrls)
Idents(seurat_obj) <- 'CARID'
vln_median(seurat_obj, features = c('IL7R', 'LEF1', 'SELL', 'IL2RA', 'GZMB', 'MKI67'),
           idents = c('CAR1', 'CAR4', 'CAR5'))
```



## Session Info

```r
sessionInfo()
```



