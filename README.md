# caRNA-seq

Author: Caleb R. Perez

This repository contains scripts used in the processing of caRNA-seq datasets generated in Perez et al., 2024 (pre-printed as of May 2024).

caRNA-seq is a scRNA-seq-based platform that allows for the simultaneous measurement of genome-wide transcriptional responses induced by a library of diverse CAR molecules, using the detection of CAR-specific barcodes to assign CAR identity to single cells. The functions in these scripts provide basic functionality for caRNA-seq analysis, including cell calling of CAR variants and conventional scRNA-seq analysis pipelines (e.g. normalization, clustering, DEG analysis, etc.). These are written to interface with the Seurat scRNA-seq analysis package (written fo v4, should be compatible with v5).

Scripts:

caRNA-seq_preprocessing.R

Provides basic dataset preprocessing functionality, including filtering on conventional QC metrics, calling cells by HTO, and calling cells by CAR BC. CAR identity is stored in a metadata column, labeled 'CARID'.

caRNA-seq_analysis.R

Provides basic analysis functionality, including normalization, clustering, DEG analysis, and geneset scoring, as well as functionality for CAR-level analysis (e.g. CAR distribution analyses).
