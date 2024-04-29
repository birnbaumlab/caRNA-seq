### CARNASEQ PIPELINE ###
# Version 4.0, edited 7/21/23
# Additions: CITEseq functionality, incorporation of MULTIseqDemux algorithm
# v4: Added ability to filter out genes (e.g. TCR VDJ genes) from count matrix or variable genes
# Author: Caleb R. Perez

'
General workflow:
seurat_obj <- create_seurat(mat_path, project_name, feature_library, +/- HTO_list) %>%
  filter_seurat(mt_thresh, unique_gene_thresh) %>%
  +/- HTO_demux() %>%
  CAR_demux()
'

### IMPORTING REQUIRED LIBRARIES ###
library(tidyverse)
library(Seurat)
library(matrixStats)
library(gtools)

### DATA PRE-PROCESSING ###
"
These functions handle dataset import and pre-processing on the basis of CAR BCs, filtering by per-cell CAR BC UMI counts, max assigned CAR BC UMI counts, and assigned CAR BC purity (% of all CAR BC UMIs assigned to max CAR BC). Functions are also provided to 
"

"
Creates seurat object with QC metrics and CARBC-related metadata added.
Uses input CARBC_df as metadata (e.g. output of process_CARBC_matrix function);
if none is given, generates a CARBC_df using given mat_path.
No filtering is performed at this stage (see filter_seurat and CAR_demux).
"
create_seurat <- function(mat_path = './filtered_feature_bc_matrix', preloaded_counts = NULL,
                          project_name = 'SeuratProject', feature_library = 'Custom', HTO_list = NULL,
                          remove_CARs = NULL, filter_genes = NULL,
                          CARBC_df = NULL, CARBC_df_file = NULL) {
  
  # Load in gene expression and CARBC/CITEseq counts if not pre-loaded
  if (is.null(preloaded_counts)) {feature_mat <- import_counts(mat_path)}
  gex_counts <- feature_mat[['Gene Expression']]
  ab_counts <- feature_mat[[feature_library]]
  
  # Filter out antibodies/CARs as desired
  ab_counts <- filter_ab_mat(ab_counts, remove_CARs)
  
  # Locate HTOs if HTO_list is given
  HTO_check <- !is.null(HTO_list)
  if (HTO_check) {
    cat('Hashtag data detected. Adding hashing antibody counts to HTO assay...\n')
    mat_list <- parse_HTO_mat(ab_counts, HTO_list)
    ab_counts <- mat_list[[1]]
    HTO_counts <- mat_list[[2]]
  }
  
  # Generate CARBC metadata if not provided
  if (is.null(CARBC_df)) {
    CARBC_matrix <- load_CARBC_matrix(preloaded_mat = ab_counts)
    CARBC_list <- process_CARBC_matrix(CARBC_matrix) #[metadata_df, counts]
    CARBC_df <-  CARBC_list[[1]]
    CARBC_counts <- CARBC_list[[2]]
  }
  
  # Filter any genes out of count matrix as desired
  if (!is.null(filter_genes)) {
    genes <- row.names(gex_counts)
    if (filter_genes == 'TCR') {
      filter_genes <- extract_TCR_genes(genes) # If desired, remove all TCR-related genes from count matrix
    }
    genes <- genes[!genes %in% filter_genes]
    gex_counts <- gex_counts[genes,]
  }
  
  # Select relevant metadata to add to Seurat object
  metadata <- CARBC_df %>% dplyr::select(Max_BC_Count, BC_Purity, Log2Count, CARID)

  # Create Seurat object, add CARBC/CITEseq counts, and HTOs if present
  seurat_obj <- CreateSeuratObject(counts = gex_counts, project = project_name)
  seurat_obj[['CARBC']] <- CreateAssayObject(counts = CARBC_counts)
  if (HTO_check) {seurat_obj[['HTO']]<- CreateAssayObject(counts = HTO_counts)}
  
  # Add CARBC metrics as metadata to Seurat object
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata)
  seurat_obj <- reorder_CARs(seurat_obj)
  CAR_list <- rownames(table(seurat_obj$CARID))
  
  # Check for citeseq data, add to Seurat object if present
  citeseq <- (length(rownames(ab_counts)) > length(CAR_list))
  if (citeseq) {
    cat('CITE-seq data detected. Adding antibody counts to ADT assay...\n')
    ab_counts <- ab_counts[!(rownames(ab_counts) %in% CAR_list), ] # remove CARBC data prior to adding
    seurat_obj[['ADT']] <- CreateAssayObject(counts = ab_counts)
  }
  
  # Add conventional mitochondrial gene count (conventional QC metric)
  seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Generate basic QC plots
  qc_plot(seurat_obj)
  
  return(seurat_obj)
}

"
Filters cells from input seurat_obj on the basis of mitochondrial gene content (percent.mt < mt_thresh) and unique genes detected (nFeature_RNA > unique_gene_thresh)
"
filter_seurat <- function(seurat_obj, mt_thresh = 10, unique_gene_thresh = 1000) {
    cat('Settings: Mitochondrial % = ', mt_thresh, ', Unique Gene Count = ', unique_gene_thresh, '.\n', sep='')
    cat('Starting cell count: ', ncol(seurat_obj), '.\n', sep = '')
    seurat_obj <- subset(seurat_obj, subset = percent.mt < mt_thresh & nFeature_RNA > unique_gene_thresh)
    cat('Filtered cell count: ', ncol(seurat_obj), '.\n', sep = '')
    cat('Filtered CAR library representation:')
    print(table(seurat_obj$CARID))
    return(seurat_obj)
  }

'
- method (str): One of MULTISeq or Purity
'
CAR_demux <- function(seurat_obj, max_filter = 3, purity_thresh = 0.9, method = 'MULTISeq',
                      filter_negatives = TRUE, filter_doublets = TRUE) {
  # Normalize CARBC counts, filter as desired
  cat('Normalizing CAR BC counts...\n')
  seurat_obj <- NormalizeData(seurat_obj, assay='CARBC', normalization.method='CLR', margin = 2)
  
  # Identify singlet CARs, filtering out doublets and negatives
  cat('Filtering CAR BC Matrix with ', method, ' algorithm', '...\n')
  cat('Starting cell count: ', ncol(seurat_obj), '.\n', sep='')
  if (method == 'MULTISeq') {
    seurat_obj <- filter_multi(seurat_obj, filter_negatives = filter_negatives, filter_doublets = filter_doublets)
  } else if (method == 'Purity') {
    seurat_obj <- filter_purity(seurat_obj, max_filter, purity_thresh)
  } else {
    cat('No cells identified. Input either MULTISeq or Purity.\n')
  }
  cat('Singlet cell count: ', ncol(seurat_obj), '.\n', sep='')
  
  cat('Final CAR library representation:')
  print(table(seurat_obj$CARID))
  
  print(CARBC_heatmap(seurat_obj))
  
  return(seurat_obj)
}

'
Parse ab-count matrix to ID and split HTO counts.
'
parse_HTO_mat <- function(counts, HTO_list) {
  HTO_counts <- counts[rownames(counts) %in% HTO_list, ]
  ab_counts <- counts[!(rownames(counts) %in% HTO_list), ]
  count_list <- list(ab_counts, HTO_counts)
  return(count_list)
}

'
Demultiplex hashed samples by HTO counts.
'
HTO_demux <- function(seurat_obj) {
  
  cat('Normalizing HTO counts...\n')
  seurat_obj <- NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR", margin = 2)
  seurat_obj <- HTODemux(seurat_obj)
  
  # Examine proportion of singlets vs. multiplets
  cat('Singlet/multiplet proportions are as follows:')
  print(table(seurat_obj$HTO_classification.global))
  cat('\n')
  print(HTOHeatmap(seurat_obj, assay = 'HTO', ncells = 1000))
  
  # Subset entire object first for singlets
  Idents(seurat_obj) <- 'HTO_classification.global'
  seurat_obj <- subset(seurat_obj, idents = 'Singlet')
  
  cat('Of singlets, hashtag proportions are as follows:')
  print(table(seurat_obj$HTO_maxID))
  cat('\n')
  
  Idents(seurat_obj) <- 'HTO_maxID'
  return(seurat_obj)
}

'
Generates basic QC scatterplots (e.g. nFeature_RNA vs. percent.mt, or nCount_RNA vs. nFeature_RNA)
'
qc_plot <- function(seurat_obj,
                    feature1 = 'nFeature_RNA',
                    feature2 = 'percent.mt',
                    feature3 = 'nCount_RNA',
                    feature4 = 'nFeature_RNA',
                    group.by = NULL) {
  pf1 <- FeatureScatter(seurat_obj, feature1 = feature1, feature2 = feature2, group.by = group.by)
  pf2 <- FeatureScatter(seurat_obj, feature1 = feature3, feature2 = feature4, group.by = group.by)
  print(pf1 + pf2)
}

qc_plot_interactive <- function(seurat_obj,
                                feature1 = 'nFeature_RNA',
                                feature2 = 'percent.mt',
                                extra_features = '',
                                yMax = NULL,
                                xMax = NULL) {
  if (!is.null(yMax) & is.null(xMax)) {
    plot <- FeatureScatter(seurat_obj, feature1 = feature1, feature2 = feature2) + ylim(NA, yMax)
  } else if (is.null(yMax) & !is.null(xMax)) {
    plot <- FeatureScatter(seurat_obj, feature1 = feature1, feature2 = feature2) + xlim(NA, xMax)
  } else if (!is.null(yMax) & !is.null(xMax)) {
    plot <- FeatureScatter(seurat_obj, feature1 = feature1, feature2 = feature2) + xlim(NA, xMax) + ylim(NA, yMax)
  } else {
    plot <- FeatureScatter(seurat_obj, feature1 = feature1, feature2 = feature2)
  }
  print(HoverLocator(plot, information = FetchData(seurat_obj, vars = c(feature1, feature2, extra_features))))
}

qc_hist <- function(seurat_obj, label_percentile = c(0.02, 0.98), defined_labels = NULL) {
  metadata <- seurat_obj@meta.data
  if (is.null(defined_labels)) {
    q1 <- quantile(metadata$nFeature_RNA, probs = label_percentile)
    q2 <- quantile(metadata$nCount_RNA, probs = label_percentile)
    q3 <- quantile(metadata$percent.mt, probs = label_percentile)
    cat('Number of gene percentiles:\n')
    print(q1)
    cat('Number of UMI percentiles:\n')
    print(q2)
    cat('Mitochondrial gene percentiles:\n')
    print(q3)
  } else {
    q1 <- defined_labels[1]
    q2 <- defined_labels[2]
    q3 <- defined_labels[3]
  }
  p1 <- ggplot(metadata, aes(x = nFeature_RNA, color = 'blue', fill = 'blue')) + 
    geom_density(alpha = 0.2) + 
    theme_classic() + xlab('Number of Genes Detected/Cell') +
    geom_vline(xintercept = q1) + NoLegend()
  p2 <- ggplot(metadata, aes(x = nCount_RNA, color = 'blue', fill = 'blue')) + 
    geom_density(alpha = 0.2) + 
    theme_classic() + xlab('Number of UMIs Detected/Cell') +
    geom_vline(xintercept = q2) + NoLegend()
  p3 <- ggplot(metadata, aes(x = percent.mt, color = 'blue', fill = 'blue')) + 
    geom_density(alpha = 0.2) + 
    theme_classic() + xlab('% mitochondrial genes') +
    geom_vline(xintercept = q3) + NoLegend()
  print(p1+p2+p3)
}

CARBC_heatmap <- function(seurat_obj, slot = 'counts', disp.max=10, print = T, 
                          subset = 2000, ...) {
  seurat_copy <- seurat_obj
  Idents(seurat_copy) <- 'CARID'
  CAR_list <- mixedsort(rownames(table(seurat_copy$CARID)))
  if (!is.null(subset)) {
    seurat_copy <- subset(seurat_copy, cells = sample(x = colnames(seurat_copy), size = subset))
  }
  plot <- DoHeatmap(seurat_copy, features = CAR_list, assay = "CARBC", slot = slot, 
                    group.by = 'CARID', angle = 90, size = 3,
                    disp.max=disp.max,...) + NoLegend()
  if (print) {print(plot)}
  return(plot)
}

filter_ab_mat <- function(ab_mat, remove_features) {
  # Remove any CARs as desired
  feature_list <- rownames(ab_mat)
  feature_list <- feature_list[!(feature_list %in% remove_features)]
  ab_mat <- ab_mat[rownames(ab_mat) %in% feature_list, ]
}



"
Parses 10X dataset to identify CAR BC matrices, and exports CAR BC matrix (filtered by per-cell UMI counts as desired), which can be used for input into Python preprocessing pipeline.
Inputs:
parent_dir (str) = Parent directory containing 10X output files, one of which is assumed to be 'filtered_feature_bc_matrix' directory
  Default value assumes filtered_feature_bc_matrix is in current directory.
output_file (str) = desired output filename or path
  Default value saves 'CAR BC Matrix.csv' file to current directory.
UMI_filter (int) = desired per-cell total UMI cutoff. All
  Default value filters out only cells with 0 detected CAR BC UMIs.
"
load_CARBC_matrix <- function(preloaded_mat = NULL,
                              parent_dir = NULL,
                              output_file = NULL,
                              library_type = 'Custom') {
  
  # If no preloaded matrix is provided, imports and isolates CARBCs
  if (is.null(preloaded_mat)) {
    cat('Importing CAR BC Matrix...\n')
    
    # If no parent_dir given, assume default CellRanger 10X directory is in current directory
    if (is.null(parent_dir)) {
      path <- './filtered_feature_bc_matrix'
    } else {
      path <- parent_dir
    }
    CARBC_matrix <- Read10X(data.dir = path)[[library_type]]
  } else {
    CARBC_matrix <- preloaded_mat
  }
  
  # Output raw CARBC file if desired
  if (!is.null(output_file)) {
    cat('CAR BC matrix saved in ', output_file, '.\n', sep='')
    write.table(CARBC_matrix, file = output_file, sep = ',', col.names = NA)
  }
  
  # Returns only CARBC matrix
  return(CARBC_matrix)
  
}

"
Filters CAR BC matrix on the basis of max BC UMI count and BC purity, where BC purity = Max BC UMI count / Total BC UMI count.
"
process_CARBC_matrix <- function(feature_mat) {
  
  # Isolate CARBCs from other potential CITEseq antibody counts
  CAR_list <- str_subset(rownames(feature_mat), pattern = 'CAR')
  CARBC_matrix <- feature_mat[rownames(feature_mat) %in% CAR_list, ]
  
  # Add CARBC metadata
  input <- t(as.matrix(CARBC_matrix))
  output <- as.data.frame(input)
  
  # Add CARBC metadata
  output['CARID'] <- colnames(output)[max.col(output, ties.method = 'first')]
  output['Max_BC_Count'] <- rowMaxs(input)
  output['Total_Count'] <- rowSums(input)
  output['BC_Purity'] <- output['Max_BC_Count'] / output['Total_Count']
  output['Log2Count'] <- log2(output['Total_Count']+1)
  #output['CARID'] <- max.col(input, ties.method = 'first')
  return(list(output, CARBC_matrix)) # metadata, counts
}

filter_purity <- function(seurat_obj,
                          max_filter = 3,
                          purity_thresh = 0.9) {
  if (purity_thresh == 0) {
    seurat_obj <- subset(seurat_obj, subset = Max_BC_Count >= max_filter)
  } else {
    seurat_obj <- subset(seurat_obj, subset = Max_BC_Count >= max_filter & BC_Purity >= purity_thresh)
  }
  return(seurat_obj)
}

filter_multi <- function(seurat_obj, filter_negatives = TRUE, filter_doublets = TRUE) {
  seurat_obj <- MULTIseqDemux(seurat_obj, assay = 'CARBC', autoThresh = TRUE)
  #seurat_obj <- MULTIseqDemux(seurat_obj, assay = 'CARBC', quantile = 0.2)
  Idents(seurat_obj) <- 'MULTI_ID'
  
  # Add global classification metadata (e.g. singlet vs. doublet vs. negative)
  MULTI_metadata <- dplyr::select(seurat_obj@meta.data, MULTI_ID)
  MULTI_metadata <- mutate(MULTI_metadata, Classification_Global = case_when(MULTI_ID == 'Negative' ~ 'Negative',
                                                                                   MULTI_ID == 'Doublet' ~ 'Doublet',
                                                                                   TRUE ~ 'Singlet'))
  MULTI_metadata <- dplyr::select(MULTI_metadata, Classification_Global)
  seurat_obj <- AddMetaData(seurat_obj, metadata = MULTI_metadata, col.name = 'Classification_Global')
  seurat_obj$Classification_Global <- factor(seurat_obj$Classification_Global, levels = c('Negative', 'Singlet', 'Doublet'))
  
  # Display and plot global classification QC
  cat('Singlet distribution is as follows:\n')
  print(table(seurat_obj$Classification_Global))
  plot_global_classifications(seurat_obj)
  
  # Filter out doublets, negatives, as determined by filter_negatives and filter_doublets
  if (filter_negatives & filter_doublets) {
    cat('Filtering out negatives and doublets...\n')
    filter_out <- c('Negative', 'Doublet')
  } else if (filter_doublets) {
    cat('Filtering out doublets only...\n')
    filter_out <- 'Doublet'
  } else if (filter_negatives) {
    cat('Filtering out negatives only...\n')
    filter_out <- 'Negative'
  } else {
    cat('Leaving negatives and doublets unfiltered...\n')
    filter_out <- NULL
  }
  if (!is.null(filter_out)) {seurat_obj <- subset(seurat_obj, idents = filter_out, invert = TRUE)}
  Idents(seurat_obj) <- 'orig.ident'
  
  # Reset CARIDs as MULTI_IDs
  seurat_obj <- reorder_CARs(seurat_obj, ident = 'MULTI_ID')
  seurat_obj$CARID <- seurat_obj$MULTI_ID
  
  return(seurat_obj)
}

'
Reorders CARs logically and removes factors with zero cells.
'
reorder_CARs <- function(seurat_obj, ident = 'CARID', defined_order = NULL) {
  if (!is.null(defined_order)) {
    CAR_list <- defined_order
  } else {
    CAR_list <- rownames(table(seurat_obj[[ident]])[table(seurat_obj[[ident]]) != 0])
    CAR_list <- mixedsort(CAR_list)
  }
  if (ident == 'CARID') {
    seurat_obj$CARID <- factor(seurat_obj$CARID, levels = CAR_list)
  } else if (ident == 'MULTI_ID') {
    seurat_obj$MULTI_ID <- factor(seurat_obj$MULTI_ID, levels = CAR_list)
  } else {
    return(CAR_list)
  }
  return(seurat_obj)
}

plot_global_classifications <- function(seurat_obj) {
  p1 <- vln_median(seurat_obj, features = 'Log2Count', group.by = 'Classification_Global') + 
    labs(y = 'Total CARBC UMIs, Log2(count+1)', title = '', x = '')
  p2 <- vln_median(seurat_obj, features = 'nCount_RNA', group.by = 'Classification_Global') + 
    labs(y = 'Total RNA UMIs Detected', title = '', x = '')
  p3 <- vln_median(seurat_obj, features = 'nFeature_RNA', group.by = 'Classification_Global') +
    labs(y = 'Unique Genes Detected', title = '', x = '')
  p4 <- vln_median(seurat_obj, features = 'percent.mt', group.by = 'Classification_Global') +
    labs(y = 'Mitochondrial Gene Composition (%)', title = '', x = '')
  print(p1|p2|p3|p4)
}

import_counts <- function(mat_path) {
  cat('Importing 10X Data...\n')
  feature_mat <- Read10X(data.dir = mat_path)
  return(feature_mat)
}

extract_TCR_genes <- function(gene_list) {
  TCR_genes <- gene_list[grepl("^TR[ABGD][CV]", gene_list)]
  return(TCR_genes)
}



