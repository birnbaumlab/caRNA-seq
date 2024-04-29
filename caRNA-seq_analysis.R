### CARNASEQ PIPELINE ###
# Version 3.0, edited 10/11/23
# v2: Added ability to filter variable genes (e.g. TCR VDJ-related genes)
# v3: Added kmeans analysis
# Slight changes to heatmap plotting functions, to allow passing of arguments to Heatmap call
# Added SCT vst.flavor = v2 argument
# Author: Caleb R. Perez

### IMPORTING REQUIRED LIBRARIES ###
library(tidyverse)
library(Seurat)
library(sctransform)
library(clustree)
library(ggrepel)
library(pheatmap)
library(DESeq2)
library(harmony)
library(EnhancedVolcano)
library(tigerstats)
library(ggpubr)
library(GSEABase)
library(GSVA)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(factoextra)
select <- dplyr::select # fixing package conflicts

### DIMENSIONAL REDUCTION AND CLUSTERING ###

'
General workflow:
normalize(seurat_obj)
pca(seurat_obj)
cluster(seurat_obj)
'

'
Normalizes gene expression and CAR BC data prior to dimensionality reduction. Gene expression is normalized by one of two methods: Log-normalization (default; workflow = "log-norm") or SCTransform (workflow = "sctransform").
Inputs:
workflow (str) = "LogNormalize" or "SCT (default)"
scale.all (bool) = Scales all genes (TRUE) or only variable genes (FALSE); default is FALSE.
scale.factor (int) = Factor to scale log-normalized UMI counts. Only relevant for workflow = "log-norm". Default is 10000.
vars.to.regress (list(str)) = Set of variables to regress out (e.g. cell cycle genes; vars.to.regress = c("S.Score", "G2M.Score")). Default is NULL, in which nothing is regressed out.
'
normalize <- function(seurat_obj,
                      workflow = 'SCT',
                      scale.all = TRUE,
                      scale.factor = 10000,
                      vars.to.regress = NULL,
                      vst.flavor = 'v2',
                      filter_variable_genes = NULL) {
  
  # Check for CITEseq data and set boolean
  citeseq <- 'ADT' %in% names(seurat_obj@assays)

  if (citeseq) {
    cat('Normalizing CITE-seq antibody counts...\n')
    seurat_obj <- NormalizeData(seurat_obj, assay='ADT', normalization.method='CLR', margin = 2)
    seurat_obj <- ScaleData(seurat_obj, assay = 'ADT', vars.to.regress = vars.to.regress)
  }
  
  cat('Normalizing gene expression counts...\n')
  seurat_obj <- NormalizeData(seurat_obj, scale.factor = scale.factor) # normalizes gene expression values to allow cell cycle scoring
  
  # Score by cell cycle
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes)
  
  if (workflow == 'LogNormalize') { # conventional normalization, scaling, etc.
    # Scale genes
    if (scale.all) {
      genes.to.scale <- rownames(seurat_obj)
    } else {
      genes.to.scale <- VariableFeatures(seurat_obj)
    }
    seurat_obj <- ScaleData(seurat_obj, features = genes.to.scale, vars.to.regress = vars.to.regress)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method='vst', nfeatures=2000)
  } else if (workflow == 'SCT') { # sctransform workflow
    seurat_obj <- SCTransform(seurat_obj, vars.to.regress = vars.to.regress, 
                              return.only.var.genes = scale.all,  # BUG should be set to !scale.all?
                              vst.flavor = vst.flavor)
  } else {
    cat('Input one of LogNormalize or SCT as workflow argument.\n')
  }
  
  # Filter variable genes if desired
  if (!is.null(filter_variable_genes)) {
    if (filter_variable_genes == 'TCR') {
      filter_variable_genes <- extract_TCR_genes(rownames(seurat_obj@assays$RNA))
    }
    cat('Filtering out the following genes from VariableFeatures:\n')
    print(filter_variable_genes)
    variable_genes <- VariableFeatures(seurat_obj)
    variable_genes <- variable_genes[!variable_genes %in% filter_variable_genes]
    VariableFeatures(seurat_obj) <- variable_genes
  }
  
  return(seurat_obj)
}

'
Calculates variance explained, used in pca function below. 
'
calculate_variance <- function(seurat_obj, reduction = 'pca') {
  cat('Variance explained by PCs...\n')
  pca <- seurat_obj[[reduction]]
  eigValues = (pca@stdev)^2  ## EigenValues
  varExplained = eigValues/sum(eigValues)
  print(cumsum(varExplained))
}

'
Calculate variance in PC space by CAR, as potential metric of CAR-intrinsic dispersity
(*writing*)
'
calculate_dispersities <- function(seurat_obj, reduction = 'pca', ident_type = 'CARID') {
  ident_list <- rownames(table(seurat_obj[[ident_type]][,1]))
  dispersity_df <- data.frame(row.names = ident_list)
  dispersity_df$dispersity <- 0
  for (ident in ident_list) {
    dispersity_df[ident,1] <- calculate_dispersity(seurat_obj, reduction, ident)
  }
  # Add dispersity to seurat obj metadata
}

calculate_dispersity <- function(seurat_obj, reduction, ident) {
  dispersity <- 0
  # Subset object
  # Calculate variance for each PC
  # Sum/weight
  return(dispersity)
}

'
Performs principal component analysis. Assumes data has been normalized and scaled appropriately,
i.e. using normalize function above.
'
pca <- function(seurat_obj, dimensionality_test = 'elbow', citeseq = F, ndims = 40, ...) {
  
  # Check for CITEseq data and set boolean
  #citeseq <- 'ADT' %in% names(seurat_obj@assays)
  
  # Run PCA on RNA counts
  cat('Performing PCA on gene expression...\n')
  seurat_obj <- RunPCA(seurat_obj, nfeatures.print=10, ...)
  calculate_variance(seurat_obj)
  
  # Display dimensionality metrics of choice
  if (dimensionality_test == 'elbow') {
    print(ElbowPlot(seurat_obj, ndims=ndims))
  } else if (dimensionality_test == 'jackstraw') {
    seurat_obj <- JackStraw(seurat_obj, num.replicate=100, dims=40)
    seurat_obj <- ScoreJackStraw(seurat_obj, dims=1:40)
    print(JackStrawPlot(seurat_obj, dims=1:40))
  }
  
  # Run PCA on ADT counts if present
  if (citeseq) {
    curr_assay <- DefaultAssay(seurat_obj)
    cat('Performing PCA on antibody counts...\n')
    DefaultAssay(seurat_obj) <- 'ADT'
    VariableFeatures(seurat_obj) <- rownames(seurat_obj[['ADT']])
    seurat_obj <- RunPCA(seurat_obj, reduction.name = 'apca', npcs=14, approx = FALSE,
                         nfeatures.print=2) #14 Ab panel
    calculate_variance(seurat_obj, reduction = 'apca')
    print(ElbowPlot(seurat_obj, reduction = 'apca', ndims = 13))
    DefaultAssay(seurat_obj) <- curr_assay
  }
  
  return(seurat_obj)
}

'
Unsupervised clustering based on input number of PCs
Leiden clustering algorithm by default
... args passed to FindClusters()
'
cluster <- function(seurat_obj,
                    reduction = 'pca',
                    pcs = 30,
                    pcs_ab = 10,
                    resolution = seq(from = 0, to = 1, by = 0.1),
                    wnn = FALSE,
                    algorithm = 4,
                    return.model = FALSE,
                    ...) {
  
  # If we wish to cluster on protein expression in addition to RNA:
  if (wnn) {
    cat('Pre-processing antibody data for weighted clustering...\n')
    seurat_obj <- FindMultiModalNeighbors(seurat_obj, reduction.list = list("pca", "apca"),
                                      dims.list = list(1:pcs, 1:pcs_ab),
                                      modality.weight.name = "RNA.weight")
    seurat_obj <- RunUMAP(seurat_obj, nn.name = 'weighted.nn', reduction.name = 'wnn.umap', reduction.key = 'wnnUMAP_', return.model = return.model)
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution, graph.name = 'wsnn', algorithm = algorithm,
                               method = 'igraph',
                               ...)
    reduction <- 'wnn.umap'
  } else {
    seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction, dims = 1:pcs)
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution, algorithm = algorithm, method = 'igraph',
                               ...)
    seurat_obj <- RunUMAP(seurat_obj, reduction = reduction, dims = 1:pcs, return.model = return.model)
    reduction <- 'umap'
  }
  
  # Output UMAP plot
  p1 <- reduction_plot(reduction, seurat_obj)
  print(p1)
  
  return(seurat_obj)
}

# Expects a single merged seurat object with metadata column that defines batches to correct
# (defined by group.by.vars)
# ... passed into normalize function
harmony <- function(merged, normalization.method, group.by.vars, ...) {
  # Normalize and pca all datasets together
  merged <- normalize(merged, workflow = normalization.method, ...)
  merged <- pca(merged)
  
  # Plot pca before harmonization
  p1 <- reduction_plot('pca', merged, ident_label = group.by.vars)
  
  # Run harmony algorithm
  if (normalization.method == "LogNormalize") {
    assay <- 'RNA'
  } else {
    assay <- 'SCT'
  }
  merged <- RunHarmony(merged, group.by.vars = group.by.vars, assay.use = assay,
                       reduction.save = 'harmony', plot_convergence = TRUE)
  
  # Plot pca post-batch correction
  p2 <- reduction_plot(reduction = 'harmony', merged, ident_label = group.by.vars)
  print(p1+p2)
  
  return(merged)
}

integrate <- function(seurat_obj_list,
                      workflow = 'LogNormalize',
                      vst.flavor = 'v2',
                      scale.all = TRUE,
                      scale.factor = 10000,
                      vars.to.regress = NULL,
                      filter_variable_genes = NULL,
                      method = 'cca',
                      integrate.adt = F, 
                      ...) {
  
  # Normalize all seurat_objs in list according to input parameters
  seurat_obj_list <- lapply(X = seurat_obj_list, FUN = function(x) {
    x <- normalize(x, workflow = workflow, scale.all = scale.all,
                   scale.factor = scale.factor, vars.to.regress = vars.to.regress,
                   vst.flavor = vst.flavor, filter_variable_genes = filter_variable_genes)
  })
  
  # Identify integration features
  if (workflow == 'SCT') {
    features <- SelectIntegrationFeatures(object.list = seurat_obj_list, nfeatures = 3000)
    seurat_obj_list <- PrepSCTIntegration(object.list = seurat_obj_list, anchor.features = features)
  } else {
    features <- SelectIntegrationFeatures(object.list = seurat_obj_list) # default of nfeatures=2000
  }

  if (method == 'rpca') {
    seurat_obj_list <- lapply(X = seurat_obj_list, FUN = RunPCA, features=features)
  }
  anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, anchor.features = features,
                                    normalization.method = workflow, reduction = method,
                                    ...)
  integrated <- IntegrateData(anchorset = anchors, normalization.method = workflow)
  DefaultAssay(integrated) <- 'integrated'
  
  if (workflow == 'LogNormalize') {
    integrated <- ScaleData(integrated)
  }
  
  
  if(integrate.adt) {
    cat('Integrating ADT assay...\n')
    integrated <- integrate_adt(seurat_obj_list, integrated, method = method)
  }
  
  return(integrated)
}

project <- function(query, reference, cluster_labels, normalization.method = 'SCT', dims = 35) {
  query <- normalize(query, workflow = normalization.method)
  anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:dims, reference.reduction = "pca")
  query <- MapQuery(anchorset = anchors, reference = reference, query = query, refdata = list(RefCluster = cluster_labels),
                   reference.reduction = 'pca', reduction.model = 'umap')
  p1 <- DimPlot(reference, reduction = 'umap', label = TRUE, repel = TRUE, group.by = cluster_labels) + 
    NoLegend() + ggtitle('Reference')
  p2 <- DimPlot(query, reduction = 'ref.umap', group.by = 'predicted.RefCluster', label = TRUE, repel = TRUE) + 
    NoLegend() + ggtitle('Projection')
  print(p1 + p2)
  return(query)
}

'
Assumes integration has been performed on RNA assay first
Returns integrated.rna object with integrated ADT data saved in ADT assay slot (overwrites raw ADT data)
This allows input into normal pca --> clustering workflow for WNN clustering.
'
integrate_adt <- function(seurat_obj_list, integrated.rna, method = 'cca') {
  for (x in seurat_obj_list) {
    DefaultAssay(x) <- 'ADT'
    VariableFeatures(x) <- rownames(x[['ADT']])
  }
  features <- SelectIntegrationFeatures(object.list = seurat_obj_list, assay = rep('ADT', each = length(list)))
  if (method == 'rpca') {
    seurat_obj_list <- lapply(X = seurat_obj_list, FUN = RunPCA, features=features)
  }
  anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, anchor.features = features,
                                    assay = rep('ADT', each = length(list)), dims = 1:13, reduction=method)
  integrated.adt <- IntegrateData(anchorset = anchors, new.assay.name = 'integrated.adt', dims = 1:13)
  DefaultAssay(integrated.adt) <- 'integrated.adt'
  integrated.adt <- ScaleData(integrated.adt)
  integrated.rna[['ADT']] <- integrated.adt[['integrated.adt']]
  return(integrated.rna)
}

### CLUSTER-DRIVEN ANALYSIS ###
'
Loops through all stored cluster resolutions, and identifies cluster markers for each set of clusters
cluster.prefix = everything prior to _res.0.X (e.g. SCT_snn for SCT_snn_res.0.1; or wsnn for WSNN_res.0.1)
'
all_cluster_markers <- function(seurat_obj, only.pos = FALSE, cluster.prefix = 'SCT_snn',
                                export.file.prefix = './cluster_markers', adj_p_thresh = 0.05) {
  resolutions <- str_subset(colnames(seurat_obj@meta.data), pattern = cluster.prefix)
  resolutions <- resolutions[resolutions != paste0(cluster.prefix, '_res.0')] # removes 1-cluster solution
  for (res in resolutions) {
    cat('Identifying cluster markers for:', res, '\n')
    Idents(seurat_obj) <- res
    filename <- paste0(export.file.prefix, '_', res, '.csv')
    cluster_markers(seurat_obj, only.pos = only.pos, export.file = filename, adj_p_thresh = adj_p_thresh)
  }
}

'
Identifies gene markers for all clusters, and outputs a CSV containing list of markers with summary statistics.
'
cluster_markers <- function(seurat_obj, only.pos = FALSE, export.file = NULL, adj_p_thresh = 0.05, ...) {
  cluster.markers <- FindAllMarkers(seurat_obj, only.pos = only.pos, ...) %>%
    filter(p_val_adj <= adj_p_thresh) %>%
    arrange(cluster, desc(avg_log2FC))
  if (!is.null(export.file)) {
    write.table(cluster.markers, export.file, sep = ",", col.names = NA)
  }
  return(cluster.markers)
}

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster, seurat_obj, grouping.var = 'orig.ident', assay = 'RNA'){
  FindConservedMarkers(seurat_obj, assay = assay,
                       ident.1 = cluster,
                       grouping.var = grouping.var) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

FindAllConservedMarkers <- function(cluster_list, seurat_obj, grouping.var = 'orig.ident', assay = 'RNA') {
  all_conserved_clusters <- map_dfr(.x = cluster_list, seurat_obj, grouping.var, assay, .f = get_conserved)
  return(all_conserved_clusters)
}

'
Identifies identity-specific markers, and outpus a CSV containing list of markers with summary statistics. Markers can be calculated relative to all other cells, or a specific other identity.
'
ident_markers <- function(seurat_obj, group.by, ident.1, ident.2 = NULL,
                          only.pos = FALSE, export.file = NULL, subset.ident = NULL,
                          adj_p_thresh = 0.05, ...) {
  markers <- FindMarkers(seurat_obj, ident.1 = ident.1, ident.2 = ident.2, group.by = group.by,
                         only.pos = only.pos, subset.ident = subset.ident, ...) %>%
    filter(p_val_adj <= adj_p_thresh) %>%
    arrange(desc(avg_log2FC))
  if (!is.null(export.file)) {
    write.table(markers, export.file, sep = ",", col.names = NA)
  }
  return(markers)
}

pathway_markers <- function(seurat_obj, assay, adj_p_thresh = 0.05, export.file = NULL, ...) {
  markers <- FindMarkers(seurat_obj, assay = assay, ...) %>% filter(p_val_adj < adj_p_thresh)
  markers$signed_log10p_adj <- log(markers$p_val_adj, 10) * sign(markers$avg_log2FC) * -1
  markers <- arrange(markers, desc(signed_log10p_adj))
  if (!is.null(export.file)) {
    write.table(markers, export.file, sep = ",", col.names = NA)
  }
  return(markers)
}

plot_marker_pVals <- function(markers, n = 10) {
  head_tail <- rbind(head(markers, n), tail(markers, n)) %>% rownames_to_column("pathway") %>%
    arrange()
  head_tail$pathway <- factor(head_tail$pathway, levels = head_tail$pathway)
  p <- ggplot(head_tail) + geom_col(aes(signed_log10p_adj, pathway), fill = "#076fa2", width = 0.6) +
    theme_bw() + labs(x='-Log10(Adj. P Value)', y='')
  print(p)
}

'
Loops through all stored clusters and calculates differentially-expressed genes across
ident groups within each cluster.
'
byCluster_ident_markers <- function(seurat_obj, group.by, ident.1, ident.2 = NULL, only.pos = FALSE, 
                                    export.file.prefix = './byCluster_markers', adj_p_thresh = 0.05) {
  clusters <- levels(seurat_obj@active.ident)
  for (cluster in clusters) {
    cat('Identifying ident markers for cluster ', cluster, '...\n')
    filename <- paste0(export.file.prefix, '_', cluster, '.csv')
    ident_markers(seurat_obj, group.by = group.by, ident.1 = ident.1, ident.2 = ident.2,
                  only.pos = only.pos, export.file = filename, subset.ident = cluster,
                  recorrect_umi=FALSE)
  }
}

'
Subset by specific idents (e.g. CARs), then identify group-specific markers. Active Ident must be set as subsets of interest.
'
CAR_stim_markers <- function(seurat_obj, subset_list = NULL, group.by = 'orig.ident', ident.1, ident.2 = NULL, only.pos = FALSE, export.file = NULL) {
  if (is.null(subset_list)) {
    subset_list <- rownames(table(seurat_obj$CARID)) # defaults to all CARs
  }
  output_list <- vector('list', length = length(subset_list))
  for (i in 1:length(subset_list)) {
    subset <- subset_list[i]
    cat('Processing subset: ', subset, '...\n', sep = '')
    markers <- FindMarkers(seurat_obj, ident.1 = ident.1, ident.2 = ident.2, group.by = group.by, subset.ident = subset,
                           only.pos = only.pos) %>%
      filter(p_val_adj <= 0.05) %>%
      arrange(desc(avg_log2FC)) %>%
      rownames_to_column(var = 'Gene')
    markers$subset <- subset
    output_list[[i]] <- markers
  }
  output <- bind_rows(output_list)
  
  if (!is.null(export.file)) {
    write.table(output, export.file, sep = ",", col.names = NA)
  }
  
  return(output)
}

'
Generate gene list based on DEG analysis using ref_object, then adds module score to object based on these lists.
'
add_DEG_score <- function(ref_object, seurat_obj, group.by, ident.1, ident.2, log2fc_thresh = 1,
                          module_name) {
  markers <- ident_markers(ref_object, group.by, ident.1, ident.2, only.pos = FALSE) %>%
    filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= log2fc_thresh) %>%
    rownames_to_column(var = 'Gene')
  up <- list(filter(markers, avg_log2FC > 0)$Gene)
  down <- list(filter(markers, avg_log2FC < 0)$Gene)
  seurat_obj <- AddModuleScore(seurat_obj, features = up, search = TRUE,
                           name = paste0(module_name, '.Up'))
  seurat_obj <- AddModuleScore(seurat_obj, features = down, search = TRUE,
                           name = paste0(module_name, '.Down'))
  return(seurat_obj)
}

'
Calculates cluster distributions and performs chi-squared enrichment analysis with respect to CAR BCs. Assumes active identity is set as cluster ID. Outputs a residuals heatmap that depicts CAR-specific enrichment in specific clusters.
... args passed to ComplexHeatmap constructor call in residuals_heatmap()
'
CAR_distribution_analysis <- function(seurat_obj, dist.file = NULL, chisq.res.file = NULL,
                                      extra_CAR_breaks = NULL, plot_all = FALSE, heatmap = TRUE, ...) {
  residuals <- ident_distribution_analysis(seurat_obj, ident2 = 'CARID', dist.file = dist.file,
                              chisq.res.file = chisq.res.file,
                              extra_CAR_breaks = extra_CAR_breaks,
                              plot_all = plot_all, heatmap = heatmap, ...)
  
  return(residuals)
}

'
Calculates cluster distributions and performs chi-squared enrichment analysis with respect to input identity. Assumes active identity is set as relevant cluster ID. Outputs a residuals heatmap that depicts ident-specific enrichment in specific clusters.
'
ident_distribution_analysis <- function(seurat_obj, ident1 = NULL, ident2, dist.file = NULL, 
                                        chisq.res.file = NULL, extra_CAR_breaks = NULL, plot_all = FALSE,
                                        heatmap=TRUE, residual_type = 'stdres', cluster = FALSE, ...) {
  if (is.null(ident1)) {
    ident_dist <- table(seurat_obj@active.ident, seurat_obj[[ident2]][,1])
  } else {
    ident_dist <- table(seurat_obj[[ident1]][,1], seurat_obj[[ident2]][,1])
  }
  if (!is.null(dist.file)) {
    write.table(ident_dist, dist.file, sep = ",", col.names = NA)
  }
  chisq <- chisq.test(ident_dist)
  if (residual_type == 'stdres') {
    residuals <- data.frame(chisq$stdres)
  } else if (residual_type == 'pearson') {
    residuals <- data.frame(chisq$residuals)
  }
  if (!is.null(chisq.res.file)) {
    write.table(residuals, dist.file, sep = ",", col.names = NA)
  }
  if (heatmap) {
    residual_heatmap(residuals_df = residuals, ident_name = ident2, extra_CAR_breaks =
                     extra_CAR_breaks, plot_all = plot_all, cluster = cluster, ...)
    #residual_heatmap <- function(residuals_df, ident_name = 'CARID', extra_CAR_breaks = NULL,
    #plot_all = FALSE, cluster = FALSE, ...)
  } 
  return(residuals)
}

# Percentage is calculated as percentage of ident1 (row) or ident2 (column)
CAR_contingency_heatmap <- function(seurat_obj, ident1 = NULL, ident2, percentage = 'column', 
                                    ...) {
  if (is.null(ident1)) {
    ident_dist <- table(seurat_obj@active.ident, seurat_obj[[ident2]][,1])
  } else {
    ident_dist <- table(seurat_obj[[ident1]][,1], seurat_obj[[ident2]][,1])
  }
  
  # Calculate row or column percentage as desired
  if (percentage == 'column') {
    mat <- colPerc(ident_dist)
    mat <- mat[-dim(mat)[1],]
  } else {
    mat <- rowPerc(ident_dist)
    mat <- mat[,-dim(mat)[2]]
  }
  
  h <- Heatmap(mat, cluster_rows = F, row_names_side = 'left', row_title = 'Cluster ID',
               rect_gp = gpar(col = "white", lwd = 2), column_dend_height = unit(1, 'in'),
               row_order = mixedsort(rownames(mat), decreasing = T), ...)
  h <- draw(h)
  return(h)
}

### PSEUDOBULKING METHODS ###

'
Calculation of pseudobulked gene expression profiles for groups defined by group.by string
'
pseudobulk <- function(seurat_obj, group.by = 'CARID', output_file = NULL, aggregate = TRUE,
                       filter_out_genes = NULL, only_keep_genes = NULL) {
  # Generate pseudobulked counts
  if (aggregate) {
    summed_counts <- AggregateExpression(seurat_obj, assays = 'RNA', group.by = group.by, slot = 'counts')
  } else {
    summed_counts <- AverageExpression(seurat_obj, assays = 'RNA', group.by = group.by, slot = 'counts')
  }
  summed_counts <- as.data.frame(summed_counts)
  colnames(summed_counts) <- sub('RNA.', '', colnames(summed_counts))
  summed_counts <- summed_counts[rowSums(summed_counts) > 1,] # filtering out 0-count genes
  
  if (!is.null(filter_out_genes)) {
    summed_counts <- summed_counts %>% filter(!row.names(summed_counts) %in% filter_out_genes)
  }
  
  if (!is.null(only_keep_genes)) {
    summed_counts <- summed_counts %>% filter(row.names(summed_counts) %in% only_keep_genes)
  }
  
  # Export as CSV if desired
  if (!is.null(output_file)) {
    write.table(summed_counts, output_file, sep = ",", col.names = NA)
  }
  
  return(summed_counts)
}

# Creates DESeq object with pseudobulked counts, using count_mat generated above
create_DESeq <- function(seurat_obj, count_mat, group.by = 'CARID') {
  ### GENERATE AND ADD RELEVANT METADATA ###
  #ident_list <- mixedsort(rownames(table(seurat_obj$CARID)))
  ident_list <- rownames(table(seurat_obj[[group.by]][,1]))
  metadata <- data.frame(Sample = ident_list, row.names = colnames(count_mat))
  
  ### CREATE DESEQ2 OBJECT ###
  deseq <- DESeqDataSetFromMatrix(countData = count_mat, colData = metadata, design = ~1)
  
  ### PRE-FILTERING GENES WITH ZERO COUNTS ###
  #deseq <- deseq[rowSums(counts(deseq)) > 1,] # already done in pseudobulk call above
  
  ### RLOG TRANSFORMATION AND PCA ###
  rld <- rlog(deseq, blind = TRUE)
  pcaData <- plotPCA(rld, intgroup = 'Sample', returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p1 <- ggplot(pcaData, aes(PC1, PC2, color=Sample)) + 
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], '% variance')) +
    ylab(paste0("PC2: ", percentVar[2], '% variance')) +
    NoLegend() + 
    theme(aspect.ratio = 1) + 
    geom_text_repel(aes(label = name))
  print(p1)
  #p1 <- plotPCA(rld, intgroup = 'Sample') + coord_fixed() + NoLegend()
  #print(p1 + geom_text_repel(aes(label = name)))
  
  ### PAIRWISE PEARSON CORRELATIONS BETWEEN SAMPLES ###
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  #pheatmap(rld_cor, annotation = metadata[,c('CARID'), drop = F]) + NoLegend()
  print(pheatmap(rld_cor))
  return(deseq)
}

km_cluster_avg <- function(seurat_obj, ident_label = NULL, assay, k, ...) {
  return (cluster_avg(seurat_obj, ident_label, assay, k, method = 'km'))
}

pam_cluster_avg <- function(seurat_obj, ident_label = NULL, assay, k, ...) {
  return (cluster_avg(seurat_obj, ident_label, assay, k, method = 'pam'))
}

multiK_cluster_avg <- function(seurat_obj, ident_label = NULL, assay, k_list, method = 'km', 
                               heatmaps = F, ...) {
  
  # Sets ident_label prefix as assay if null
  if (is.null(ident_label)) {
    ident_label <- assay
  }
  
  # Create ident_list using prefix and list of k values
  ident_list <- paste0(ident_label, '_k', k_list)
  
  # Loop through different k values, cluster, and save corresponding idents
  i <- 1
  for (ident in ident_list) {
    k <- k_list[i]
    cat('Clustering with k =', k, '...\n')
    seurat_obj <- cluster_avg(seurat_obj, ident_label = ident, assay = assay,
                              k = k, method = method, heatmap = heatmaps, ...)
    cat('CAR classes saved in ', ident, '.\n')
    i <- i + 1
  }
  
  return(seurat_obj)
}

# K-means clustering on pseudobulked average expression profiles using cluster::kmeans()
# Prints heatmap of average expression profiles by CAR after clustering, split by resulting km clusters
# Adds metadata information to seurat_obj with km cluster labels if ident_label passed as string,
# else, returns cluster vector
cluster_avg <- function(seurat_obj, ident_label = NULL, assay, k, method = 'km',
                        heatmap = T, return_heatmap = F, row_names = F, ...) {
  scaled_avg <- compute_pb_avg(seurat_obj, assay = assay, group.by = 'CARID', 
                               seurat_scale = F, transpose = T)
  set.seed(1)
  if (method == 'km') {
    km <- kmeans(scaled_avg, centers = k, iter.max = 100, nstart = 100)
    clusters <- km$cluster
  } else if (method == 'pam') {
    pam <- cluster::pam(scaled_avg, k = k, nstart = 100)
    clusters <- pam$clustering
  }
  
  name <- paste0(assay, ' Z Scores')
  h <- Heatmap(t(scaled_avg), show_row_names = row_names, name = name, column_split = clusters, ...)
  h <- draw(h)
  
  # Draw heatmap as desired
  if (heatmap) { # will do this once import_clustering function is fixed
    #name <- paste0(assay, ' Z Scores')
    #h <- Heatmap(t(scaled_avg), show_row_names = F, name = name, column_split = clusters, ...)
    #h <- draw(h)
  }
  
  # Save ident as desired
  if (!is.null(ident_label)) {
    seurat_obj <- import_clustering_heatmap(seurat_obj, h, ident_label)
    #seurat_obj <- import_clustering(seurat_obj, clusters, ident_label)
    return(seurat_obj)
  } else {
    if (return_heatmap) {
      return(h)
    } else {
      return (clusters)
    }
  }
}

# K-means clustering on pseudobulked average expression profiles using ComplexHeatmap
# Prints heatmap of average expression profiles by CAR after clustering, split by resulting km clusters
# Adds metadata information to seurat_obj with km cluster labels if ident_label passed as string
km_cluster_avg_heatmap <- function(seurat_obj, ident_label = NULL, assay, k, ...) {
  name <- paste0(assay, ' Z Scores')
  set.seed(1)
  h <- DoAverageHeatmap(seurat_obj, group.by = 'CARID', assays = assay, features = rownames(seurat_obj[[assay]]),
                        show_row_names = F, column_dend_height = unit(1, "in"), name = name,
                        clustering_method_columns = 'average', column_km = k, column_km_repeats = 100)
  if (!is.null(ident_label)) {
    seurat_obj <- import_clustering_heatmap(seurat_obj, h, ident_label)
    return(seurat_obj)
  }
}

# Sweep across k values for k-means clustering on pseudobulked CARs to determine optimal k
# Uses factoextra package
# ... passed into compute_pb_avg function
km_sweep <- function(seurat_obj, assay, k.max = 10, ...) {
  
  df <- compute_pb_avg(seurat_obj, assay = assay, transpose = T, ...)
  
  # Elbow method
  p1 <- fviz_nbclust(df, kmeans, method = "wss", k.max = k.max)
  
  # Silhouette method
  p2 <- fviz_nbclust(df, kmeans, method = "silhouette", k.max = k.max)
  
  # Gap statistic
  p3 <- fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", k.max = k.max,
                     nboot = 500, verbose = F)
  
  print(p1+p2+p3)
}

# Calculate and plot a consensus matrix given a set of cluster vectors
# Still writing...
consensus_matrix <- function(cluster_lists, k) {
  # Create an empty adjacency matrix
  N <- length(cluster_lists[[1]])
  adjacency_matrix <- matrix(0, nrow = N, ncol = N)
  
  # Iterate through each cluster membership list
  for (i in 1:length(cluster_lists)) {
    cluster_list <- cluster_lists[[i]]
    
    # Ensure order is the same across lists, get cluster values
    #cluster_list <- cluster_list[mixedsort(names(cluster_list))]
    #cluster_list <- unname(cluster_list)
    
    # Increment the count for sample pairs belonging to the same cluster
    for (cluster_id in 1:k) {
      samples_in_cluster <- which(cluster_list == cluster_id)
      if (length(samples_in_cluster) > 1) {
        sample_pairs <- combn(samples_in_cluster, 2)
        adjacency_matrix[sample_pairs] <- adjacency_matrix[sample_pairs] + 1
      }
    }
    
    return(adjacency_matrix)
  }
  
  # Divide the adjacency matrix by M to get the fraction of time each sample clusters with every other sample
  adjacency_matrix <- adjacency_matrix / length(cluster_lists)
  
  return(adjacency_matrix)
}

### METADATA ADDITION ###

# Labels cells as CD4 or CD8 T cells based on cluster membership, antibody expression,
# and/or RNA expression
# If clusters are given, they must be in a named list labeled as CD4 and CD8 cluster idents,
# and seurat_obj default identity should be set to these clusters
ID_CD4_CD8 <- function(seurat_obj, clusters = NULL,
                       filter_ADT = F, filter_RNA = T, column_name = 'CD4_8') {
  
  # Initialize cell lists with all cells
  CD4.cells <- colnames(seurat_obj)
  CD8.cells <- colnames(seurat_obj)
  
  # ID by cluster first if passed in
  if (!is.null(clusters)) {
    cat('Labeling by KNN cluster...\n')
    CD4.cells.clust <- WhichCells(seurat_obj, idents = unlist(clusters['CD4'], use.names = F))
    CD8.cells.clust <- WhichCells(seurat_obj, idents = unlist(clusters['CD8'], use.names = F))
    CD4.cells <- Reduce(intersect, list(CD4.cells, CD4.cells.clust))
    CD8.cells <- Reduce(intersect, list(CD8.cells, CD8.cells.clust))
  }
  
  default_assay <- DefaultAssay(seurat_obj)
  
  # Then by RNA
  if (filter_RNA) {
    cat('Labeling by RNA expression...\n')
    DefaultAssay(seurat_obj) <- 'RNA'
    CD4.cells.RNA <- WhichCells(seurat_obj, expression = CD4 > 0 & CD8A == 0)
    CD8.cells.RNA <- WhichCells(seurat_obj, expression = (CD8A > 0 & CD4 == 0) | (CD8B > 0 & CD4 == 0))
    CD4.cells <- Reduce(intersect, list(CD4.cells, CD4.cells.RNA))
    CD8.cells <- Reduce(intersect, list(CD8.cells, CD8.cells.RNA))
  }
  
  # Then by ADT
  if (filter_ADT) {
    cat('Labeling by ADT expression...\n')
    DefaultAssay(seurat_obj) <- 'ADT'
    CD4.cells.ADT <- WhichCells(seurat_obj, expression = CD4.1 > 0 & CD8 == 0)
    CD8.cells.ADT <- WhichCells(seurat_obj, expression = (CD8 > 0 & CD4.1 == 0))
    CD4.cells <- Reduce(intersect, list(CD4.cells, CD4.cells.ADT))
    CD8.cells <- Reduce(intersect, list(CD8.cells, CD8.cells.ADT))
  }
  
  DefaultAssay(seurat_obj) <- default_assay
  metadata <- data.frame(row.names = colnames(seurat_obj),
                         x = rep('Unassigned', length(colnames(seurat_obj))))
  metadata[CD4.cells,1] <- 'CD4'
  metadata[CD8.cells,1] <- 'CD8'
  seurat_obj <- AddMetaData(seurat_obj, metadata, col.name = column_name)
  return(seurat_obj)
}


'
Set of functions to add metadata information, including CAR groupings or geneset-level scores.
'

add_CAR_level_metadata <- function(seurat_obj, metadata_file) {
  metadata <- read.csv(file=metadata_file)
  colnames(metadata)[1] <- 'CARID'
  current_metadata <- seurat_obj@meta.data %>% dplyr::select(CARID) %>% rownames_to_column('CellBC')
  
  merged <- left_join(current_metadata, metadata, by = 'CARID') %>%
    column_to_rownames('CellBC') %>% dplyr::select(-CARID)
  
  seurat_obj <- AddMetaData(seurat_obj, merged, col.name = colnames(merged))
}

add_geneset_scores <- function(seurat_obj, genesets, geneset_type = 'csv', 
                               analysis_mode = 'UCell', assay_name, ...) {
  if (geneset_type == 'csv') {
    genesets <- read.csv(genesets, fill = FALSE)
  } else if (geneset_type == 'gmt') {
    genesets <- getGmt(genesets, geneIdType = SymbolIdentifier()) %>% geneIds()
    names(genesets) <- gsub("_", ".", names(genesets))
  } else {
    cat('Input one of csv or gmt as geneset_type.\n')
  }
  if (analysis_mode == 'gsva') {
    gene_expression <- as.matrix(seurat_obj@assays$RNA@data)
    gsva <- gsva(gene_expression, genesets)
    seurat_obj[[assay_name]] <- CreateAssayObject(data = gsva)
  } else if (analysis_mode == 'module.score' | analysis_mode == 'UCell') {
    seurat_obj <- AddScore(seurat_obj = seurat_obj, method = analysis_mode, 
                           assay_name = assay_name, genesets = genesets, ...)
  } else {
    cat('Input one of gsva, module.score, or UCell as analysis_mode argument.\n')
  }
  return(seurat_obj)
}

# AddModuleScore using UCell or Seurat methods
# Stored in a new assay slot as defined by assay_name, rather than metadata slots for downstream processing
AddScore <- function(seurat_obj, method, assay_name, genesets, ...) {
  seurat_copy <- seurat_obj
  if (method == 'module.score') { # Untested
    seurat_copy <- AddModuleScore(seurat_copy, features = genesets, name = assay_name, search = TRUE)
    feature_names <- paste0(assay_name, 1:length(names(genesets)))
  } else {
    seurat_copy <- AddModuleScore_UCell(seurat_copy, features = genesets, name = NULL, ...)
    feature_names <- names(genesets)
  }
  metadata <- seurat_copy@meta.data %>% dplyr::select(any_of(feature_names)) %>% t()
  seurat_obj[[assay_name]] <- CreateAssayObject(counts = metadata)
  rm(seurat_copy)
  return(seurat_obj)
}
# Import by-CAR clustering results into Seurat object as metadata
# NEEDS FIXING TO CONVERT TO WORK WITH CLUSTER VECTOR FROM KMEANS/PAM/ETC
import_clustering <- function(seurat_obj, cluster_vector, column_name) {
  

  
  # Get column labels
  #col_labels <- colnames(heatmap@ht_list[[1]]@matrix)
  
  # Turn list of indices into list of labels
  cluster_list <- lapply(clusters, function(indices) col_labels[indices])
  
  # Add to Seurat object via add_CAR_groups below
  seurat_obj <- add_CAR_groups(seurat_obj, CAR_subsets = cluster_list, ident_name = column_name)
  return(seurat_obj)
}

# Import by-CAR clustering results into Seurat object as metadata
import_clustering_csv <- function(seurat_obj, cluster.csv, column_name) {
  
  # Import cluster IDs from cluster.csv, and process into cluster_list (input to add_CAR_groups)
  cluster_list <- list()
  cluster_IDs <- read_csv(file = cluster.csv, show_col_types = FALSE)
  for (i in 1:ncol(cluster_IDs)) {
    CARs <- pull(drop_na(cluster_IDs[,i]))
    current_cluster <- colnames(cluster_IDs)[i]
    cluster_list[[current_cluster]] <- CARs
  }

  # Add to Seurat object via add_CAR_groups below
  seurat_obj <- add_CAR_groups(seurat_obj, CAR_subsets = cluster_list, ident_name = column_name)
  return(seurat_obj)
}

# Retrieves cluster memberships from input heatmap
# Assumes k-means clustering was performed on columns of heatmap
import_clustering_heatmap <- function(seurat_obj, heatmap, column_name) {
  
  # Extract clusters and associated column label indices
  clusters <- column_order(heatmap)
  clusters <- clusters[order(names(clusters))] # sort
  
  # Get column labels
  col_labels <- colnames(heatmap@ht_list[[1]]@matrix)
  
  # Turn list of indices into list of labels
  cluster_list <- lapply(clusters, function(indices) col_labels[indices])
  
  # Add Metadata
  seurat_obj <- add_CAR_groups(seurat_obj, CAR_subsets = cluster_list, ident_name = column_name)
  return(seurat_obj)
}

'
Add group-level metadata to Seurat object based on CAR groupings defined in input groups.csv file.
Can be combined with ident_markers with group.by = ident_name to identify group-specific markers.
CAR_subsets should be formatted as a list of lists, in which each individual list represents
an EXCLUSIVE group of CARs (i.e. one CAR cannot be represented in multiple lists)
'
add_CAR_groups <- function(seurat_obj, CAR_subsets, ident_name) {
  metadata <- data.frame(row.names = colnames(seurat_obj),
                         x = rep('None', length(colnames(seurat_obj))))
  Idents(seurat_obj) <- 'CARID'
  subset_names <- names(CAR_subsets)
  i <- 1
  for (CAR_subset in CAR_subsets) {
    cells <- WhichCells(seurat_obj, idents = CAR_subset)
    metadata[cells,1] <- subset_names[i]
    i <- i+1
  }
  seurat_obj <- AddMetaData(seurat_obj, metadata, col.name = ident_name)
  return(seurat_obj)
}

# same as above, but with option to use other metadata columns to define groupings
add_ident_groups <- function(seurat_obj, ident_subsets, group.by, new_ident_name) {
  metadata <- data.frame(row.names = colnames(seurat_obj),
                         x = rep('None', length(colnames(seurat_obj))))
  Idents(seurat_obj) <- group.by
  subset_names <- names(ident_subsets)
  i <- 1
  for (ident_subset in ident_subsets) {
    cells <- WhichCells(seurat_obj, idents = ident_subset)
    metadata[cells,1] <- subset_names[i]
    i <- i+1
  }
  seurat_obj <- AddMetaData(seurat_obj, metadata, col.name = new_ident_name)
  return(seurat_obj)
}

'
Generates groupings based on individual ICD identities (i.e. groups CARs by presence or absence of given ICD)
Outputs groups.csv file for input into add_CAR_groups in order to add grouping metadata to Seurat object.
If keep_union is set to TRUE, CARs with combinations of given ICDs will be retained as a separate grouping; otherwise, they will be removed.
'
generate_ICD_groupings <- function(CAR_ICD_file, ICD_list, output_file = NULL, keep_union = TRUE) {
   CAR_table <- read.csv(file = CAR_ICD_file, header = TRUE, stringsAsFactors = FALSE) %>% 
     data.frame(row.names=1)
   
   output_list <- list()

   for (ICD_name in ICD_list) {
     CAR_subset <- rownames(filter(CAR_table, if_any(everything(), ~ grepl(ICD_name, .x, fixed = TRUE))))
     output_list[[ICD_name]] <- CAR_subset
   }
   
   # calculate union
   union <- Reduce(intersect, output_list)
   
   # remove elements based on union
   for (ICD_name in ICD_list) {
     filtered <- output_list[[ICD_name]]
     filtered <- filtered[!(filtered %in% union)]
     output_list[[ICD_name]] <- filtered
   }
   
   if (keep_union) {
     output_list[['Multiple']] <- union
   }
   
   return(output_list)
}

# Computes average pseudobulked by group.by variable
# Can be scaled by feature, and returned transposed as desired
# Uses R scaling by default, but can also use seurat ScaleData() function, which sets
# a max scaled value -- useful when genes are features to limit effect of outliers
compute_pb_avg <- function(seurat_obj, assay, group.by = 'CARID', 
                           features = NULL, slot = 'data',
                           transpose = F, seurat_scale = F, scale = T) {
  
  # Use Seurat to average by ident
  avg <- AverageExpression(seurat_obj, group.by = group.by, slot = slot, assays = assay,
                           features = features, return.seurat = seurat_scale)
  if (scale) {
    if (seurat_scale) {
      avg <- ScaleData(avg)
      Z <- GetAssayData(avg, slot = 'scale.data')
    } else {
      avg <- avg[[1]]
      Z <- t(scale(t(avg)))
    }
  } else {
    Z <- avg[[1]]
  }
  
  if (transpose) {
    Z <- t(Z)
  }
  return(Z)
}

'
Average expression heatmap
'
DoAverageHeatmap <- function(seurat_obj, assays, group.by, features = NULL, slot = 'data',
                             return.Z.scores = F, transpose = F, gene_features = F, 
                             scale = T, ...) {
  # avg <- AverageExpression(seurat_obj, group.by = group.by, slot = slot, return.seurat = TRUE, assays = assays, features=features)
  # avg <- ScaleData(avg)
  # scaled_avg <- GetAssayData(avg, slot = 'scale.data')
  # if (transpose) {
  #   scaled_avg <- t(scaled_avg)
  # }
  scaled_avg <- compute_pb_avg(seurat_obj, assay = assays, group.by = group.by, features = features,
                               slot = slot, transpose = transpose, seurat_scale = gene_features,
                               scale = scale)
  
  p1 <- Heatmap(scaled_avg, row_names_side = 'left', ...)
  # Defaults are euclidean dist, complete linkage
  p1 <- draw(p1)
  if (return.Z.scores) {
    return(scaled_avg)
  }
  return(p1)
}



### PLOTTING & VISUALIZATION

reduction_plot <- function(reduction,
                           seurat_obj,
                           ident_label = NULL,
                           export_file = NULL,
                           width = 5,
                           height = 5,
                           legend = FALSE,
                           axes = TRUE,
                           title = FALSE) {
  if (legend) {
    p1 <- DimPlot(seurat_obj, reduction = reduction, pt.size=0.25, group.by = ident_label)
  } else {
    p1 <- DimPlot(seurat_obj, reduction = reduction, pt.size=0.25, group.by = ident_label, label=TRUE,
                  label.box = TRUE, repel = F, label.size = 6) + NoLegend()
  }
  
  if (!title) {
    p1 <- p1 + labs(title=NULL)
  }
  if (!axes) {
    p1 <- p1 + NoAxes()
  }
  print(p1)
  if (!is.null(export_file)) {
    pdf(file = export_file, width = width, height = height, pointsize = 12, bg = 'transparent')
    print(p1)
    dev.off()
  }
  return(p1)
}

'
PCA plot
'
pcaPlot <- function(seurat_obj,
                    ident_label = NULL,
                    export_file = NULL,
                    width = 5,
                    height = 5,
                    legend = FALSE,
                    axes = TRUE,
                    title = FALSE) {
  reduction_plot('pca', seurat_obj, ident_label, export_file, width, height, legend, axes, title)
}
'
UMAP plot
'
umapPlot <- function(seurat_obj,
                 ident_label = NULL,
                 export_file = NULL,
                 width = 5,
                 height = 5,
                 legend = FALSE,
                 axes = TRUE,
                 title = FALSE) {
  reduction_plot('umap', seurat_obj, ident_label, export_file, width, height, legend, axes, title)
}
'
UMAP plot
'
wnnPlot <- function(seurat_obj,
                     ident_label = NULL,
                     export_file = NULL,
                     width = 5,
                     height = 5,
                     legend = FALSE,
                     axes = TRUE,
                     title = FALSE) {
  reduction_plot('wnn.umap', seurat_obj, ident_label, export_file, width, height, legend, axes, title)
}

find_feature_order <- function(seurat_obj, feature, assay, group.by, slot = 'data') {
  avg <- AverageExpression(seurat_obj, features = feature, assays = assay, group.by = group.by)
  order <- avg[[assay]] %>% t() %>% data.frame() %>% arrange(desc(.)) %>% rownames()
  return(order)
}

ordered_vln <- function(seurat_obj, feature, assay, group.by, stat = 'mean', slot = 'data', ...) {
  copy <- seurat_obj
  order <- find_feature_order(copy, feature, assay, group.by, slot = slot)
  copy <- reorder_CARs(copy, ident = 'CARID', defined_order = order)
  if (stat == 'median') {
    p <- vln_median(copy, features = feature, group.by = group.by, assay = assay, ...)
  } else {
    p <- vln_mean(copy, features = feature, group.by = group.by, assay = assay, ...)
  }
  return(p)
}

#chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

'
Heatmap of standardized residuals calculated based on CAR distribution within clusters
'
residual_heatmap <- function(residuals_df, ident_name = 'CARID', extra_CAR_breaks = NULL,
                             plot_all = FALSE, cluster = FALSE, legend_title = 'Pearson Residuals',
                             ...) {
  if (cluster) {
    res_matrix <- residuals_df %>% pivot_wider(names_from = Var2, values_from = Freq) %>% 
      column_to_rownames(var = 'Var1') %>% as.matrix()
    p1 <- Heatmap(res_matrix, cluster_rows = F, row_names_side = 'left', row_title = 'Cluster ID', 
            name = legend_title, rect_gp = gpar(col = "white", lwd = 2),
            column_dend_height = unit(1, "in"),
            row_order = mixedsort(rownames(res_matrix), decreasing = T), ...)
  } else {
    if (ident_name == 'CARID') {
      breaks <- c('CAR5', 'CAR10', 'CAR15', 'CAR20', 'CAR25', 'CAR30', 'CAR35')
      if (!is.null(extra_CAR_breaks)) {
        breaks <- append(breaks, extra_CAR_breaks)
      } else if (plot_all) {
        levels(residuals_df$Var2) <- sub('CAR', '', levels(residuals_df$Var2))
        breaks <- levels(residuals_df$Var2)
      }
      p1 <- ggplot(residuals_df, aes(x = Var2, y = Var1, fill = Freq)) + 
        geom_tile(colour = 'white', size = 0.25) + labs(x='CAR ID', y='Cluster ID', fill = 'Residuals') + 
        scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
        scale_x_discrete(expand = c(0,0), breaks=breaks) + 
        theme_grey(base_size = 14) + 
        theme(plot.background=element_blank(),panel.border=element_blank())
    } else {
      p1 <- ggplot(residuals_df, aes(x = Var2, y = Var1, fill = Freq)) + 
        geom_tile(colour = 'white', size = 0.25) + labs(x=ident_name, y='Cluster ID', fill = 'Residuals') + 
        scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
        theme_grey(base_size = 14) + 
        theme(plot.background=element_blank(),panel.border=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
    }
  }

  print(p1)
  return(p1)
}

'
Exports multiple contour plots
'
contourPlots <- function(seurat_obj,
                         ident_list,
                         xHigh,
                         xLow,
                         yHigh,
                         yLow,
                         prefix = '',
                         ident_type = 'CARID',
                         parent_dir = '.',
                         reduction = 'umap',
                         axes = FALSE) {
  for (ident in ident_list) {
    filename <- paste0(parent_dir,'/',prefix,ident,'.pdf')
    pdf(file = filename, width = 3, height = 3, pointsize = 12)
    p1 <- contourPlot(seurat_obj, ident, ident_type, reduction, axes = axes)
    print(p1 & xlim(xLow, xHigh) & ylim(yLow, yHigh))
    dev.off()
  }
}

'
Single contour plot of CAR identity
'
contourPlot <- function(seurat_obj,
                        ident,
                        ident_type = 'CARID',
                        reduction = 'umap',
                        axes = FALSE) {
  
  Idents(seurat_obj) <- ident_type
  cells <- WhichCells(seurat_obj, idents = ident)
  subset_obj <- subset(seurat_obj, cells = cells)
  UMAP <- data.frame(Embeddings(subset_obj, reduction = reduction))
  key1 <- paste0(Key(subset_obj)[reduction],'1')
  key2 <- paste0(Key(subset_obj)[reduction],'2')
  if(axes) {
    p1 <- ggplot(UMAP, aes_string(x = key1, y = key2)) + stat_density2d(aes(color = stat(nlevel))) +
      scale_colour_distiller(palette = "RdYlBu") + theme_classic() + NoLegend()
  } else {
    p1 <- ggplot(UMAP, aes_string(x = key1, y = key2)) + stat_density2d(aes(color = stat(nlevel))) +
      scale_colour_distiller(palette = "RdYlBu")  + theme_void() + NoLegend()
  }

  return(p1)
}

'
Exports multiple CAR BC feature plots
'
CAR_feature_plots <- function(seurat_obj,
                              ident_list,
                              prefix = '',
                              parent_dir = '.',
                              reduction = 'umap') {
  for (ident in ident_list) {
    filename <- paste0(parent_dir,'/',prefix,ident,'.pdf')
    pdf(file = filename, width = 3, height = 3, pointsize = 9)
    p1 <- CAR_feature_plot(seurat_obj, ident, reduction = reduction)
    print(p1)
    dev.off()
  }
}

'
Single feature plot
'
CAR_feature_plot <- function(seurat_obj, ident, reduction = 'umap') {
  p1 <- FeaturePlot(seurat_obj, features = ident, reduction = reduction) + NoAxes() + NoLegend()
  return(p1)
}

'Violin plot with median line'
vln_median <- function(seurat_obj, features, group.by = 'CARID', slot = 'data', 
                       medLine_size = 6, ...) {
  #p1 <- VlnPlot(seurat_obj, features = features, pt.size=0, group.by = group.by) + NoLegend() + 
  #  stat_summary(fun = median, geom='point', size = 6, colour = "black", shape = 95)
  p1 <- VlnPlot(seurat_obj, features = features, pt.size=0, group.by = group.by, slot = slot, ...) + NoLegend() 
  p1 <- p1 & stat_summary(fun = median, geom='point', size = medLine_size, colour = "black", shape = 95)
  return(p1)
}

'Violin plot with median line'
vln_mean <- function(seurat_obj, features, group.by = 'CARID', slot = 'data', ...) {
  #p1 <- VlnPlot(seurat_obj, features = features, pt.size=0, group.by = group.by) + NoLegend() + 
  #  stat_summary(fun = median, geom='point', size = 6, colour = "black", shape = 95)
  p1 <- VlnPlot(seurat_obj, features = features, pt.size=0, group.by = group.by, slot = slot, ...) + NoLegend() 
  p1 <- p1 & stat_summary(fun = mean, geom='point', size = 6, colour = "black", shape = 95)
  return(p1)
}

StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

# Functions:
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)),
          plot.title = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}
#### END STACKED VLNPLOT ###

volcano_plot <- function(markers = NULL, seurat_obj = NULL, labeled_genes = NULL, 
                         FCcutoff = 1, adj_pCutoff = 0.05, ...) {
  if (is.null(markers)) {
    markers <- FindMarkers(seurat_obj, logfc.threshold = 0, ...)
  }
  if (is.null(labeled_genes)) {
    labeled_genes <- markers %>% filter(p_val_adj < adj_pCutoff & abs(avg_log2FC) > FCcutoff) %>% rownames()
  }
  bonferroniCorrection <- calculateBonferroni(markers)
  p1 <- EnhancedVolcano(markers, lab = rownames(markers), selectLab = labeled_genes, drawConnectors = TRUE,
                  x = 'avg_log2FC', y = 'p_val', cutoffLineType = 'blank', FCcutoff = FCcutoff, 
                  pCutoff = adj_pCutoff/bonferroniCorrection, title = NULL, legendPosition = 'none',
                  subtitle = NULL, caption = NULL, col = c('gray', 'gray', 'gray', 'red3'),
                  axisLabSize = 7, colAlpha = 1)
  print(p1)
  return(markers)
}

calculateBonferroni <- function(markers) {
  markers <- filter(markers, p_val_adj != 0 & p_val_adj != 1)
  return(markers[1,'p_val_adj'] / markers[1,'p_val'])
}

# Find topN cluster markers
find_topN <- function(seurat_obj, marker_file = NULL, marker_obj = NULL, n = 5) {
  if (!is.null(marker_file)) {
    markers <- read_csv(file = marker_file) %>% column_to_rownames(var='...1')
  } else if (!is.null(marker_obj)) {
    markers <- marker_obj %>% filter(p_val_adj <= 0.05) %>% arrange(cluster, desc(avg_log2FC))
  } else {
    markers <- cluster_markers(seurat_obj, only.pos = T, export.file = NULL, adj_p_thresh = 0.05)
  }
  topn <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = n)
  topn_genes <- unique(topn$gene)
  return(topn_genes)
}

# Plotting cluster markers
top_n_DotPlot <- function(seurat_obj, marker_file = NULL, marker_obj = NULL, n = 5, rotate = F, ...) {
  topn_genes <- find_topN(seurat_obj, marker_file, marker_obj, n)
  p <- DotPlot(seurat_obj, features = topn_genes, ...)
  if (rotate) {
    p <- p + coord_flip()
  }
  print(p)
  return(p)
}


### K-MEANS CLASS ANALYSIS

# Define column orders by k-means clusterings
get_km_cluster_list <- function(seurat_obj, ident_name) {
  
  # Define contingency table of CAR x class membership
  tab <- table(seurat_obj[['CARID']][,1], seurat_obj[[ident_name]][,1])
  
  # Get the row names for each column with non-zero values
  non_zero_row_names <- apply(tab != 0, MARGIN = 2, function(col) row.names(tab)[col])
  
  # Create a named list from the results
  result_list <- as.list(non_zero_row_names)
  
  # Set the names of the list elements to the column names
  names(result_list) <- colnames(tab)
  
  return(result_list)
}

# Define column splits using ordered cluster_list
get_column_splits <- function(result_list) {
  col_splits <- c()
  k <- 1
  for (x in result_list) {
    current_k <- paste0('k', k)
    split <- rep(current_k, length(x))
    col_splits <- c(col_splits, split)
    k <- k + 1
  }
  
  return(col_splits)
}

# Get split_df putting above two functions together, and sorting by CAR ID
get_splits_df <- function(seurat_obj, ident_name) {
  result_list <- get_km_cluster_list(seurat_obj, ident_name)
  col_splits <- get_column_splits(result_list)
  col_names <- unlist(result_list, use.names = F)
  df <- data.frame(CAR = col_names, splits = col_splits) %>% dplyr::slice(mixedorder(CAR))
  return(df)
}

# helper function for computing pseudobulked data in correct orientation
get_pb_data <- function(seurat_obj, assay) {
  return(compute_pb_avg(seurat_obj, assay = assay, transpose = T))
}

# Calculate kmeans helper function
calculate_km <- function(seurat_obj, data = NULL, k, assay) {
  if (is.null(data)) {
    data <- compute_pb_avg(seurat_obj, assay = assay, transpose = T)
  }
  set.seed(1)
  return(kmeans(data, centers = k, iter.max = 100, nstart = 100))
}

# Calculating distances between centroids and candidates -- which is most representative of class?
select_centroids <- function(seurat_obj, k, assay = 'SCENIC', exclude_CARs = NULL) {
  data <- get_pb_data(seurat_obj, assay)
  km.result <- calculate_km(data = data, k = k, assay = assay)
  centroids <- km.result$centers
  clusters <- km.result$cluster
  if (!is.null(exclude_CARs)) {
    clusters <- clusters[!(names(clusters) %in% exclude_CARs)]
  }
  selected_indices <- unique(clusters)
  for (i in selected_indices) {
    cluster_data <- data[clusters == i, ]
    centroid <- centroids[i, ]
    if (dim(cluster_data)[1] < 2) {
      cat('One CAR in cluster. Automatically assigning to centroid...')
      cat('Centroid #', i, 'is ', rownames(cluster_data)[1])
    } else {
      distances <- sqrt(rowSums((cluster_data - centroid)^2))
      selected_indices[i] <- which.min(distances)
      cat('Closest to centroid #', i, ' is...')
      print(distances[selected_indices[i]])
    }
  }
}

# Maximizing distance between CAR candidates
select_maxDist <- function(seurat_obj, assay, k, exclude_CARs = NULL) {
  
  # Calculate km
  data <- get_pb_data(seurat_obj, assay)
  km.result <- calculate_km(data = data, k = k, assay = assay )
  clusters <- km.result$cluster
  
  if (!is.null(exclude_CARs)) {
    clusters <- clusters[!(names(clusters) %in% exclude_CARs)]
  }
  
  # make vectors containing each cluster
  list <- list()
  for (i in 1:k) {list[[i]] <- names(clusters)[clusters == i]}
  list <- Filter(length, list) # removes zero-CAR clusters
  
  # Find all combinations
  dist_vect <- c()
  combos <- expand.grid(list, stringsAsFactors = F)
  for (i in 1:nrow(combos)) {
    currCombo <- unlist(combos[i,], use.names = F)
    currData <- data[currCombo,]
    total_dist <- Rfast::total.dist(currData)
    dist_vect[i] <- total_dist
  }
  combos$Dist <- dist_vect
  combos <- combos %>% arrange(desc(Dist))
  return(combos)
}

############### STILL WRITING ###############

'
Screen different factors (such as resolutions) and output DimPlots for each
'
reduction_plots <- function(seurat_obj, group.by_list, ...) {
  plots <- lapply(colnames(data2), plot_data_column, data = data2)
}

'
Balloon plot of CAR-cluster residuals.
'
balloonPlot <- function(seurat_obj, ident_name, export_file = NULL, width = 6, height = 3, 
                        base_size = 20, percentage = 'Cluster', ...) {
  # Calculate residuals and ident contingency table
  cluster_distribution <- table(seurat_obj@active.ident, seurat_obj[[ident_name]][,1])
  residuals <- ident_distribution_analysis(seurat_obj, ident2 = ident_name, heatmap=F, ...)
  
  # Calculate %s for plotting
  if (percentage == 'Cluster') {
    cluster_freq <- rowPerc(cluster_distribution)
    cluster_freq <- data.frame(cluster_freq[,-ncol(cluster_freq)])
  } else {
    cluster_freq <- colPerc(cluster_distribution)
    cluster_freq <- data.frame(cluster_freq[-nrow(cluster_freq),])
    cluster_freq <- rownames_to_column(cluster_freq, var = 'Var1') %>%
      pivot_longer(cols = !Var1) %>% arrange(name) %>% setNames(c('Var1', 'Var2', 'Freq'))
  }

  cluster_freq$Res <- residuals$Freq
  p <- ggballoonplot(cluster_freq, x = 'Var1', y = 'Var2', size = 'Freq', fill = 'Res', 
                size.range = c(1, 25)) + scale_fill_gradient2(low='blue', high='red') + 
    guides(size = FALSE, fill=guide_colorbar(title = 'Std. Residuals', frame.colour='black',
                                             barheight=10, barwidth=0.75))
  p <- p+theme_bw(base_size=base_size)+labs(x='', y='')
  print(p)
  
  if (!is.null(export_file)) {
    ggsave(filename = export_file, plot=p, dev = 'eps',
           width = width, height = height, units = 'in', dpi = 300)
  }
  return(p)
}

'
Apply one of various batch correction methods. Options include: integration or harmony.
Replaces normalize, pca, and cluster
'
batch_correct <- function(seurat_obj, method, normalization.method, ...) {
  if (method == 'harmony') {
    seurat_obj <- harmony(seurat_obj, normalization.method = normalization.method, ...)
  } else if (method == 'integration') {
    seurat_obj <- integrate(seurat_obj, workflow = normalization.method, ...)
  }
  return(seurat_obj)
}



