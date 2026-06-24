# =============================================================================
# Hierarchical Clustering-Based Sample QC and Outlier Removal Pipeline
# =============================================================================
# Purpose:
#   Performs library-size-based quality control on pre-filtered circRNA BSJ
#   count matrices. For each dataset and filter level, samples are clustered
#   within each biological condition using hierarchical clustering on library
#   sizes, allowing visual identification of outlier samples. Manually curated
#   outliers are then removed, and the cleaned matrices are saved for downstream
#   simulation.
#
# Workflow:
#   1. Load pre-filtered count matrices (autofilter / min1 / min5) and metadata
#   2. Build DGEList objects
#   3. Cluster samples by library size within each condition
#   4. Inspect clustering plots and summary tables
#   5. Drop manually identified outlier samples
#   6. Intersect non-outlier samples across filter levels for consistency
#   7. Save cleaned count matrices and metadata for simulation input
#
# Inputs:
#   - Pre-filtered count matrices (CSV, genes × samples)
#   - Sample metadata CSVs (tab-separated, must contain Sample_name and Source)
#
# Outputs (written to for_simulation/{dataset}/):
#   - {dataset}_metadata.csv            : Metadata for retained samples
#   - {dataset}_{filter}_counts.csv     : Cleaned count matrix per filter level
#

# =============================================================================
# Input: Pre-filtered count matrix paths
#   Each dataset entry is a list of three paths in order:
#   [1] autofilter, [2] min1, [3] min5
# =============================================================================

datasets <- list(
  `BC` = list("/path/to/BC_filtered/automated_filtered_BC.csv","/path/to/BC_filtered/min1_filtered_BC.csv","/path/to/BC_filtered/min5_filtered_BC.csv"),
  `EBC1` = list("/path/to/EBC1_filtered/automatic_filtered_EBC.csv","/path/to/EBC1_filtered/min1_filtered_EBC.csv","/path/to/EBC1_filtered/min5_filtered_EBC.csv"),
  `EBC2` = list("/path/to/EBC2_filtered/automatic_filtered_EBC2.csv","/path/to/EBC2_filtered/min1_filtered_EBC2.csv","/path/to/filtered_datasets/EBC2_filtered/min5_filtered_EBC2.csv"),
  `HCC-PBMC` = list("/path/to/HCC_PBMC_filtered/automated_filtered_HCC_PBMC.csv","/path/to/HCC_PBMC_filtered/min1_filtered_HCC_PBMC.csv","/path/to/HCC_PBMC_filtered/min5_filtered_HCC_PBMC.csv"),
  `HCC-tissue` = list("/path/to/HCC_Tissue_filtered/automated_filtered_HCC_Tissue.csv","/path/to/HCC_Tissue_filtered/min1_filtered_HCC_Tissue.csv", "/path/to/HCC_Tissue_filtered/min5_filtered_HCC_Tissue.csv")
)

# =============================================================================
# Input: Sample metadata paths (tab-separated)
#   Required columns: Sample_name (matches count matrix column names),
#                     Source (biological condition / group label)
# =============================================================================
metadata <- list(
  `BC` = "/path/to/HC_filtered_final_datasets/metadata/PRJNA553624_metadata.csv",
  `EBC1` = "/path/to/Metadata_for_HC.csv",
  `EBC2` = "/path/to/HC_filtered_final_datasets/metadata/2025-11-25_EBC2.csv",
  `HCC-PBMC` = "/path/to/HC_filtered_final_datasets/metadata/PBMC_metadata_cleaned.csv",
  `HCC-tissue` = "/path/to/HC_filtered_final_datasets/metadata/PRJNA716508_metadata_cleaned.csv"
)

# =============================================================================
# Function: hc_filtered
# =============================================================================
# Clusters samples by library size within each biological condition using
# hierarchical clustering (Euclidean distance, complete linkage), and returns
# a diagnostic plot and summary table for manual outlier inspection.
#
# Parameters:
#   dgelist_filtered : A DGEList object with counts and sample metadata.
#                      Rownames of $samples must be valid sample IDs.
#   metadata         : data.frame or data.table of sample-level metadata.
#   nclust           : Maximum number of clusters per condition (default 8).
#                      Reduced automatically if a condition has fewer samples.
#   Run              : Column name in metadata containing sample IDs that match
#                      DGEList rownames (default "Run").
#   Sample           : Column name in metadata containing the group/condition
#                      label (default "Source").
#   filter           : Filter level label, used in plot subtitle (e.g. "min1").
#   dataset          : Dataset label, used in plot subtitle (e.g. "HCC-PBMC").
#
# Returns a named list:
#   matched_metadata        : Metadata rows reordered to match DGEList sample order
#   sample_libsize_clust_dt : data.table with SampleID, lib_size, Condition,
#                             and LibGroup (cluster assignment) per sample
#   plot                    : ggplot2 boxplot of library sizes faceted by condition
#   hc_summary_dt           : Per-cluster summary (N, min, median, max lib size)
#   lib_size_sample_summary : data.frame of raw library sizes with sample IDs
# =============================================================================
hc_filtered <- function(dgelist_filtered,
                        metadata,
                        nclust = 8,
                        Run = "Run",
                        Sample = "Source",
                        filter, dataset) {

  stopifnot(inherits(dgelist_filtered, "DGEList"))
  stopifnot(is.data.frame(metadata) || data.table::is.data.table(metadata))

  # ---- 1) Extract library sizes from DGEList ----
  lib_sizes <- dgelist_filtered$samples$lib.size
  sample_ids <- rownames(dgelist_filtered$samples)
  if (is.null(sample_ids)) {
    stop("dgelist_filtered$samples must have rownames corresponding to sample IDs (e.g., Run).")
  }
  names(lib_sizes) <- sample_ids

  # ---- 2) Match metadata rows to DGEList sample order ----
  if (!Run %in% colnames(metadata)) stop(sprintf("Column '%s' not found in metadata.", Run))
  if (!Sample %in% colnames(metadata)) stop(sprintf("Column '%s' not found in metadata.", Sample))

  matched_metadata <- metadata[match(sample_ids, metadata[[Run]]), , drop = FALSE]
  if (anyNA(matched_metadata[[Run]])) {
    missing <- sample_ids[is.na(matched_metadata[[Run]])]
    stop("These DGEList samples were not found in metadata[[Run]]: ",
         paste(missing, collapse = ", "))
  }

  # sanity check
  stopifnot(all(matched_metadata[[Run]] == sample_ids))

  # ---- 3) Build base data.table (one row per sample) ----
  lib_size_dt <- data.table::data.table(
    SampleID = matched_metadata[[Run]],
    Condition = matched_metadata[[Sample]],
    lib_size = as.numeric(lib_sizes)
  )

  # ---- 4) Hierarchical clustering within each condition ----
  # Euclidean distance on raw library sizes, complete linkage.
  # k is capped at nclust but reduced automatically for small conditions.
  clust_dt <- lib_size_dt[, {
    x <- lib_size
    k <- min(nclust, length(x))
    hc <- stats::hclust(stats::dist(x, method = "euclidean"), method = "complete")
    .(SampleID = SampleID,
      lib_size = lib_size,
      Condition.x = Condition,
      LibGroup = as.integer(stats::cutree(hc, k = k)))
  }, by = Condition]


  print(clust_dt)

  # ---- 5) Diagnostic plot: library size distribution per cluster per condition ----
  plot <- ggplot2::ggplot(
    clust_dt,
    ggplot2::aes(x = factor(LibGroup), y = lib_size)
  ) +
    ggplot2::geom_boxplot(varwidth = TRUE, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, shape = 1) +
    ggplot2::facet_wrap(~Condition.x, scales = "free_x") +
    ggplot2::labs(
      title ="Library Size Clustering (within each condition)",
      subtitle = paste("Dataset: ", dataset, " – filter:", filter),
      x = "Library Size Cluster",
      y = "Library Size (BSJ counts)"
    ) +
    ggplot2::theme_minimal()

  # ---- 6) Per-cluster summary statistics ----
  summary_dt <- clust_dt[, .(
    Nsamples = .N,
    MinLibSize = min(lib_size),
    MedianLibSize = median(lib_size),
    MaxLibSize = max(lib_size)
  ), by = .(Condition.x, LibGroup)][order(Condition.x, LibGroup)]
  
  # Simple per-sample library size table (used for IQR-based inspection)
  lib_size_sample_summary <- data.frame(libsizes = dgelist_filtered$samples$lib.size,
                                        sample_ids = rownames(dgelist_filtered$samples))

  return(list(
    matched_metadata = matched_metadata,
    sample_libsize_clust_dt = clust_dt,
    plot = plot,
    hc_summary_dt = summary_dt,
    lib_size_sample_summary = lib_size_sample_summary
    ))
}

# =============================================================================
# Assign filter-level names to each dataset's path list
#   Order must match the three paths provided per dataset above:
#   [1] autofilter, [2] min1, [3] min5
# =============================================================================
filter_names <- c("autofilter", "min1", "min5")

datasets <- lapply(datasets, function(x) {
  stopifnot(length(x) == length(filter_names))
  setNames(x, filter_names)
})


all_results <- lapply(names(datasets), function(ds) {
  
  message("Processing dataset: ", ds)

  # Load metadata
  meta <- read.csv(
    metadata[[ds]],
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  ds_filters <- datasets[[ds]]
  
  filter_results <- lapply(names(ds_filters), function(filt) {
    
    message("  Filter: ", filt)

    # Load count matrix
    mtx <- read.csv(
      ds_filters[[filt]],
      sep = ',',
      header = TRUE,
      row.names = 1,
      check.names = FALSE
    )

    # Verify all count matrix columns are present in metadata, then align rows
    stopifnot(all(colnames(mtx) %in% meta$Sample_name))
    meta <- meta[match(colnames(mtx), meta$Sample_name), ]
    stopifnot(all(meta$Sample_name == colnames(mtx)))

    # Build DGEList (group used for display / downstream consistency)
    dge <- DGEList(counts = mtx, group=meta$Source)

    # Run library-size clustering QC
    hc_res <- hc_filtered(
      dgelist_filtered = dge,
      metadata         = meta,
      nclust           = 8,
      Run              = "Sample_name",
      Sample           = "Source",
      filter           = filt,
      dataset          = ds
    )
    
    list(
      dge      = dge,
      hc       = hc_res,
      n_genes  = nrow(mtx),
      lib_size = dge$samples$lib.size
    )
  })
  
  setNames(filter_results, names(ds_filters))
})

names(all_results) <- names(datasets)

View(all_results)

# =============================================================================
# Step 2: Manually curated outlier sample IDs to remove
#   Samples identified as outliers from hc_filtered plots (Example: extreme library sizes isolated in their own cluster). Empty vectors = no outliers.
# =============================================================================
drop_samples <- list(
  BC = c("SRR11600338"),
  EBC1 = c(),
  EBC2 = c("Healthy7_Run2"),
  `HCC-PBMC` = c("C1", "C3", "C4"),
  `HCC-tissue` = c("SRR14027947", "SRR14027941")
)

# =============================================================================
# Function: filter_and_save_dataset
# =============================================================================
# Removes outlier samples from all filter-level DGEList objects for one dataset,
# then takes the intersection of retained samples across filter levels to ensure
# all output matrices have a consistent set of samples. Saves cleaned count
# matrices and metadata to disk.
#
# Parameters:
#   res     : Named list of filter-level results for one dataset
#             (output of the all_results[[ds]] lapply block).
#             Each element must contain $dge (DGEList) and $hc$matched_metadata.
#   ds_name : Dataset name (must match a key in drop_samples).
#   outdir  : Output directory path; created if it does not exist.
#
# Outputs written to outdir:
#   {ds_name}_metadata.csv          : Filtered metadata (retained samples only)
#   {ds_name}_{filter}_counts.csv   : Cleaned count matrix per filter level
#
# Returns:
#   Named list (one entry per filter) each containing $dge (cleaned DGEList)
# =============================================================================
filter_and_save_dataset <- function(res, ds_name, outdir) {
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Determine retained samples per filter ----
  final_samples_list <- lapply(names(res), function(filt) {
    dge <- res[[filt]]$dge
    
    setdiff(
      colnames(dge$counts),
      drop_samples[[ds_name]]
    )
  })
  
  # Take intersection so all filter-level matrices share the same sample set.
  # This prevents downstream tools from encountering mismatched sample sets.
  final_samples <- Reduce(intersect, final_samples_list)
  
 # ---- Save metadata once (same sample set applies to all filter levels) ----
  meta <- res[[1]]$hc$matched_metadata
  meta_filt <- meta[meta$Sample_name %in% final_samples, ]
  write.csv(
    meta_filt,
    file = file.path(outdir, paste0(ds_name, "_metadata.csv")),
    row.names = FALSE
  )
  
  # ---- Subset and save each filter-level count matrix ----
  filtered <- lapply(names(res), function(filt) {
    dge <- res[[filt]]$dge
    dge_filt <- dge[, final_samples, keep.lib.sizes = FALSE]
    
    write.csv(
      dge_filt$counts,
      file = file.path(outdir, paste0(ds_name, "_", filt, "_counts.csv"))
    )
    
    list(
      dge = dge_filt
    )
  })
  
  names(filtered) <- names(res)
  return(filtered)
}

# =============================================================================
# Step 3: Apply outlier removal and save cleaned matrices for all datasets
#   Output root: for_simulation/{dataset}/
# =============================================================================
filtered_results <- lapply(names(all_results), function(ds) {
  message("Filtering and saving dataset: ", ds)
  filter_and_save_dataset(
    res     = all_results[[ds]],
    ds_name = ds,
    outdir  = file.path("/media/meteor/FatDawg/Benchmark_Paper/HC_dataset_filtering/for_simulation", ds)
  )
})

names(filtered_results) <- names(all_results)
