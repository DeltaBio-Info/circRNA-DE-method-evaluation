# =============================================================================
# RNA-seq Simulation Pipeline using SPsimSeq
# =============================================================================

# Purpose:
#  This script simulates the circRNA datasets utilized in the study from real observed count matrices using the SPsimSeq framework. For each input dataset and filtering strategy, two different differential expression (DE) scenarios were generated
#  DE0 (no DE-Signal) and DE10-Signal (10% of the circRNAs are DE).

# Inputs:
#   - Count matrices in CSV file format: genes x samples, one per filtering strategy per dataset
#   - Metadata CSV file format: sample-level information including group labels

# Outputs:
#  - counts.csv: Simulated count matrix
#  - TP.csv: True positive (DE genes) annotations (rowData)
#  - design.csv: Simulated sample design/group assignment (colData)

library(SPsimSeq)
library(edgeR)
library(SingleCellExperiment)
library(SimSeq)

# Set global seed for full reproducability
set.seed(42)

# =============================================================================
# Main Function: simulate_dataset
# =============================================================================
# Simulates n.sim count datasets from a single real count matrix using SPsimSeq and then saves results to disk.
# Skips if all output files exist.

# Parameters:
# count_matrix_path: Path to the real count matrix (gene x sample) in CSV file format
# metadata_path: Path to the sample metadata CSV (must contain Sample_name and Source columns)
# output_path: Directory where the output CSVs will be written
# DE_fraction: Proportion of genes to simulated as differetially expressed (for DE0 = 0 and for DE10-Signal = 0.1 (10% DE)
# dataset_name: Short identified for the source dataset (used in output filenames)
# filter_name: Label for the filtering strategy applied to counts (used in output filenames)
#
#
#
# Output files per replicate (naming: {dataset}_{filter}_DE{frac}_rep{i}_*):
#   counts.csv: Simulated count matrix
#   TP.csv: Gene-level metadata including DE status (true positives)
#   design.csv: Sample-level metadata with group assignments


simulate_dataset <- function(count_matrix_path,
                             metadata_path,
                             output_path,
                             DE_fraction,
                             dataset_name,
                             filter_name) {
  
  #### 1. Load input data ####
  raw_counts <- read.csv(count_matrix_path, sep = ",", header = TRUE, row.names = 1)
  study_design <- read.csv(metadata_path, sep = ",", header = TRUE, stringsAsFactors = TRUE)

  # Ensure metadata rows align with count matrix columns (important to ensure sample label consistency)
  stopifnot(all(colnames(raw_counts) %in% study_design$Sample_name))
  study_design <- study_design[match(colnames(raw_counts), study_design$Sample_name), ]
  
  #### 2. Build EdgeR DGEList ####
  dge <- DGEList(counts = raw_counts, group = study_design$Source)
  filtered_counts <- dge$counts
  
  #### 3. Set simulation scale based on dataset size ####
  #  - Small datasets (<12 samples): simulate 16 samples, 30 replicates
  #  - Larger datasets: simulate same number as original (even), 50 replicates
  
  n_orig <- ncol(filtered_counts)
  if (n_orig <= 12) {
    tot.samples <- 16    # small datasets get 16 samples (8 per group)
    n.sim <- 30
  } else {
    tot.samples <- 2 * floor(n_orig / 2) # size for bigger datasets, must be able to be divided by 2
    n.sim <- 50
  }

  #### 4. Skip early if all outputs already exist ####
  expected_files <- unlist(
    lapply(seq_len(n.sim), function(i) {
      prefix <- paste(
        dataset_name,
        filter_name,
        paste0("DE", DE_fraction),
        paste0("rep", i),
        sep = "_"
      )
      c(
        file.path(output_path, paste0(prefix, "_counts.csv")),
        file.path(output_path, paste0(prefix, "_TP.csv")),
        file.path(output_path, paste0(prefix, "_design.csv"))
      )
    })
  )
  
  if (all(file.exists(expected_files))) {
    message("  -> Outputs already exist, skipping simulation.")
    return(invisible(NULL))
  }

  #### 5. Validate and configure group structure (must be exactly 2 factors: Healthy and Cancer) ####
  grp <- factor(study_design$Source)
  stopifnot(nlevels(grp) == 2)
  group.config <- rep(0.5, 2)
  
  message("Simulating ", n.sim, " datasets for ", dataset_name, " (", filter_name,
          "), DE_fraction = ", DE_fraction, ", total_samples = ", tot.samples)
   
  #### 6. Run SPsimSeq simulation ####
  #
  # Parameters:
  #  n.genes:                Simulate a minimum of 1000 genes (or fewer if the input has less)
  #  pDE:                    Fraction of genes that are truly DE
  #  lfc.thrld:              Minimum |log2 fold-change| for a gene to be called DE (0.5)
  #  t.thrld:                Minimum t-statistics threshold for DE gene selection (2.5)
  #  genewiseCor:            Preserve gene-wise correlation structure from real data
  #  log.CPM.transform:      Use logCPM space for fitting the simulation model
  #  variable.lib.size:      Keep library sizes fixed across simulated samples
  sim_data <- SPsimSeq(
    n.sim = n.sim,
    s.data = filtered_counts,
    group = grp,
    group.config = group.config,
    n.genes = min(1000, nrow(filtered_counts)),
    tot.samples = tot.samples,
    pDE = DE_fraction,
    lfc.thrld = 0.5,
    t.thrld = 2.5,
    genewiseCor = TRUE,
    log.CPM.transform = TRUE,
    variable.lib.size = FALSE,
    result.format = "list",
    #return.details = TRUE,
    verbose = TRUE
  )
  sim_data
  cat("Generated", length(sim_data), "datasets.\n")

  #### 7. Save each simulated replicate to disk ####
  # File naming pattern: {dataset_name}_{filter_name}_DE{fraction}_rep{i}_{type}.csv
  # Skips individual replicates that are already saves
  for (i in seq_along(sim_data)) {
    sim <- sim_data[[i]]
    prefix <- paste(dataset_name, filter_name, paste0("DE", DE_fraction), paste0("rep", i), sep = "_")
    
    message("  Saving replicate ", i, "...")
    counts_file <- file.path(output_path, paste0(prefix, "_counts.csv"))
    tp_file     <- file.path(output_path, paste0(prefix, "_TP.csv"))
    design_file <- file.path(output_path, paste0(prefix, "_design.csv"))
    
    if (file.exists(counts_file) &&
        file.exists(tp_file) &&
        file.exists(design_file)) {
      message("  Replicate ", i, " already exists, skipping.")
      next
    }
    write.csv(sim$counts, counts_file)
    write.csv(sim$rowData, tp_file)
    write.csv(sim$colData, design_file, row.names = FALSE)
  }
  
  message("Done with ", dataset_name, " (", filter_name, "), DE_fraction = ", DE_fraction, "\n")
}

# =============================================================================
# Main: Iterate over all datasets, filtering strategies, and DE fractions
# =============================================================================
# --- Directory configuration ---
base_outdir <- "/path/to/simulated_datasets/BSJ_only"
input_folder <- "/path/to/hc_filtered_datasets/for_simulation/BSJ_only"

# DE fractions to simulate: 0% DE (null) and 10% DE (signal)
DE_settings <- c(0, 0.1)

# Loop over each dataset subdirectory
datasets <- list.dirs(input_folder, recursive = FALSE)

for (ds_path in datasets) {
  ds_name <- basename(ds_path)
  message("Dataset: ", ds_name)

  # Metadata file follows naming convention: {dataset_name}_metadata.csv
  metadata_path <- file.path(ds_path, paste0(ds_name, "_metadata.csv"))

  # Find all filtered count matrices in the dataset folder
  count_files <- list.files(ds_path, pattern = "_counts\\.csv$", full.names = TRUE)
  
  for (cnt in count_files) {
    # Extract filter label from filename: {dataset_name}_{filter_name}_counts.csv
    filter_name <- sub("\\.csv$", "", sub(paste0(ds_name, "_|_counts"), "", basename(cnt)))
    message("  Filter: ", filter_name)
    
    for (DEf in DE_settings) {
      # Create output directory: base_outdir/{dataset}/{filter}/DE_{fraction}/
      outdir <- file.path(base_outdir, ds_name, filter_name, paste0("DE_", DEf))
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      
      simulate_dataset(
        count_matrix_path = cnt,
        metadata_path = metadata_path,
        output_path = outdir,
        DE_fraction = DEf,
        dataset_name = ds_name,
        filter_name = filter_name
      )
    }
  }
}

