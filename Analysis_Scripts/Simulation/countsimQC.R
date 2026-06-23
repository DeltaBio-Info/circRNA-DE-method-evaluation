# =============================================================================
# countsimQC Report Generation Pipeline
# =============================================================================
# Purpose:
#   Generates countsimQC quality-control reports that compare simulated RNA-seq
#   count matrices against their original (real) count matrix. For each
#   dataset across three different filtering strategies and DE fraction combinations, an HTML report and
#   accompanying plots are produced, allowing for visual and statistical assessment
#   of how faithfully SPsimSeq reproduced the distributional properties of the real data.
#
# Directory structure expected:
#  original/inputs/{dataset}/{dataset}_{filter}_counts.csv
#  simulated_inputs/{dataset}/{filter}_counts/{DE_folder}/*_counts.csv
#
# Outputs (per dataset, filter, and DE combination)
#  base_outdir/{dataset}/{filter}/{DE_folder}/countsimQC_report.html
#  base_outdir/{dataset}/{filter}/{DE_folder}/ generated plots

library(countsimQC)

#### Directory configuration ####
base_outdir      <- "/path/to/countsimQC_reports/BSJ_only"
original_inputs  <- "/path/to/HC_dataset_filtering/for_simulation/BSJ_only"
simulated_inputs <- "/path/to/Simulated_datasets/BSJ_only"


# =============================================================================
# Main: Iterate over datasets → filters → DE fractions
# =============================================================================
# Each subdirectory of original_inputs is treated as one dataset
datasets <- list.dirs(original_inputs, recursive = FALSE, full.names = FALSE)

for (dataset in datasets) {
  
  message("Processing dataset: ", dataset)
  
  orig_dir <- file.path(original_inputs, dataset)
  sim_base <- file.path(simulated_inputs, dataset)

  # Skip datasets that have no simulated output yet
  if (!dir.exists(sim_base)) next
  
  # ---------------------------------------------------------------------------
  # Part 1: Filtering strategies
  #   Subdirectories of sim_base follow the pattern {filter_name}_counts/
  #   (Example: autofilter_counts/, min1_counts/)
  # ---------------------------------------------------------------------------
  filter_dirs <- list.dirs(sim_base, recursive = FALSE, full.names = FALSE)
  
  for (filter_dir in filter_dirs) {
    print(filter_dir)

    # Strip trailing "_counts" suffix to get the clean filter label
    filter_name <- sub("_counts$", "", filter_dir)
    sim_filter_dir <- file.path(sim_base, filter_dir)
    
  # ---------------------------------------------------------------------------
  # Part 2: DE fraction folders
  #   Each filter directory contains subfolder for DE setting: DE_0/, DE_0.1/
  # ---------------------------------------------------------------------------
    de_dirs <- list.dirs(sim_filter_dir, recursive = FALSE, full.names = FALSE)
    
    for (de_dir in de_dirs) {
      
      message("  Filter: ", filter_name, " | ", de_dir)

      # Create output directory for this specific combination
      out_dir <- file.path(
        base_outdir,
        dataset,
        filter_name,
        de_dir
      )
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      # -----------------------------------------------------------------------
      # Load original (real) count matrix
      #   Expected filename: {dataset}_{filter_name}_counts.csv
      #   Rows = genes, columns = samples
      # -----------------------------------------------------------------------
      orig_file <- file.path(
        orig_dir,
        paste0(dataset, "_", filter_name, "_counts.csv")
      )
      
      if (!file.exists(orig_file)) {
        warning("Missing original file: ", orig_file)
        next
      }
      
      original_matrix <- as.matrix(
        read.csv(orig_file, row.names = 1, check.names = FALSE)
      )
      
      # -----------------------------------------------------------------------
      # Load all simulated count matrices for this filter and DE combination
      #   Each file is one SPsimSeq replicate (*_counts.csv)
      # -----------------------------------------------------------------------
      sim_dir <- file.path(sim_filter_dir, de_dir)
      
      sim_files <- list.files(
        sim_dir,
        pattern = "_counts\\.csv$",
        full.names = TRUE
      )
      
      if (length(sim_files) == 0) {
        warning("No simulated files in ", sim_dir)
        next
      }

      # Read each replicate into a named list of matrices
      count_list <- lapply(sim_files, function(f) {
        as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
      })

      # Name each list entry by its replicate identifier
      names(count_list) <- sub("_counts\\.csv$", "", basename(sim_files))

      # Append the real data as "original" so countsimQC can compare against it
      count_list[["original"]] <- original_matrix

      # -----------------------------------------------------------------------
      # Generate countsimQC report
      #
      # ddsList    : Named list of count matrices, the simulated replicates and original. countsimQC treats the entry named "original" as the reference for comparison.
      # outputFile : Fixed report filename written into outputDir
      # savePlots  : Export individual plots as files aloongside the HTML
      # dpi        : Set resolution for the saved plot images
      # description: Free-text label embdedded in the report header
      # -----------------------------------------------------------------------
      countsimQCReport(
        ddsList     = count_list,
        outputFile = "countsimQC_report.html",
        outputDir  = out_dir,
        savePlots  = TRUE,
        description = paste(
          "Dataset:", dataset,
          "| Filter:", filter_name,
          "|", de_dir
        ),
        dpi = 300
      )
    }
  }
}
