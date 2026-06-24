#### Pipeline to run circMeta functions circJuncDE and circCLRDE ####

# Authors: Erda Qorri


############################################################
## Import libraries
############################################################

library(circMeta)
library(dplyr)

# =============================================================================
# Function 1: clean_ciri3_files
# =============================================================================
# Cleans CIRI3 output files and converts them into a format compatible with circMeta.
#
# Parameters:
#     input_dir            :            Directory containing CIRI3 output files.
#     output_dir           :            Directory where cleaned files will be written.
#     pattern              :            File pattern used to identify CIRI3 files.
#     overwrite            :            Logical indicating whether existing output
#                                       files should be replaced.
#
#
# Returns:
#        cleaned_files     :            Character vector containing paths to all
#                                       cleaned CIRI3 files.
# =============================================================================

clean_ciri3_files <- function(input_dir,
                              output_dir,
                              pattern = "\\.ciri$",
                              overwrite = FALSE) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Search for input files matching the defined pattern
  files <- list.files(
    input_dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE
  )

  files <- files[!grepl("BSJ_Matrix|FSJ_Matrix", files)]

  if (length(files) == 0) {
    stop("No CIRI files found in: ", input_dir)
  }

  canonical_chr <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  canonical_no_chr <- c(as.character(1:22), "X", "Y", "M", "MT")

  cleaned_files <- character(length(files))

  for (i in seq_along(files)) {

    f <- files[i]

    message("Cleaning: ", basename(f))

    # Read the CIRI files
    df <- read.delim(
      f,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    required_cols <- c("chr", "circRNA_start", "circRNA_end")
    missing_cols <- setdiff(required_cols, names(df))

    if (length(missing_cols) > 0) {
      stop(
        "Missing required columns in file: ", basename(f),
        "\nMissing columns: ", paste(missing_cols, collapse = ", ")
      )
    }

    ### File Incompatibility 1 ###
    # remove Score column if present: CIRI3 output files have an additional Score column which is not compatible with circMeta
    df <- df[, names(df) != "Score", drop = FALSE]

    ### File Incompatibility 2 ###
    # fix chromosome naming
    df$chr <- as.character(df$chr)

    idx <- df$chr %in% canonical_no_chr
    df$chr[idx] <- paste0("chr", df$chr[idx])

    # due to different genome source change chrMT to chrM otherwise throws an error
    df$chr[df$chr == "chrMT"] <- "chrM"

    # rebuild circRNA_ID
    df$circRNA_ID <- paste0(
      df$chr,
      ":",
      df$circRNA_start,
      "|",
      df$circRNA_end
    )

    # keep only canonical chromosomes: scaffolds and contigs were removed as they trigger an error from circMeta
    df <- df[df$chr %in% canonical_chr, , drop = FALSE]

    out_file <- file.path(output_dir, basename(f))

    if (file.exists(out_file) && !overwrite) {
      stop(
        "Output file already exists: ", out_file,
        "\nUse overwrite = TRUE if you want to replace it."
      )
    }

    write.table(
      df,
      file = out_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )

    cleaned_files[i] <- out_file
  }

  message("Cleaned files written to: ", output_dir)

  return(cleaned_files)
}


# =============================================================================
# Function 2: run_circmeta_juncDE
# =============================================================================
# Runs circMeta::circJuncDE to identify differentially expressed circRNAs
# directly from CIRI output files.
#
# Parameters:
#     ciri_dir             :            Directory containing cleaned CIRI files.
#     metadata_path        :            Path to metadata file.
#     output_dir           :            Directory where results will be written.
#     dataset_name         :            Dataset label used for output naming.
#     sample_col           :            Metadata column containing sample IDs.
#     condition_col        :            Metadata column containing biological
#                                       conditions.
#     case_label           :            Label corresponding to the case group.
#     control_label        :            Label corresponding to the control group.
#     file_pattern         :            Pattern used to identify CIRI files.
#     id_source            :            Method used to extract sample IDs
#                                       ("filename" or "folder").
#     file_prefix          :            Prefix removed from filenames when
#                                       extracting sample IDs.
#     file_suffix          :            Suffix removed from filenames when
#                                       extracting sample IDs.
#     de_method            :            Differential expression method passed to
#                                       circJuncDE.
#     fdr_cutoff           :            FDR threshold used to identify significant
#                                       circRNAs.
#
#
# Returns:
#        result_granges    :            GRanges object returned by circJuncDE.
#        result_df         :            Complete differential expression results.
#        significant_df    :            Results passing the specified FDR cutoff.
#        sig_df_with_id    :            Significant results with circRNA_ID added.
#        design            :            Binary design vector used by circJuncDE.
#        ordering_check    :            Table verifying sample ordering and
#                                       condition assignment.
#        files             :            Vector of CIRI files included in the
#                                       analysis.
# =============================================================================

run_circmeta_juncDE <- function(ciri_dir,
                                metadata_path,
                                output_dir,
                                dataset_name,
                                sample_col = "SampleID",
                                condition_col = "Condition",
                                case_label = "EBC",
                                control_label = "Healthy",
                                file_pattern = "\\.ciri$",
                                id_source = c("filename", "folder"),
                                file_prefix = "CIRI3_",
                                file_suffix = ".ciri",
                                de_method = "pois.ztest",
                                fdr_cutoff = 0.05) {

  id_source <- match.arg(id_source)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ##########################################################
  ## Read metadata
  ##########################################################

  metadata <- read.csv(
    metadata_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!sample_col %in% names(metadata)) {
    stop("Sample column not found in metadata: ", sample_col)
  }

  if (!condition_col %in% names(metadata)) {
    stop("Condition column not found in metadata: ", condition_col)
  }

  ##########################################################
  ## List CIRI files
  ##########################################################

  ciri_files <- list.files(
    ciri_dir,
    pattern = file_pattern,
    full.names = TRUE,
    recursive = TRUE
  )

  ciri_files <- ciri_files[!grepl("BSJ_Matrix|FSJ_Matrix", ciri_files)]

  if (length(ciri_files) == 0) {
    stop("No CIRI files found in: ", ciri_dir)
  }

  ciri_files <- sort(ciri_files)

  ##########################################################
  ## Extract sample IDs
  ##########################################################

  if (id_source == "filename") {

    sample_ids <- basename(ciri_files)

    if (!is.null(file_prefix) && file_prefix != "") {
      sample_ids <- sub(paste0("^", file_prefix), "", sample_ids)
    }

    if (!is.null(file_suffix) && file_suffix != "") {
      sample_ids <- sub(paste0("\\", file_suffix, "$"), "", sample_ids)
    }

  } else if (id_source == "folder") {

    sample_ids <- basename(dirname(ciri_files))
  }

  ##########################################################
  ## Match metadata to file order
  ##########################################################

  metadata_ordered <- metadata[
    match(sample_ids, metadata[[sample_col]]),
    ,
    drop = FALSE
  ]

  missing_idx <- is.na(metadata_ordered[[sample_col]])

  if (any(missing_idx)) {
    missing_samples <- sample_ids[missing_idx]

    stop(
      "Sanity Check: Some CIRI files were not found in metadata:\n",
      paste(missing_samples, collapse = "\n")
    )
  }

  if (!all(sample_ids == metadata_ordered[[sample_col]])) {
    stop("Sanity Check: File order and metadata order do not match.")
  }

  ##########################################################
  ## Check condition labels
  ##########################################################

  valid_conditions <- c(case_label, control_label)

  unexpected_conditions <- setdiff(
    unique(metadata_ordered[[condition_col]]),
    valid_conditions
  )

  if (length(unexpected_conditions) > 0) {
    stop(
      "Unexpected condition labels found: ",
      paste(unexpected_conditions, collapse = ", "),
      "\nExpected only: ",
      paste(valid_conditions, collapse = ", ")
    )
  }

  ##########################################################
  ## Create design vector
  ##########################################################

  design <- ifelse(
    metadata_ordered[[condition_col]] == case_label,
    1,
    0
  )

  names(design) <- sample_ids

  ordering_check <- data.frame(
    file = basename(ciri_files),
    sample_id = sample_ids,
    condition = metadata_ordered[[condition_col]],
    design = design,
    stringsAsFactors = FALSE
  )

  message("\nDataset: ", dataset_name)
  message("DE method: ", de_method)
  message("Case label: ", case_label, " = 1")
  message("Control label: ", control_label, " = 0")
  message("\nDesign table:")
  print(table(ordering_check$condition, ordering_check$design))

  write.csv(
    ordering_check,
    file.path(output_dir, paste0(dataset_name, "_sample_ordering_check.csv")),
    row.names = FALSE
  )

  ##########################################################
  ## Run circMeta
  ##########################################################

  res <- circJuncDE(
    files = ciri_files,
    designs = design,
    circ.method = "CIRI",
    DE.method = de_method
  )

  res_df <- as.data.frame(res)

  if (!"fdr" %in% names(res_df)) {
    stop("The result does not contain an 'fdr' column.")
  }

  sig_df <- res_df %>%
    filter(fdr < fdr_cutoff)

  # Extract coordinates from results
  results_coords <- data.frame(
    chr = as.character(res_df$seqnames),
    start = res_df$start,
    end = res_df$end,
    stringsAsFactors = FALSE
  )

  message("Preview of result coordinates (first 6 rows):")
  message(paste(capture.output(head(results_coords)), collapse = "\n"))

  # Read cleaned CIRI files to get circRNA_IDs
  all_circrna_ids <- NULL
  for (f in ciri_files) {
    df <- read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    all_circrna_ids <- rbind(all_circrna_ids, df[, c("chr", "circRNA_start", "circRNA_end", "circRNA_ID")])
  }

  # Remove duplicates (keep unique circRNAs)
  all_circrna_ids <- all_circrna_ids[!duplicated(all_circrna_ids[, 1:3]), ]

  # Merge results with circRNA_IDs
  sig_df_with_id <- merge(
    sig_df,
    all_circrna_ids,
    by.x = c("seqnames", "start", "end"),
    by.y = c("chr", "circRNA_start", "circRNA_end"),
    all.x = TRUE
  )


  ##########################################################
  ## Save outputs
  ##########################################################

  all_out <- file.path(
    output_dir,
    paste0(dataset_name, "_circMeta_", de_method, "_all_results.csv")
  )

  sig_out <- file.path(
    output_dir,
    paste0(dataset_name, "_circMeta_", de_method, "_FDR", fdr_cutoff, "_significant.csv")
  )

  sig_out_with_id <- file.path(
    output_dir,
    paste0(dataset_name, "_circMeta_", de_method, "_FDR", fdr_cutoff, "_significant_with_circRNA_ID.csv")
  )

  write.csv(res_df, all_out, row.names = TRUE)
  write.csv(sig_df, sig_out, row.names = TRUE)
  write.csv(sig_df_with_id, sig_out_with_id, row.names = TRUE)

  message("\nFinished circJuncDE.")
  message("All results: ", all_out)
  message("Significant results: ", sig_out)
  message("Significant results with circRNA_ID: ", sig_out_with_id)
  message("Number significant at FDR < ", fdr_cutoff, ": ", nrow(sig_df))

  return(list(
    result_granges = res,
    result_df = res_df,
    significant_df = sig_df,
    sig_df_with_id = sig_df_with_id,
    design = design,
    ordering_check = ordering_check,
    files = ciri_files
  ))
}


# =============================================================================
# Function 3: run_circmeta_CLRDE
# =============================================================================
# Runs circMeta::circCLRDE on the output generated by circJuncDE and combines
# CLR-based differential expression statistics with the original circRNA
# annotation.
#
# Parameters:
#     junc_res             :            GRanges object returned by circJuncDE.
#     output_dir           :            Directory where results will be written.
#     dataset_name         :            Dataset label used for output naming.
#     de_method            :            Differential expression method used by
#                                       circCLRDE ("wald", "lr", "fisher",
#                                       or "chisq").
#     is_shrink            :            Logical indicating whether shrinkage
#                                       estimation should be applied.
#     is_equalrho          :            Logical indicating whether equal rho
#                                       should be assumed.
#     is_peudo             :            Logical indicating whether pseudocounts
#                                       should be added.
#     fdr_cutoff           :            FDR threshold used to identify significant
#                                       circRNAs.
#
#
# Returns:
#        result                    :   Complete merged results containing
#                                       circRNA annotation, junction-level
#                                       statistics and CLR statistics.
#        significant_df            :   Results passing the specified FDR cutoff.
#        clr_input                 :   CLR input matrices generated by getCLR
#                                       (x1, n1, x2, n2).
#        method                    :   Differential expression method used.
#        is_shrink                 :   Shrinkage setting used during analysis.
#        is_equalrho               :   Equal rho setting used during analysis.
#        is_peudo                  :   Pseudocount setting used during analysis.
#        runtime                   :   Execution time of circCLRDE.
#        all_results_path          :   Path to the complete results file.
#        significant_results_path  :   Path to the significant results file.
# =============================================================================

#### NOTE: is.equalrho and is.peudo are kept at default

run_circmeta_CLRDE <- function(junc_res,
                               output_dir,
                               dataset_name,
                               de_method = "wald",
                               is_shrink = TRUE,
                               is_equalrho = FALSE,
                               is_peudo = TRUE,
                               fdr_cutoff = 0.05) {

  de_method <- match.arg(de_method, choices = c("wald", "lr", "fisher", "chisq"))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ##########################################################
  ## Extract CLR matrices
  ##########################################################

  tmp <- getCLR(junc_res)

  # Validate dimensions
  stopifnot(all(dim(tmp$x1) == dim(tmp$n1)))
  stopifnot(all(dim(tmp$x2) == dim(tmp$n2)))

  # Validate row name alignment within groups
  stopifnot(all(rownames(tmp$x1) == rownames(tmp$n1)))
  stopifnot(all(rownames(tmp$x2) == rownames(tmp$n2)))

  # Validate row name alignment across groups
  stopifnot(all(rownames(tmp$x1) == rownames(tmp$x2)))

  ##########################################################
  ## Extract annotation from GRanges - BEFORE running CLR
  ## This is the key step to preserve circRNA identity
  ## since circCLRDE works on matrices and loses annotation
  ##########################################################

  res_junc_ann <- as.data.frame(junc_res) %>%
    dplyr::rename(
      lfc_junc  = lfc,
      pval_junc = pval,
      fdr_junc  = fdr
    )

  # Critical check: row order must match between annotation and CLR matrices
  if (!all(rownames(tmp$x1) == rownames(res_junc_ann))) {
    stop("Row order mismatch between CLR matrices and GRanges annotation. Cannot safely merge.")
  }

  message("Row order validated: CLR matrices and GRanges annotation are aligned.")

  ##########################################################
  ## Run circCLRDE
  ##########################################################

  message("Running circCLRDE...")
  message("  DE method:   ", de_method)
  message("  Shrinkage:   ", is_shrink)
  message("  Equal rho:   ", is_equalrho)
  message("  Pseudocount: ", is_peudo)

  runtime <- system.time({
    res <- circCLRDE(
      x1         = tmp$x1,
      n1         = tmp$n1,
      x2         = tmp$x2,
      n2         = tmp$n2,
      DE.method  = de_method,
      is.shrink  = is_shrink,
      is.equalrho = is_equalrho,
      is.peudo   = is_peudo
    )
  })

  ##########################################################
  ## Prepare CLR results dataframe
  ##########################################################

  res_clr_df <- as.data.frame(res)

  # Handle FDR column - fisher and chisq return pval only
  if (!"fdr" %in% names(res_clr_df)) {
    if (de_method %in% c("fisher", "chisq")) {
      warning("fisher/chisq methods do not return FDR. Applying BH correction to p-values.")
      res_clr_df$fdr <- p.adjust(res_clr_df$pval, method = "BH")
    } else {
      stop("No 'fdr' column found in circCLRDE output.")
    }
  }

  # Rename to avoid column collisions when merging with junction results
  res_clr_df <- res_clr_df %>%
    dplyr::rename(
      stat_clr = stat,
      pval_clr = pval,
      fdr_clr  = fdr
    )

  ##########################################################
  ## Merge annotation + CLR results
  ## Row order is preserved so cbind is safe
  ##########################################################

  final <- cbind(res_junc_ann, res_clr_df)

  # Add circRNA_ID as explicit identifier
  final$circRNA_ID <- paste0(
    final$seqnames, ":",
    final$start, "|",
    final$end
  )

  # Reorder columns so circRNA_ID is first
  final <- final %>%
    dplyr::select(circRNA_ID, everything())

  ##########################################################
  ## Filter significant results
  ##########################################################

  sig_df <- final[final$fdr_clr < fdr_cutoff, , drop = FALSE]

  ##########################################################
  ## Output file naming
  ##########################################################

  method_tag <- paste0(
    de_method,
    "_shrink_",    is_shrink,
    "_equalrho_",  is_equalrho,
    "_peudo_",     is_peudo
  )

  all_out <- file.path(
    output_dir,
    paste0(dataset_name, "_circMeta_circCLRDE_", method_tag, "_all_results.csv")
  )

  sig_out <- file.path(
    output_dir,
    paste0(dataset_name, "_circMeta_circCLRDE_", method_tag, "_FDR", fdr_cutoff, "_significant.csv")
  )

  write.csv(final,   all_out, row.names = TRUE)
  write.csv(sig_df,  sig_out, row.names = TRUE)

  ##########################################################
  ## Summary messages
  ##########################################################

  message("\nFinished circCLRDE")
  message("Method:          ", de_method)
  message("Shrinkage:       ", is_shrink)
  message("Equal rho:       ", is_equalrho)
  message("Pseudocount:     ", is_peudo)
  message("Total circRNAs:  ", nrow(final))
  message("Significant at FDR < ", fdr_cutoff, ": ", nrow(sig_df))
  message("All results:     ", all_out)
  message("Significant:     ", sig_out)

  ##########################################################
  ## Return
  ##########################################################

  return(list(
    result          = final,        # full merged results with annotation
    significant_df  = sig_df,       # FDR filtered
    clr_input       = tmp,          # raw x1, n1, x2, n2 matrices
    method          = de_method,
    is_shrink       = is_shrink,
    is_equalrho     = is_equalrho,
    is_peudo        = is_peudo,
    runtime         = runtime,
    all_results_path       = all_out,
    significant_results_path = sig_out
  ))
}
