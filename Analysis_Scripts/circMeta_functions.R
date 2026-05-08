############################################################
## circMeta circJuncDE reusable pipeline
############################################################

library(circMeta)
library(dplyr)

############################################################ NOTE: I validated that the circRNA ID backtrack method I used works and the circRNA IDs are fine ############################################################

############################################################
## 1. Clean CIRI3 files
############################################################

# input_dir: directory containing CIRI files to clean
# output_dir: directory where cleaned files will be written
# pattern: regex pattern to identify files (defaults to files ending in .ciri)
# overwrite: boolean flag (defaults to FALSE) to prevent accident overwriting of the original files


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


############################################################
## 2. Run circJuncDE
############################################################

# ciri_dir: directory containing cleaned CIRI files
# metadata_path: csv file with sample metadata
# output_dir: directory where to save the results (does not need to exist, the script will create it)
# dataset_name: name for the dataset being analyzed (will be used to name the output files)
# sample_col: column name in the metadata identifying samples (default: Sample_ID)
# condition_col: column name for experimental conditons (default: Condition)
# case_label: label for case or disease group (Default: EBC)
# control_label: label for control group (default: healthy)
# file_pattern: regex to identify CIRI files (default: .ciri)
# id_source: how to extract sample IDs, from filename or folder name (default: both options available)
# file_prefix: prefix to remove from filenames (default: "CIRI3_")
# file_suffix: suffix to remove from filenames (default: ".ciri")

######## important params ########
# de_method: statistical method for DE analysis (default: pois.ztest" - Poisson Z-test)
# fdr_cutoff: FDR threshold significance (default: 0.05)


runtime_log_path <- "/media/alexandria-kouri/Data_8TB/Manuscript_Revisions/circRNA_DE_revision/circMeta/CIRI3_outs/runtime_log_BSJ_Junc_only.csv"

log_runtime <- function(dataset, method, fdr_cutoff, runtime, log_path) {
  entry <- data.frame(
    Dataset = dataset,
    DE_Method = method,
    FDR_Cutoff = fdr_cutoff,
    Runtime_seconds = as.numeric(runtime["elapsed"]),
    stringsAsFactors = FALSE
  )

  if (file.exists(log_path)) {
    existing <- read.csv(log_path, stringsAsFactors = FALSE)
    entry <- rbind(existing, entry)
  }

  write.csv(entry, log_path, row.names = FALSE)
  message("Runtime logged: ", dataset, " | ", method, " | ", round(as.numeric(runtime["elapsed"]), 2), "s")
}

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

### Returns list ###
# result_granges: Raw GRanges object from circJuncDE (has full annotation information)
# result_df: results as dataframe
# significant_df: significant results only (pass the FDR filer)
# design: the design vector used
# ordering_check: sample order verification table
# files: CIRI files analyzed

############################################################
## 3. Example run for your EBC dataset
############################################################

cleaned_ebc_files <- clean_ciri3_files(
  input_dir = "/media/alexandria-kouri/Data_8TB/Manuscript_Revisions/circRNA_DE_revision/circMeta/CIRI3_outs/BC/",
  output_dir = "/media/alexandria-kouri/Data_8TB/Manuscript_Revisions/circRNA_DE_revision/circMeta/CIRI3_outs/BC/BC_circMeta_cleaned",
  pattern = "\\.ciri$",
  overwrite = TRUE
)

runtime <- system.time({
ebc_res <- run_circmeta_juncDE(
  ciri_dir = "/media/alexandria-kouri/Data_8TB/Manuscript_Revisions/circRNA_DE_revision/circMeta/CIRI3_outs/BC/BC_circMeta_cleaned",
  metadata_path = "/media/alexandria-kouri/Data_8TB/Manuscript_Revisions/circRNA_DE_revision/circMeta/CIRI3_outs/BC/BC_metadata.csv",
  output_dir = "/media/alexandria-kouri/Data_8TB/Manuscript_Revisions/circRNA_DE_revision/circMeta/CIRI3_outs/BC/BC_circMETA_DE_outs",
  dataset_name = "BC",
  sample_col = "Sample_ID",
  condition_col = "Condition",
  case_label = "BC",
  control_label = "Healthy",
  id_source = "filename",
  file_prefix = "CIRI3_",
  file_suffix = ".ciri",
  de_method = "pois.ztest", # pois.ztest
  fdr_cutoff = 0.05
)
})

log_runtime("BC", "pois.ztest", 0.05, runtime, runtime_log_path)

ebc_res$result_granges
ebc_res$files
ebc_res$ordering_check

ebc_res$sig_df_with_id %>% View()





######################### Run circmetaCLRDE #########################
#### NOTE: is.equalrho and is.peudo are kept at default

run_circmeta_CLRDE <- function(junc_res,
                               output_dir,
                               dataset_name,
                               de_method = "wald",
                               is_shrink = TRUE,
                               is_equalrho = FALSE,
                               is_peudo = TRUE,
                               fdr_cutoff = 0.05) {

  # there are different options for de_method offered by circMeta
  de_method <- match.arg(de_method, choices = c("wald", "lr", "fisher", "chisq"))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # CLR extraction as given in the circMeta repository
  tmp <- getCLR(junc_res)

  # Sanity checks that both junction and total reads have the same dim and rownames (circRNA identifers) are aligned between them
  stopifnot(all(dim(tmp$x1) == dim(tmp$n1)))
  stopifnot(all(dim(tmp$x2) == dim(tmp$n2)))
  stopifnot(all(rownames(tmp$x1) == rownames(tmp$n1)))
  stopifnot(all(rownames(tmp$x2) == rownames(tmp$n2)))
  stopifnot(all(rownames(tmp$x1) == rownames(tmp$x2)))

  runtime <- system.time({
    res <- circCLRDE(
      x1 = tmp$x1,
      n1 = tmp$n1,
      x2 = tmp$x2,
      n2 = tmp$n2,
      DE.method = de_method,
      is.shrink = is_shrink,
      is.equalrho = is_equalrho,
      is.peudo = is_peudo
    )
  })

  res_df <- as.data.frame(res)

  if (!"fdr" %in% names(res_df)) {
    stop("No 'fdr' column found in circCLRDE output.")
  }

  # these two methods do not return fdr so the chunk of code above will throw an error
  if(de_method %in% c("fisher", "chisq") && !"fdr" %in% names(res_df)) {
    warning("fisher/chisq methods do not return FDR. Using raw p-value instead.")
    sig_df <- res_df[res_df$pval < fdr_cutoff, , drop = FALSE]
  }

  sig_df <- res_df[res_df$fdr < fdr_cutoff, , drop = FALSE]

  method_tag <- paste0(de_method, "_shrink_", is_shrink)

  all_out <- file.path(
    output_dir,
    paste0(dataset_name, "_circMeta_circCLRDE_", method_tag, "_all_results.csv")
  )

  sig_out <- file.path(
    output_dir,
    paste0(dataset_name, "_circMeta_circCLRDE_", method_tag, "_FDR", fdr_cutoff, "_significant.csv")
  )

  write.csv(res_df, all_out, row.names = TRUE)
  write.csv(sig_df, sig_out, row.names = TRUE)

  message("Finished circCLRDE")
  message("Method: ", de_method)
  message("Shrinkage: ", is_shrink)
  message("Rho: ", is_equalrho)
  message("Pseudo Count Added: ", is_peudo)
  message("Total circRNAs: ", nrow(res_df))
  message("Significant at FDR < ", fdr_cutoff, ": ", nrow(sig_df))

  return(list(
    result = res,
    result_df = res_df,
    significant_df = sig_df,
    clr_input = tmp,
    method = de_method,
    is_shrink = is_shrink,
    runtime = runtime,
    all_results_path = all_out,
    significant_results_path = sig_out
  ))
}

bc_clr_wald <- run_circmeta_CLRDE(
  junc_res = ebc_res$result_granges,
  output_dir = "/media/alexandria-kouri/Data_8TB/Manuscript_Revisions/circRNA_DE_revision/circMeta/CIRI3_outs/BC/BC_circMETA_DE_outs",
  dataset_name = "BC",
  de_method = "wald",
  is_shrink = TRUE,
  is_equalrho = FALSE,
  is_peudo = TRUE,
  fdr_cutoff = 0.05
)
bc_clr_wald$

bc_clr_wald$result %>% View()



################### Separate circRNA ID backtracking function ###################
### circMeta loses the circRNA ID information but it can be backtracked to the original results file ###
### res_junc (output of circJuncDE) -> getCLR() -> circCLRDE() (works only on matrices and not GRanges)
### the circRNA ID annotation information exists in the res_junc GRanges object, the main output of circJuncDE
### Important that the circRNA ID mapping happens at the non-FDR filtered matrix

# Step 1: Extract annotation information from res_junc
res_junc_ann <- as.data.frame(ebc_res$result_granges)
View(res_junc_ann)
head(res_junc_ann)

res_junc_ann <- res_junc_ann %>% rename(lfc = "lfc_junc",
                                    pval = "pval_junc",
                                    fdr = "fdr_junc")

# Step 2: Ensure that the row order of the x1 matches that of res_junc_ann
all(rownames(tmp$x1) == rownames(res_junc_ann)) # TRUE

# Step 3: Attach to the circCLRDE results
res_clr_df <- as.data.frame(bc_clr_wald$result)
res_clr_df <- res_clr_df %>% rename(stat = "stat_clr",
                                    pval = "pval_clr",
                                    fdr = "fdr_clr")

# Step 4: Merge both dataframes together
final_merged <- cbind(res_junc_ann, res_clr_df)

final_merged$circRNA_ID <- paste0(
  final_merged$seqnames, ":",
  final_merged$start, "|",
  final_merged$end
)

# Step 5: Validate that the rownames match
final_merged <- res_junc_ann[match(rownames(tmp$x1), rownames(res_junc_ann)), ]
final <- cbind(final_merged, res_clr_df)


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

bc_clr_wald <- run_circmeta_CLRDE(
  junc_res = ebc_res$result_granges,
  output_dir = "/media/alexandria-kouri/Data_8TB/Manuscript_Revisions/circRNA_DE_revision/circMeta/CIRI3_outs/BC/BC_circMETA_DE_outs",
  dataset_name = "BC",
  de_method = "wald",
  is_shrink = TRUE,
  is_equalrho = FALSE,
  is_peudo = TRUE,
  fdr_cutoff = 0.05
)

bc_clr_wald$significant_df %>% View()



ebc_res$result_granges %>% head()
bc_clr_wald$clr_input$x1 %>% head()
bc_clr_wald$clr_input$x2 %>% head()

bc_clr_wald$clr_input$x1 -> x1
bc_clr_wald$clr_input$x2 -> x2
cbind(x1, x2) -> x12
View(x12)

head(x2)
bc_clr_wald$clr_input$n1 -> n1
bc_clr_wald$clr_input$n2 -> n2
cbind(n1, n2)  -> n12
View(bc_clr_wald$clr_input)

bc_clr_wald$significant_df %>% head()
bc_clr_wald$result %>% View()

######################### Testing #########################

ebc_res$result_granges -> res

head(res)

ciri_dir <- "/media/alexandria-kouri/Data_8TB/Manuscript_Revisions/circRNA_DE_revision/EBC1/EBC1_clean/"  # Change this to your actual path

# Get list of CIRI files
ciri_files <- list.files(
  ciri_dir,
  pattern = "\\.ciri$",
  full.names = TRUE,
  recursive = TRUE
)

ciri_files <- ciri_files[!grepl("BSJ_Matrix|FSJ_Matrix", ciri_files)]

head(ciri_files)  # Check which files were found
length(ciri_files)  # How many files?

# Step 4: Read files and collect circRNA_IDs
all_circrna_ids <- NULL

for (f in ciri_files) {
  message("Reading: ", basename(f))  # Show progress

  df <- read.delim(
    f,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  head(df)  # Look at what you read

  # Extract just these 4 columns
  subset_df <- df[, c("chr", "circRNA_start", "circRNA_end", "circRNA_ID")]

  head(subset_df)  # Check what was extracted

  # Append to growing data frame
  all_circrna_ids <- rbind(all_circrna_ids, subset_df)
}

head(all_circrna_ids)
dim(all_circrna_ids)  # How many total circRNAs collected?

# Step 5: Remove duplicates
# If the same circRNA appears in multiple samples, keep only one

all_circrna_ids_unique <- all_circrna_ids[!duplicated(all_circrna_ids[, 1:3]), ]

dim(all_circrna_ids)          # Original count
dim(all_circrna_ids_unique)   # After removing duplicates

# Step 6: Merge your results with circRNA_IDs
# This is the KEY step - it matches coordinates to find circRNA IDs

# First look at what you're merging FROM
head(res[, c("seqnames", "start", "end")])

# Look at what you're merging TO
head(all_circrna_ids_unique)

# Now do the merge
res_with_id <- merge(
  res,
  all_circrna_ids_unique,
  by.x = c("seqnames", "start", "end"),  # Column names in res
  by.y = c("chr", "circRNA_start", "circRNA_end"),  # Column names in all_circrna_ids_unique
  all.x = TRUE  # Keep all rows from res, even if no match
)

head(res_with_id)
dim(res_with_id)  # Should still have same number of rows as res

# Check for any NAs in circRNA_ID (failed matches)
sum(is.na(res_with_id$circRNA_ID))  # How many didn't match?

# Step 7: Check the results
head(res_with_id)

# Look specifically at circRNA_IDs
head(res_with_id$circRNA_ID)

# Filter to significant ones if you want
fdr_cutoff <- 0.05
sig_res_with_id <- res_with_id[res_with_id$fdr < fdr_cutoff, ]

head(sig_res_with_id)
dim(sig_res_with_id)

circClass(luad_files[1:2],circ.method=c('CIRI'),gene=NULL,gexon=NULL)
