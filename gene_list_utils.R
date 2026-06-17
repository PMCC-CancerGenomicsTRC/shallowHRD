##### GENE LIST UTILITIES #####
# This file provides functions to load, validate, and process custom gene lists
# for copy-number analysis in shallowHRD

#' Load and validate a custom gene list from file
#'
#' @param file_path Path to the gene list file (TSV or CSV)
#' @param verbose Print validation messages
#'
#' @return A data frame with columns: gene, chr, position
#'
load_gene_list <- function(file_path, verbose = TRUE) {
  if (!file.exists(file_path)) {
    stop("Gene list file not found: ", file_path)
  }
  
  # Try to detect delimiter
  first_line <- readLines(file_path, n = 1)
  if (grepl("\t", first_line)) {
    gene_df <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
  } else if (grepl(",", first_line)) {
    gene_df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  } else {
    stop("Could not detect delimiter (tab or comma) in gene list file")
  }
  
  if (verbose) {
    cat("Loaded gene list with", nrow(gene_df), "entries\n")
    cat("Columns:", paste(colnames(gene_df), collapse = ", "), "\n")
  }
  
  return(gene_df)
}


#' Validate and normalize a gene list
#'
#' @param gene_df Data frame with gene annotations
#' @param verbose Print validation details
#'
#' @return Validated and normalized data frame with columns: gene, chr, position
#'
validate_and_normalize_gene_list <- function(gene_df, verbose = TRUE) {
  # Check required columns
  has_gene <- "gene" %in% colnames(gene_df)
  has_chr <- "chr" %in% colnames(gene_df)
  has_position <- "position" %in% colnames(gene_df)
  has_start_end <- ("start" %in% colnames(gene_df)) && ("end" %in% colnames(gene_df))
  
  if (!has_gene) {
    stop("Gene list must have a 'gene' column")
  }
  if (!has_chr) {
    stop("Gene list must have a 'chr' column")
  }
  if (!has_position && !has_start_end) {
    stop("Gene list must have either 'position' column or both 'start' and 'end' columns")
  }
  
  # If start/end provided, derive position as midpoint
  if (has_start_end && !has_position) {
    if (verbose) {
      cat("Deriving representative position as midpoint of start/end\n")
    }
    gene_df$position <- as.integer((as.numeric(gene_df$start) + as.numeric(gene_df$end)) / 2)
  }
  
  # Normalize chromosome format (remove 'chr' prefix if present, convert to numeric)
  gene_df$chr <- as.character(gene_df$chr)
  gene_df$chr <- sub("^chr", "", gene_df$chr, ignore.case = TRUE)
  gene_df$chr <- as.numeric(gene_df$chr)
  
  # Validate chromosome range for hg19 (1-22, X=23, Y=24)
  valid_chrs <- c(1:22, 23, 24)
  invalid_chrs <- gene_df$chr[!gene_df$chr %in% valid_chrs]
  if (length(invalid_chrs) > 0) {
    stop("Invalid chromosome values: ", paste(unique(invalid_chrs), collapse = ", "),
         "\nValid values for hg19: 1-22, 23 (X), 24 (Y)")
  }
  
  # Validate position coordinates
  gene_df$position <- as.numeric(gene_df$position)
  if (any(is.na(gene_df$position))) {
    stop("Some gene positions could not be converted to numeric values")
  }
  if (any(gene_df$position < 1)) {
    stop("Gene positions must be >= 1")
  }
  
  # Check for duplicates
  dup_genes <- gene_df$gene[duplicated(gene_df$gene)]
  if (length(dup_genes) > 0) {
    warning("Duplicate gene entries found: ", paste(unique(dup_genes), collapse = ", "))
  }
  
  # Keep only required columns
  gene_df <- gene_df[, c("gene", "chr", "position")]
  
  if (verbose) {
    cat("Validation passed for", nrow(gene_df), "genes\n")
    cat("Chromosome range: ", min(gene_df$chr), "-", max(gene_df$chr), "\n")
    cat("Position range: ", min(gene_df$position), "-", max(gene_df$position), "\n")
  }
  
  return(gene_df)
}


#' Extract CN and ratio metrics for a single gene
#'
#' This function mimics the repeated per-gene code in the original script.
#' It extracts both initial (per-window) and final (segmented) metrics for a gene.
#'
#' @param gene_name Name of the gene
#' @param chr Chromosome number
#' @param position Representative coordinate for the gene
#' @param B_initial Initial ratio file (all windows)
#' @param C_final Final segmented ratio file
#' @param THR Threshold for CN calculation
#'
#' @return A list with elements: gene, chr, start, ratio_point_initial,
#'         CN_to_baseline_point_initial, ratio_segment_initial,
#'         CN_to_baseline_segment_initial, ratio_segment_final,
#'         CN_to_baseline_segment_final
#'
extract_gene_metrics <- function(gene_name, chr, position, B_initial, C_final, THR) {
  # Find the closest window in the initial ratio file for this gene's chromosome
  B_chr <- B_initial[B_initial$chr == chr, ]
  
  if (nrow(B_chr) == 0) {
    warning("No windows found for ", gene_name, " on chr", chr)
    return(NULL)
  }
  
  # Find closest window using Closest function (or equivalent)
  closest_idx <- Closest(B_chr$start, position)[1]
  
  if (is.na(closest_idx)) {
    warning("Could not find closest window for ", gene_name)
    return(NULL)
  }
  
  higlight_initial <- B_chr[closest_idx, ]
  
  # Find closest segment in final segmentation
  C_chr <- C_final[C_final$chr == chr, ]
  
  if (nrow(C_chr) == 0) {
    warning("No segments found for ", gene_name, " on chr", chr)
    return(NULL)
  }
  
  closest_final_idx <- Closest(C_chr$start, position)[1]
  
  if (is.na(closest_final_idx)) {
    warning("Could not find closest segment for ", gene_name)
    return(NULL)
  }
  
  higlight_final <- C_chr[closest_final_idx, ]
  
  # Calculate CN metrics
  CN_baseline_initial_point <- round(higlight_initial$ratio / THR, 3)
  CN_baseline_initial_segment <- round(higlight_initial$ratio_median / THR, 3)
  CN_baseline_final_segment <- round(higlight_final$ratio_median / THR, 3)
  
  # Return as list (will be converted to data frame row later)
  result <- list(
    gene = gene_name,
    chr = chr,
    start = higlight_initial$start,
    ratio_point_initial = higlight_initial$ratio,
    CN_to_baseline_point_initial = CN_baseline_initial_point,
    ratio_segment_initial = higlight_initial$ratio_median,
    CN_to_baseline_segment_initial = CN_baseline_initial_segment,
    ratio_segment_final = higlight_final$ratio_median,
    CN_to_baseline_segment_final = CN_baseline_final_segment
  )
  
  return(result)
}


#' Build gene summary table from custom gene list
#'
#' @param gene_list Data frame with columns: gene, chr, position
#' @param B_initial Initial ratio file (all windows)
#' @param C_final Final segmented ratio file
#' @param THR Threshold for CN calculation
#' @param verbose Print progress
#'
#' @return Data frame with one row per gene and amplification/deletion metrics
#'
build_gene_summary_table <- function(gene_list, B_initial, C_final, THR, verbose = TRUE) {
  if (verbose) {
    cat("Building gene summary table for", nrow(gene_list), "genes...\n")
  }
  
  results_list <- list()
  
  for (i in seq_len(nrow(gene_list))) {
    gene_row <- gene_list[i, ]
    
    if (verbose && i %% 10 == 0) {
      cat("  Processing gene", i, "of", nrow(gene_list), "\n")
    }
    
    metrics <- extract_gene_metrics(
      gene_name = gene_row$gene,
      chr = gene_row$chr,
      position = gene_row$position,
      B_initial = B_initial,
      C_final = C_final,
      THR = THR
    )
    
    if (!is.null(metrics)) {
      results_list[[i]] <- metrics
    }
  }
  
  # Convert list of lists to data frame
  results_df <- do.call(rbind, lapply(results_list, as.data.frame, stringsAsFactors = FALSE))
  
  if (is.null(results_df) || nrow(results_df) == 0) {
    stop("No gene metrics could be extracted")
  }
  
  # Ensure column order
  col_order <- c("gene", "chr", "start", "ratio_point_initial",
                 "CN_to_baseline_point_initial", "ratio_segment_initial",
                 "CN_to_baseline_segment_initial", "ratio_segment_final",
                 "CN_to_baseline_segment_final")
  results_df <- results_df[, col_order]
  
  # Convert numeric columns to appropriate types
  results_df$chr <- as.integer(results_df$chr)
  results_df$start <- as.integer(results_df$start)
  results_df$ratio_point_initial <- as.numeric(results_df$ratio_point_initial)
  results_df$CN_to_baseline_point_initial <- as.numeric(results_df$CN_to_baseline_point_initial)
  results_df$ratio_segment_initial <- as.numeric(results_df$ratio_segment_initial)
  results_df$CN_to_baseline_segment_initial <- as.numeric(results_df$CN_to_baseline_segment_initial)
  results_df$ratio_segment_final <- as.numeric(results_df$ratio_segment_final)
  results_df$CN_to_baseline_segment_final <- as.numeric(results_df$CN_to_baseline_segment_final)
  
  if (verbose) {
    cat("Successfully extracted metrics for", nrow(results_df), "genes\n")
  }
  
  return(results_df)
}


#' Get default (built-in) gene list for hg19
#'
#' This reproduces the hard-coded gene list from the original script.
#' Returns a data frame suitable for validate_and_normalize_gene_list().
#'
#' @return Data frame with columns: gene, chr, position
#'
get_default_gene_list <- function() {
  gene_data <- data.frame(
    gene = c("DOCK7", "HCN3", "KLHL12", "RBBP5", "TSC22D2", "PIK3CA", "ANKRD17", "CDKN2AIP",
             "RTN4IP1", "AIM1", "INTS10", "PPP2R2A", "BRF2", "ZNF703", "MYC", "CD274_PDL1",
             "MTAP", "SEPHS1", "ZMIZ1", "WAPAL", "PTEN", "CSTF3", "BBS1", "CTTN", "CCND1",
             "ATG16L2", "INTS4", "CCDC77", "FOXM1", "YEATS4", "MDM2", "BRCA2", "C13orf23",
             "TRIM13", "SDCCAG1", "SNAP23", "IGF1R", "CYB5B", "P53", "ELAC2", "MAP2K4",
             "ERAL1", "NF1", "ERBB2", "BRCA1", "PHB", "SUPT4H1", "RAD51C", "GALK1", "AKAP8",
             "BRD4", "PIK3R2", "CCNE1", "NOSIP", "C20orf111", "ZNF217", "TSHZ2", "SAMD10", "PCNT"),
    chr = c(1, 1, 1, 1, 3, 3, 4, 4, 6, 6, 8, 8, 8, 8, 8, 9, 9, 10, 10, 10, 10, 11, 11, 11,
            11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 14, 15, 15, 16, 17, 17, 17, 17, 17, 17,
            17, 17, 17, 17, 17, 19, 19, 19, 19, 19, 20, 20, 20, 20, 21),
    position = c(63037227, 155253447, 202878988, 205073188, 150155147, 178912013, 74031804,
                 184368003, 107048506, 106914242, 19692253, 26189610, 37704083, 37555419,
                 128751439, 5460548, 21834858, 13374861, 80828723, 88238281, 89677535,
                 33144578, 66289588, 70263658, 69462583, 72533068, 77647740, 525161, 2976593,
                 69769087, 69223209, 32932025, 39598114, 50581891, 50284154, 42806544,
                 99349764, 69479345, 7581274, 12908137, 11985670, 27185057, 29563320,
                 37864629, 41236847, 47486829, 56426497, 56791453, 73754412, 15477397,
                 15394840, 18272658, 30309059, 50071269, 42831995, 52196994, 51850383,
                 62608232, 47804876),
    stringsAsFactors = FALSE
  )
  
  return(gene_data)
}


#' Determine gene list to use (custom or default)
#'
#' @param custom_gene_list_path Path to custom gene list file, or NULL to use default
#' @param verbose Print messages about which list is being used
#'
#' @return Data frame with validated gene list (columns: gene, chr, position)
#'
get_gene_list <- function(custom_gene_list_path = NULL, verbose = TRUE) {
  if (!is.null(custom_gene_list_path)) {
    if (verbose) {
      cat("Loading custom gene list from:", custom_gene_list_path, "\n")
    }
    gene_df <- load_gene_list(custom_gene_list_path, verbose = verbose)
  } else {
    if (verbose) {
      cat("Using built-in default gene list\n")
    }
    gene_df <- get_default_gene_list()
  }
  
  # Validate and normalize
  gene_df <- validate_and_normalize_gene_list(gene_df, verbose = verbose)
  
  return(gene_df)
}
