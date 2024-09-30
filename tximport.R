#!/usr/bin/env Rscript

# Suppress package startup messages
suppressPackageStartupMessages({
  library(tximport)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# Function for error handling
handle_error <- function(e) {
  message("Error: ", e$message)
  print(sessionInfo())
  quit(status = 1)
}

# Main function
main <- function() {
  tryCatch({
    message("R version: ", R.version.string)
    message("tximport version: ", packageVersion("tximport"))

    # Parse command line arguments
    args <- commandArgs(trailingOnly=TRUE)
    if (length(args) != 3) {
      stop("Usage: Rscript tximport.R <tx2gene_file> <salmon_dir> <output_prefix>")
    }
    tx2gene_file <- args[1]
    salmon_dir <- args[2]
    output_prefix <- args[3]

    # Read tx2gene file
    message("Reading tx2gene file...")
    tx2gene <- read_tsv(tx2gene_file, col_names = c("tx", "gene"), show_col_types = FALSE)

    # Print tx2gene info for debugging
    message("TX2GENE structure:")
    print(head(tx2gene))
    print(dim(tx2gene))

    # Get Salmon quant files
    message("Searching for Salmon quant files...")
    quant_files <- list.files(salmon_dir, pattern = "quant.sf$", recursive = TRUE, full.names = TRUE)
    if (length(quant_files) == 0) {
      stop("No quant.sf files found in the specified directory")
    }
    names(quant_files) <- basename(dirname(quant_files))

    # Print quant files for debugging
    message("Quant files found:")
    print(quant_files)

    # Function to clean transcript IDs
    clean_ids <- function(ids) {
        str_split(ids, "\\|", n = 2) %>% sapply("[", 1)
    }

    # Read the first quant file to get transcript IDs
    first_quant <- read_tsv(quant_files[1], show_col_types = FALSE)
    salmon_tx_ids <- clean_ids(first_quant$Name)

    message("Number of transcripts before filtering: ", nrow(tx2gene))
    message("Number of transcripts in Salmon output: ", length(salmon_tx_ids))

    # Filter tx2gene to only include transcripts present in Salmon output
    tx2gene_filtered <- tx2gene %>%
      filter(tx %in% salmon_tx_ids)

    message("Number of transcripts after filtering: ", nrow(tx2gene_filtered))

    missing_tx <- setdiff(tx2gene$tx, salmon_tx_ids)
    message("Number of transcripts in tx2gene but not in Salmon output: ", length(missing_tx))
    message("First few missing transcripts: ", paste(head(missing_tx), collapse=", "))

    message("TPM values for first few transcripts in first sample before processing:")
    print(head(first_quant$TPM))

    # Run Tximport with different settings
    message("Running tximport...")
    txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene_filtered, txOut = TRUE, ignoreTxVersion = TRUE)

    message("TPM values for first few transcripts in first sample after processing:")
    print(head(txi$abundance[,1]))

    message("Summarizing to gene level...")
    txi.s <- summarizeToGene(txi, tx2gene = tx2gene_filtered, countsFromAbundance = "scaledTPM")
    txi.ls <- summarizeToGene(txi, tx2gene = tx2gene_filtered, countsFromAbundance = "lengthScaledTPM")
    txi.dtu <- txi  # Already at transcript level

    gi <- summarizeToGene(txi, tx2gene_filtered)
    gi.s <- txi.s
    gi.ls <- txi.ls

    # Save results
    message("Saving results...")
    saveRDS(txi, file = paste0(output_prefix, ".txi.rds"))
    saveRDS(txi.s, file = paste0(output_prefix, ".txi.s.rds"))
    saveRDS(txi.ls, file = paste0(output_prefix, ".txi.ls.rds"))
    saveRDS(txi.dtu, file = paste0(output_prefix, ".txi.dtu.rds"))
    saveRDS(gi, file = paste0(output_prefix, ".gi.rds"))
    saveRDS(gi.s, file = paste0(output_prefix, ".gi.s.rds"))
    saveRDS(gi.ls, file = paste0(output_prefix, ".gi.ls.rds"))

    write_tsv(as_tibble(gi$abundance, rownames = "gene_id"), file = paste0(output_prefix, ".gene_tpm.tsv"))
    write_tsv(as_tibble(gi$counts, rownames = "gene_id"), file = paste0(output_prefix, ".gene_counts.tsv"))
    write_tsv(as_tibble(gi.s$abundance, rownames = "gene_id"), file = paste0(output_prefix, ".gene_tpm_scaled.tsv"))
    write_tsv(as_tibble(gi.s$counts, rownames = "gene_id"), file = paste0(output_prefix, ".gene_counts_scaled.tsv"))
    write_tsv(as_tibble(gi.ls$abundance, rownames = "gene_id"), file = paste0(output_prefix, ".gene_tpm_length_scaled.tsv"))
    write_tsv(as_tibble(gi.ls$counts, rownames = "gene_id"), file = paste0(output_prefix, ".gene_counts_length_scaled.tsv"))
    write_tsv(as_tibble(txi$abundance, rownames = "transcript_id"), file = paste0(output_prefix, ".transcript_tpm.tsv"))
    write_tsv(as_tibble(txi$counts, rownames = "transcript_id"), file = paste0(output_prefix, ".transcript_counts.tsv"))
    write_tsv(as_tibble(txi.s$abundance, rownames = "transcript_id"), file = paste0(output_prefix, ".transcript_tpm_scaled.tsv"))
    write_tsv(as_tibble(txi.s$counts, rownames = "transcript_id"), file = paste0(output_prefix, ".transcript_counts_scaled.tsv"))
    write_tsv(as_tibble(txi.ls$abundance, rownames = "transcript_id"), file = paste0(output_prefix, ".transcript_tpm_length_scaled.tsv"))
    write_tsv(as_tibble(txi.ls$counts, rownames = "transcript_id"), file = paste0(output_prefix, ".transcript_counts_length_scaled.tsv"))
    write_tsv(as_tibble(txi.dtu$abundance, rownames = "transcript_id"), file = paste0(output_prefix, ".transcript_tpm_dtu_scaled.tsv"))
    write_tsv(as_tibble(txi.dtu$counts, rownames = "transcript_id"), file = paste0(output_prefix, ".transcript_counts_dtu_scaled.tsv"))

    write_tsv(tx2gene_filtered, "tximport.tx2gene.tsv")
    write.table(txi$abundance, "suppa_tpm.txt", sep = "\t", quote = FALSE, row.names = TRUE)

    message("Tximport completed successfully")
  }, error = handle_error)
}
