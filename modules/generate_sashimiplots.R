#!/usr/bin/env Rscript

library(optparse)

# Define las opciones que se necesitan para hacer el sashimiplot
option_list = list(
  make_option(c("-d", "--rMATSdir"), type="character", default=NULL, help="Directory of rMATS results", metavar="character"),
  make_option(c("-e", "--eventTypes"), type="character", default="SE,A3SS,A5SS,MXE,RI", help="Event types separated by commas", metavar="character"),
  make_option(c("-x", "--bam1"), type="character", default=NULL, help="File with list of BAMs for condition 1", metavar="character"),
  make_option(c("-y", "--bam2"), type="character", default=NULL, help="File with list of BAMs for condition 2", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="sashimi_plots", help="Output directory", metavar="character"),
  make_option(c("-l", "--label1"), type="character", default="Condition1", help="Label for condition 1", metavar="character"),
  make_option(c("-m", "--label2"), type="character", default="Condition2", help="Label for condition 2", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Filtra por el tipo de evento
filter_events_by_type <- function(file, event_type) {
  if (file.exists(file)) {
    events <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    filtered_events <- events[events$eventType == event_type,]
    return(filtered_events)
  } else {
    cat("Warning: The file", file, "does not exist.\n")
    return(data.frame())  # Return an empty data frame if the file does not exist
  }
}

# Se leen los archivos bam
bam1_files <- readLines(opt$bam1)
bam2_files <- readLines(opt$bam2)

# Se genera el output
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# Se procesa el evento
event_types <- strsplit(opt$eventTypes, ",")[[1]]

# Se genera los archivos en el directorio
rmats_files <- list.files(path = opt$rMATSdir, pattern = "\\.txt$", full.names = TRUE)

for (event_type in event_types) {
  # Se busca el evento
  event_files <- grep(paste0(event_type, ".MATS.(JC|JCEC).txt$"), rmats_files, value = TRUE)

  for (event_file in event_files) {
    if (file.exists(event_file)) {
      filtered_events <- filter_events_by_type(event_file, event_type)

      cat("Filtered events for", event_type, ":", nrow(filtered_events), "\n")

      # Se genera el sashimi con el evemto filtrado
      for (i in 1:nrow(filtered_events)) {
        event <- filtered_events[i,]
        out_dir <- file.path(opt$outdir, event_type, event$geneSymbol)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        # Se realiza el comando del sashimi
        cmd <- paste("rmats2sashimiplot",
                     "--b1", shQuote(paste(bam1_files, collapse=",")),
                     "--b2", shQuote(paste(bam2_files, collapse=",")),
                     "--event-type", event_type,
                     "-e", shQuote(event_file),
                     "--l1", shQuote(opt$label1),
                     "--l2", shQuote(opt$label2),
                     "--exon_s", "1",
                     "--intron_s", "5",
                     "-o", shQuote(out_dir))

        # Se informa de la ejecucion
        cat("Executing command for", event$geneSymbol, ":", cmd, "\n")
        system(cmd)

        # Verify if sashimiplots were generated
        if (length(list.files(out_dir)) > 0) {
          cat("Sashimiplot generated for", event$geneSymbol, "in", out_dir, "\n")
        } else {
          cat("Error: No sashimiplot generated for", event$geneSymbol, "\n")
        }
      }
    } else {
      cat("Warning: The file", event_file, "does not exist.\n")
    }
  }
}

cat("Sashimiplot generation process completed.\n")
