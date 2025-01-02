#!/usr/bin/env Rscript

# Cargar bibliotecas necesarias
suppressPackageStartupMessages({
    library(tximport)
    library(readr)
    library(dplyr)
    library(stringr)
    library(tibble)
})

# Función para manejo de errores
handle_error <- function(e) {
    message("Error: ", e$message)
    print(sessionInfo())
    quit(status = 1)
}

# Función para verificar archivo
check_file_exists <- function(file_path) {
    if (!file.exists(file_path)) {
        stop(sprintf("Archivo no encontrado: %s", file_path))
    }
}

# Función principal
main <- function() {
    tryCatch({
        # Verificar argumentos
        args <- commandArgs(trailingOnly=TRUE)
        if (length(args) != 3) {
            stop("Uso: tximport.R <tx2gene> <salmon_dir> <output_prefix>")
        }

        tx2gene_file <- args[1]
        salmon_dir <- args[2]
        output_prefix <- args[3]

        # Verificar archivos
        check_file_exists(tx2gene_file)
        check_file_exists(salmon_dir)

        # Leer archivo tx2gene
        message("Leyendo archivo tx2gene...")
        tx2gene <- read_tsv(tx2gene_file, 
                           col_names = c("tx", "gene_id", "gene_name"),
                           show_col_types = FALSE)

        # Guardar una copia de tx2gene con versiones intactas
        tx2gene_with_versions <- tx2gene

        # Remover números de versión solo para el procesamiento interno
        tx2gene$tx <- sub("\\.\\d+$", "", tx2gene$tx)
        
        message("Primeras líneas del archivo tx2gene procesado:")
        print(head(tx2gene))

        # Obtener archivos de Salmon
        fns <- list.files(salmon_dir, pattern = "quant.sf$", 
                         recursive = TRUE, full.names = TRUE)
        if (length(fns) == 0) {
            stop("No se encontraron archivos quant.sf en: ", salmon_dir)
        }
        names(fns) <- basename(dirname(fns))

        # Procesar cada archivo de Salmon
        message("Procesando archivos de Salmon...")
        for (file in fns) {
            message("Procesando archivo: ", file)
            quant_data <- read_tsv(file, show_col_types = FALSE)
            # Limpiar el ID y guardar versión original
            original_names <- quant_data$Name
            quant_data$Name <- sub("\\|.*$", "", quant_data$Name)
            quant_data$Name <- sub("\\.\\d+$", "", quant_data$Name)
            write_tsv(quant_data, file)
        }

        # Crear tx2gene simplificado para tximport
        tx2gene_simple <- tx2gene[, c("tx", "gene_id")]
        message("Primeras líneas de tx2gene_simple:")
        print(head(tx2gene_simple))

        # Ejecutar tximport
        message("Ejecutando tximport...")
        txi <- tximport(fns, type = "salmon", tx2gene = tx2gene_simple,
                       txOut = TRUE, ignoreTxVersion = TRUE)

        # Generar versiones sumarizadas
        message("Generando versiones sumarizadas...")
        txi.s <- summarizeToGene(txi, tx2gene = tx2gene_simple,
                                countsFromAbundance = "scaledTPM")
        txi.ls <- summarizeToGene(txi, tx2gene = tx2gene_simple,
                                 countsFromAbundance = "lengthScaledTPM")
        txi.dtu <- txi

        # Generar objetos a nivel de genes
        gi <- summarizeToGene(txi, tx2gene_simple)
        gi.s <- txi.s
        gi.ls <- txi.ls

        # Guardar resultados
        message("Guardando resultados...")
        
        # Guardar archivos RDS
        saveRDS(txi, file = paste0(output_prefix, ".txi.rds"))
        saveRDS(txi.s, file = paste0(output_prefix, ".txi.s.rds"))
        saveRDS(txi.ls, file = paste0(output_prefix, ".txi.ls.rds"))
        saveRDS(txi.dtu, file = paste0(output_prefix, ".txi.dtu.rds"))
        saveRDS(gi, file = paste0(output_prefix, ".gi.rds"))
        saveRDS(gi.s, file = paste0(output_prefix, ".gi.s.rds"))
        saveRDS(gi.ls, file = paste0(output_prefix, ".gi.ls.rds"))

        # Función para agregar columnas de ID y escribir archivos
        write_with_ids <- function(data, file_path, id_col = "gene_id") {
            df <- data.frame(id = rownames(data), data)
            colnames(df)[1] <- id_col
            write_tsv(df, file_path)
        }

        # Guardar archivos a nivel de genes
        write_with_ids(gi$abundance, paste0(output_prefix, ".gene_tpm.tsv"))
        write_with_ids(gi$counts, paste0(output_prefix, ".gene_counts.tsv"))
        write_with_ids(gi.s$abundance, paste0(output_prefix, ".gene_tpm_scaled.tsv"))
        write_with_ids(gi.s$counts, paste0(output_prefix, ".gene_counts_scaled.tsv"))
        write_with_ids(gi.ls$abundance, paste0(output_prefix, ".gene_tpm_length_scaled.tsv"))
        write_with_ids(gi.ls$counts, paste0(output_prefix, ".gene_counts_length_scaled.tsv"))

        # Función para guardar archivos de transcriptos con información de genes
        write_transcript_data <- function(data, tx2gene, file_path) {
            tx_data <- data.frame(
                tx = rownames(data),
                tx2gene[match(rownames(data), tx2gene$tx), ],
                data
            )
            write_tsv(tx_data, file_path)
        }

        # Guardar archivos a nivel de transcriptos
        write_transcript_data(txi$abundance, tx2gene_with_versions, 
                            paste0(output_prefix, ".transcript_tpm.tsv"))
        write_transcript_data(txi$counts, tx2gene_with_versions,
                            paste0(output_prefix, ".transcript_counts.tsv"))
        write_transcript_data(txi.s$abundance, tx2gene_with_versions,
                            paste0(output_prefix, ".transcript_tpm_scaled.tsv"))
        write_transcript_data(txi.s$counts, tx2gene_with_versions,
                            paste0(output_prefix, ".transcript_counts_scaled.tsv"))
        write_transcript_data(txi.ls$abundance, tx2gene_with_versions,
                            paste0(output_prefix, ".transcript_tpm_length_scaled.tsv"))
        write_transcript_data(txi.ls$counts, tx2gene_with_versions,
                            paste0(output_prefix, ".transcript_counts_length_scaled.tsv"))
        write_transcript_data(txi.dtu$abundance, tx2gene_with_versions,
                            paste0(output_prefix, ".transcript_tpm_dtu_scaled.tsv"))
        write_transcript_data(txi.dtu$counts, tx2gene_with_versions,
                            paste0(output_prefix, ".transcript_counts_dtu_scaled.tsv"))

        # Guardar tx2gene
        write_tsv(tx2gene_with_versions, "tximport.tx2gene.tsv")

        # Preparar y guardar TPM para SUPPA
        message("Preparando archivo TPM para SUPPA...")

        # Crear matriz TPM con IDs de transcriptos
        tpm_data <- txi$abundance

        # Para asegurar IDs únicos, usar el mapeo original de tx2gene_with_versions
        transcript_ids <- tx2gene_with_versions$tx
        names(transcript_ids) <- sub("\\.\\d+$", "", tx2gene_with_versions$tx)

        # Mapear los IDs manteniendo las versiones
        new_rownames <- transcript_ids[rownames(tpm_data)]

        # Eliminar filas con nombres NA y duplicados
        valid_rows <- !is.na(new_rownames) & !duplicated(new_rownames)
        tpm_data <- tpm_data[valid_rows, ]
        new_rownames <- new_rownames[valid_rows]

        # Asignar los nuevos nombres de fila
        rownames(tpm_data) <- new_rownames

        # Verificar y limpiar datos TPM
        tpm_data[is.na(tpm_data)] <- 0
        tpm_data <- round(tpm_data, digits = 6)

        # Guardar archivo TPM para SUPPA
        write.table(tpm_data, "suppa_tpm.txt", 
                sep = "\t", 
                quote = FALSE, 
                row.names = TRUE,
                col.names = TRUE)

        # Verificar archivo generado
        message("Verificando archivo SUPPA TPM...")
        tpm_check <- read.table("suppa_tpm.txt", header = TRUE, sep = "\t", check.names = FALSE)
        message("Dimensiones del archivo TPM: ", paste(dim(tpm_check), collapse = " x "))
        message("Primeras líneas del archivo TPM para SUPPA:")
        print(head(tpm_check))

        message("Tximport completado exitosamente")
        
    }, error = handle_error)
}

# Ejecutar función principal
main()

# Imprimir información de la sesión
sessionInfo()
