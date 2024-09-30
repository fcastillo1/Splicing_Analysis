#!/usr/bin/env Rscript

# Cargar librerías necesarias
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(ggrepel)
    library(EnhancedVolcano)
    library(dplyr)
    library(stringr)
})

# Función para registrar mensajes
log_message <- function(message) {
    cat(paste0(Sys.time(), ": ", message, "\n"))
}

# Función para leer archivos delimitados de manera flexible
read_delim_flexible <- function(file, header = TRUE, row.names = NULL, check.names = FALSE) {
    ext <- tolower(tail(strsplit(basename(file), split = "\\.")[[1]], 1))
    separator <- if (ext %in% c("tsv", "txt")) "\t" else if (ext == "csv") "," else stop(paste("Unknown file extension:", ext))
    read.delim(file, sep = separator, header = header, row.names = row.names, check.names = check.names)
}

# Parsear argumentos de línea de comandos
args <- commandArgs(trailingOnly = TRUE)
params <- list()
for (i in seq(1, length(args), 2)) {
    params[[sub("^--", "", args[i])]] <- args[i + 1]
}

log_message("Parámetros utilizados:")
print(params)

# Verificar archivo GTF
log_message(paste("GTF file:", params$gtf_file))
if (is.null(params$gtf_file) || !file.exists(params$gtf_file)) {
    stop(paste("El archivo GTF no existe o no se ha especificado:", params$gtf_file))
}

# Verificar otros archivos de entrada
required_files <- c("count_file", "sample_file")
for (file in required_files) {
    if (!file.exists(params[[file]])) {
        stop(paste("El archivo", file, "no existe:", params[[file]]))
    }
}

tryCatch({
    # Leer archivos de entrada
    log_message("Leyendo archivos de entrada")
    counts <- read_delim_flexible(params$count_file)
    samples <- read_delim_flexible(params$sample_file)

    # Manejar la columna gene_id si existe
    if ("gene_id" %in% colnames(counts)) {
        rownames(counts) <- counts$gene_id
        counts <- counts[, -which(colnames(counts) == "gene_id")]
    }

    # Asegurarse de que las muestras en counts coincidan con las de samples
    counts <- counts[, samples$sample_id]

    log_message(paste("Dimensiones de la matriz de conteos:", paste(dim(counts), collapse="x")))
    log_message(paste("Número de muestras:", nrow(samples)))

    # Parsear el contraste
    contrast_parts <- strsplit(params$contrast, ",")[[1]]
    if (length(contrast_parts) != 3) {
        stop("El formato del contraste debe ser: variable,referencia,objetivo")
    }
    contrast_variable <- contrast_parts[1]
    reference_level <- contrast_parts[2]
    target_level <- contrast_parts[3]

    log_message(paste("Contraste:", paste(contrast_parts, collapse=" vs ")))

    # Asegurar que la condición sea un factor con el nivel de referencia correcto
    samples[[contrast_variable]] <- factor(samples[[contrast_variable]], levels = c(reference_level, target_level))

    # Crear DESeqDataSet
    log_message("Creando DESeqDataSet")
    dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                  colData = samples,
                                  design = as.formula(paste("~", contrast_variable)))

    # Filtrar genes de baja expresión
    log_message("Filtrando genes de baja expresión")
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]

    # Ejecutar DESeq2
    log_message("Ejecutando DESeq2")
    dds <- DESeq(dds)

    # Obtener resultados
    log_message("Obteniendo resultados")
    res <- results(dds, contrast = c(contrast_variable, target_level, reference_level))

    if (params$shrink_lfc == "true") {
        res <- lfcShrink(dds, contrast = c(contrast_variable, target_level, reference_level), res = res, type = "ashr")
    }

    # Función mejorada para crear mapeo de ID a nombre
    create_id_to_name_map <- function(gtf_file, batch_size = 1000000) {
        log_message("Creando mapeo de ID a nombre desde el archivo GTF")
        
        id_to_name <- list()
        con <- file(gtf_file, "r")
        batch <- readLines(con, n = batch_size)
        
        while(length(batch) > 0) {
            gene_lines <- grep("^[^#].+gene\\t", batch, value = TRUE)
            
            if(length(gene_lines) > 0) {
                gene_ids <- sapply(strsplit(gene_lines, "\t"), function(x) {
                    attr <- x[9]
                    gsub(".*gene_id \"([^\"]+)\".*", "\\1", attr)
                })
                
                gene_names <- sapply(strsplit(gene_lines, "\t"), function(x) {
                    attr <- x[9]
                    gsub(".*gene_name \"([^\"]+)\".*", "\\1", attr)
                })
                
                id_to_name <- c(id_to_name, setNames(gene_names, gene_ids))
            }
            
            batch <- readLines(con, n = batch_size)
        }
        
        close(con)
        
        id_to_name <- unlist(id_to_name)
        log_message(paste("Se encontraron", length(id_to_name), "mapeos de ID a nombre"))
        
        return(id_to_name)
    }

    # Crear mapeo de ID a nombre
    id_to_name_map <- create_id_to_name_map(params$gtf_file)

    # Verificar el mapeo
    if (!is.null(id_to_name_map)) {
        log_message(paste("Número total de mapeos ID a nombre:", length(id_to_name_map)))
        log_message(paste("Primeros 5 mapeos:", paste(head(names(id_to_name_map)), head(id_to_name_map), sep="=", collapse=", ")))
        
        # Verificar cuántos IDs en los resultados tienen un nombre de gen correspondiente
        matching_ids <- sum(rownames(res) %in% names(id_to_name_map))
        log_message(paste("Número de IDs en los resultados que tienen un nombre de gen correspondiente:", matching_ids))
    } else {
        log_message("No se pudo crear el mapeo de ID a nombre. Se usarán los IDs originales.")
    }

    # Función para actualizar los resultados con nombres de genes
    update_results_with_gene_names <- function(res, id_to_name) {
        log_message("Actualizando resultados de DESeq2 con nombres de genes")
        
        res_df <- as.data.frame(res)
        res_df$gene_name <- id_to_name[rownames(res_df)]
        res_df$gene_name[is.na(res_df$gene_name)] <- rownames(res_df)[is.na(res_df$gene_name)]
        
        return(res_df)
    }

    # Actualizar resultados con nombres de genes
    res_with_names <- update_results_with_gene_names(res, id_to_name_map)

    # Guardar resultados
    log_message("Guardando resultados")
    write.csv(res_with_names, file = paste0(params$output_prefix, "_results.csv"))

    # PCA Plot
    log_message("Generando PCA Plot")
    vsd <- varianceStabilizingTransformation(dds)
    pcaData <- plotPCA(vsd, intgroup = c(contrast_variable), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))

    pdf(paste0(params$output_prefix, "_pca_plot.pdf"), width = 18, height = 10)
    print(ggplot(pcaData, aes(x = PC1, y = PC2, color = !!sym(contrast_variable), shape = !!sym(contrast_variable))) +
        geom_point(size = 3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle("Principal Component Analysis") +
        theme_bw(base_size = 10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
              legend.title = element_text(face = "bold", size = 14),
              legend.text = element_text(size = 12),
              axis.title = element_text(face = "bold", size = 16),
              axis.text = element_text(size = 12),
              legend.position = "right") +
        scale_color_brewer(palette = "Set1") +
        theme_minimal() + 
        scale_shape_manual(values = c(16, 17, 15)) +
        coord_fixed(ratio = 1) +
        stat_ellipse(aes(color = !!sym(contrast_variable)), type = "t", level = 0.95, linetype = 2))
    dev.off()
    log_message("PCA Plot generado con éxito")

    # Heatmap
    log_message("Generando Heatmap")
    topGenes <- head(order(res$padj), 30)
    mat <- assay(vsd)[topGenes, ]
    mat <- t(scale(t(mat)))
    df <- as.data.frame(colData(dds)[, contrast_variable, drop = FALSE])

    rownames(mat) <- res_with_names$gene_name[topGenes]

    pdf(paste0(params$output_prefix, "_heatmap.pdf"), width = 12, height = 14)
    pheatmap(mat, 
             annotation_col = df,
             main = "Top 30 Differentially Expressed Genes",
             show_rownames = TRUE,
             fontsize_row = 10,
             fontsize_col = 10,
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
    dev.off()
    log_message("Heatmap generado con éxito")

    # MA Plot
    log_message("Generando MA Plot")
    pdf(paste0(params$output_prefix, "_ma_plot.pdf"), width = 12, height = 10)
    plotMA(res, ylim = c(-2, 2))
    dev.off()
    log_message("MA Plot generado con éxito")

    # Volcano Plot
    # Preparación de datos para los gráficos de volcán
    volcano_data <- as.data.frame(res)
    volcano_data$gene_name <- res_with_names$gene_name
    volcano_data$diffexpressed <- "Not Significant"
    volcano_data$diffexpressed[volcano_data$log2FoldChange > 1 & volcano_data$padj < 0.05] <- "Up"
    volcano_data$diffexpressed[volcano_data$log2FoldChange < -1 & volcano_data$padj < 0.05] <- "Down"

    # Seleccionar los top genes más significativos
    top_genes <- volcano_data %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>%
      head(15)

    # Gráfico de volcán con etiquetas
    log_message("Generando Volcano Plot con etiquetas")
    pdf(paste0(params$output_prefix, "_volcano_plot_labeled.pdf"), width = 12, height = 10)
    print(EnhancedVolcano(volcano_data,
        lab = volcano_data$gene_name,
        x = 'log2FoldChange',
        y = 'padj',
        title = paste('Differential Expression Analysis:', target_level, 'vs', reference_level),
        subtitle = paste("Total =", nrow(volcano_data), "variables"),
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = 1.5,
        labSize = 4.0,
        col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
        colAlpha = 0.5,
        legendPosition = 'right',
        legendLabSize = 12,
        legendIconSize = 4.0,
        drawConnectors = TRUE,
        widthConnectors = 0.5,
        selectLab = top_genes$gene_name))
    dev.off()
    log_message("Volcano Plot con etiquetas generado con éxito")

    # Gráfico de volcán sin etiquetas
    log_message("Generando Volcano Plot sin etiquetas")
    pdf(paste0(params$output_prefix, "_volcano_plot_no_labels.pdf"), width = 10, height = 8)
    print(ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed)) +
      geom_point(size = 1.5, alpha = 0.7) +
      scale_color_manual(values = c("Down" = "#00BFC4", "Not Significant" = "grey", "Up" = "#F8766D"),
                         labels = c(paste("Down regulated (", sum(volcano_data$diffexpressed == "Down"), ")", sep=""),
                                    "Not significant",
                                    paste("Up regulated (", sum(volcano_data$diffexpressed == "Up"), ")", sep=""))) +
      theme_minimal() +
      theme(legend.position = "top") +
      xlab("log2(Fold Change)") +
      ylab("-log10(Adjusted p-value)") +
      ggtitle("Volcano Plot", subtitle = paste("Total =", nrow(volcano_data), "variables")) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed"))
    dev.off()
    log_message("Volcano Plot sin etiquetas generado con éxito")

    # Dispersion Plot
    log_message("Generando Dispersion Plot")
    pdf(paste0(params$output_prefix, "_dispersion_plot.pdf"), width = 10, height = 8)
    plotDispEsts(dds)
    dev.off()
    log_message("Dispersion Plot generado con éxito")

    # Guardar conteos normalizados
    log_message("Guardando conteos normalizados")
    normalized_counts <- counts(dds, normalized = TRUE)
    write.csv(normalized_counts, 
              file = paste0(params$output_prefix, "_normalized_counts.csv"))

    # Guardar factores de tamaño
    log_message("Guardando factores de tamaño")
    # Guardar factores de tamaño
    log_message("Guardando factores de tamaño")
    size_factors <- sizeFactors(dds)
    write.csv(data.frame(sample = names(size_factors), size_factor = size_factors),
              file = paste0(params$output_prefix, "_sizefactors.csv"), 
              row.names = FALSE)

    # Guardar objeto dds
    log_message("Guardando objeto dds")
    saveRDS(dds, file = paste0(params$output_prefix, "_dds.rds"))

    # Log
    log_message("Generando archivo de log")
    sink(paste0(params$output_prefix, "_R_sessionInfo.log"))
    print(sessionInfo())
    sink()

    log_message("Análisis DESeq2 completado")

}, error = function(e) {
    log_message(paste("Error en el análisis principal:", e$message))
    print(traceback())
    quit(status = 1)
})

# Si se llega a este punto, el script se ha ejecutado con éxito
quit(status = 0)
