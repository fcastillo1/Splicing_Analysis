#!/usr/bin/env Rscript

# Cargar las bibliotecas necesarias
suppressPackageStartupMessages({
    library(optparse)
    library(gprofiler2)
    library(ggplot2)
    library(htmlwidgets)
    library(yaml)
})

# Definir las opciones de línea de comandos
option_list <- list(
    make_option(c("--de_file"), type="character", default=NULL, help="Archivo de resultados de DESeq2"),
    make_option(c("--organism"), type="character", default="hsapiens", help="Organismo para gProfiler"),
    make_option(c("--correction_method"), type="character", default="g_SCS", help="Método de corrección"),
    make_option(c("--significance_threshold"), type="double", default=0.05, help="Umbral de significancia"),
    make_option(c("--padj_threshold"), type="double", default=0.05, help="Umbral de p-valor ajustado"),
    make_option(c("--lfc_threshold"), type="double", default=0.58, help="Umbral de log2 fold change"),
    make_option(c("--output_prefix"), type="character", default="gprofiler", help="Prefijo para archivos de salida")
)

# Parsear los argumentos de línea de comandos
opt <- parse_args(OptionParser(option_list=option_list))

# Función para manejar errores
handle_error <- function(e) {
    cat("Error:", conditionMessage(e), "\n")
    quit(status = 1)
}

# Usar tryCatch para manejar errores
tryCatch({
    # Leer el archivo de resultados de DESeq2
    deg_data <- read.csv(opt$de_file, header=TRUE, row.names=1)

    # Generar el background
    background_genes <- unique(gsub('"', '', deg_data$gene_name))
    write.table(background_genes, file="background_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

    # Filtrar los genes DE
    de_genes <- background_genes[abs(deg_data$log2FoldChange) > opt$lfc_threshold & deg_data$padj < opt$padj_threshold]

    # Imprimir información de depuración
    cat("Número total de genes en el background:", length(background_genes), "\n")
    cat("Número de genes DE:", length(de_genes), "\n")
    cat("Primeros genes DE:", head(de_genes), "\n")

    # Ejecutar gProfiler
    gost_results <- gost(query = de_genes,
                         organism = opt$organism,
                         correction_method = opt$correction_method,
                         user_threshold = opt$significance_threshold,
                         domain_scope = "custom",
                         custom_bg = background_genes,
                         sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF", "MIRNA", "HPA", "CORUM"))

    # Procesar y guardar resultados
    if (!is.null(gost_results$result) && nrow(gost_results$result) > 0) {
        # Guardar todas las rutas enriquecidas
        write.table(gost_results$result, 
                    file=paste0(opt$output_prefix, ".gprofiler2.all_enriched_pathways.tsv"), 
                    sep="\t", 
                    row.names=FALSE, 
                    quote=FALSE)
        
        # Guardar resultados como archivo RDS
        saveRDS(gost_results, file=paste0(opt$output_prefix, ".gprofiler2.gost_results.rds"))
        
        # Crear y guardar el gráfico interactivo
        saveWidget(gostplot(gost_results, capped = FALSE, interactive = TRUE), 
                   file=paste0(opt$output_prefix, ".gprofiler2.gostplot.html"))
        
        # Crear y guardar el gráfico PNG
        png(paste0(opt$output_prefix, ".gprofiler2.gostplot.png"), width=1000, height=800)
        print(gostplot(gost_results, capped = FALSE, interactive = FALSE))
        dev.off()
        
        # Crear gráficos para cada fuente
        for (source in unique(gost_results$result$source)) {
            subset_results <- gost_results$result[gost_results$result$source == source, ]
            
            # Guardar resultados de la fuente
            write.table(subset_results,
                        file=paste0(opt$output_prefix, ".gprofiler2.", source, ".sub_enriched_pathways.tsv"),
                        sep="\t",
                        row.names=FALSE,
                        quote=FALSE)
            
            # Crear gráfico de barras
            p <- ggplot(subset_results, aes(x=reorder(term_name, -log10(p_value)), y=-log10(p_value))) +
                 geom_bar(stat="identity") +
                 coord_flip() +
                 labs(title=paste("Enriched", source, "pathways"), x="", y="-log10(p-value)") +
                 theme_minimal() +
                 theme(axis.text.y = element_text(size = 8))
            
            # Guardar gráfico de barras
            ggsave(paste0(opt$output_prefix, ".gprofiler2.", source, ".sub_enriched_pathways.png"), 
                   p, width=10, height=8)
        }
    } else {
        cat("No se encontraron resultados significativos o el resultado está vacío.\n")
        # Crear archivos vacíos para satisfacer los requisitos de salida de Nextflow
        file.create(paste0(opt$output_prefix, ".gprofiler2.all_enriched_pathways.tsv"))
        file.create(paste0(opt$output_prefix, ".gprofiler2.gost_results.rds"))
        file.create(paste0(opt$output_prefix, ".gprofiler2.gostplot.html"))
        file.create(paste0(opt$output_prefix, ".gprofiler2.gostplot.png"))
    }

    # Guardar información de la sesión
    writeLines(capture.output(sessionInfo()), con=paste0(opt$output_prefix, "_R_sessionInfo.log"))

    # Guardar versiones del software
    versions <- list(
        R_version = R.version.string,
        gprofiler2_version = as.character(packageVersion("gprofiler2")),
        ggplot2_version = as.character(packageVersion("ggplot2")),
        htmlwidgets_version = as.character(packageVersion("htmlwidgets"))
    )
    write_yaml(versions, file="versions.yml")

}, error = handle_error)
