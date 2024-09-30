process GPROFILER2_GOST {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/gprofiler2", mode: 'copy'


    input:
    tuple val(meta), path(de_file)
    tuple val(meta), path(background_file)

    output:
    tuple val(meta), path("*.gprofiler2.all_enriched_pathways.tsv"), emit: all_enrich
    tuple val(meta), path("*.gprofiler2.gost_results.rds"), emit: rds
    tuple val(meta), path("*.gprofiler2.gostplot.png"), optional: true, emit: plot_png
    tuple val(meta), path("*.gprofiler2.gostplot.html"), optional: true, emit: plot_html
    tuple val(meta), path("*R_sessionInfo.log"), emit: session_info
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(gprofiler2)
    library(ggplot2)
    library(htmlwidgets)

    # Leer el archivo de resultados de DESeq2
    deg_data <- read.csv("${de_file}", header=TRUE, row.names=1)

    # Leer el archivo de background
    background_genes <- readLines("${background_file}")

    # Filtrar los genes DE
    de_genes <- deg_data\$gene_name[abs(deg_data\$log2FoldChange) > 0.58 & deg_data\$padj < 0.05]

    # Imprimir información de depuración
    cat("Número total de genes en el background:", length(background_genes), "\\n")
    cat("Número de genes DE:", length(de_genes), "\\n")
    cat("Primeros genes DE:", head(de_genes), "\\n")

    # Ejecutar gProfiler
    gost_results <- gost(query = de_genes,
                         organism = "hsapiens",
                         correction_method = "g_SCS",
                         user_threshold = 0.05,
                         domain_scope = "custom",
                         custom_bg = background_genes,
                         sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF", "MIRNA", "HPA", "CORUM"))

    # Guardar resultados como archivo RDS
    saveRDS(gost_results, file="${prefix}.gprofiler2.gost_results.rds")

    # Verificar si hay resultados y guardarlos
    if (!is.null(gost_results\$result) && nrow(gost_results\$result) > 0) {
        # Convertir columnas de lista a cadenas
        gost_results\$result <- as.data.frame(lapply(gost_results\$result, function(x) {
            if (is.list(x)) paste(unlist(x), collapse = ";") else x
        }))

        # Guardar todas las rutas enriquecidas
        write.table(gost_results\$result, 
                    file="${prefix}.gprofiler2.all_enriched_pathways.tsv", 
                    sep="\\t", 
                    row.names=FALSE, 
                    quote=FALSE)

        # Crear y guardar el gráfico interactivo
        saveWidget(gostplot(gost_results, capped = FALSE, interactive = TRUE), 
                   file="${prefix}.gprofiler2.gostplot.html")

        # Crear y guardar el gráfico PNG
        png("${prefix}.gprofiler2.gostplot.png", width=1000, height=800)
        print(gostplot(gost_results, capped = FALSE, interactive = FALSE))
        dev.off()
    } else {
        cat("No se encontraron resultados significativos en gProfiler.\\n")
        # Crear un archivo vacío para satisfacer la salida
        file.create("${prefix}.gprofiler2.all_enriched_pathways.tsv")
    }

    # Guardar información de la sesión
    writeLines(capture.output(sessionInfo()), con="${prefix}_R_sessionInfo.log")

    # Guardar versiones del software
    writeLines(
        c('"${task.process}":',
          paste0('    gprofiler2: ', packageVersion('gprofiler2')),
          paste0('    r-base: ', R.version.string),
          paste0('    r-ggplot2: ', packageVersion('ggplot2'))),
        "versions.yml"
    )
    """
}
