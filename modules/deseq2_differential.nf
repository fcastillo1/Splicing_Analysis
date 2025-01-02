// Inicio del proceso de deseq2
process DESEQ2_DIFFERENTIAL {
    tag "DESeq2_analysis"
    label 'process_medium'
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path(counts)
    path(samplesheet)
    val(contrast_variable)
    val(reference_level)
    val(target_level)
    path(gtf)

    output:
    path "deseq2_results.csv", emit: results
    path "normalized_counts.csv", emit: normalized_counts
    path "*.pdf", emit: plots
    path "*.rds", emit: rds
    path "versions.yml", emit: versions

    // Script de funcionamiento
    script:
    """
    #!/usr/bin/env Rscript

    library(DESeq2)
    library(tximport)
    library(readr)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(rtracklayer)
    library(EnhancedVolcano)

    # 1. Leer datos
    counts_matrix <- read.delim("salmon.merged.gene_counts_length_scaled.tsv",
                              row.names=1,
                              check.names=FALSE)

    sample_info <- read.csv("${samplesheet}", stringsAsFactors=FALSE)

    # 2. Asegurar que el orden de las muestras coincide
    sample_info <- sample_info[match(colnames(counts_matrix), sample_info\$sample_id),]

    # 3. Crear objeto DESeq2
    dds <- DESeqDataSetFromMatrix(
        countData = round(counts_matrix),
        colData = sample_info,
        design = formula(paste0("~", "${contrast_variable}"))
    )

    # 4. Definir el nivel de referencia
    dds\$${contrast_variable} <- relevel(factor(dds\$${contrast_variable}),
                                       ref = "${reference_level}")

    # 5. Filtrar genes con baja expresión
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]

    # 6. Ejecutar DESeq2
    dds <- DESeq(dds)

    # 7. Obtener resultados
    res <- results(dds,
                  contrast=c("${contrast_variable}","${target_level}","${reference_level}"),
                  alpha=${params.alpha})

    # 8. Añadir información de genes
    gtf_data <- import("${gtf}")
    gene_info <- data.frame(
        gene_id = gtf_data[gtf_data\$type == "gene"]\$gene_id,
        gene_name = gtf_data[gtf_data\$type == "gene"]\$gene_name,
        stringsAsFactors = FALSE
    )

    # 9. Procesar y guardar resultados
    res_df <- as.data.frame(res)
    res_df\$gene_id <- rownames(res_df)
    res_df <- merge(res_df, gene_info, by="gene_id", all.x=TRUE)

    # Añadir columnas para el volcano plot
    res_df\$significant <- ifelse(res_df\$padj < ${params.alpha} & abs(res_df\$log2FoldChange) > ${params.lfc_threshold}, 
                                ifelse(res_df\$log2FoldChange > 0, "Up", "Down"), "NS")

    write.csv(res_df, "deseq2_results.csv", row.names=FALSE)

    # 10. Generar y guardar visualizaciones
    # MA Plot
    pdf("ma_plot.pdf", width=10, height=8)
    plotMA(res, main="MA Plot")
    dev.off()

    # PCA Plot
    vsd <- vst(dds, blind=FALSE)
    pdf("pca_plot.pdf", width=10, height=8)
    plotPCA(vsd, intgroup="${contrast_variable}")
    dev.off()

    # Volcano Plot
    pdf("volcano_plot.pdf", width=12, height=10)
    volcano_plot <- EnhancedVolcano(res_df,
        lab = res_df\$gene_name,
        x = 'log2FoldChange',
        y = 'padj',
        title = paste('Differential Expression:', "${target_level}", 'vs', "${reference_level}"),
        pCutoff = ${params.alpha},
        FCcutoff = ${params.lfc_threshold},
        pointSize = 3.0,
        labSize = 4.0,
        col = c('grey', 'blue', 'blue', 'red'),
        colAlpha = 0.4,
        legendPosition = 'right',
        legendLabSize = 12,
        legendIconSize = 4.0,
        drawConnectors = TRUE,
        widthConnectors = 0.75
    )
    print(volcano_plot)
    dev.off()

    # Heatmap de los top genes
    top_genes <- head(order(res\$padj), 50)
    mat <- assay(vsd)[top_genes, ]
    annotation_col <- data.frame(
        Condition = dds\$${contrast_variable},
        row.names = colnames(mat)
    )
    pdf("heatmap.pdf", width=12, height=10)
    pheatmap(mat,
            annotation_col = annotation_col,
            scale = "row",
            main = "Top 50 Differentially Expressed Genes",
            fontsize_row = 8,
            fontsize_col = 10)
    dev.off()

    # 11. Guardar conteos normalizados
    normalized_counts <- counts(dds, normalized=TRUE)
    write.csv(normalized_counts, "normalized_counts.csv")

    # 12. Guardar objeto DESeq2
    saveRDS(dds, "dds.rds")

    # Generar resumen de resultados
    sig_genes <- sum(res_df\$padj < ${params.alpha} & !is.na(res_df\$padj))
    up_genes <- sum(res_df\$significant == "Up", na.rm=TRUE)
    down_genes <- sum(res_df\$significant == "Down", na.rm=TRUE)

    writeLines(
        c(
            "Resumen del análisis:",
            paste("Total genes analizados:", nrow(res_df)),
            paste("Genes significativos (padj <", ${params.alpha}, "):", sig_genes),
            paste("Genes up-regulados:", up_genes),
            paste("Genes down-regulados:", down_genes)
        ),
        "analysis_summary.txt"
    )

    # Crear archivo de versiones
    writeLines(
        c(
            '"${task.process}":',
            paste('    r-base:', R.Version()\$version.string),
            paste('    bioconductor-deseq2:', packageVersion('DESeq2')),
            paste('    enhancedvolcano:', packageVersion('EnhancedVolcano'))
        ),
        'versions.yml'
    )
    """
}
