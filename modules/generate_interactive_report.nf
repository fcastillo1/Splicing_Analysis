process GENERATE_INTERACTIVE_REPORT {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    tuple val(meta), path(deseq2_results)
    tuple val(meta), path(normalized_counts)
    path gtf_file
    path samplesheet
    tuple val(meta), path(gprofiler_results)
    tuple val(meta), path(gprofiler_plot_png)
    
    path rmd_template

    output:
    path "DESeq2_Interactive_Report.html"

    script:
    """
    # Renombrar los archivos de entrada si es necesario
    [ "$deseq2_results" != "deseq2_results.csv" ] && cp $deseq2_results deseq2_results.csv || true
    [ "$normalized_counts" != "normalized_counts.csv" ] && cp $normalized_counts normalized_counts.csv || true
    [ "$gtf_file" != "annotation.gtf" ] && cp $gtf_file annotation.gtf || true
    [ "$samplesheet" != "samplesheet.csv" ] && cp $samplesheet samplesheet.csv || true
    [ "$gprofiler_results" != "deseq2.gprofiler2.all_enriched_pathways.tsv" ] && cp $gprofiler_results deseq2.gprofiler2.all_enriched_pathways.tsv || true
    [ "$gprofiler_plot_png" != "deseq2.gprofiler2.gostplot.png" ] && cp $gprofiler_plot_png deseq2.gprofiler2.gostplot.png || true
    [ "$rmd_template" != "DESeq2_Interactive_Report.Rmd" ] && cp $rmd_template DESeq2_Interactive_Report.Rmd || true

    # Establecer la variable de entorno RESULTS_DIR
    export RESULTS_DIR="${params.outdir}"

    # Listar archivos en el directorio actual para depuración
    ls -l

    # Renderizar el informe
    Rscript -e "rmarkdown::render('DESeq2_Interactive_Report.Rmd', output_file = 'DESeq2_Interactive_Report.html')"
    """
}
