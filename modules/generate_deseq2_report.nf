process GENERATE_DESEQ2_REPORT {
    tag "interactive_report"
    publishDir "${params.outdir}/deseq2/report", mode: 'copy'
    
    // Usar tu imagen Docker personalizada
    container 'rnaseq_pipeline:latest'
    
    input:
    path 'deseq2_results.csv'
    path 'normalized_counts.csv'
    path 'samplesheet.csv'
    path 'DESeq2_Interactive_Report.Rmd'
    
    output:
    path "deseq2_interactive_report.html", emit: report
    path "versions.yml", emit: versions
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Renderizar reporte
    rmarkdown::render('DESeq2_Interactive_Report.Rmd',
        output_file = 'deseq2_interactive_report.html',
        params = list(
            deseq2_results = "deseq2_results.csv",
            normalized_counts = "normalized_counts.csv",
            samplesheet = "samplesheet.csv"
        )
    )
    
    # Generar versions.yml
    writeLines(
        c(
            '"GENERATE_DESEQ2_REPORT":',
            paste('    r-base:', R.Version()\$version.string),
            paste('    r-rmarkdown:', packageVersion('rmarkdown')),
            paste('    r-flexdashboard:', packageVersion('flexdashboard')),
            paste('    r-deseq2:', packageVersion('DESeq2')),
            paste('    r-tidyverse:', packageVersion('tidyverse')),
            paste('    r-plotly:', packageVersion('plotly'))
        ),
        'versions.yml'
    )
    """
}
