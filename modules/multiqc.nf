process multiqc_run {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"
    path "versions.yml", emit: versions

    script:
    """
    multiqc .
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed -e "s/multiqc, version //g")
    END_VERSIONS
    """
}
