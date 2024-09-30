process fastqc_raw {
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    label 'fastqc'

    input:
    tuple val(sample_id), path(reads), val(condition)

    output:
    tuple val(sample_id), path("*_fastqc.{zip,html}"), emit: fastqc_output

    script:
    """
    fastqc ${reads} --outdir . --threads ${task.cpus}
    """
}
