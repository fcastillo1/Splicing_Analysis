// Inicio de proceso de fastqc para los datos despues de hacer trimming
process fastqc_trimmed {
    publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy'
    label 'fastqc'

    input:
    tuple val(sample_id), path(trimmed_reads), val(condition)

    output:
    tuple val(sample_id), path("*_fastqc.{zip,html}"), emit: fastqc_output
    

    // Script de funcionamiento
    script:
    """
    fastqc ${trimmed_reads} --outdir . --threads ${task.cpus}
    
    """
}
