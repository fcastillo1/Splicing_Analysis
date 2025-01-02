// Inicio de proceso de fastqc_raw para los datos cuando estan crudos
process fastqc_raw {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads), val(condition)
    
    output:
    tuple val(sample_id), path("*_fastqc.{zip,html}"), emit: fastqc_output
    
    // Script de funcionamiento
    script:
    """
    #!/bin/bash
    fastqc ${reads.join(' ')} --outdir . --threads $task.cpus
    """
}
