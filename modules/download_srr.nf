// incio del proceso de descarga
process download_srr {
    publishDir "${params.outdir}/downloaded", mode: 'copy'

    input:
    tuple val(srr_id), val(condition)

    output:
    tuple val(srr_id), path("${srr_id}_1.fastq.gz"), val(condition)

// Script de funcionamiento
    script:
    """
    # Descarga y Procesamiento
    prefetch ${srr_id}
    fastq-dump --gzip --split-3 ${srr_id}
    
    if [ -f ${srr_id}_1.fastq.gz ]; then
        echo "Single-end read detected"
    elif [ -f ${srr_id}.fastq.gz ]; then
        mv ${srr_id}.fastq.gz ${srr_id}_1.fastq.gz
    else
        echo "Error: No fastq file generated" >&2
        exit 1
    fi
    """
}
