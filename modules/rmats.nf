process RMATS_ORIGINAL_GTF {
    cache false
    publishDir "${params.outdir}/rMATS_out/original_gtf", mode: 'copy'
    label 'rmats'

    input:
    path original_gtf
    path bam1
    path bam2

    output:
    path "*", emit: results

    script:
    """
    # Ejecucion de rMATS
    echo "Executing rMATS with original GTF"
    rmats.py --b1 ${bam1} \
        --b2 ${bam2} \
        --gtf ${original_gtf} \
        -t paired \
        --readLength ${params.readlength} \
        --nthread ${task.cpus} \
        --od ./ \
        --tmp ./tmp_original \
        --libType fr-firststrand \
        --novelSS \
        --mil ${params.mil} \
        --mel ${params.mel} \
        --cstat ${params.cstat ?: '0.0001'} \
        --variable-read-length \
        --allow-clipping
    """
}
