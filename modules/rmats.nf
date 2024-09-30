process RMATS_ORIGINAL_GTF {
   // cache false
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
    echo "Executing rMATS with original GTF"
    rmats.py --b1 ${bam1} \
             --b2 ${bam2} \
             --gtf ${original_gtf} \
             --od ./ \
             --tmp ./tmp_original \
             -t paired \
             --readLength ${params.readlength} \
             --nthread ${task.cpus} \
             --libType fr-firststrand \
             --novelSS \
             --mil ${params.mil} \
             --mel ${params.mel}
    """
}

process RMATS_MERGED_GTF {
    publishDir "${params.outdir}/rMATS_out/merged_gtf", mode: 'copy'
    label 'rmats'

    input:
    path merged_gtf
    path bam1
    path bam2

    output:
    path "*", emit: results

    script:
    """
    echo "Executing rMATS with merged GTF"
    echo "GTF file: ${merged_gtf}"
    echo "BAM1 file: ${bam1}"
    echo "BAM2 file: ${bam2}"
    
    rmats.py --b1 ${bam1} \
             --b2 ${bam2} \
             --gtf ${merged_gtf} \
             --od ./ \
             --tmp ./tmp_merged \
             -t paired \
             --readLength ${params.readlength} \
             --nthread ${task.cpus} \
             --libType fr-firststrand \
             --novelSS \
             --mil ${params.mil} \
             --mel ${params.mel}

    echo "rMATS execution completed. Output files:"
    ls -l
    """
}
