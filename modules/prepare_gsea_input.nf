process PREPARE_GSEA_INPUT {
    tag "$meta.id"
    publishDir "${params.outdir}/gsea/input", mode: 'copy'

    input:
    tuple val(meta), path(deseq_results)

    output:
    tuple val(meta), path("${meta.id}.rnk"), emit: rnk
    
    script:
    """
    awk -F',' 'NR>1 && \$8!="NA" && \$4!="NA" {print \$8"\t"\$4}' ${deseq_results} > ${meta.id}.rnk
    """
}
