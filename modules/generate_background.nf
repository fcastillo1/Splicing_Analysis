process GENERATE_BACKGROUND {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(deseq2_results)

    output:
    tuple val(meta), path("${meta.id}_background.txt"), emit: background

    script:
    """
    # Genera los datos de background para deseq2
    awk -F ',' '{gsub(/"/, "", \$7); print \$7}' ${deseq2_results} > ${meta.id}_background.txt
    """
}
