process SAMTOOLS_SORT_INDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/samtools", mode: 'copy'
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam), val(condition)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), val(condition), emit: sorted_bam_and_index
    path "versions.yml", emit: versions

    script:
    """
    # Orden e index de los archivos BAM
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam $bam
    samtools index ${sample_id}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1s/^.*samtools //p')
    END_VERSIONS
    """
}
