process SALMON_QUANT {
    tag "$sample_id"
    publishDir "${params.outdir}/salmon_quant", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(condition)
    path index

    output:
    tuple val(sample_id), path("${sample_id}"), val(condition), emit: results
    path "versions.yml", emit: versions

    script:
    def single_end = reads instanceof Path || reads.size() == 1

    if (single_end) {
        """
        salmon quant \\
            -i $index \\
            -l ${params.lib_type} \\
            -r ${reads} \\
            --validateMappings \\
            -o $sample_id \\
            -p $task.cpus \\
            --seqBias \\
            --gcBias \\
            --numBootstraps 100

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(salmon --version | sed 's/salmon //g')
        END_VERSIONS
        """
    } else {
        def read1 = reads[0]
        def read2 = reads[1]

        """
        # Ejecucion de cuantificacion
        salmon quant \\
            -i $index \\
            -l ${params.lib_type} \\
            -1 $read1 \\
            -2 $read2 \\
            --validateMappings \\
            -o $sample_id \\
            -p $task.cpus \\
            --seqBias \\
            --gcBias \\
            --numBootstraps 100

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(salmon --version | sed 's/salmon //g')
        END_VERSIONS
        """
    }
}
